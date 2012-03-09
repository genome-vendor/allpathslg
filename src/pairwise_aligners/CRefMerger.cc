///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "Fastavector.h"
#include "PrintAlignment.h"
#include "Superb.h"
#include "feudal/BinaryStream.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "pairwise_aligners/CRefMerger.h"
#include "pairwise_aligners/SmithWatBandedA.h"

/**
 * class CRefMerger
 * Constructor
 */
CRefMerger::CRefMerger( const int min_overlap,
			const int swband_ratio,
			const vecbvec &ref ) :
  min_overlap_ ( min_overlap ),
  swband_ratio_ ( swband_ratio ),
  ref_ ( ref )
{ }

/**
 * class CRefMerger
 * Constructor
 */
CRefMerger::CRefMerger( const int min_overlap,
			const int swband_ratio,
			const vecbvec &ref,
			const vecbvec &in_contigs,
			const vec<look_align> &in_aligns ) :
  min_overlap_ ( min_overlap ),
  swband_ratio_ ( swband_ratio ),
  ref_ ( ref )
{
  this->Setup( in_contigs, in_aligns );
}
  
/**
 * class CRefMerger
 * Setup
 */
void CRefMerger::Setup( const vecbvec &in_contigs,
			const vec<look_align> &in_aligns )
{
  contigs_.reserve( in_aligns.size( ) );
  aligns_.reserve( in_aligns.size( ) );
  
  // Build contigs and aligns (one align per contig).
  for (size_t idx=0; idx<in_aligns.size( ); idx++) {
    look_align al = in_aligns[idx];
    ForceAssert( al.IsProper( ) );
    ForceAssert( al.Fw1( ) );
    
    int cid_orig = al.query_id;
    al.query_id = contigs_.size( );
    contigs_.push_back( in_contigs[cid_orig] );
    aligns_.push_back( al );
  }
  
  // Sort aligns by start on target.
  sort( aligns_.begin( ), aligns_.end( ), order_lookalign_TargetBegin( ) );
}

/**
 * class CRefMerger
 * Merge
 */
void CRefMerger::Merge( ostream *log )
{
  if ( log )
    *log << Date( ) << ": iter 0 (N50 = "
	 << this->ContigsN50( ) << ")"
	 << endl;
  
  int iter = 0;
  while ( 1 ) {
    int n_merges = this->MergeIteration( );
    iter++;
    if ( log )
      *log << Date( ) << ": iter "
	   << iter << " ("
	   << n_merges << " merges, N50 = "
	   << this->ContigsN50( ) << ")"
	   << endl;
    if ( n_merges < 1 ) break;
  }

}

/**
 * class CRefMerger
 * Save
 */
void CRefMerger::Save( const String &out_head,
		       const int min_clen,
		       ostream *log )
{
  String contigsa_file = out_head + ".contigs.fasta";
  String contigsb_file = out_head + ".contigs.fastb";
  String supers_file = out_head + ".superb";
  String aligns_file = out_head + ".qltout";
  
  this->CompactifyData( min_clen );

  if ( log ) *log << Date( ) << ": saving merged contigs (fasta)" << endl;
  Ofstream( cg_out, contigsa_file );
  for (size_t ii=0; ii<contigs_.size( ); ii++) {
    fastavector cgfasta( contigs_[ii] );
    cgfasta.Print( cg_out, "contig_" + ToString( ii ) );
  }
  cg_out.close();
  
  if ( log ) *log << Date( ) << ": saving merged contigs (fastb)" << endl;
  contigs_.WriteAll( contigsb_file );

  if ( log ) *log << Date( ) << ": saving (trivial) supers" << endl;
  vec<superb> supers( contigs_.size( ) );
  for (int ii=0; ii<(int)supers.size( ); ii++) {
    supers[ii].SetNtigs( 1 );
    supers[ii].SetTig( 0, ii );
    supers[ii].SetLen( 0, contigs_[ii].size( ) );
  }
  WriteSuperbs( supers_file, supers );

  if ( log ) *log << Date( ) << ": saving aligns on reference" << endl;
  WriteLookAligns( aligns_file, aligns_ );
}

/**
 * class CRefMerger
 * ContigsN50
 */
int CRefMerger::ContigsN50( ) const
{
  vec<int> lens;
  lens.reserve( contigs_.size( ) );
  for (int ii=0; ii<(int)contigs_.size( ); ii++)
    if ( contigs_[ii].size( ) > 0 )
      lens.push_back( (int)contigs_[ii].size( ) );
  sort( lens.begin( ), lens.end( ) );
  return N50( lens );
}

/**
 * class CRefMerger
 * MergeIteration
 * private
 */
int CRefMerger::MergeIteration( )
{
  typedef pair<double,int> weight;
  typedef pair<int,int> twins;

  // List and sort pairs of contigs with an implied overlap >= min_overlap_.
  vec< pair<weight,twins> > candidates;
  for (size_t ii=0; ii<aligns_.size( ); ii++) {
    for (size_t jj=ii+1; jj<aligns_.size( ); jj++) {
      if ( this->ImpliedOverlap( ii, jj ) < min_overlap_ ) break;
      double al_err = Max( aligns_[ii].ErrorRate( ), aligns_[jj].ErrorRate( ) );
      int len_ii = aligns_[ii].a.Pos2( ) - aligns_[ii].a.pos2( );
      int len_jj = aligns_[jj].a.Pos2( ) - aligns_[jj].a.pos2( );
      int al_len = Min( len_ii, len_jj );
      weight c_weight = make_pair( al_err, - al_len );
      twins c_twins = make_pair( ii, jj );
      candidates.push_back( make_pair( c_weight, c_twins ) );
    }
  }
  sort( candidates.begin( ), candidates.end( ) );
  
  // Merge pairs of overlapping contigs.
  int n_merges = 0;
  vec<bool> tagged( aligns_.size( ), false );
  for (size_t candid=0; candid<candidates.size( ); candid++) {
    const twins &selected = candidates[candid].second;
    int ii = selected.first;
    int jj = selected.second;
    if ( tagged[ii] || aligns_[ii].query_length < 1 ) continue;
    if ( tagged[jj] || aligns_[jj].query_length < 1 ) continue;
    if ( this->MergeContigs( ii, jj ) ) {
      tagged[ii] = true;
      tagged[jj] = true;
      n_merges++;
    }
  }

  // Compactify data, if needed.
  if ( n_merges > 0 ) this->CompactifyData( );

  // Done.
  return n_merges;
}

/**
 * class CRefMerger
 * MergeContigs
 * private
 */
bool CRefMerger::MergeContigs( const int idx1, const int idx2 )
{
  // Insufficient implied overlap.
  int overlap = ImpliedOverlap( idx1, idx2 );
  if ( overlap < min_overlap_ ) return false;

  // Second contig appears fully embedded, skip it. 
  const look_align &al1 = aligns_[idx1];
  const look_align &al2 = aligns_[idx2];
  if ( al1.a.pos2( ) <= al2.a.pos2( ) && al2.a.Pos2( ) <= al1.a.Pos2( ) )
    return false;

  // Align contigs.
  int offset = this->PosOn1( al2.a.pos2( ), al1.a );
  int band = Max( 1, overlap / swband_ratio_ );
  
  align al;
  int dummy;
  const bvec &b1 = contigs_[al1.query_id];
  const bvec &b2 = contigs_[al2.query_id];
  SmithWatBandedA2<unsigned int>( b1, b2, offset, band, al, dummy );
  if ( ! this->IsPerfect( b1, b2, al ) ) return false;
  
  // Merge contig from al2 onto contig from al1, and zero contig from al2.
  bvec frag1( b1, 0, al.pos1( ) );
  contigs_[al1.query_id] = Cat( frag1, b2 );
  contigs_[al2.query_id].resize( 0 );
  
  // Merge al1 and al2 onto al1, reset al2 by setting its query length to 0.
  align nAl1;
  int nOff = - ( al1.a.pos2( ) - al1.a.pos1( ) );
  int nBand = 0;
  for (int ii=0; ii<al1.a.Nblocks( ); ii++)
    nBand += Abs( al1.a.Gaps( ii ) );
  for (int ii=0; ii<al2.a.Nblocks( ); ii++)
    nBand += Abs( al2.a.Gaps( ii ) );
  
  const bvec &nB1 = contigs_[al1.query_id];
  const bvec &nB2 = ref_[al1.target_id];
  SmithWatBandedA2<unsigned int>( nB1, nB2, nOff, nBand, nAl1, dummy );
  
  #ifdef VERBOSE
  cout << "ALIGN between " << al1.query_id
       << " and " << al2.query_id << "  (offset = " << offset << ")\n"
       << " al1:   " << this->BriefAlignInfo( al1.a ) << "\n"
       << " al2:   " << this->BriefAlignInfo( al2.a ) << "\n"
       << " al1+2: " << this->BriefAlignInfo( nAl1 ) << "\n";
  PrintVisualAlignment( True, cout, nB1, nB2, nAl1 );
  #endif
  
  aligns_[idx1].query_length = (unsigned int)contigs_[al1.query_id].size( );
  aligns_[idx1].a = nAl1;
  aligns_[idx2].query_length = 0;
  
  // Align found.
  ForceAssert( this->IsValid( idx1, &cout ) );
  ForceAssert( this->IsValid( idx2, &cout ) );
  return true;
}
  
/**
 * class CRefMerger
 * CompactifyData
 * private
 */
void CRefMerger::CompactifyData( const int min_clen )
{
  vec<look_align> keeper_als;
  vecbvec keeper_cgs;
  keeper_als.reserve( aligns_.size( ) );
  keeper_cgs.reserve( contigs_.size( ) );

  // Keep only valid aligns (query_length > 0), and contig >= min_clen.
  for (int ii=0; ii<aligns_.isize( ); ii++) {
    if ( (int)aligns_[ii].query_length < Max( 1, min_clen ) ) {
      contigs_[ aligns_[ii].query_id ].resize( 0 );
      continue;
    }

    look_align al = aligns_[ii];
    int orig_id = al.query_id;
    al.query_id = (unsigned int)keeper_cgs.size( );
    keeper_als.push_back( al );
    keeper_cgs.push_back( contigs_[orig_id] );
  }
  swap( keeper_als, aligns_ );
  swap( keeper_cgs, contigs_ );
  
  // Internal consistency check.
  int n_errors = 0;
  for (int ii=0; ii<aligns_.isize( ); ii++)
    if ( ! this->IsValid( ii, &cout ) )
      n_errors++;
  ForceAssert( n_errors == 0 );
}

/**
 * CRefMerger
 * PosOn1
 * private
 *
 * NB: there is a PosOn1 member method for class align, but it returns
 * an error if pos2 corresponds to a gap. Here, if pos2 corresponds to
 * a gap, we slide off the gap, until the first base on 1.
 */
int CRefMerger::PosOn1( int pos2, const align &al ) const
{
  if ( pos2 < al.pos2( ) || pos2 > al.Pos2( ) ) return -1;
  
  int p1 = al.pos1( );
  int p2 = al.pos2( );
  for (int jj=0; jj<al.Nblocks( ); jj++) {
    int gap = al.Gaps( jj );
    int len = al.Lengths( jj );
    if ( gap > 0 ) {
      if ( p2 <= pos2 && p2 + gap >= pos2 ) return p1 + gap;
      p2 += gap;
    }
    if ( gap < 0 ) p1 += -gap; 
    for (int kk=0; kk<len; kk++) {
      if ( p2 == pos2 ) return p1;
      if ( p2 > pos2 ) return -1;
      p1++;
      p2++;
    }
  }

  return -1;
}

/**
 * CRefMerger
 * PosOn2
 * private
 *
 * NB (same as PosOn1): there is a PosOn2 member method for class
 * align, but it returns an error if pos1 corresponds to a gap. Here,
 * if pos1 corresponds to a gap, we slide off the gap, until the first
 * base on 2.
 */
int CRefMerger::PosOn2( int pos1, const align &al ) const
{
  if ( pos1 < al.pos1( ) || pos1 > al.Pos1( ) ) return -1;
  
  int p1 = al.pos1( );
  int p2 = al.pos2( );
  for (int jj=0; jj<al.Nblocks( ); jj++) {
    int gap = al.Gaps( jj );
    int len = al.Lengths( jj );
    if ( gap > 0 ) p2 += gap;
    if ( gap < 0 ) {
      if ( p1 <= pos1 && p1 - gap >= pos1 ) return p2 - gap;
      p1 += -gap;
    }
    for (int kk=0; kk<len; kk++) {
      if ( p1 == pos1 ) return p2;
      if ( p1 > pos1 ) return -1;
      p1++;
      p2++;
    }
  }
  
  return -1;
}

/**
 * class CRefMerger
 * IsPerfect
 * private
 */
bool CRefMerger::IsPerfect( const bvec &b1,
			    const bvec &b2,
			    const align &al ) const
{
  return al.Errors( b1, b2 ) < 1;
}

/**
 * class CRefMerger
 * ImpliedOverlap
 * private
 */
int CRefMerger::ImpliedOverlap( const int idx1, const int idx2 ) const
{
  const look_align &al1 = aligns_[idx1];
  const look_align &al2 = aligns_[idx2];
  if ( al1.target_id != al2.target_id ) return 0;
  
  return al1.a.Pos2( ) - al2.a.pos2( );
}

/**
 * class CRefMerger
 * IsValid
 * private
 */
bool CRefMerger::IsValid( const int idx, ostream *log ) const
{
  bool is_valid = true;

  // Contig has been merged.
  const look_align &al = aligns_[idx];
  const int qid = al.query_id;
  if ( contigs_[qid].size( ) < 1 || al.query_length < 1 ) {
    if ( ! ( contigs_[qid].size( ) < 1 && al.query_length < 1 ) ) {
      if ( log )
	*log << "al." << idx
	     << ": query_length is " << al.query_length
	     << ", but contig_[" << qid
	     << "] has length " << contigs_[qid].size( )
	     << " (improperly deleted align?)\n";
      return false;
    }
    return true;
  }
  
  // query_length and size of contig do not match.
  if ( (int)contigs_[qid].size( ) != (int)al.query_length ) {
    is_valid = false;
    if ( log )
      *log << "al." << idx
	   << ": query_length is " << al.query_length
	   << ", but contig_[" << qid
	   << "] has length " << contigs_[qid].size( )
	   << "\n";
  }
  
  // Invalid align.
  int pos1 = al.a.pos1( );
  if ( pos1 < 0 || pos1 > (int)al.query_length ) {
    is_valid = false;
    if ( log )
      *log << "al." << idx
	   << ": pos1 is " << pos1
	   << ", but query has length " << al.query_length
	   << "\n";
  }
  
  int Pos1 = al.a.Pos1( );
  if ( Pos1 < 0 || Pos1 > (int)al.query_length ) {
    is_valid = false;
    if ( log )
      *log << "al." << idx
	   << ": Pos1 is " << Pos1
	   << ", but query has length " << al.query_length
	   << "\n";
  }
  
  int pos2 = al.a.pos2( );
  if ( pos2 < 0 || pos2 > (int)al.target_length ) {
    is_valid = false;
    if ( log )
      *log << "al." << idx
	   << ": pos2 is " << pos2
	   << ", but target has length " << al.target_length
	   << "\n";
  }
  
  int Pos2 = al.a.Pos2( );
  if ( Pos2 < 0 || Pos2 > (int)al.target_length ) {
    is_valid = false;
    if ( log )
      *log << "al." << idx
	   << ": Pos2 is " << Pos2
	   << ", but target has length " << al.target_length
	   << "\n";
  }
  
  // Done.
  return is_valid;
}

/**
 * class CRefMerger
 * BriefAlignInfo
 * private
 */
String CRefMerger::BriefAlignInfo( const align &al ) const
{
  String str_al
    = "[p1, " + ToString( al.pos1( ) ) + "]  "
    + "[p2, " + ToString( al.pos2( ) ) + "]  ";
  for (int ii=0; ii<al.Nblocks( ); ii++)
    str_al
      += "[" + ToString( al.Gaps( ii ) )
      +  ", " + ToString( al.Lengths( ii ) )
      +  "]  ";
  str_al
    += "[P1, " + ToString( al.Pos1( ) ) + "]  "
    +  "[P2, " + ToString( al.Pos2( ) ) + "]";
  
  return str_al;
}

