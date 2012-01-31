///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "VecUtilities.h"
#include "efasta/EfastaTools.h"
#include "paths/AssemblyCleanupTools.h"


Bool scompare( const superb& s1, const superb& s2 ){
  int l1 = s1.FullLength();
  int l2 = s2.FullLength();
  return (l1 > l2);
}


Assembly::Assembly( const String in_superb_file, const String in_efasta_file ){
  // Loading scaffolds
  cout << Date( ) << ": loading superb file" << endl;
  ReadSuperbs( in_superb_file, scaffolds );
  
  // reading contig information
  LoadEfastaIntoStrings(in_efasta_file, efastas);
  tigMap.resize( efastas.size() );
  for ( int tid = 0; tid < efastas.isize(); tid++ )
    tigMap[tid] = ToString(tid);

  scaffMap.resize( scaffolds.size() );
  for ( int sid = 0; sid < scaffolds.isize(); sid++ )
    scaffMap[sid] = ToString(sid);
}


Assembly::Assembly( const vec<superb>& scaffoldsIn, const vec<efasta>& efastasIn ){

  scaffolds = scaffoldsIn;
  efastas   = efastasIn;
  tigMap.resize( efastas.size() );
  for ( int tid = 0; tid < efastas.isize(); tid++ )
    tigMap[tid] = ToString(tid);

  scaffMap.resize( scaffolds.size() );
  for ( int sid = 0; sid < scaffolds.isize(); sid++ )
    scaffMap[sid] = ToString(sid);
}


size_t Assembly::scaffoldsTotLen() const{
  size_t len = 0;
  for ( size_t is = 0; is < scaffolds.size(); is++ )
    len += scaffolds[is].FullLength();
  return len;
}

size_t Assembly::scaffoldsRedLen() const{
  size_t len = 0;
  for ( size_t is = 0; is < scaffolds.size(); is++ )
    len += scaffolds[is].ReducedLength();
  return len;
}

size_t Assembly::scaffoldsNtigs() const{
  size_t ntigs = 0;
  for ( size_t is = 0; is < scaffolds.size(); is++ )
    ntigs += scaffolds[is].Ntigs();
  return ntigs;
}

// check integrity of scafolds and contigs data: contig size in superb == contig size in efasta,
//  each contig used once and only once
void Assembly::check_integrity() const{
  
  cout << Date() << ": checking integrity" << endl;
  vec<int> tigs_used( efastas.size(), 0);
  for ( size_t i = 0; i < efastas.size( ); i++ )
  {    vec<String> s(1);
       s[0] = efastas[i];
       ValidateEfastaRecord(s);    }
  for ( size_t si = 0; si < scaffolds.size(); si++ ){
    const superb & s = scaffolds[si];
    ForceAssertGt( s.Ntigs(), 0 );
    for ( int tpos = 0; tpos < s.Ntigs(); tpos++ ){
      size_t tid = s.Tig(tpos);
      ForceAssertLt( tid, efastas.size() );
      if ( efastas[tid].Length1() != s.Len(tpos) ){
	PRINT5( si, tpos, tid, s.Len(tpos), efastas[tid].Length1() );
	ForceAssertEq( efastas[tid].Length1(), s.Len(tpos) );
      }
      tigs_used[tid]++;
    }
  }
  vec<size_t> unused_tigs, overused_tigs;
  for ( size_t tid = 0; tid < tigs_used.size(); tid++ ){
    if ( tigs_used[tid] == 0 )
      unused_tigs.push_back( tid );
    else if ( tigs_used[tid] > 1 )
      overused_tigs.push_back(tid);
    
  }
  
  if ( unused_tigs.size() > 0 || overused_tigs.size() > 0 ){
    
    if ( unused_tigs.size() > 0 ){
      int max_un_size = efastas.at( unused_tigs[0] ).Length1();
      for ( size_t i = 0; i < unused_tigs.size(); i++ )
	if (  efastas.at( unused_tigs[i] ).Length1() > max_un_size )
	  max_un_size = efastas.at( unused_tigs[i] ).Length1();

      cout << "maximum size of unused contig is : " << max_un_size << endl;
    }

    PRINT2( unused_tigs.size(), overused_tigs.size() );
    ForceAssert( unused_tigs.size() == 0 && overused_tigs.size() == 0 );
  }
  
  return;
}


void Assembly::remove_small_scaffolds(const int min_scaffold_size) {
  cout << Date() << " removing small scaffolds: " << endl;
  cout << "initial number of scaffolds = " << scaffolds.size() << endl;
  for ( int si = 0; si < scaffolds.isize(); si++ )
    if ( scaffolds.at(si).ReducedLength() < min_scaffold_size ){
      scaffolds.erase( scaffolds.begin() + si );
      scaffMap.erase( scaffMap.begin() + si );
      si--;      
    }
  for ( int si = 0; si < scaffolds.isize(); si++ )
    ForceAssertGe( scaffolds.at(si).ReducedLength(), min_scaffold_size );
  cout << "final number of scaffolds = " << scaffolds.size() << endl;
}

void Assembly::remove_contigs( const vec<Bool>& to_remove )
{
  vec<int> offsets( efastas.size(), 0 );

  int offset = 0;
  for ( size_t tid = 0; tid < efastas.size(); tid++ ){
    if ( to_remove[tid] ){
      offsets[tid] = -1;
      offset++;
    }
    else{ offsets[tid] = offset; }
  }
    
  ForceAssertEq( offsets.size(), efastas.size() );
  for ( int tid = 0; tid < offsets.isize(); tid++ ){
    if ( offsets[tid] > 0 ){
      efastas[ tid - offsets[tid] ] = efastas[tid];
      tigMap[tid - offsets[tid] ] = tigMap[tid];
    }      
  }
  efastas.resize( efastas.size() - offset );
  tigMap.resize( tigMap.size() - offset );
  ForceAssertEq( efastas.size(), tigMap.size() );
  
  for ( size_t si = 0; si < scaffolds.size(); si++ ){
    for ( int tpos = 0; tpos < scaffolds[si].Ntigs(); tpos++ ){
      int tid = scaffolds[si].Tig(tpos);
      if ( offsets[tid] >= 0 ){
	int newtid = tid - offsets[tid];
	ForceAssertGe( newtid, 0 );
	scaffolds[si].SetTig( tpos, newtid );
      }
      else{
	scaffolds[si].RemoveTigByPos( tpos );
	tpos--;
      }
    }
  }
}

   
void Assembly::remove_small_contigs( const int min_size_solo, const int min_size_in ){
  cout << Date() << " removing small contigs and renumbering" << endl;
  vec<int> offsets( efastas.size(), 0 );
  vec<Bool> tigsToRemove( efastas.size(), False);
  for ( size_t si = 0; si < scaffolds.size(); si++ ){
    if ( scaffolds[si].Ntigs() == 1 ){
      if ( scaffolds[si].Len(0) < min_size_solo )
	tigsToRemove[ scaffolds[si].Tig(0) ] = True;
    }
    else{
      for ( int tpos = 0; tpos < scaffolds[si].Ntigs(); tpos++ ){
	if ( scaffolds[si].Len(tpos) < min_size_in )
	  tigsToRemove[ scaffolds[si].Tig(tpos) ] = True;
      }
    }
  }
  remove_contigs(tigsToRemove);
}

void Assembly::remove_unused_contigs(){
  cout << Date() << ": removing unused contigs and renumbering" << endl;
   vec<int> tigs_used( efastas.size(), 0);
  for ( size_t si = 0; si < scaffolds.size(); si++ ){
    const superb & s = scaffolds[si];
    ForceAssertGt( s.Ntigs(), 0 );
    for ( int tpos = 0; tpos < s.Ntigs(); tpos++ ){
      size_t tid = s.Tig(tpos);
      ForceAssertLt( tid, efastas.size() );
      if ( efastas[tid].Length1() != s.Len(tpos) ){
	PRINT5( si, tpos, tid, s.Len(tpos), efastas[tid].Length1() );
	ForceAssertEq( efastas[tid].Length1(), s.Len(tpos) );
      }
      tigs_used[tid]++;
    }
  }

  size_t unusedCt = 0;
  vec<int> offsets( efastas.size(), 0 );
  int offset = 0;
  for ( size_t tid = 0; tid < efastas.size(); tid++ ){
    if ( ! tigs_used[tid] ){
      offsets[tid] = -1;
      offset++;
      unusedCt++;
    }
    else{ offsets[tid] = offset; }
  }
  cout << Date( ) << ": found " << unusedCt << " unused contigs, removing" << endl;
  ForceAssertEq( offsets.size(), efastas.size() );
  for ( int tid = 0; tid < offsets.isize(); tid++ ){
    if ( offsets[tid] > 0 ){
      efastas[ tid - offsets[tid] ] = efastas[tid];
      tigMap[tid - offsets[tid] ] = tigMap[tid];
    }      
  }
  efastas.resize( efastas.size() - offset );
  tigMap.resize( tigMap.size() - offset );
  ForceAssertEq( efastas.size(), tigMap.size() );
  
  cout << Date() << ": updating scaffolds tig ids" << endl;
  for ( size_t si = 0; si < scaffolds.size(); si++ ){
    for ( int tpos = 0; tpos < scaffolds[si].Ntigs(); tpos++ ){
      int tid = scaffolds[si].Tig(tpos);
      if ( offsets[tid] >= 0 ){
	int newtid = tid - offsets[tid];
	ForceAssertGe( newtid, 0 );
	scaffolds[si].SetTig( tpos, newtid );
      }
      else{
	scaffolds[si].RemoveTigByPos( tpos );
	tpos--;
      }
    }
  }
  cout << Date() << ": done with removing unused contigs" << endl;
}


void Assembly::reorder(){
  // Sorting scaffolds according to size and renumbering contigs according
  //   to sequential appearance in scaffolds
  cout << Date() << " sorting scaffolds" << endl;
  SortSync( scaffolds, scaffMap, scompare );
  renumber();
}



// renumber all the contigs sequentially according to the scaffold
void Assembly::renumber(){
  cout << Date() << " renumbering contigs for ordered scaffolds" << endl;
  vec<String> otigMap;  
  vec<superb> oscaffolds = scaffolds;
  vec<efasta> oefastas;
  int cTid = -1;
  for ( size_t is = 0; is < scaffolds.size(); is++ ){
    for ( int tpos = 0; tpos < scaffolds[is].Ntigs(); tpos++ ){
      cTid++;
      int oTid = scaffolds[is].Tig(tpos);
      oscaffolds.at(is).SetTig( tpos, cTid );
      oefastas.push_back( efastas.at(oTid) );
      otigMap.push_back( tigMap.at(oTid) );
    }
  }
  efastas.resize(0);
  tigMap.resize(0);
  size_t newScaffoldsTotLen = 0, newScaffoldsRedLen = 0;
  for ( size_t is = 0; is < scaffolds.size(); is++ ){
    newScaffoldsTotLen += scaffolds[is].FullLength();
    newScaffoldsRedLen += scaffolds[is].ReducedLength();
  }
  scaffolds = oscaffolds; 
  oscaffolds.resize(0);
  efastas = oefastas; 
  oefastas.resize(0);
  tigMap = otigMap;
}

void Assembly::dedup() {
  // Remove duplicate scaffolds...currently only handles case of singleton contigs
  // which are duplciates fw or rc. --bruce 8 Jun 2011
  cout << Date() << " removing duplicate scaffolds: " << endl;
  cout << "initial number of scaffolds = " << scaffolds.size() << endl;
  int removed_count = 0;
  int removed_size = 0;
  for ( int si = 0; si < scaffolds.isize(); si++ ) {
    if (scaffolds[si].Ntigs() != 1) continue;
    for ( int sj = si + 1; sj < scaffolds.isize(); sj++ ) {
      if (scaffolds[sj].Ntigs() != 1) continue;
      efasta &tigi = efastas[scaffolds[si].Tig(0)];
      efasta &tigj = efastas[scaffolds[sj].Tig(0)];
      if (tigi.size() != tigj.size()) continue;
      basevector ibases, jbases;
      tigi.FlattenTo(ibases);
      tigj.FlattenTo(jbases);
      bool rc = False;
      if (ibases != jbases) {
	jbases.ReverseComplement();
	rc = True;
	if (ibases != jbases) continue;
      }
      /*
      cout << "scaffold " << sj << " duplicate of " << si
	   << " (length " << scaffolds[si].FullLength()
	   << (rc ? " rc" : "") << ")" << endl;
      */
      scaffolds.erase( scaffolds.begin() + sj );
      scaffMap.erase( scaffMap.begin() + sj );
      sj--;
      removed_count++;
      removed_size += ibases.size();
    }
  }
  cout << "removed " << removed_count << " duplicate scaffolds"
       << " (" << removed_size << " bases)" << endl;
  cout << "final number of scaffolds = " << scaffolds.size() << endl;
  remove_unused_contigs();
}


void Assembly::Write( const String head_out ) const {

  // writing output
  cout << Date() << ": writing output files" << endl;
  WriteSuperbs( head_out + ".superb", scaffolds );
  WriteSummary( head_out + ".summary", scaffolds );


  
  Ofstream( efout, head_out + ".contigs.efasta" );
  for ( size_t id = 0; id < efastas.size(); id++ )
    efastas[id].Print(efout, "contig_" + ToString(id) );
}

void Assembly::WriteExtra( const String head_out ) const{

  vec<fastavector> fastas(efastas.size());
  for ( size_t id = 0; id < efastas.size(); id++ )
    efastas[id].FlattenTo( fastas[id] );
  Ofstream( fout, head_out + ".contigs.fasta" );
  for ( size_t id = 0; id < fastas.size(); id++ )
    fastas[id].Print(fout, "contig_" + ToString(id) );

  {
    vecfastavector vec_tigs;
    for ( size_t i = 0; i < fastas.size( ); i++ )
      vec_tigs.push_back_reserve( fastas[i] );
    vec_tigs.WriteAll( head_out + ".contigs.vecfasta" ); 
  }

  vecbasevector bases( efastas.size() );
  for ( size_t id = 0; id < efastas.size(); id++ )
    efastas[id].FlattenTo( bases[id] );
  
  bases.WriteAll( head_out + ".contigs.fastb" );
  
  Ofstream( cmout, head_out + ".contigs.mapping" );
  for ( size_t id = 0; id < tigMap.size(); id++ )
    cmout << ToString(id) + " from " + ToString( tigMap[id] ) << "\n";

  Ofstream( smout, head_out + ".superb.mapping" );
  for ( size_t is = 0; is < scaffMap.size(); is++ )
    smout << ToString(is) + " from " + ToString( scaffMap[is] ) << "\n";

  WriteScaffoldedEFasta( head_out + ".assembly.efasta", efastas, scaffolds );
  WriteScaffoldedFasta( head_out + ".assembly.fasta", fastas, scaffolds );

}
