///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#include "TokenizeString.h"
#include "FastIfstream.h"
#include "fastg/FastgTools.h"
#include "efasta/EfastaTools.h"

#define Err(message)                                      \
{    cout << message << "\n" << endl;   \
     TracebackThisProcess( );    }




String fastg_meta::GetVersion(){ return version_; }

String fastg_meta::GetFileHeader( const String assembly_name = "" ){
    return "# begin:FASTG:version=" + version_ + 
      ":assembly_name=\"" + assembly_name + "\";";
}
 
String fastg_meta::GetFileFooter(){
  return "# end:FASTG;";
}


// helper functions

Bool areConsecutive( const vec<int>& slengths ){
  Bool Consecutive = True;
  for ( size_t li = 1; li < slengths.size(); li++ )
    if ( slengths[li] > slengths[li-1] + 1 ){
      Consecutive = False;
      break;
    }
  return Consecutive;
}


Bool areTandem( const vec<String>& alts, String& record ){
  record.clear();
  if ( alts.size() < 2 || 
       ( alts.size() == 2 && (alts[0] == "" || alts[1] == "") )   )
    return False;

  vec<int> slengths( alts.size() );
  for ( size_t ia = 0; ia < alts.size(); ia++ )
    slengths[ia] = alts[ia].isize();
  Sort( slengths );
  
  Bool isTandem = False;
  
  int minPositSize = slengths.front() == 0 ? slengths.at(1) : slengths.at(0);
  for ( int tsize = 1; tsize <= minPositSize; tsize++ ){
    Bool isFitting = True;
    for ( size_t il = 0; il < slengths.size(); il++ ){
      if ( slengths[il] != 0 && slengths[il] % tsize != 0 ){
	isFitting = False;
	break;
      }
    }
    if ( isFitting ){
      Bool FoundDifferent = False;
      String firstTandem = alts[0].size() > 0 ? alts[0].substr(0,tsize) : alts[1].substr(0,tsize);
      for ( size_t ia = 0; ia < alts.size(); ia++ ){
	if ( alts[ia].size() == 0 ) continue;
	for ( size_t is = 0; is * tsize < alts[ia].size(); is += tsize )
	  if ( alts[ia].substr(is, tsize) != firstTandem ){
	    FoundDifferent = True;
	    break;
	  }
	if ( FoundDifferent ) break;
      }

      if (! FoundDifferent ){
	isTandem = True;
	vec<int> tlengths( slengths.size() );
	for ( size_t il = 0; il < slengths.size(); il++ )
	  tlengths[il] = slengths[il]/tsize;
	
	record += "tandem";
	if ( areConsecutive( tlengths ) ){
	  record += ":size=(";
	  record += ToString( alts.front().isize()/tsize ) + "," +
	    ToString(slengths.front()/tsize ) + ".." + ToString(slengths.back()/tsize ) + ")|";
	}else{
	  isTandem = False;
	  break;
	}
	record += firstTandem + "]";
	break;
      }    
    }
  }

  return isTandem;
}




basefastg::basefastg( ) { }

basefastg::basefastg( const String& s ) : String(s) { }

basefastg::basefastg( const efasta& ef ) {
  (*this).clear();

  vec< vec<String> > blocks;
  ef.GetBlocks( blocks );
  
  for ( size_t ib = 0; ib < blocks.size(); ib++ ){
    if ( blocks[ib].size() == 1 ){
      *this += blocks[ib][0];
    }
    else if ( blocks[ib].size() > 1 ){
      int csize = blocks[ib][0].size();
      *this += "[" + ToString( csize ) + ":";    

      String entry = "";
      if ( areTandem( blocks[ib], entry ) ){
	*this += entry;
      }else{
	*this += "alt|";
	for ( size_t ia = 0; ia < blocks[ib].size(); ia++ ){
	  *this += blocks[ib][ia];
	  if ( ia + 1 < blocks[ib].size() ) *this += ",";
	}
	*this += "]";
      }	
    }else
      Err("unknown efasta format");
    
  }
  ForceAssertEq( (*this).Length1(), ef.Length1() );
  ForceAssertEq( (*this).MinLength(), ef.MinLength() );
  ForceAssertEq( (*this).MaxLength(), ef.MaxLength() );
}


basefastg::basefastg( const int& sep, const int& dev ){
  (*this).clear();
  *this += 
    "[" + ToString(sep) + ":gap:size=(" + ToString(sep) + "," + ToString( sep - 3 * dev ) + ".." + 
    ToString( sep + 3 * dev ) + ")]";
}

basefastg::basefastg( const superb& s, const vec<efasta>& econtigs ){
  (*this).clear();
  for ( int it = 0; it < s.Ntigs(); it++) {
    int tigId = s.Tig( it );
    (*this) += basefastg( econtigs[tigId] );
    if ( it < s.Ntigs() -1 )
      (*this) += basefastg( s.Gap(it), s.Dev(it) );
  } 
}

basefastg::basefastg( const superb& s, const vec<basefastg>& fcontigs_bases ){
  (*this).clear();
  for ( int it = 0; it < s.Ntigs(); it++) {
    int tigId = s.Tig( it );
    (*this) += fcontigs_bases[tigId];
    if ( it < s.Ntigs() -1 )
      (*this) += basefastg( s.Gap(it), s.Dev(it) );
  } 
}

basefastg::basefastg( const superb& s, const vec<recfastg>& fcontigs ){
  (*this).clear();
  for ( int it = 0; it < s.Ntigs(); it++) {
    int tigId = s.Tig( it );
    (*this) += fcontigs[tigId].bases();
    if ( it < s.Ntigs() -1 )
      (*this) += basefastg( s.Gap(it), s.Dev(it) );
  } 
}


void basefastg::AsSections( vec<String>& sections ) const {
  sections.clear();
  int nSections = 0;
  for ( int p = 0; p < (*this).isize(); p++ )
    if ( (*this)[p] == '[' ) nSections++;
  sections.reserve(nSections);
  String section;
  for ( int p = 0; p < (*this).isize(); p++ ){
    if ( (*this)[p] == '[' ){
      if ( section.nonempty() ){
	sections.push_back( section );
	section = "";
      }
      section += '[';
    }else if ( (*this)[p] == ']' ){
      section += ']';
      sections.push_back( section );
      section = "";
    }else if ( p == (*this).isize() -1 ){
      section += (*this)[p];
      sections.push_back( section );
      section = "";
    } else
      section += (*this)[p];
  }
  
  String sumb = "";
  int sumlen = 0;
  for ( size_t i = 0; i < sections.size(); i++ ){
    sumlen += sections[i].isize();
    sumb += sections[i];
    ForceAssertGt( sections[i].size(), 0u );
  }

  ForceAssertEq( (*this).isize(), sumlen );
  ForceAssertEq( *this, sumb );
}

int basefastg::Length1() const {
  int length1 = 0;
  vec<String> sections;
  AsSections( sections );
  for ( size_t i = 0; i < sections.size(); i++ ){
    if ( ! sections[i].Contains("[") )
      length1 += sections[i].isize();
    else
      length1 += atoi( sections[i].After("[").Before(":").c_str() );
  }
  return length1;
}

int basefastg::MaxLength() const {
  int maxLen = 0;
  vec<String> sections;
  AsSections( sections );
  for ( size_t is = 0; is < sections.size(); is++ ){
    if ( ! sections[is].Contains("[") )
      maxLen += sections[is].isize();
    else if ( sections[is].Contains("gap") ) {
      maxLen += atoi( sections[is].After("..").Before(")").c_str() );
    }else if ( sections[is].Contains("tandem") ){
      String srepeat = sections[is].After("|").Before("]");
      int maxRep = atoi( sections[is].After("..").Before(")").c_str() );
      maxLen += srepeat.isize() * maxRep;
    }else if ( sections[is].Contains("alt") ) {
      vec<char> separators; separators.push_back(',');
      vec<String> svalues;
      TokenizeStrictly( sections[is].After("|").Before("]"), separators, svalues );
      vec<int> ivalues( svalues.size() );
      for ( size_t k = 0; k < svalues.size(); k++ )
	ivalues[k] = svalues[k].isize();
      maxLen += Max(ivalues);
    }else
      Err("unknown format");
    
  }
  return maxLen;
}


int basefastg::MinLength() const {
  int minLen = 0;
  vec<String> sections;
  AsSections( sections );
  for ( size_t is = 0; is < sections.size(); is++ ){
    if ( ! sections[is].Contains("[") )
      minLen += sections[is].isize();
    else if ( sections[is].Contains("gap") ) {
      minLen += atoi( sections[is].Before("..").After(",").c_str() );
    }else if ( sections[is].Contains("tandem") ){
      String srepeat = sections[is].After("|").Before("]");
      int minRep = atoi( sections[is].Before("..").After(",").c_str() );
      minLen += srepeat.isize() * minRep;
    }else if ( sections[is].Contains("alt") ) {
      vec<char> separators; separators.push_back(',');
      vec<String> svalues;
      TokenizeStrictly( sections[is].After("|").Before("]"), separators, svalues );
      vec<int> ivalues( svalues.size() );
      for ( size_t k = 0; k < svalues.size(); k++ )
	ivalues[k] = svalues[k].isize();
      minLen += Min(ivalues);
    }else
      Err("unknown format");
    
  }
  return minLen;
}


Bool basefastg::IsGapless() const {
  return ! (*this).Contains("gap");
}


headfastg::headfastg( ) { }

headfastg::headfastg( const String& id ) : String(id) { }

headfastg::headfastg( const String& id, const vec<String>& next_ids ){
  *this = id;
  if ( next_ids.size() > 0 ){
    *this += ":";
    for ( int i = 0; i < next_ids.isize(); i++ ){
      *this += next_ids[i];
      if ( i < next_ids.isize() -1 )
	*this += ",";
    } 
  }
}

headfastg::headfastg( const String& id, const vec<int>& next_ids ){
  *this = id;
  if ( next_ids.size() > 0 ){
    *this += ":";
    for ( int i = 0; i < next_ids.isize(); i++ ){
      *this += ToString(next_ids[i]);
      if ( i < next_ids.isize() -1 )
	*this += ",";
    } 
  }
  
}

recfastg::recfastg( ){ }

recfastg::recfastg( const headfastg& header, const basefastg& bases){ 
  header_ = header;
  bases_  = bases;
}

void recfastg::Set( const headfastg& header, const basefastg& bases ){
  header_ = header;
  bases_  = bases;
}

int recfastg::Length1() const {
  return bases_.Length1();
}

int recfastg::MaxLength() const {
  return bases_.MaxLength();
}

int recfastg::MinLength() const {
  return bases_.MinLength();
}

Bool recfastg::IsGapless() const{
  return bases_.IsGapless();
}

const basefastg& recfastg::bases() const { return bases_; }
const headfastg& recfastg::header() const { return header_; }

void recfastg::Print( ostream& out ) const {
  String id = header_;
  String totrim =">,";
  id.TrimInPlace( totrim.c_str() );
  out << ">" << id << ";";

  for ( size_type i = 0; i < bases_.size(); ++i )
    {
      if ( !(i % 80) ) out << '\n';
      out << bases_[i];
    }
  out << '\n';
}

void recfastg::Print( ostream& out, const String& id ) const {
  out << ">" << id << ";";
  for ( size_type i = 0; i < bases_.size(); ++i )
    {
      if ( !(i % 80) ) out << '\n';
      out << bases_[i];
    }
  out << '\n';
}

void recfastg::AsScaffold( superb& s, vec<recfastg>& fcontigs, int& lastTid ) const {
  vec<String> sections;
  bases_.AsSections( sections );
  fcontigs.clear();
  s.SetNtigs( sections.size() ); // will be trimmed later
  basefastg tbases;
  for ( int i = 0; i < sections.isize(); i++ ){
    String& section = sections[i];
    if ( ! section.Contains("gap") ){
      tbases += section;
    }else{ 
      ForceAssertNe( i, 0 ); 
      ForceAssertNe( i, sections.isize() );
      if ( ! tbases.IsGapless() )
	cout << tbases << endl;
      ForceAssert( tbases.IsGapless() );
      lastTid++;
      fcontigs.push_back( recfastg( ToString(lastTid), tbases ) );
      s.SetTig( fcontigs.isize() -1, lastTid );
      s.SetLen( fcontigs.isize() -1, tbases.Length1() );
      tbases.clear();

      String gapdef = section.Between( "gap:size=(", ")" );
      int gapsize = atoi( gapdef.Before(",").c_str() );
      gapdef = gapdef.After(",");
      int min = atoi( gapdef.Before("..").c_str() );
      int max = atoi( gapdef.After("..").c_str() );
      s.SetGap( fcontigs.isize() -1, gapsize );
      s.SetDev( fcontigs.isize() -1, round( double(max - min + 1.0) / 6.0 ) );
    }
  }
  if ( tbases.nonempty() ){
    ForceAssert( tbases.IsGapless() );
    lastTid++;
    fcontigs.push_back( recfastg( ToString(lastTid), tbases ) );
    s.SetTig( fcontigs.isize() -1, lastTid );
    s.SetLen( fcontigs.isize() -1, tbases.Length1() );
    tbases.clear();
  }
  s.SetNtigs( fcontigs.isize() );
}


void recfastg::AsFasta( fastavector& fa ) const {
  efasta fe;
  AsEfasta( fe );
  fe.FlattenTo( fa );
}

void recfastg::AsFastb( basevector& fb ) const {
  fastavector fa;
  (*this).AsFasta( fa );
  fb = fa.ToBasevector();
}

void recfastg::AsEfasta( efasta& fe ) const {
  fe.clear();
  vec<char> separators; separators.push_back(',');
  vec<String> sections;
  int minGapSize = 0, maxGapSize = 0;
  bases_.AsSections( sections );
  for ( size_t is = 0; is < sections.size(); is++ ){

    if ( ! sections[is].Contains("[") )
      fe += sections[is];
    else{

      if ( sections[is].After("[").Contains("[") || sections[is].After("]").Contains("]") )
	Err("unknown formmat containing double square brackets");

      if ( sections[is].Contains(":gap") ){
	int gapsize = atoi( sections[is].After("[").Before(":").c_str() );
	for ( int i = 1; i <= gapsize; i++ ) 
	  fe += "N";
	minGapSize += atoi( sections[is].Before("..").After(",").c_str() );
	maxGapSize += atoi( sections[is].After("..").Before(")").c_str() );
      }else if ( sections[is].Contains(":alt") ) {
	int alen = atoi( sections[is].After("[").Before(":").c_str() );
	vec<String> svalues;
	TokenizeStrictly( sections[is].After("|").Before("]"), separators, svalues );
	fe += "{";
	ForceAssertEq( alen, svalues.front().isize() );
	for ( int isv = 0; isv < svalues.isize(); isv++ ){
	  fe += svalues[isv];
	  if ( isv < svalues.isize() -1 )
	    fe += ",";
	}
	fe += "}";
      }else if ( sections[is].Contains(":tandem") ) {
	int alen = atoi( sections[is].After("[").Before(":").c_str() );
	int nrepeat1 = atoi( sections[is].After("size=(").Before(",").c_str() );
	String srepeat = sections[is].After("|").Before("]");
	ForceAssertEq( alen, nrepeat1 * srepeat.isize() );

	vec<int> nreps; nreps.push_back( nrepeat1 );
	String srange = sections[is].After("size=(").Before(")");
	if ( srange.Contains("..") ){
	  srange = srange.After(",");
	  int min = atoi( srange.Before("..").c_str() );
	  int max = atoi( srange.After("..").c_str() );
	  ForceAssertLe( min, max );
	  for ( int nr = min; nr <= max; nr++ )
	    if ( nr != nrepeat1 )
	      nreps.push_back(nr);
	}else{
	  Err("wrong tandem format in: " + sections[is] ); 
	}

	fe += "{";
	for ( size_t inr = 0; inr < nreps.size(); inr++ ){
	  int nr = nreps[inr];
	  String scumul = "";
	  for ( int i = 1; i <= nr; i++ )
	    scumul += srepeat;
	  fe += scumul;
	  if ( inr +1 < nreps.size() )
	    fe += ",";
	}
	fe += "}";
      }else{
	Err("Err. converting to efasta: unknown format in " + sections[is] );
      }
    }
  }
 
  ForceAssertEq( bases_.Length1(), fe.Length1() );
  ForceAssertEq( bases_.MinLength(), fe.MinLength() + minGapSize );
  ForceAssertEq( bases_.MaxLength(), fe.MaxLength() + maxGapSize );
  
}


Bool recfastg::ReadRecord( ifstream& in ){
  header_.clear();
  bases_.clear();
  headfastg header;
  basefastg bases;
  String line;
  char c;

  if ( in.eof() || in.fail() ) return False;

  while( in.peek() != '>' ){
    in.get(c);
    if ( in.fail() || in.eof() )
      return False;
  }
  in.get(c);
  ForceAssertEq(c,'>');
  
  while( in.peek() != ';' ){
    in.get(c);
    if ( in.fail() || in.eof() )
      return False;
    if ( c != '\n' && c != ' ' && c != '\t' ) header += c;
  }
  in.get(c);
  header += c;
  ForceAssertEq(c,';');
  //PRINT(header);

  while( in.peek() != '>' && in.peek() != '#' ){
    in.get(c);
    if ( in.fail() || in.eof() )
      break;
    if ( c != '\n' && c != ' ' && c != '\t' ) bases += c;
  }
  if ( bases.Contains("[alt") )
    Err("Error: wrong format, contains \"[alt\"");
  //PRINT(bases);
  Set( header, bases );
  return True;
}



void LoadFastg( const String& fn, vec<recfastg>& records ){
  ifstream in( fn.c_str() );
  recfastg record;

  while( record.ReadRecord( in ) )
    records.push_back( record );
  
  cout << Date() << ": " << fn << " loaded" << endl;
}



void WriteFastg( const String& fn, const vec<recfastg>& records ){ 
  ofstream out( fn.c_str() );

  fastg_meta FGM;
  out << FGM.GetFileHeader( fn ) << "\n";

  for ( size_t i = 0; i < records.size(); i++ )
    records[i].Print( out );

  out << FGM.GetFileFooter() << "\n";
}
