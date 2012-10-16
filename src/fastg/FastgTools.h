///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef FASTG_TOOLS_H
#define FASTG_TOOLS_H


#include "CoreTools.h"
#include "efasta/EfastaTools.h"


// define version

namespace FastgVersion
{ 
  static const String fastg_current_version = "0.10"; 
}


// class for reporting meta data and formatted headers and footers

class fastg_meta{

 public:
  fastg_meta( ) : version_( FastgVersion::fastg_current_version ) { }
    
  String GetVersion();
  String GetFileHeader(const String assembly_name);
  String GetFileFooter();
  
 private:
  const String version_;

};


class recfastg;

// Class basefastg represents bases record in fastg string.

class basefastg : public String {
   
 public:
  
  basefastg( );
  basefastg( const String& s );
  basefastg( const efasta& ef );
  basefastg( const int& sep, const int& dev );
  basefastg( const superb& s, const vec<efasta>& econtigs);
  basefastg( const superb& s, const vec<basefastg>& fcontigs);
  basefastg( const superb& s, const vec<recfastg>& fcontigs);
  
  void AsSections( vec<String>& sections ) const;

  int Length1() const;
  int MinLength() const;
  int MaxLength() const;
  Bool IsGapless() const;
};


class headfastg : public String {
   
 public:
  
  headfastg( );
  headfastg( const String& id );
  headfastg( const String& id, const vec<String>& next_ids );
  headfastg( const String& id, const vec<int>& next_ids );

};

class recfastg {

 private:
  headfastg header_;
  basefastg bases_;

 public:
  recfastg();
  recfastg( const headfastg& header, const basefastg& bases); 


  void Set( const headfastg& header, const basefastg& bases);

  int Length1() const;
  int MinLength() const;
  int MaxLength() const;
  Bool IsGapless() const;

  Bool ReadRecord( ifstream& in );

  const basefastg& bases() const;
  const headfastg& header() const;

  // Prints in a following format: "><header_>;\n" followed by the full bases
  // sequence and ambiguity information; breaks
  // long sequences nicely into 80-character lines

  void Print( ostream& out ) const;
  
  void Print( ostream& out, const String& id ) const;

  void AsScaffold( superb& s, vec<recfastg>& fcontigs, int& lastTid ) const;

  void AsFasta( fastavector& fa ) const;
  void AsFastb( basevector& fb ) const;
  void AsEfasta( efasta& fe ) const;

};


void LoadFastg( const String& fn, vec<recfastg>& records );

void WriteFastg( const String& fn, const vec<recfastg>& records );




#endif
