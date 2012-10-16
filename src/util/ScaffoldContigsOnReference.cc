/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Basevector.h"
#include "Fastavector.h"
#include "FetchReads.h"
#include "Superb.h"
#include "lookup/LookAlign.h"
#include "paths/SaveScaffoldGraph.h"
#include "util/RunCommand.h"
#include "util/ScaffoldContigsOnRef.h"
// MakeDepend: dependency EvalScaffolds

/**
 * ScaffoldContigsOnReference
 *
 * Align contigs to a reference, and generate scaffolds.
 *
 * HEAD_REF: head of reference genome
 * MAX_GAP: max gap allowed between adjacent contigs
 * MAX_OVERLAP: max overlap allowed between adjacent contigs
 * MIN_GAP_DEV: min threshold for gap dev (usually set to gap_size / 4 )
 * FORCE: do not use cached aligns
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( ASSEMBLY_IN );
  CommandArgument_String( ASSEMBLY_OUT );
  CommandArgument_String( HEAD_REF );
  CommandArgument_Int_OrDefault( MAX_GAP, 2000 );
  CommandArgument_Int_OrDefault( MAX_OVERLAP, 1000 );
  CommandArgument_Int_OrDefault( MIN_GAP_DEV, 20 );
  CommandArgument_Bool_OrDefault( FORCE, False );
  EndCommandArguments;

  // Dir and file names.
  String tmp_dir = ASSEMBLY_OUT + ".tmp";

  String contigs_file = ASSEMBLY_IN + ".contigs.fasta";
  String supers_file = ASSEMBLY_IN + ".superb";
  String lookup_file = HEAD_REF + ".lookup";
  String aligns_file = tmp_dir + "/aligns.qlt";

  Mkpath( tmp_dir );

  // Load.
  vec<look_align> aligns;
  if ( FORCE || ! IsRegularFile( aligns_file ) ) {
    String comm
      = "EvalScaffolds LOOKUP=" + lookup_file
      + " SCAFFOLDS=" + ASSEMBLY_IN
      + " OUT_DIR=" + tmp_dir;
    RunCommand( comm );
  }
  LoadLookAligns( aligns_file, aligns );
  
  vec<fastavector> contigs;
  LoadFromFastaFile( contigs_file, contigs );
  
  cout << Date( ) << ": loading supers" << endl;
  vec<superb> initial_supers;
  ReadSuperbs( supers_file, initial_supers );
  
  // Scaffold contigs, and save (logging within).
  ScaffoldContigsOnRef( initial_supers, contigs, aligns, ASSEMBLY_OUT,
			MAX_GAP, MAX_OVERLAP, MIN_GAP_DEV, &cout );

  cout << Date( ) << ": ScaffoldContigsOnReference done" << endl;
  
}
