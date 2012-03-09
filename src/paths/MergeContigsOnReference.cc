/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include <omp.h>
#include "CoreTools.h"
#include "MainTools.h"
#include "Basevector.h"
#include "lookup/LookAlign.h"
#include "pairwise_aligners/CRefMerger.h"
#include "paths/reporting/ReftigUtils.h"
#include "paths/SaveScaffoldGraph.h"
#include "util/RunCommand.h"
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

/**
 * MergeContigsOnReference
 *
 * Align contigs to a reference, and merge consecutive contigs, if
 * they overlap by >= MIN_OVERLAP bases (perfect matches only). If a
 * repetitive contig is involved in multiple overlaps, all are taken
 * into account. For example, if unique sequences A and B both align
 * repetitive sequence R, then we merge both A + R -> A', and B + R ->
 * B'.
 *
 * K: needed by the aligner (GetAlignsFast)
 * CONTIGS: full path name of input fastb
 * REF_HEAD: it loads <REF_HEAD>.{fastb,lookup}
 * OUT_DIR: full path name of output dir
 * OUT_HEAD: head name in OUT_DIR for output files
 * MIN_CLEN: only save contigs >= MIN_CLEN
 * MIN_OVERLAP: minimum (perfect) overlap required for merging
 * SWBAND_RATIO: sw band, defined as ( overlap / SWBAND_RATIO )
 * NUM_THREADS: use all available if 0
 * FW_ONLY: discard contigs that align rc on reference
 * FORCE: do not load cached aligns, regenerate them
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_Int( K );
  CommandArgument_String( CONTIGS );
  CommandArgument_String( REF_HEAD );
  CommandArgument_String( OUT_DIR );
  CommandArgument_String_OrDefault( OUT_HEAD, "merged" );
  CommandArgument_Int_OrDefault( MIN_CLEN, 1000 );
  CommandArgument_Int_OrDefault( MIN_OVERLAP, 12 );
  CommandArgument_Int_OrDefault( SWBAND_RATIO, 6 );
  CommandArgument_Int_OrDefault( NUM_THREADS, 0);
  CommandArgument_Bool_OrDefault( FW_ONLY, True );
  CommandArgument_Bool_OrDefault( FORCE, False );
  EndCommandArguments;

  // Thread control.
  NUM_THREADS = configNumThreads( NUM_THREADS );
  omp_set_num_threads( NUM_THREADS );

  // Dir and file names.
  String alignsFile = OUT_DIR + "/aligns.qlt";
  String outHead = OUT_DIR + "/" + OUT_HEAD;
  String lookupFile = REF_HEAD + ".lookup";
  String targetFile = REF_HEAD + ".fastb";
  
  // Needed.
  vec<String> needed;
  needed.push_back( CONTIGS );
  needed.push_back( lookupFile );
  needed.push_back( targetFile );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;
  
  Mkpath( OUT_DIR );
  
  // Get and filter aligns.
  vec<look_align> aligns;
  GetAlignsFast( K, CONTIGS, lookupFile, alignsFile, aligns, !FORCE, OUT_DIR );
  {
    vec<look_align> select;
    select.reserve( aligns.size( ) );
    for (size_t ii=0; ii<aligns.size( ); ii++) {
      if ( aligns[ii].Rc1( ) ) continue;
      if ( ! aligns[ii].IsProper( ) ) continue;
      select.push_back( aligns[ii] );
    }
    swap( select, aligns );
  }
  
  // Load.
  cout << Date( ) << ": loading contigs" << endl;
  vecbvec contigs( CONTIGS );

  cout << Date( ) << ": loading reference" << endl;
  vecbvec targets( targetFile );
  
  // Merge contigs.
  CRefMerger merger( MIN_OVERLAP, SWBAND_RATIO, targets, contigs, aligns );
  merger.Merge( &cout );
  merger.Save( outHead, MIN_CLEN, &cout );
  
  // Done.
  cout << Date( ) << ": MergeContigsOnReference done" << endl;
  
}
