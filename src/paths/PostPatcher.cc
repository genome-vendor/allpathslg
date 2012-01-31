///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// PostPatcher.  Close gaps in scaffolds, bootstrapping off alignments generated by
// UnipathPatcher.  Note that the arguments to this program need to match those
// passed to UnipathPatcher.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

// MakeDepend: dependency PatcherCottage
// MakeDepend: dependency PostPatcherBuildJoins
// MakeDepend: dependency QueryLookupTable

#include <omp.h>

#include "Basevector.h"
#include "Fastavector.h"
#include "FastIfstream.h"
#include "FetchReads.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "ParseSet.h"
#include "Qualvector.h"
#include "Superb.h"
#include "VecTemplate.h"
#include "pairwise_aligners/AlignTwoBasevectors.h"
#include "paths/AssemblyEdit.h"
#include "paths/HyperFastavector.h"
#include "paths/PostPatcherBridgeGap.h"
#include "paths/RegapSupers.h"
#include "paths/ScaffoldsUtils.h"
#include "paths/UnipathFixerTools.h"
#include "util/SearchFastb2Core.h"


int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String(SUBDIR);
     CommandArgument_Int(K);
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
	   "Number of threads to use (use all available processors if set to 0)");
     CommandArgument_String_OrDefault(SCAFFOLDS_IN, "linear_scaffolds");
     CommandArgument_String_OrDefault(SCAFFOLDS_OUT, SCAFFOLDS_IN + ".patched");
     CommandArgument_String_OrDefault(HEAD, "all_reads");
     CommandArgument_String_OrDefault(JUMP_READS, "jump_reads_filt_cpd");
     CommandArgument_String_OrDefault(FRAG_READS, "frag_reads_filt_cpd");
     CommandArgument_String_OrDefault(UNIPATH_PATCH_DIR, "unipath_patch");
     CommandArgument_String_OrDefault(POST_PATCH_DIR, "post_patch");
     CommandArgument_Bool_OrDefault(CHECKPOINT, False);
     CommandArgument_Bool_OrDefault(MERGE_OVERLAPPING_CONTIGS, True);
     CommandArgument_String_OrDefault_Doc(LOG, "",
          "{option_1,...,option_n} where option_i is one of\n"
          "ALL: all options\n"
          "SOME: some options (ATTEMPTED_JOINS + ALIGNS + ALIGNS_ALL)\n"
          "ATTEMPTED_JOINS: print attempted joins\n"
          "ASSEMBLY: log minimal assembly details\n"
          "ASSEMBLYi: log assembly details, level i, i in {1,2,3}\n"
          "ALIGN_READS: alignments of reads used\n"
          "ALIGNS: alignments for successful joins\n"
          "ALIGNS_ALL: alignments for attempted joins\n"
          "CORRECT: details of error correction during gap closing\n"
          "ACTION: log calls and returns during joining\n"
          "PATCHES: print sequences of patches\n");
     CommandArgument_Bool_OrDefault(INSTALL, True);
     CommandArgument_Bool_OrDefault(WRITE, True);
     CommandArgument_String_OrDefault_Doc(JOINS, "",
          "If specified, do only these joins.");
     CommandArgument_String_OrDefault_Doc(LR, "",
          "a list of pairs \"(u1 u2)\" of contigs to attempt to join "
          "(rather than pairs that are identified for joining)");
     CommandArgument_Bool_OrDefault(MOC_STRINGENT, False);
     CommandArgument_LongLong_OrDefault(MAX_MEMORY_GB,0);
     CommandArgument_Double_OrDefault_Doc(BRIDGEGAP_MAX_SIGMA, 5.0, "If positive, "
          "max standard deviations a bridge may vary from estimated gap size");
     EndCommandArguments;

     SetMaxMemory(MAX_MEMORY_GB<<30);

     // Start.

     double clock = WallClockTime( );

     // Thread control

     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads( NUM_THREADS );

     // Sanity check some arguments.

     ForceAssert( !SCAFFOLDS_IN.Contains( "/" ) );
     ForceAssert( !SCAFFOLDS_OUT.Contains( "/" ) );
     if ( !INSTALL && WRITE ) 
     {    cout << "!INSTALL && WRITE doesn't make sense" << endl;
          cout << "Abort." << endl;
          exit(1);    }

     // Heuristic constants.

     const int MAX_PATHS = 100;
     const int MAX_PATH_ITERATIONS = 100000;
     const int MAX_EDGES = 2000;
     const int MAX_READS = 1000;
     const int MAX_JOINS = 40;
     const int MIN_OVERLAP_END = 25;
     const int prox = 30;

     // Parse arguments.

     vec< pair<int,int> > LR_to_process;
     {    vec<String> LRs;
          ParseStringSet( LR, LRs );
          for ( int i = 0; i < LRs.isize( ); i++ )
          {    int L = LRs[i].Between( "(", " " ).Int( );
               int R = LRs[i].Between( " ", ")" ).Int( );
               LR_to_process.push( L, R );    }
          Sort(LR_to_process);    }
     vec<int> joins;
     ParseIntSet( JOINS, joins );

     // Define logging options.

     vec<String> log_options, allowed_log_options;
     ParseStringSet( LOG, log_options );
     allowed_log_options.push_back( "ALL", "SOME", "ACTION", "ATTEMPTED_JOINS" );
     allowed_log_options.push_back( "ASSEMBLY", "CORRECT", "PATCHES" );
     allowed_log_options.push_back( "ASSEMBLY1", "ASSEMBLY2", "ASSEMBLY3" );
     allowed_log_options.push_back( "ALIGN_READS", "ALIGNS", "ALIGNS_ALL" );
     Sort(allowed_log_options);
     if ( !BinSubset( log_options, allowed_log_options ) )
     {    cout << "Illegal log option." << endl << "Abort." << endl;
          exit(1);    }
     #define REQUESTED(x) Member( log_options, String(x) )
     Bool log_all = REQUESTED( "ALL" );
     Bool log_some = log_all || REQUESTED( "SOME" );
     Bool log_attempted_joins = log_some || REQUESTED( "ATTEMPTED_JOINS" );
     Bool log_align_reads = log_all || REQUESTED( "ALIGN_READS" );
     Bool log_aligns = log_some || REQUESTED( "ALIGNS" );
     Bool log_aligns_all = log_some || REQUESTED( "ALIGNS_ALL" );
     Bool log_action = log_all || REQUESTED( "ACTION" );
     Bool log_patches = log_all || REQUESTED( "PATCHES" );

     // Define directories.

     String data_dir = PRE + "/" + DATA, run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
     String tmp_dir = sub_dir + "/tmp";
     Mkdir777( sub_dir + "/" + POST_PATCH_DIR );
     Mkpath( tmp_dir );
     cout << Date( ) << ": " << run_dir << endl;
     String ch_head 
          = sub_dir + "/" + POST_PATCH_DIR + "/PostPatcher." + SCAFFOLDS_IN + ".";
     Mkdir777( run_dir + "/post_patch" );

     // Check for existence of files.

     String unibases_file = run_dir + "/" + HEAD + ".unibases.k" + ToString(K);
     String supers_file = sub_dir + "/" + SCAFFOLDS_IN + ".superb";
     String tigsa_file = sub_dir + "/" + SCAFFOLDS_IN + ".contigs.fasta";
     String reads_file = run_dir + "/" + FRAG_READS + ".fastb";
     String quals_file = run_dir + "/" + FRAG_READS + ".qualb";
     String jreads_file = run_dir + "/" + JUMP_READS + ".fastb";
     String jquals_file = run_dir + "/" + JUMP_READS + ".qualb";
     String pairs_file = run_dir + "/" + FRAG_READS + ".pairs";
     String jpairs_file = run_dir + "/" + JUMP_READS + ".pairs";
     vec<String> required;
     required.push_back( unibases_file, supers_file );
     required.push_back( tigsa_file, reads_file, quals_file, jreads_file );
     required.push_back( jquals_file, pairs_file, jpairs_file );
     for ( int i = 0; i < required.isize( ); i++ )
     {    if ( !IsRegularFile( required[i] ) )
          {    cout << "Can't find " << required[i] << endl;
               cout << "Abort." << endl;
               exit(1);    }    }
     String results_file = sub_dir + "/" + POST_PATCH_DIR + "/PostPatcher." 
          + SCAFFOLDS_OUT + ".RESULTS";
     if ( IsRegularFile(results_file) && IsOlder( results_file, supers_file ) )
     {    cout << "Results file older than input scaffolds file, not OK." << endl;
          cout << "Abort." << endl;
          exit(1);    }

     // Test for ready to join.

     String JOINDATA_file = ch_head + "JOINDATA", JDLEN_file = ch_head + "JDLEN";
     Bool ready_to_join = CHECKPOINT && IsRegularFile(JOINDATA_file)
          && IsRegularFile(JDLEN_file);

     // Compute read lengths.

     vec<uint16_t> read_len, jread_len;
     if ( !ready_to_join )
     {    cout << Date( ) << ": computing fragment read lengths" << endl;
          String READLEN_file 
               = sub_dir + "/" + POST_PATCH_DIR + "/PostPatcher.stable.READLEN";
          if ( CHECKPOINT && IsRegularFile(READLEN_file) )
               BinaryRead3( READLEN_file, read_len );
          else
          {    vecbasevector reads(reads_file);
               read_len.resize( reads.size( ) );
               for ( size_t i = 0; i < reads.size( ); i++ )
               {    ForceAssertLt( reads[i].size( ), 65536u );
                    read_len[i] = reads[i].size( );    }
               if (CHECKPOINT) BinaryWrite3( READLEN_file, read_len );    }
          cout << Date( ) << ": computing jump read lengths" << endl;
          String JREADLEN_file 
               = run_dir + "/" + POST_PATCH_DIR + "/PostPatcher.stable.JREADLEN";
          if ( CHECKPOINT && IsRegularFile(JREADLEN_file) )
               BinaryRead3( JREADLEN_file, jread_len );
          else
          {    vecbasevector jreads(jreads_file);
               jread_len.resize( jreads.size( ) );
               for ( size_t i = 0; i < jreads.size( ); i++ )
               {    ForceAssertLt( jreads[i].size( ), 65536u );
                    jread_len[i] = jreads[i].size( );    }
               if (CHECKPOINT) BinaryWrite3( JREADLEN_file, jread_len );    }    }

     // Load scaffolds.

     cout << Date( ) << ": loading scaffolds" << endl;
     vec<superb> scaffolds;
     ReadSuperbs( supers_file, scaffolds );
     int nscaff = scaffolds.size( );

     // Load the contigs.
     
     cout << Date( ) << ": loading contigs" << endl;
     String TIGS_file = ch_head + "TIGS";
     vecbasevector tigs;
     if ( CHECKPOINT && IsRegularFile(TIGS_file) ) tigs.ReadAll(TIGS_file);
     else
     {    FetchReads( tigs, 0, tigsa_file );
          tigs.WriteAll(TIGS_file);    }
     int ntigs = tigs.size( );
     vec<fastavector> tigsa;
     if ( !( ready_to_join && !INSTALL ) ) LoadFromFastaFile( tigsa_file, tigsa );

     // Map the unibases to the contigs.  There is a potential problem with the way 
     // we're doing this: we require that each unibase maps perfectly.

     cout << Date( ) << ": mapping unibases to contigs" << endl;
     String UALIGNS_file = ch_head + "UALIGNS";
     cout << ch_head << endl;
     cout << UALIGNS_file << endl;
     vec< triple<int64_t,int64_t,int> > UALIGNS;
     if ( !ready_to_join )
     {    if ( CHECKPOINT && IsRegularFile(UALIGNS_file) )
               BinaryRead3( UALIGNS_file, UALIGNS );
          else
          {    const int max_placements = 1;
               SearchFastb2( unibases_file, TIGS_file, K, &UALIGNS, 0, max_placements );
               if (CHECKPOINT) BinaryWrite3( UALIGNS_file, UALIGNS );    }    }

     // Load segments generated by UnipathPatcher.

     cout << Date( ) << ": loading segments" << endl;
     vec<segalign> SEGS, JSEGS;
     String SEGS_file = run_dir + "/" + UNIPATH_PATCH_DIR + "/UnipathPatcher.SEGS";
     String JSEGS_file = run_dir + "/" + UNIPATH_PATCH_DIR + "/UnipathPatcher.JSEGS";
     String RALIGNS_file = ch_head + "RALIGNS", JRALIGNS_file = ch_head + "JRALIGNS";
     if ( !ready_to_join )
     {    if ( !CHECKPOINT || !IsRegularFile(RALIGNS_file) )
               BinaryRead3( SEGS_file, SEGS );
          if ( !CHECKPOINT || !IsRegularFile(JRALIGNS_file) )
               BinaryRead3( JSEGS_file, JSEGS );    }

     // Load the unibases.

     vecbasevector unibases;
     if ( !ready_to_join ) unibases.ReadAll(unibases_file);

     // Index the unibase alignments.

     cout << Date( ) << ": indexing unibase alignments" << endl;
     vec<size_t> U_START;
     if ( !ready_to_join )
     {    U_START.resize( unibases.size( ) + 1 );
          size_t POS = 0;
          for ( int64_t u = 0; u <= (int64_t) unibases.size( ); u++ )
          {    while( POS < UALIGNS.size( ) && UALIGNS[POS].first < u ) ++POS;
               U_START[u] = POS;    }    }

     // Index scaffolds.

     cout << Date( ) << ": indexing scaffolds" << endl;
     vec<int> to_scaffold, to_scaffold_pos;
     if ( !ready_to_join )
     {    to_scaffold.resize(ntigs, -1), to_scaffold_pos.resize(ntigs, -1);
          for ( int i = 0; i < scaffolds.isize( ); i++ ) 
          {    const superb& s = scaffolds[i];
               for ( int j = 0; j < s.Ntigs( ); j++ ) 
               {    ForceAssertEq( to_scaffold[ s.Tig(j) ], -1 );
                    to_scaffold[ s.Tig(j) ] = i;
                    to_scaffold_pos[ s.Tig(j) ] = j;    }    }    }

     // Map the reads to the contigs.  Note that we assume reads are fw on the
     // unibases.

     cout << Date( ) << ": start process of mapping reads to contigs" << endl;
     vec< triple<int64_t,int,int> > RALIGNS, JRALIGNS; // (rid, tig, pos)
     if ( !ready_to_join )
     {    if ( CHECKPOINT && IsRegularFile(RALIGNS_file) )
               BinaryRead3( RALIGNS_file, RALIGNS );
          else
          {    cout << Date( ) << ": mapping fragment reads" << endl;
               for ( size_t i = 0; i < SEGS.size( ); i++ )
               {    const segalign& a = SEGS[i];
                    int64_t rid = a.rid; 
                    int u = a.u, rpos = a.rpos, upos = a.upos;
                    if ( U_START[u] < U_START[u+1] )
                    {    int tig = UALIGNS[ U_START[u] ].second;
                         int tpos = UALIGNS[ U_START[u] ].third;
                         int read_start_on_tig;
                         if ( tpos >= 0 ) read_start_on_tig = tpos + upos - rpos;
                         else
                         {    read_start_on_tig 
                                   = -tpos-1 + unibases[u].isize( ) - upos
                                   - read_len[rid] + rpos;    }
                         read_start_on_tig 
                              = Max( 0, read_start_on_tig ); // don't like!
                         if ( tpos >= 0 ) RALIGNS.push(rid, tig, read_start_on_tig);
                         else RALIGNS.push(rid, tig, -read_start_on_tig-1);    }    }
               ParallelUniqueSort(RALIGNS);
               cout << Date( ) << ": found " << RALIGNS.size( )
                    << " alignments of fragment reads" << endl;
               if (CHECKPOINT) BinaryWrite3( RALIGNS_file, RALIGNS );    }
          if ( CHECKPOINT && IsRegularFile(JRALIGNS_file) )
               BinaryRead3( JRALIGNS_file, JRALIGNS );
          else
          {    cout << Date( ) << ": mapping jump reads" << endl;
               for ( size_t i = 0; i < JSEGS.size( ); i++ )
               {    const segalign& a = JSEGS[i];
                    int64_t rid = a.rid; 
                    int u = a.u, rpos = a.rpos, upos = a.upos;
                    if ( U_START[u] < U_START[u+1] )
                    {    int tig = UALIGNS[ U_START[u] ].second;
                         int tpos = UALIGNS[ U_START[u] ].third;
                         int read_start_on_tig;
                         if ( tpos >= 0 ) read_start_on_tig = tpos + upos - rpos;
                         else
                         {    read_start_on_tig 
                                   = -tpos-1 + unibases[u].isize( ) - upos
                                   - jread_len[rid] + rpos;    }
                         read_start_on_tig 
                              = Max( 0, read_start_on_tig ); // don't like!
                         if ( tpos >= 0 ) JRALIGNS.push(rid, tig, read_start_on_tig);
                         else JRALIGNS.push(rid, tig, -read_start_on_tig-1);   }   }
               ParallelUniqueSort(JRALIGNS);
               if (CHECKPOINT) BinaryWrite3( JRALIGNS_file, JRALIGNS );    }
          Destroy(SEGS), Destroy(JSEGS), Destroy(UALIGNS), Destroy(U_START); 
          Destroy(unibases);    }

     // Set up data structures to track reads that might be in or near gaps.
     // We should but do not at present allow for the possibility that only one
     // read is in or near the gap; both are always included.

     vec< vec< vec<opair> > > ingap, jingap;

     // Find stuff in or near gaps.

     if ( !ready_to_join )
     {
     String INGAP_file = ch_head + "INGAP", JINGAP_file = ch_head + "JINGAP";
     if ( CHECKPOINT && IsRegularFile(INGAP_file) && IsRegularFile(JINGAP_file) )
     {    BinaryReader::readFile(INGAP_file.c_str(),&ingap);
          BinaryReader::readFile(JINGAP_file.c_str(),&jingap);    }
     else
     {    ingap.resize(nscaff), jingap.resize(nscaff);
          for ( int i = 0; i < nscaff; i++ )
          {    ingap[i].resize( scaffolds[i].Ngaps( ) );
               jingap[i].resize( scaffolds[i].Ngaps( ) );    }

          // Find stuff in or near gaps.

          cout << Date( ) << ": loading pairs" << endl;
          PairsManager pairs(pairs_file), jpairs(jpairs_file);
          pairs.makeCache( ), jpairs.makeCache( );
          int64_t nreads = MastervecFileObjectCount(reads_file);
          int64_t njreads = MastervecFileObjectCount(jreads_file);
          {    cout << Date( ) << ": indexing fragment read to tig alignments" 
                    << endl;
               vec<size_t> R_START(nreads+1);
               {    size_t POS = 0;
                    for ( int64_t u = 0; u <= nreads; u++ )
                    {    while( POS < RALIGNS.size( ) && RALIGNS[POS].first < u ) 
                              ++POS;
                         R_START[u] = POS;    }    }
               cout << Date( ) << ": looking for fragment reads near gaps" << endl;
               for ( size_t pi = 0; pi < pairs.nPairs( ); pi++ )
               {    
                    // First look for reads placed near gaps.
     
                    int64_t id1 = pairs.ID1(pi), id2 = pairs.ID2(pi);
                    int lib = pairs.libraryID(pi);
                    int sep = pairs.getLibrarySep(lib);
                    int dev = pairs.getLibrarySD(lib);
                    for ( int idpass = 1; idpass <= 2; idpass++ )
                    {    int64_t id = ( idpass == 1 ? id1 : id2 );
                         int64_t idp = ( idpass == 1 ? id2 : id1 );
     
                         // Suppose that read id is placed essentially uniquely.
                         // (No, now we allow multiple placements.  We at least need
                         // to allow placements on both sides of a gap.)
     
                         // if ( R_START[id+1] - R_START[id] == 1 )
                         for ( size_t l = R_START[id]; l < R_START[id+1]; l++ )
                         {    int tig = RALIGNS[l].second, pos = RALIGNS[l].third;
                              int xpos = ( pos >= 0 ? pos : -pos-1 );
                              int s = to_scaffold[tig], sp = to_scaffold_pos[tig];
                              if ( s < 0 ) continue;
                              const superb& S = scaffolds[s];
                              int tiglen = tigsa[tig].size( );
     
                              // Is read id placed in or near a gap?

                              if ( xpos <= prox && sp > 0 ) 
                              {    if ( pos >= 0 )
                                        ingap[s][sp-1].push( id, idp, sep, dev );
                                   else 
                                   {    ingap[s][sp-1].push( 
                                             idp, id, sep, dev );    }    }
                              if ( xpos + read_len[id] >= tiglen - prox 
                                   && sp < S.Ngaps( ) )
                              {    if ( pos >= 0 )
                                        ingap[s][sp].push( id, idp, sep, dev );
                                   else ingap[s][sp].push( idp, id, sep, dev );    }

                              // Is read idp unplaced?  Where might it belong?
                              // For now we don't bother with this case, since our
                              // fragment pairs are short and at present we always
                              // uses pairs in the gap assemblies.
     
                                    }    }    }    }
          {    cout << Date( ) << ": indexing jump read to tig alignments" << endl;
               vec<size_t> JR_START(njreads+1);
               {    size_t POS = 0;
                    for ( int64_t u = 0; u <= njreads; u++ )
                    {    while( POS < JRALIGNS.size( ) && JRALIGNS[POS].first < u ) 
                              ++POS;
                         JR_START[u] = POS;    }    }
               cout << Date( ) << ": looking for jump reads near gaps" << endl;
               for ( size_t pi = 0; pi < jpairs.nPairs( ); pi++ )
               {    
                    // First look for reads placed near gaps.
     
                    int64_t id1 = jpairs.ID1(pi), id2 = jpairs.ID2(pi);
                    for ( int idpass = 1; idpass <= 2; idpass++ )
                    {    int64_t id = ( idpass == 1 ? id1 : id2 );
                         int64_t idp = ( idpass == 1 ? id2 : id1 );

                         // Suppose that read id is uniquely placed.
                         // (No, now we allow multiple placements.  We at least need
                         // to allow placements on both sides of a gap.)
     
                         // if ( JR_START[id+1] - JR_START[id] == 1 )
                         for ( size_t l = JR_START[id]; l < JR_START[id+1]; l++ )
                         {    int tig = JRALIGNS[l].second, pos = JRALIGNS[l].third;
                              int xpos = ( pos >= 0 ? pos : -pos-1 );
                              int s = to_scaffold[tig], sp = to_scaffold_pos[tig];
                              if ( s < 0 ) continue;
                              const superb& S = scaffolds[s];
                              int tiglen = tigsa[tig].size( );
                              int sep = jpairs.getLibrarySep( jpairs.libraryID(pi) );
                              int mean = sep + jread_len[id] + jread_len[idp];
                              int dev = jpairs.getLibrarySD( jpairs.libraryID(pi) );
     
                              // Is read id placed in or near a gap?
     
                              if ( xpos <= prox && sp > 0 ) 
                              {    if ( pos >= 0 )
                                        jingap[s][sp-1].push( id, idp, sep, dev );
                                   else 
                                   {    jingap[s][sp-1].push( 
                                             idp, id, sep, dev );    }    }
                              if ( xpos + jread_len[id] >= tiglen - prox 
                                   && sp < S.Ngaps( ) )
                              {    if ( pos >= 0 )
                                        jingap[s][sp].push( id, idp, sep, dev );
                                   else jingap[s][sp].push( idp, id, sep, dev );    }
     
                              // Is read idp unplaced?  Where might it belong?
                              // We assume that the pair is an outie.  But then
                              // aren't the pairs reversed?  What we're doing here
                              // may be nonsensical.
     
                              if ( JR_START[idp+1] - JR_START[idp] == 0 )
                              {    const int fudge = 300;
                                   const int dev_mult = 3;
                                   if ( pos < 0 )
                                   {    int startp = xpos + mean
                                             - jread_len[idp] - S.Len(sp);
                                        for ( int r = sp; r < S.Ngaps( ); r++ )
                                        {    int low = -dev_mult * S.Dev(r) - fudge;
                                             int high = S.Gap(r) 
                                                  + dev_mult * S.Dev(r) + fudge;
                                             if ( startp < low ) break;
                                             if ( startp <= high )
                                             {    if ( pos >= 0 )
                                                  {    jingap[s][r].push( id, idp, 
                                                            sep, dev );    }
                                                  else 
                                                  {    jingap[s][r].push( idp, id, 
                                                            sep, dev );    }    }
                                             startp -= 
                                                  S.Gap(r) - S.Len(r+1);    }    }
                                   else
                                   {    int startp = pos + jread_len[id] - mean;
                                        for ( int r = sp; r >= 1; r-- )
                                        {    int low = -dev_mult * S.Dev(r-1)
                                                  - S.Gap(r-1) - fudge;
                                             int high 
                                                  = dev_mult * S.Dev(r-1) + fudge;
                                             if ( startp > high ) break;
                                             if ( startp >= low )
                                             {    if ( pos >= 0 )
                                                  {    jingap[s][r-1].push( id, idp, 
                                                            sep, dev );    }
                                                  else 
                                                  {    jingap[s][r-1].push( idp, id, 
                                                            sep, dev );    }    }
                                             startp += 
                                                  S.Gap(r-1) + S.Len(r-1);    }    }
                                                  }    }    }    }    }
          BinaryWriter::writeFile(INGAP_file.c_str(),ingap);
          BinaryWriter::writeFile(JINGAP_file.c_str(),jingap); }
     Destroy(RALIGNS), Destroy(JRALIGNS), Destroy(read_len), Destroy(jread_len);
     Destroy(to_scaffold), Destroy(to_scaffold_pos);
     }

     // Temporarily dump data for scaffolds.

     if ( !ready_to_join && log_all )
     {    int s = 0;
          const superb& S = scaffolds[s];
          cout << "\nSCAFFOLD " << s << endl;
          for ( int i = 0; i < S.Ntigs( ); i++ )
          {    cout << "\nCONTIG " << S.Tig(i) << "\n";
               if ( i == S.Ngaps( ) ) break;
               cout << "\nfragment pairs in gap:\n";
               for ( int j = 0; j < ingap[s][i].isize( ); j++ )
               {    const opair& p = ingap[s][i][j];
                    cout << p.id1 << " " << p.id2 << " " << p.sep << " +/- " 
                         << p.dev << "\n";    }
               cout << "\njump pairs in gap:\n";
               for ( int j = 0; j < jingap[s][i].isize( ); j++ )
               {    const opair& p = jingap[s][i][j];
                    cout << p.id1 << " " << p.id2 << " " << p.sep << " +/- " 
                         << p.dev << "\n";    }    }
          cout << "\n";    }

     // Build join data.

     vec<size_t> join_data_offsets;
     if ( !ready_to_join )
     {    SystemSucceed( "PostPatcherBuildJoins" + ARGC(PRE) + ARGC(DATA) + ARGC(RUN)
               + ARGC(SUBDIR) + ARGC(K) + ARGC(SCAFFOLDS_IN) + ARGC(JUMP_READS)
               + ARGC(FRAG_READS) + ARGC(POST_PATCH_DIR) + ARG(NH, True) );    }
     Bool have_results = CHECKPOINT && IsRegularFile( run_dir + "/" + POST_PATCH_DIR 
          + "/PostPatcher." + SCAFFOLDS_OUT + ".RESULTS" );
     if ( !have_results )
     {    BinaryReader::readFile( JDLEN_file.c_str(), &join_data_offsets );
          cout << Date( ) << ": join data created, memory usage = "
               << ToStringAddCommas( MemUsageBytes( ) ) << endl;
          // Note: something is wrong at this point.  On bushbaby, we observe 135 GB
          // memory usage, whereas the total amount accounted for by undestroyed
          // data structures is 3 GB.  Perhaps Destroy is not working.
          // (Should be different now.)
          }

     // Try to make joins. 

     vec<int> gap_to_scaffold, gap_to_scaffold_pos;
     int ngaps = 0;
     for ( int i = 0; i < scaffolds.isize( ); i++ )
     {    for ( int j = 0; j < scaffolds[i].Ngaps( ); j++ )
          {    gap_to_scaffold.push_back(i);
               gap_to_scaffold_pos.push_back(j);    
               ngaps++;    }    }
     vec< pair<int,int> > Xjoiners;
     for ( int i = 0; i < ngaps; i++ )
     {    const superb& S = scaffolds[ gap_to_scaffold[i] ];
          int m1 = S.Tig( gap_to_scaffold_pos[i] );
          int m2 = S.Tig( gap_to_scaffold_pos[i] + 1 );
          Xjoiners.push( m1, m2 );    }
     bvec3 new_stuff;
     vec< vec<int> > start_stop;
     vec<String> reports;
     double join_clock = WallClockTime( );
     int attempted_joins, npatches;
     vec<Bool> joined, perfect;
     if ( !have_results )
     {    MakeJoins( Xjoiners, tigs, "PostPatcher", K, data_dir, run_dir, POST_PATCH_DIR, 
               attempted_joins, joined, npatches, perfect, new_stuff, start_stop, 
               reports, NUM_THREADS, joins, LOG, log_attempted_joins, log_action, 
               MAX_PATHS, MAX_PATH_ITERATIONS, MAX_EDGES, MAX_READS, MAX_JOINS,
               MIN_OVERLAP_END, JOINDATA_file, join_data_offsets,
               LR_to_process );    }

     // Save join results.

     if ( !have_results )
     {    BinaryWriter w(results_file.c_str());
          w.write(new_stuff);
          w.write(start_stop);
          w.close(); }
     else
     {    BinaryReader r(results_file.c_str());
          r.read( &new_stuff );
          r.read( &start_stop );    }

     // Print reports.

     if ( !have_results )
     {    Bool have_report = False;
          for ( int i = 0; i < reports.isize( ); i++ )
               if ( reports[i].size( ) > 0 ) have_report = True;
          if (have_report) cout << "\n";
          for ( int i = 0; i < reports.isize( ); i++ )
               cout << reports[i];
          if (have_report) cout << "\n";
          cout << Date( ) << ": " << TimeSince(join_clock) << " used in joining" 
               << endl;    }

     // Install the patches.  For a given gap, if there is more than one patch, we
     // require that the patches have only mismatches between each other, and the
     // total mismatch fraction is at most 5%.  This code is quadratic in the 
     // scaffold size, since each edit causes the entire scaffold to be edited.  We 
     // go backwards through the gaps so that we go backwards through each super, 
     // and thus the edits don't step on each other.

     vec<assembly_edit> edits;
     if (INSTALL)
     {    for ( int i = ngaps-1; i >= 0; i-- )
          {    int nbridges = new_stuff[i].size( ) / 6;
               if ( nbridges == 0 ) continue;

               int s = gap_to_scaffold[i], p = gap_to_scaffold_pos[i];
               superb& S = scaffolds[s];
	       int m1 = S.Tig(p), m2 = S.Tig(p+1);

	       vec<int> gaps, starts, stops;
               vec<basevector> bridges;
               for ( int j = 0; j < nbridges; j++ ) 
               {    bridges.push_back( new_stuff[i][6*j] );
		    int L_start = start_stop[i][2*j]; 
		    int R_stop = start_stop[i][2*j+1];
		    int gap = L_start + bridges[j].isize( ) 
		         - R_stop - (int) tigsa[m1].size( );
		    starts.push_back( L_start );
		    stops.push_back( R_stop );
		    gaps.push_back(gap);    }
               if (log_patches)
               {    cout << "\nscaffold " << s << ", from contig "
                         << m1 << " to " << m2 << ":\n";
                    cout << "predicted gap = " << S.Gap(p) << " +/- " << S.Dev(p)
                         << "\n";
                    cout << "list of bridged gaps [devs off by] =";
                    for ( int j = 0; j < nbridges; j++ )
                    {    int actual = gaps[j];
                         double offby 
                              = double( actual - S.Gap(p) ) / double( S.Dev(p) );
                         cout << " " << actual << "[" << std::fixed 
                              << setprecision(2) << offby << "]";    }
                    cout << "\npatches:\n";
                    for ( int j = 0; j < nbridges; j++ )
                    {    bridges[j].Print( cout, "bridge_" + ToString(j+1)
                              + ":" + ToString( start_stop[i][2*j] ) + "-"
                              + ToString( start_stop[i][2*j+1] )
			      + "s" + ToString(bridges[j].size())
					   + "g" + ToString(gaps[j]));    }    }
               UniqueSort(starts), UniqueSort(stops);
               if ( !starts.solo( ) || !stops.solo( ) ) 
               {    if (log_patches)
                    {    cout << "Rejecting patches because more than one start "
                              << "or stop." << endl;    }
                    continue;    }

               // Minimize edit and save.

               int start1 = starts[0], stop2 = stops[0];
               while( start1 < (int) tigsa[m1].size( ) )
               {    Bool agree = True;
                    for ( int j = 0; j < bridges.isize( ); j++ )
                    {    if ( bridges[j].size( ) == 0
                              || tigsa[m1][start1] != as_base( bridges[j][0] ) )
                         {    agree = False;    }    }
                    if ( !agree ) break;
                    start1++;
                    for ( int j = 0; j < bridges.isize( ); j++ )
                    {    bridges[j].SetToSubOf( 
                              bridges[j], 1, bridges[j].isize( ) - 1 );    }    }
               while( stop2 > 0 )
               {    Bool agree = True;
                    for ( int j = 0; j < bridges.isize( ); j++ )
                    {    if ( bridges[j].size( ) == 0
                              || tigsa[m2][stop2-1] 
                              != as_base( bridges[j][ bridges[j].isize( ) - 1 ] ) )
                         {    agree = False;    }    }
                    if ( !agree ) break;
                    stop2--;
                    for ( int j = 0; j < bridges.isize( ); j++ )
                         bridges[j].resize( bridges[j].isize( ) - 1 );    }
               edits.push( assembly_edit::GAP_CLOSER, m1, start1, m2, stop2, 
                    bridges );    }    }

     // Merge overlapping contigs.

     if (MERGE_OVERLAPPING_CONTIGS) {
       cout << Date( ) << ": merging overlapping contigs" << endl;
       int merge_count = 0;
       for (int s = 0; s < nscaff; ++s) {
	 superb &S = scaffolds[s];
	 for (int p = S.Ngaps() - 1; p >= 0; --p) {
	   int gap = S.Gap(p);
	   if (gap >= 0) continue;

	   int dev = S.Dev(p);
	   int t1 = S.Tig(p);
	   int t2 = S.Tig(p+1);
	   
	   int min_overlap = -gap - BRIDGEGAP_MAX_SIGMA * dev;
           if ( MOC_STRINGENT && min_overlap < 0 ) continue;
	   int max_overlap = -gap + BRIDGEGAP_MAX_SIGMA * dev;
  
	   min_overlap = max(min_overlap, 1);

	   bitvector c1_amb, c2_amb;
	   bvec c1 = tigsa[t1].ToBasevector(&c1_amb);
	   bvec c2 = tigsa[t2].ToBasevector(&c2_amb);

	   align a;
	   int rc;
	   int overlap = AlignTwoBasevectors(c1, c2, a, min_overlap, max_overlap, 0.05, NULL, rc);

	   if (overlap > 0) {
	     /* Debugging
	     PRINT2(s,p);
	     PRINT4(t1, t2, min_overlap, max_overlap);
	     if (overlap >= min_overlap && overlap <= max_overlap) cout << "GOOD ";
	     else cout << "BAD ";
	     PRINT(overlap);
	     PRINT4(a.pos1(), a.Pos1(), a.pos2(), a.Pos2());
	     */

             if (log_patches)
             {    cout << "merging contigs of size " << c1.size( ) << " and "
                       << c2.size( ) << " along overlap of size " << overlap
                       << "\n";    }

             vec<basevector> bridges;
             basevector b;
             bridges.push_back(b);
             edits.push( assembly_edit::GAP_CLOSER, t1, a.pos1( ), t2, 0, bridges );

	     ++merge_count;
	   }
	 }
       }
       cout << Date( ) << ": " << merge_count << " contigs merged" << endl;
     }

     // Write edits.

     if (WRITE)
     {    String outputFile = sub_dir + "/" + SCAFFOLDS_OUT + ".edits";
          BinaryWriter::writeFile( outputFile.c_str( ), edits );    }

     // Finish up.

     if ( !have_results )
     {    cout << "\nSUMMARY:\n";
          cout << attempted_joins << " joins attempted" << endl;
          cout << Sum(joined) << " joins made, yielding " << npatches
               << " patches" << endl;
          if (log_aligns)
          {    cout << Sum(joined) - Sum(perfect) << " imperfect joins made" 
                    << endl;    }    }
     cout << "\n" << Date( ) << ": done, time used = "
          << TimeSince(clock) << endl;    }
