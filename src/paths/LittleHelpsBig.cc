/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// LittleHelpsBig.  Use unipaths from little K1 to enlarge unipaths from big K2.

#include "Basevector.h"
#include "MainTools.h"
#include "paths/GetNexts.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerPath.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/Unipath.h"
#include "paths/UnibaseUtils.h"


static inline 
String Tag(String S = "LHB") { return Date() + " (" + S + "): "; } 



int main(int argc, char *argv[])
{
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_String(IN_HEAD1);
  CommandArgument_String(IN_HEAD2);
  CommandArgument_Int(K1);
  CommandArgument_Int(K2);
  CommandArgument_String(OUT_HEAD);
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  EndCommandArguments;
  
  RunTime();  

  // Thread control
   
  NUM_THREADS = configNumThreads(NUM_THREADS);
 
  // Load unibases.
  cout << Tag() << "Loading unibases" << endl;
  
  String KS1 = ToString(K1), KS2 = ToString(K2);
  BaseVecVec unibases1(IN_HEAD1 + ".unibases.k" + KS1);
  BaseVecVec unibases2(IN_HEAD2 + ".unibases.k" + KS2);
  
  const size_t nuni1 = unibases1.size();
  const size_t nuni2 = unibases2.size();

  cout << Tag() << "n_unibases(K=" << KS1 << ")= " << nuni1 << endl;
  cout << Tag() << "n_unibases(K=" << KS2 << ")= " << nuni2 << endl;
  
  // Make output directory, if necessary.
  if (OUT_HEAD.Contains("/")) Mkpath(OUT_HEAD.RevBefore("/"));


  BaseVecVec new_reads = unibases2;
  {  // Scoping out intermediate data structures
  
    // Build unipaths2 and find terminal unibases.
    // Terminal unibases, or "leaves", are unibases that have only a single
    // adjacency on the unipath adjacency graph.  Leaves that are anchored at
    // their start are called "sinks", while leaves anchored at their end are
    // called "sources".
    vec<longlong> sinks;
    vec<longlong> sources;
    
    { // scoping bracket
      cout << Tag() << "Loading unipaths" << endl;

      vecKmerPath unipaths2(IN_HEAD2 + ".unipaths.k" + KS2);
      digraph A2;
      BinaryRead(IN_HEAD2 + ".unipath_adjgraph.k" + KS2, A2);
      HyperKmerPath h2;
      BuildUnipathAdjacencyHyperKmerPath(K2, A2, unipaths2, h2);
    
      // Find terminal unibases.
      cout << Tag() << "Finding sorces and sinks" << endl;
      for (longlong v = 0; v < h2.N(); v++) {
	if (h2.To(v).empty() && h2.From(v).solo())
	  sources.push_back(h2.EdgeObjectIndexByIndexFrom(v, 0));
	if (h2.From(v).empty() && h2.To(v).solo())
	  sinks.push_back(h2.EdgeObjectIndexByIndexTo(v, 0));
      }
    }
    const size_t n_sinks = sinks.size();
    const size_t n_sources = sources.size();

    cout << Tag() << "n_sinks=   " << n_sinks << endl;
    cout << Tag() << "n_sources= " << n_sources << endl;

    
    cout << Tag() << "Creating local data structures" << endl;
    
    // Determine which unibases1 can follow which.  This is not quite right
    // because it assumes that K1-1 overlap is enough for adjacency.
  
    // Build paths, using both sets of unibases and K=K1.
    BaseVecVec unibases = unibases1;
    unibases.Append(unibases2);
    vecKmerPath paths;
    ReadsToPathsCoreY(unibases, K1, paths, 
		       OUT_HEAD + ".LittleHelpsBig.unibases", NUM_THREADS);
  
    // Form paths database for for unibases1.
    vecKmerPath paths1;
    vec<tagged_rpint> pathsdb1;
    for (size_t i = 0; i < nuni1; i++)
      paths1.push_back_reserve(paths[i]);
    CreateDatabase(paths1, pathsdb1);
  
    // Pre-processing step:
    // For each sink, and each source, determine whether or not an alignment is
    // possible using that unibase.  We do this by examining the unibases1.
    vec<bool>   sinks_good(n_sinks,   true);
    vec<bool> sources_good(n_sources, true);
  
    for (size_t j1 = 0; j1 < n_sinks; j1++) {
      // Find the last kmer in this unipath (which is from unipaths2).
      size_t sink = sinks[j1] + nuni1;
      kmer_id_t kmer = paths[sink].Stop();
      
      // This kmer must appear EXACTLY once in the pathsdb1.
      vec<longlong> rpints;
      Contains(pathsdb1, kmer, rpints);
      if (!rpints.solo()) sinks_good[j1] = false;
    }
  
    for (size_t j2 = 0; j2 < n_sources; j2++) {
      // Find the first kmer in this unipath.
      size_t source = sources[j2] + nuni1;
      kmer_id_t kmer = paths[source].Start();
      
      // This kmer must appear EXACTLY once in the pathsdb1.
      vec<longlong> rpints;
      Contains(pathsdb1, kmer, rpints);
      if (!rpints.solo()) sources_good[j2] = false;
    }
  
  
  
    // Create a BaseVecVec containing the last K2 bases of each of the unibases
    // in sinks, and the first K2 bases of each of the unibases in sources.
    // Then, path these bases using K1, and create a pathsdb.  Now we can look up
    // a K1-mer in the pathsdb to quickly determine whether or not an overlap of
    // size > K1 is possible.
    cout << Tag() << "Building leafs." << endl;
    BaseVecVec leaf_bases;
    const size_t n_total_leaves = n_sources + n_sinks;
    leaf_bases.Reserve(n_total_leaves + ((K2-1) * n_total_leaves) / 16, n_total_leaves);
    for (size_t j1 = 0; j1 < n_sinks; j1++) {
      const BaseVec & sink = unibases[sinks[j1] + nuni1];
      leaf_bases.push_back(BaseVec(sink, sink.size() - (K2-1), K2-1));
    }
    for (size_t j2 = 0; j2 < n_sources; j2++) {
      const BaseVec & source = unibases[sources[j2] + nuni1];
      leaf_bases.push_back(BaseVec(source, 0, K2-1));
    }

    cout << Tag() << "Pathing." << endl;
    vecKmerPath leaf_paths;
    vecKmerPath leaf_pathsrc;
    vec<tagged_rpint> leaf_pathsdb;
    ReadsToPathsCoreY(leaf_bases, K1, leaf_paths, leaf_pathsrc, leaf_pathsdb);
  


    // Now, gather a set of eligible pairs of sink and source unibases.
    // A pair of sink and source unibases (call them S1 and S2) is considered
    // "eligible" iff the last K1-mer of S1 appears in S2 *and* the first K1-mer
    // of S2 appears in S1.  If this criterion is not met, these unibases cannot
    // possibly overlap by K1 or more bases.
    cout << Tag() << "Building list of eligible pairs." << endl;
    vec< pair<size_t, size_t> > eligibles;
  

    // Loop over all sink unipaths, and examine the last kmer in each.
    cout << Tag() << "Loop over sinks." << endl;
    for (size_t j1 = 0; j1 < n_sinks; dots_pct(j1++, n_sinks)) {
      kmer_id_t kmer = leaf_paths[j1].Stop();
      
      // Find all other appearances of this kmer in leaf_pathsdb.
      vec<longlong> places;
      Contains(leaf_pathsdb, kmer, places);
      const size_t n_places = places.size();

      for (size_t ip = 0; ip < n_places; ip++) {
	const longlong id = leaf_pathsdb[places[ip]].PathId();
	if (id >= 0 &&    // ignore RC alignments
            id >= longlong(n_sinks))  // ignore alignments to other sink unipaths
          eligibles.push_back(make_pair(j1, size_t(id - n_sinks)));
      }
    }
    


    // Loop over all source unipaths, and examine the first kmer in each.
    cout << Tag() << "Loop over sources." << endl;
    for (size_t j2 = 0; j2 < n_sinks; dots_pct(j2++, n_sinks)) {
      kmer_id_t kmer = leaf_paths[n_sinks + j2].Start();
      
      // Find all other appearances of this kmer in leaf_pathsdb.
      vec<longlong> places;
      Contains(leaf_pathsdb, kmer, places);
      const size_t n_places = places.size();
      
      for (size_t ip = 0; ip < n_places; ip++) {
	const longlong id = leaf_pathsdb[places[ip]].PathId();
	if (id >= 0 &&  // ignore RC alignments
            id < longlong(n_sinks)) // ignore alignments to other source unipaths
          eligibles.push_back(make_pair(size_t(id), j2));
      }
    }
  
    // Every eligible sink/source pair should now appear TWICE in eligibles.
    // Sort the list of eligibles, and keep track of how many times each s/s
    // pair appeared - later, we'll drop all s/s pairs that appeared only once.
    cout << Tag() << "Unique sort." << endl;
    vec<unsigned> eligible_counts;
    UniqueSortAndCount(eligibles, eligible_counts);
  
  
  
  
    // MAIN ALGORITHM
    // For each eligible sink/source pair, try to bridge.

    const size_t n_eligibles = eligibles.size();

    cout << Tag() << "Looking for bridges" << endl;
    cout << Tag() << "n_eligibles = " << n_eligibles << endl;

    double timestamp = WallClockTime();
  
  
    BaseVecVec merged_paths;
    size_t n_merges = 0;
    
    // Loop over sink/source pairs.
    for (size_t ie = 0; ie < n_eligibles; ie++) {
      size_t j1 = eligibles[ie].first;
      size_t j2 = eligibles[ie].second;
      
      if (sinks_good[j1] &&
          sources_good[j2] &&
          eligible_counts[ie] != 1) {
        
        size_t sink = sinks[j1] + nuni1;
        size_t source = sources[j2] + nuni1;
        kmer_id_t x1 = paths[sink].Stop();
        kmer_id_t x2 = paths[source].Start();
        
        // Finally, check whether there is a direct overlap between
        // the paths.  This is an expensive function call, so we only
        // perform it if we have to.
        size_t overlap = LargestOverlap(unibases[sink],
                                        unibases[source],
                                        K2-1, K1);
        if (overlap > 0) {
          
          // We've found an overlap!  Merge the unibases together and add the
          // resulting BaseVec to the merged_paths.
          BaseVec merge = unibases[sink];
          merge.resize(merge.size() - overlap);
          merge = Cat(merge, unibases[source]);
          merged_paths.push_back_reserve(merge);
          n_merges++;
        }
      }
      dots_pct(ie, n_eligibles);
    }
  
    cout << "Number of merges made: " << n_merges << endl;
    cout << "Elapsed time: " << (WallClockTime() - timestamp) << " s" << endl;
  
    // Combine old unibases and patches together

    // GetNexts() needs vec<vec<in> >. This could be a problem.
    vec< vec<int> > nexts2(nuni2);
    GetNexts(K2, unibases2, nexts2);
    new_reads.Append(merged_paths);
    for (size_t id1 = 0; id1 < nuni2; id1++) {
      for (size_t j = 0; j < nexts2[id1].size(); j++) {
	size_t id2 = nexts2[id1][j];
	BaseVec b = unibases2[id1];
	b.resize(b.size() + 1);
	b.Set(b.size() - 1, unibases2[id2][K2-1]);
	new_reads.push_back_reserve(b);
      }
    }

  } // end scoping of intermediate data structures

  // Now build the new unipaths.
  cout << Tag() << "Building new unipaths" << endl;

  vecKmerPath new_paths;
  vecKmerPath new_pathsrc;
  vecKmerPath new_unipaths;
  vec<tagged_rpint> new_pathsdb;
  vec<tagged_rpint> new_unipathsdb;
  digraph new_adj_graph;
  ReadsToPathsCoreY(new_reads, K2, new_paths, new_pathsrc, new_pathsdb,
                     OUT_HEAD + ".LittleHelpsBig.pseudo_reads", NUM_THREADS);
  Unipath(new_paths, new_pathsrc, new_pathsdb, new_unipaths, new_unipathsdb);
  KmerBaseBroker new_kbb(K2, new_paths, new_pathsrc, new_pathsdb, new_reads);
  BuildUnipathAdjacencyGraph(new_paths, new_pathsrc, new_pathsdb, 
                             new_unipaths, new_unipathsdb, new_adj_graph);

  // Write output files.
  cout << Tag() << "Writing output files" << endl;
  new_unipaths.WriteAll(OUT_HEAD + ".unipaths.k" + KS2);
  BinaryWrite3(OUT_HEAD + ".unipathsdb.k" + KS2, new_unipathsdb);
  BinaryWrite(OUT_HEAD + ".unipath_adjgraph.k" + KS2, new_adj_graph);
  BaseVecVec new_unibases;
  new_unibases.reserve(new_unipaths.size());
  for (size_t i = 0; i < new_unipaths.size(); i++)
    new_unibases.push_back(new_kbb.Seq(new_unipaths[i]));
  new_unibases.WriteAll(OUT_HEAD + ".unibases.k" + KS2);
  
  cout << Tag() << "Done!" << endl;
  return 0;
}
