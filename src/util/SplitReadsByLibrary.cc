///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

const char *DOC =
  "Splits a fastb file containing paired reads into sets of one or two fastb files, "
  "a set for each library, containing the 1st reads (A) and the 2nd reads (B) "
  "(in corresponding order).";

#include "MainTools.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "PairsManager.h"
#include "feudal/IncrementalWriter.h"


int main( int argc, char *argv[] )
{
  RunTime( );
  BeginCommandArguments;
  CommandDoc(DOC);
  CommandArgument_String_Doc(READS_IN,
    "Fastb file with associated .pairs or pairto file to split.");
  CommandArgument_String_OrDefault_Doc(READS_OUT, "",
    "Output fastb files are READS_OUT.<libname>.<A|B|AB>.fastb");
  CommandArgument_Bool_OrDefault_Doc(QUALS, False,
    "Split associated qualb files.");
  CommandArgument_Bool_OrDefault_Doc(SPLIT_PAIRS, False,
    "Split into A and B reads or a single interleaved AB read fastb.");
  CommandArgument_Bool_OrDefault_Doc(WRITE_PAIRS, True,
    "When SPLIT_PAIRS=False, should we write a pairto file?");
  CommandArgument_StringSet_OrDefault_Doc(LIBRARIES,"",
    "Libraries to extract - default is all libraries");
  EndCommandArguments;

  if (READS_OUT == "")
    READS_OUT = READS_IN;

  // Strip .fastb from filenames

  READS_IN = READS_IN.SafeBefore(".fastb");
  READS_OUT = READS_OUT.SafeBefore(".fastb");

  size_t nreads = MastervecFileObjectCount(READS_IN + ".fastb");

  // Load pairing information

  cout << Date() << " Loading pairing information..." << endl;
  PairsManager pairs;
  if (IsRegularFile(READS_IN + ".pairs"))
    pairs.Read(READS_IN + ".pairs");
  else
    pairs.ReadFromPairtoFile(READS_IN + ".pairto", nreads);
  ForceAssertEq(pairs.nReads(), nreads);
  size_t npairs = pairs.nPairs();
  size_t nlibs = pairs.nLibraries();

  cout << Date() << " Found " << nlibs << " libraries" << endl;
  cout << Date() << " Found " << npairs * 2 << " paired reads" << endl;

  // Determine which libraries to extract

  vec<bool> libs_to_extract(nlibs, false);
  if (LIBRARIES.empty()) 
    libs_to_extract.assign(nlibs, true);
  else {
    int unknown = 0;
    for (size_t i = 0; i < LIBRARIES.size(); i++) {
      int pos = Position(pairs.getLibraryNames(), LIBRARIES[i]);
      if (pos != -1)
	libs_to_extract[pos] = true;
      else {
	cout << " Unknown library name: " << LIBRARIES[i] << endl;
	unknown++;
      }
    }
    if (unknown)  {
      cout << "ERROR: Unable to find " + ToString(unknown) + " libraries." << endl;
      exit(1);
    }
  }
    
  // New pairings (if not splitting by A and B read)

  vec<PairsManager> pairs_out(nlibs);
  vec<longlong> lib_read_count(nlibs,0);

  // Prepare incremental writers for each library

  vec<IncrementalWriter<bvec>* > readsA(nlibs), readsBB(nlibs);
  vec<IncrementalWriter<qvec>* > qualsA(nlibs), qualsBB(nlibs);
  vec<IncrementalWriter<bvec>* >& readsB = (SPLIT_PAIRS ? readsBB : readsA);
  vec<IncrementalWriter<qvec>* >& qualsB = (SPLIT_PAIRS ? qualsBB : qualsA);

  for (size_t i = 0; i < nlibs; ++i) {
    String libName = pairs.getLibraryName(i);
    if (libs_to_extract[i])
      cout << Date() << " Extracting library: " << libName << endl;
    else {
      cout << Date() << " Skipping library  : " << libName << endl;
      continue;
    }

    if (SPLIT_PAIRS) {
      String readsNameA = READS_OUT + "." + libName + ".A.fastb";
      readsA[i] = new IncrementalWriter<bvec>(readsNameA.c_str());
      String readsNameB = READS_OUT + "." + libName + ".B.fastb";
      readsB[i] = new IncrementalWriter<bvec>(readsNameB.c_str());

      if (QUALS) {
	String qualsNameA = READS_OUT + "." + libName + ".A.qualb";
	qualsA[i] = new IncrementalWriter<qvec>(qualsNameA.c_str());
	String qualsNameB = READS_OUT + "." + libName + ".B.qualb";
	qualsB[i] = new IncrementalWriter<qvec>(qualsNameB.c_str());
      }
    } else {
      String readsNameAB = READS_OUT + "." + libName + ".AB.fastb";
      readsA[i] = new IncrementalWriter<bvec>(readsNameAB.c_str());

      if (QUALS) {
	String qualsNameAB = READS_OUT + "." + libName + ".AB.qualb";
	qualsA[i] = new IncrementalWriter<qvec>(qualsNameAB.c_str());
      }
    }
   
  }

  // Split reads by library

  {
    cout << Date() << " Loading reads..." << endl;
    vecbasevector reads(READS_IN + ".fastb");

    cout << Date() << " Splitting reads by library" << endl;
    for (size_t pairId = 0; pairId < npairs; pairId++ ) {
      int library = pairs.libraryID(pairId);
      if (!libs_to_extract[library]) continue;
      readsA[library]->add(reads[pairs.ID1(pairId)]);
      readsB[library]->add(reads[pairs.ID2(pairId)]);
    }
  }


  // Split quals by library

  if (QUALS) {
    vecqualvector quals;
    cout << Date() << " Loading quals..." << endl;
    quals.ReadAll(READS_IN + ".qualb");
    ForceAssertEq(nreads, quals.size());

    cout << Date() << " Splitting quals by library" << endl;
    for (size_t pairId = 0; pairId < npairs; pairId++ ) {
      int library = pairs.libraryID(pairId);
      if (!libs_to_extract[library]) continue;
      qualsA[library]->add(quals[pairs.ID1(pairId)]);
      qualsB[library]->add(quals[pairs.ID2(pairId)]);
    }
  }

  // Prepare new pairs files

  if (SPLIT_PAIRS == false && WRITE_PAIRS == True) {
    cout << Date() << " Splitting pairing info by library" << endl;
    for (size_t pairId = 0; pairId < npairs; pairId++ ) {
      int library = pairs.libraryID(pairId);
      if (!libs_to_extract[library]) continue;
      longlong id1 = lib_read_count[library]++;
      longlong id2 = lib_read_count[library]++;
      pairs_out[library].addPair(id1, id2, pairs.sep(pairId), pairs.sd(pairId), pairs.libraryName(pairId), true);
    }

    for (size_t i = 0; i < nlibs; ++i) {
      if (!libs_to_extract[i]) continue;
      String pairsNameAB = READS_OUT + "." + pairs.getLibraryName(i) + ".AB.pairs";
      pairs_out[i].Write(pairsNameAB);
    }
  }

  // Close open files

  for (size_t i = 0; i < nlibs; ++i) {
    if (!libs_to_extract[i]) continue;
    readsA[i]->close();
    if (SPLIT_PAIRS) readsB[i]->close();
    if (QUALS) {
      qualsA[i]->close();
      if (SPLIT_PAIRS) qualsB[i]->close();
    }
  }

  // Write library statistics

  ofstream libstats;
  OpenOfstream( libstats, READS_OUT + ".lib_stats");
  pairs.printLibraryStats( libstats );
  libstats.close();
  
  cout << Date() << " Finished" << endl;
  
}
