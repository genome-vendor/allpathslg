///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

const char *DOC =
  "A simple utility to print out source information for read IDs.  Will "
  "iterate back through as much historical read-renumbering information it "
  "has (from modules which use ReadTracker).";

/* ReadTrack
 *
 * Input files required:
 *
 * <READS>.readtrack
 * 
 * Other Arguments:
 * 
 * READ_IDS
 * A set of read IDs in braces (e.g., "{2345,123132,3454363}").
 *
 * Bruce Walker
 * 16 Dec 09
 ******************************************************************************/

#include "MainTools.h"
#include "util/ReadTracker.h"
#include "util/RunCommand.h"

// Return filename portion of path (after last slash)
String basename(const String path)
{
  size_t n = path.rfind("/");
  return path.substr(n+1);
}


int main( int argc, char **argv )
{
  RunTime( );

  BeginCommandArguments;
  CommandDoc(DOC);
  CommandArgument_String(READS);
  CommandArgument_String_OrDefault_Doc(READ_IDS, "",
    "Read IDs to track {id1,id2,id3,...} or track all by default.");
  EndCommandArguments;

  ReadTracker rt, rtnext;
  rt.Load(READS);

  vec<longlong> read_ids;
  uint32_t n_track;
  
  bool track_all = (READ_IDS == "");
  if (track_all) {
    n_track = rt.size();
  } else {
    ParseLongLongSet(READ_IDS, read_ids, False);
    n_track = read_ids.size();
  }

  for (size_t i = 0; i < n_track; ++i) {
    uint32_t r = (track_all ? i : read_ids[i]);
    ForceAssertLt((uint32_t)r, rt.size());
    String source = rt.GetReadSource(r);
    uint32_t index = rt.GetReadIndex(r);
    cout << basename(READS) << ":" << r << " " << basename(source) << ":" << index;
    while (!track_all) {
      rtnext.Load(source);
      if (rtnext.size() == 0) break;
      ForceAssertLt((uint32_t)index, rtnext.size());
      source = rtnext.GetReadSource(index);
      index = rtnext.GetReadIndex(index);
      cout << " " << basename(source) << ":" << index;
    }
    cout << endl;
  }


  
}


