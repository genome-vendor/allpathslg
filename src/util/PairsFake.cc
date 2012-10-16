///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//  
//  Generates a fake .pairs file for unpaired reads.
//

#include "MainTools.h"
#include "Basevector.h"
#include "PairsManager.h"


int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String_Doc(HEAD, "looks for base vectors in '<HEAD>.fastb' and generated '<HEAD>.pairs'");
  EndCommandArguments;

  cout << "Creating pairs file '" << HEAD << ".pairs'." << endl << endl;
  PairsManager pairs(MastervecFileObjectCount(HEAD + ".fastb"));
  pairs.Write(HEAD + ".pairs");

}
