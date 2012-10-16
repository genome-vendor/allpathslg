///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file ValidateReadLengths.cc
 * \author tsharpe
 * \date Mar 21, 2012
 *
 * \brief
 */

#include "MainTools.h"
#include "Basevector.h"
#include "feudal/FeudalFileReader.h"
#include <algorithm>
#include <cstring>
#include <iostream>

int main( int argc, char** argv )
{
    RunTime( );
    BeginCommandArguments;
    CommandArgument_String( NAME );
    EndCommandArguments;

    bool err = false;
    FeudalFileReader fRdr((NAME+".fastb").c_str());
    FeudalFileReader qRdr((NAME+".qualb").c_str());
    if ( fRdr.getNElements() != qRdr.getNElements() )
    {
        std::cout << "There are " << fRdr.getNElements()
                    << " reads in the fastb file, and " << qRdr.getNElements()
                    << " reads in the qualb file" << std::endl;
        err = true;
    }
    size_t nnn = std::min(fRdr.getNElements(),qRdr.getNElements());
    typedef bvec::size_type bsiz_t;
    for ( size_t idx = 0; idx < nnn; ++idx )
    {
        bsiz_t nCalls = 0;
        memcpy(&nCalls,fRdr.getFixedData(idx,sizeof(bsiz_t)),sizeof(bsiz_t));
        size_t nQuals = qRdr.getDataLen(idx);
        if ( nCalls != nQuals )
        {
            std::cout << "For read " << idx << " there are " << nCalls
                        << " calls in the fastb file, and " << nQuals
                        << " quals in the qualb file" << std::endl;
            err = true;
        }
    }
    if ( !err )
        std::cout << "Everything is fine." << std::endl;
}
