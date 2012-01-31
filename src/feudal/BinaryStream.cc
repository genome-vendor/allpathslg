///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file BinaryStream.cc
 * \author tsharpe
 * \date Aug 13, 2009
 *
 * \brief
 */
#include "feudal/BinaryStream.h"
#include "system/ErrNo.h"
#include "system/Exit.h"
#include <iostream>

void BinaryWriter::fail( char const* operation )
{
    ErrNo err;
    std::cout << "BinaryWriter failed to " << operation << " the file "
              << mFilename << err << std::endl;
    CRD::exit(1);
}

void BinaryWriter::write( void const* buf, size_t len )
{
    do
    {
        using std::min;
        ssize_t result = ::write(mFD, buf, min(len,1ul<<30));
        if ( result > 0 )
        {
            len -= result;
            buf = reinterpret_cast<char const*>(buf) + result;
        }
        else if ( result != -1 || errno != EINTR )
            fail("write");
    }
    while ( len );
}

void BinaryReader::readLoop( char* buf, size_t len )
{
    size_t remain = 0;
    while ( len )
    {
        remain = fillBuf();
        if ( !remain )
        {
            std::cout
               << "BinaryReader attempted to read past the end of file "
               << mFR.getFilename() << std::endl;
            CRD::exit(1);
        }
        if ( remain > len )
            remain = len;
        memcpy(buf, mpBuf, remain);
        buf += remain;
        len -= remain;
    }
    mpBuf += remain;
}

void BinaryReader::testToken()
{
    MagicToken tok;
    if ( !read(&tok).isValid() )
    {
        std::cout << "Reading binary file " << mFR.getFilename()
                  << " failed: Initial token is invalid." << std::endl;
        CRD::exit(1);
    }
}
