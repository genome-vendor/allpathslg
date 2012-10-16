///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file KMerContext.h
 * \author tsharpe
 * \date Aug 14, 2012
 *
 * \brief
 */
#ifndef KMERS_KMERCONTEXT_H_
#define KMERS_KMERCONTEXT_H_

#include "dna/Bases.h"
#include "system/Assert.h"


// keeps track of a kmer's predecessor and successor base codes
class KMerContext
{
public:
    KMerContext() : mVal(0) {}
    KMerContext( unsigned char predCode, unsigned char succCode )
    : mVal(pred2Val(predCode)|succ2Val(succCode)) {}

    KMerContext& operator|=( KMerContext const& context )
    { mVal |= context.mVal; return *this; }

    // bit 0 is set if A is a predecessor, bit1=C, bit2=G, bit3=T
    unsigned char getPredecessors() const { return mVal >> 4; }

    bool isPredecessor( unsigned char predCode ) const
    { return pred2Val(predCode) & mVal; }

    void setPredecessor( unsigned char predCode )
    { *this |= finalContext(predCode); }

    void removePredecessor( unsigned char predCode )
    { mVal &= ~pred2Val(predCode); }

    unsigned getPredecessorCount() const
    { return sideCount(getPredecessors()); }

    // bit 0 is set if A is a successor, bit1=C, bit2=G, bit3=T
    unsigned char getSuccessors() const { return mVal & 0xfu; }

    bool isSuccessor( unsigned char succCode ) const
    { return succ2Val(succCode) & mVal; }

    void setSuccessor( unsigned char succCode )
    { *this |= initialContext(succCode); }

    void removeSuccessor( unsigned char succCode )
    { mVal &= ~succ2Val(succCode); }

    unsigned getSuccessorCount() const
    { return sideCount(getSuccessors()); }

    KMerContext rc() const
    { return KMerContext(gRCVals[mVal]); }

    static KMerContext initialContext( unsigned char succCode )
    { return KMerContext(succ2Val(succCode)); }

    static KMerContext finalContext( unsigned char predCode )
    { return KMerContext(pred2Val(predCode)); }

private:
    explicit KMerContext( unsigned char val ) : mVal(val) {}

    static unsigned char pred2Val( unsigned char predCode )
    { return Base::val2Bits(predCode) << 4; }

    static unsigned char succ2Val( unsigned char succCode )
    { return Base::val2Bits(succCode); }

    static unsigned sideCount( unsigned char bits )
    { AssertLe(bits,0xfu); return gSideCounts[bits&0xful]; }

    unsigned char mVal;
    static unsigned char gRCVals[];
    static unsigned gSideCounts[];
};

#endif /* KMERS_KMERCONTEXT_H_ */
