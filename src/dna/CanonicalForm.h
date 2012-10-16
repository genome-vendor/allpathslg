///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file CanonicalForm.h
 * \author tsharpe
 * \date Aug 30, 2012
 *
 * \brief
 */
#ifndef DNA_CANONICALFORM_H_
#define DNA_CANONICALFORM_H_

#include <iterator>
#include <ostream>

namespace CanonicalForm
{

enum Value { FWD, REV, PALINDROME };

/// Get CanonicalForm when K is a compile-time constant (like for kmers).
/// The Itr type must have base codes as its value_type.
/// You'll get a compiler error if Itr isn't a random access iterator.
template <class Itr, unsigned K>
inline Value getForm( Itr beg, unsigned v=K )
{ if ( K&1 ) return beg[K/2] & 2 ? REV : FWD;
  Itr end(beg+K);
  while ( beg != end )
  { typename Itr::value_type f = *beg;
    typename Itr::value_type r = *--end ^ 3;
    if ( f < r ) return FWD;
    if ( r < f ) return REV;
    ++beg; }
  return PALINDROME; }

template <class Itr, unsigned K>
inline bool isFwd( Itr beg, unsigned v=K )
{ return getForm<Itr,K>(beg) == FWD; }

template <class Itr, unsigned K>
inline bool isRev( Itr beg, unsigned v=K )
{ return getForm<Itr,K>(beg) == REV; }

template <class Itr, unsigned K>
inline bool isPalindrome( Itr beg, unsigned v=K )
{ return (K&1) ? false : getForm<Itr,K>(beg) == PALINDROME; }

/// Get CanonicalForm when the length can only be determined at run-time
/// (like for a bvec).
/// The Itr type must have base codes as its value_type.
/// You'll get a compiler error if Itr isn't a random access iterator.
template <class Itr>
inline Value getForm( Itr beg, Itr end )
{ using std::distance; size_t len = distance(beg,end);
  if ( len & 1 ) return beg[len/2] & 2 ? REV : FWD;
  while ( beg != end )
  { typename Itr::value_type f = *beg;
    typename Itr::value_type r = *--end ^ 3;
    if ( f < r ) return FWD;
    if ( r < f ) return REV;
    ++beg; }
  return PALINDROME; }

inline Value complement( Value form )
{ switch ( form )
  { case FWD: form = REV; break;
    case REV: form = FWD; break;
    case PALINDROME: break; }
  return form; }

inline std::ostream& operator<<( std::ostream& os, CanonicalForm::Value form )
{ return os << "+-|"[form]; }

}
#endif /* DNA_CANONICALFORM_H_ */
