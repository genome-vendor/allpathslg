///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ARRAY_H
#define ARRAY_H

#include "CoreTools.h"
#include "Vec.h"

// ==============================================================================
// A fixed-size simple c-style 2D array that allocates faster than vec< vec<T> >
// The elements are uninitialized.
// 
// No copy or assignments are allowed yet.
//
// Usage: 
//      RecArray<int> r(3,4);
//      r[1][2] = 1;
//      cout << "nrow=" << r.size1() << " ncol=" << r.size2() << end;
// ==============================================================================
template <typename T>
class RecArray {
public:
    // ctor and dtor
    RecArray( int m, int n ) : m_(m), n_(n)  {
        pt_ = new T* [m];
        for( int i = 0; i < m ; i++ )
            pt_[i] = new T [n];
    }
    ~RecArray() { 
        for( int i = 0; i < m_ ; i++ ) {
            delete [] pt_[i];
        }
        delete [] pt_;
    }
    // array operations
    T*& operator[] ( int i ) { return pt_[i]; }
    int size1() const { return m_; }
    int size2() const { return n_; }

private:
    // Disable copy and assign constructors
    RecArray( const RecArray& );
    RecArray& operator=( const RecArray&); 
    T ** pt_;
    int m_;
    int n_;
};




// ==================================================================
// Smooth an 1-dimensional array of double values using an Gaussian. 
// 
// ( The input array is either a vec<double> for normal array,  or 
// a map<int, double> for sparse array )

// =======================================================

template <typename T>
void SmoothArrayGaussian ( const T& distr, int delta,
         T & distr_s )
{
     distr_s = distr;
     SmoothArrayGaussian( distr_s, delta ); // Automatic dispatcher 
}
// With normal array input as a vec<int>
void SmoothArrayGaussian ( vec<double> & distr, int delta );
// Here the array is sparse, and is represented using std::map<int, double>
void SmoothArrayGaussian ( map<int, double> & distr, int delta );

// Return normalized Gaussain function in an array of length ( 2 * width_n_delta * delta + 1 )
// with center  at index  width_n_delta * delta ,  and width_n_delta * delta on boths sides
vec<double> GaussianKernal( int delta, int width_n_delta = 4 );

#endif
