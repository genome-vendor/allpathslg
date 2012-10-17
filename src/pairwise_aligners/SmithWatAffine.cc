///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// SmithWatAffine( S, T )
//
// Not optimized

#include "Basevector.h"
#include "math/Functions.h"
#include "math/Array.h"

#include "Alignment.h"
#include "ShortVector.h"
#include "pairwise_aligners/SmithWatAffine.h"


#include <omp.h>
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

unsigned int SmithWatAffine( const basevector& S, const basevector& T, 
			     alignment& a,
			     bool penalize_left_gap,
			     bool penalize_right_gap,
                             const int mismatch_penalty,
                             const int gap_open_penalty,
                             const int gap_extend_penalty )
{
     ForceAssert( penalize_left_gap );
     ForceAssert( penalize_right_gap );

     ForceAssertGt( S.size(), 0u );
     ForceAssertGt( T.size(), 0u );

     unsigned int n = S.size( ), N = T.size( );

     //     ForceAssertLe( n, N );

     avector<char> s, t;
     s.resize(n);
     for ( unsigned int i = 0; i < n; i++ )
          s(i) = S[i];
     t.resize(N);
     for ( unsigned int i = 0; i < N; i++ )
          t(i) = T[i];

     const int Infinity = 100000000;
     int best_score = Infinity;
     vec< vec<unsigned int> > score_x;
     vec< vec<unsigned int> > score_y;
     vec< vec<unsigned int> > score_z;

     vec< vec<unsigned char> > x_from;
     vec< vec<unsigned char> > y_from;
     vec< vec<unsigned char> > z_from;

     score_x.resize( n+1 );
     score_y.resize( n+1 );
     score_z.resize( n+1 );
     x_from.resize( n+1 );
     y_from.resize( n+1 );
     z_from.resize( n+1 );

     for ( unsigned int i = 0; i <= n; ++i )
     {
         score_x[i].resize( N+1 );
         score_y[i].resize( N+1 );
         score_z[i].resize( N+1 );
         x_from[i].resize( N+1 );
         y_from[i].resize( N+1 );
         z_from[i].resize( N+1 );
     }
     
     score_x[0][0] = 0;
     score_y[0][0] = Infinity;
     score_z[0][0] = Infinity;
     x_from[0][0] = 's';
     y_from[0][0] = 's';
     z_from[0][0] = 's';

     for ( unsigned int i = 1; i <= n; i++ )
     {    score_x[i][0] = Infinity;
	  score_y[i][0] = Infinity;
	  score_z[i][0] = gap_open_penalty + gap_extend_penalty * i;
	  x_from[i][0] = 's';
	  y_from[i][0] = 's';
	  z_from[i][0] = 's';   } 

     for ( unsigned int j = 1; j <= N; j++) 
       {  score_x[0][j] = Infinity;
          score_y[0][j] = gap_open_penalty + gap_extend_penalty * j;
	  score_z[0][j] = Infinity; 
	  x_from[0][j] = 's';
	  y_from[0][j] = 's';
	  z_from[0][j] = 's';   } 

     for ( unsigned int i = 1; i <= n; i++ )
     {   for ( unsigned int j = 1; j <= N; j++ )
	  {    unsigned int x_x = score_x[i-1][j-1] + mismatch_penalty * ( s(i-1) != t(j-1) );
	       unsigned int x_y = score_y[i-1][j-1] + mismatch_penalty * ( s(i-1) != t(j-1) );
	       unsigned int x_z = score_z[i-1][j-1] + mismatch_penalty * ( s(i-1) != t(j-1) );
	       unsigned int y_x = score_x[i][j-1] + gap_open_penalty;
	       unsigned int y_y = score_y[i][j-1] + gap_extend_penalty;
	       unsigned int y_z = Infinity; //score_z[i][j-1] + gap_open_penalty;
	       unsigned int z_x = score_x[i-1][j] + gap_open_penalty;
	       unsigned int z_y = Infinity; //score_y[i-1][j] + gap_open_penalty;
	       unsigned int z_z = score_z[i-1][j] + gap_extend_penalty;

	       score_x[i][j] = Min( Min( x_x, x_y ), x_z ); 		               
	       score_y[i][j] = Min( Min( y_x, y_y ), y_z );
	       score_z[i][j] = Min( Min( z_x, z_y ), z_z );    

	       if ( x_x <= x_y )
	       {    if ( x_x <= x_z ) x_from[i][j] = 'x';
	            else x_from[i][j] =  'z';    }
	       else
	       {    if ( x_y <= x_z ) x_from[i][j] = 'y';
	            else x_from[i][j] =  'z';    }

	       if ( y_x <= y_y )
	       {    if ( y_x <= y_z ) y_from[i][j] = 'x';
	            else y_from[i][j] =  'z';    }
	       else
	       {    if ( y_y <= y_z ) y_from[i][j] = 'y';
	            else y_from[i][j] =  'z';    }

	       if ( z_x <= z_y )
	       {    if ( z_x <= z_z ) z_from[i][j] = 'x';
	            else z_from[i][j] =  'z';    }
	       else
	       {    if ( z_y <= z_z ) z_from[i][j] = 'y';
	            else z_from[i][j] =  'z';    }    }    }

     best_score = Min( score_x[n][N], Min( score_y[n][N], score_z[n][N] ) );

     vec< vec<unsigned char> > *from;
     if ( score_x[n][N] <= score_y[n][N] )
     {    if ( score_x[n][N] <= score_z[n][N] ) from = &x_from;
          else from = &z_from;    }
     else
     {    if ( score_y[n][N] <= score_z[n][N] ) from = &y_from;
          else from = &z_from;    }
     
     int i = int(n);
     int j = int(N);
     int lcount = 0, g1count = 0, g2count = 0;
     int last_length = 0;
     avector<int> gaps(0), lengths(0);
     while(1)
     {  
          unsigned char dir = (*from)[i][j];
	  //cout << dir;
          if ( from == &x_from )
          {    if ( g1count > 0 )
               {    if ( last_length > 0 )
		    {    gaps.Prepend( g1count );
		         lengths.Prepend( last_length );    }
	            g1count = 0;    }
               if ( g2count > 0 )
               {    if ( last_length > 0 )
	            {    gaps.Prepend( -g2count );
                         lengths.Prepend( last_length );    }
                    g2count = 0;    }
               ++lcount;
               --i;
               --j;   }
          else if ( from == &z_from )  // gap on long sequence
          {    if ( lcount > 0 )
               {    last_length = lcount;
                    lcount = 0;    }
               ForceAssert( g1count == 0 );
               ++g2count;
               --i;    }
          else                           // gap on short sequence
          {    if ( lcount > 0 )
               {    last_length = lcount;
                    lcount = 0;    }
               ForceAssert( g2count == 0 );
               ++g1count;
               --j;    }

	  if ( dir == 'x') from = &x_from;
	  else if ( dir == 'y') from = &y_from;
	  else from = &z_from;

	  if( (*from)[i][j] == 's' ) break;

    }

     //cout << "\n";

     if ( g1count != 0 ) gaps.Prepend( g1count );
     else if ( g2count != 0 ) gaps.Prepend( -g2count );
     else gaps.Prepend(0);

     lengths.Prepend( lcount );

     int pos1 = i;
     int pos2 = j;

     if ( gaps(0) < 0 )
     {   pos2 -= gaps(0);
         gaps(0) = 0;    }

     if ( gaps(0) > 0 )
     {   pos1 += gaps(0);
         gaps(0) = 0;    }

     int errors = best_score;
     a = alignment( pos1, pos2, errors, gaps, lengths );

     return best_score;    }


unsigned int SmithWatAffineParallel(const basevector & S, const basevector & T,
                            alignment & a,
                            bool penalize_left_gap,
                            bool penalize_right_gap,
                            const int mismatch_penalty,
                            const int gap_open_penalty,
                            const int gap_extend_penalty
                            )
{
    //cout <<  Date( ) << ": SmithWatAffineParallel start" << endl;
    //double clock = WallClockTime( );
    ForceAssert(penalize_left_gap);
    ForceAssert(penalize_right_gap);

    ForceAssertGt(S.size(), 0u);
    ForceAssertGt(T.size(), 0u);

    unsigned int n = S.size(), N = T.size();

    //     ForceAssertLe( n, N );

    avector < char >s, t;
    s.resize(n);
    for (unsigned int i = 0; i < n; i++)
        s(i) = S[i];
    t.resize(N);
    for (unsigned int i = 0; i < N; i++)
        t(i) = T[i];

    const int Infinity = 100000000;
    int best_score = Infinity;

    RecArray<int> score_x( n+1, N+1 );
    RecArray<int> score_y( n+1, N+1 );
    RecArray<int> score_z( n+1, N+1 );

    RecArray<unsigned char> x_from( n+1, N+1 );
    RecArray<unsigned char> y_from( n+1, N+1 );
    RecArray<unsigned char> z_from( n+1, N+1 );

    score_x[0][0] = 0;
    score_y[0][0] = Infinity;
    score_z[0][0] = Infinity;
    x_from[0][0] = 's';
    y_from[0][0] = 's';
    z_from[0][0] = 's';

    for (unsigned int i = 1; i <= n; i++) {
        score_x[i][0] = Infinity;
        score_y[i][0] = Infinity;
        score_z[i][0] = gap_open_penalty + gap_extend_penalty * i;
        x_from[i][0] = 's';
        y_from[i][0] = 's';
        z_from[i][0] = 's';
    }

    for (unsigned int j = 1; j <= N; j++) {
        score_x[0][j] = Infinity;
        score_y[0][j] = gap_open_penalty + gap_extend_penalty * j;
        score_z[0][j] = Infinity;
        x_from[0][j] = 's';
        y_from[0][j] = 's';
        z_from[0][j] = 's';
    }


    // parallelize the smith-waterman calculation.
    // The idea is to divide the n+1 by N+1 matrix into many (BxB) blocks. 
    // Every block will only depend on the block on it's left, up, and up-left
    // blocks, and all the block(i,j) with same i + j can be calculated parallelized.
    //
    // block size = 200 is chosen heuristically, which seems most efficient
    //
    const int block_size = 200;
    int n_block_s = n / block_size;
    if ( n_block_s * block_size < (int)n ) n_block_s += 1;
    int n_block_t = N / block_size;
    if ( n_block_t * block_size < (int)n ) n_block_t += 1;
    //cout <<  Date( ) << ": start the main loop" << endl;
    //cout <<  Date( ) << ": time used = " << TimeSince(clock) << endl;
    for ( int level = 0; level < n_block_s + n_block_t; ++level ) {
    // number of blocks in the diagonal that are independent to each other
    // start and stop of the block index
    // i_block_s in [0, n_block_s)
    // i_block_t in [0, n_block_t)
    // where i_block_s + i_block_t = level

    // in other words i_block_s is in ( level - n_block_t, level ] && [0, n_block_s)
    
    #pragma omp parallel for 
    for ( int i_block_s = max( 0, level - n_block_t + 1 ); 
          i_block_s < min( n_block_s, level + 1 ); i_block_s++ ) {
        int i_block_t = level - i_block_s;
        ForceAssertGe( i_block_t , 0 );
        ForceAssertLt( i_block_t , n_block_t );

        int istart = i_block_s * block_size + 1;
        int istop = min( istart + block_size, (int)n + 1 );
        int jstart = i_block_t * block_size + 1;
        int jstop = min( jstart + block_size, (int)N + 1 );


    for ( int i = istart; i < istop; i++) {
        for ( int j = jstart; j < jstop; j++) {
            unsigned int x_x =
                score_x[i - 1][j - 1] + mismatch_penalty * (s(i - 1) !=
                                                            t(j - 1));
            unsigned int x_y =
                score_y[i - 1][j - 1] + mismatch_penalty * (s(i - 1) !=
                                                            t(j - 1));
            unsigned int x_z =
                score_z[i - 1][j - 1] + mismatch_penalty * (s(i - 1) !=
                                                            t(j - 1));
            unsigned int y_x = score_x[i][j - 1] + gap_open_penalty;
            unsigned int y_y = score_y[i][j - 1] + gap_extend_penalty;
            unsigned int y_z = Infinity;        //score_z[i][j-1] + gap_open_penalty;
            unsigned int z_x = score_x[i - 1][j] + gap_open_penalty;
            unsigned int z_y = Infinity;        //score_y[i-1][j] + gap_open_penalty;
            unsigned int z_z = score_z[i - 1][j] + gap_extend_penalty;

            score_x[i][j] = Min(Min(x_x, x_y), x_z);
            score_y[i][j] = Min(Min(y_x, y_y), y_z);
            score_z[i][j] = Min(Min(z_x, z_y), z_z);

            if (x_x <= x_y) {
                if (x_x <= x_z)
                    x_from[i][j] = 'x';
                else
                    x_from[i][j] = 'z';
            }
            else {
                if (x_y <= x_z)
                    x_from[i][j] = 'y';
                else
                    x_from[i][j] = 'z';
            }

            if (y_x <= y_y) {
                if (y_x <= y_z)
                    y_from[i][j] = 'x';
                else
                    y_from[i][j] = 'z';
            }
            else {
                if (y_y <= y_z)
                    y_from[i][j] = 'y';
                else
                    y_from[i][j] = 'z';
            }

            if (z_x <= z_y) {
                if (z_x <= z_z)
                    z_from[i][j] = 'x';
                else
                    z_from[i][j] = 'z';
            }
            else {
                if (z_y <= z_z)
                    z_from[i][j] = 'y';
                else
                    z_from[i][j] = 'z';
            }
        }
    }



    // end parallelization 
    }
    }

    best_score = Min(score_x[n][N], Min(score_y[n][N], score_z[n][N]));

    //vec < vec < unsigned char > >*from;
    RecArray<unsigned char> *from;
    if (score_x[n][N] <= score_y[n][N]) {
        if (score_x[n][N] <= score_z[n][N])
            from = &x_from;
        else
            from = &z_from;
    }
    else {
        if (score_y[n][N] <= score_z[n][N])
            from = &y_from;
        else
            from = &z_from;
    }

    int i = int (n);
    int j = int (N);
    int lcount = 0, g1count = 0, g2count = 0;
    int last_length = 0;
    avector < int >gaps(0), lengths(0);
    while (1) {
        unsigned char dir = (*from)[i][j];
        //cout << dir;
        if (from == &x_from) {
            if (g1count > 0) {
                if (last_length > 0) {
                    gaps.Prepend(g1count);
                    lengths.Prepend(last_length);
                }
                g1count = 0;
            }
            if (g2count > 0) {
                if (last_length > 0) {
                    gaps.Prepend(-g2count);
                    lengths.Prepend(last_length);
                }
                g2count = 0;
            }
            ++lcount;
            --i;
            --j;
        }
        else if (from == &z_from)       // gap on long sequence
        {
            if (lcount > 0) {
                last_length = lcount;
                lcount = 0;
            }
            ForceAssert(g1count == 0);
            ++g2count;
            --i;
        }
        else                    // gap on short sequence
        {
            if (lcount > 0) {
                last_length = lcount;
                lcount = 0;
            }
            ForceAssert(g2count == 0);
            ++g1count;
            --j;
        }

        if (dir == 'x')
            from = &x_from;
        else if (dir == 'y')
            from = &y_from;
        else
            from = &z_from;

        if ((*from)[i][j] == 's')
            break;

    }

    //cout << "\n";

    if (g1count != 0)
        gaps.Prepend(g1count);
    else if (g2count != 0)
        gaps.Prepend(-g2count);
    else
        gaps.Prepend(0);

    lengths.Prepend(lcount);

    int pos1 = i;
    int pos2 = j;

    if (gaps(0) < 0) {
        pos2 -= gaps(0);
        gaps(0) = 0;
    }

    if (gaps(0) > 0) {
        pos1 += gaps(0);
        gaps(0) = 0;
    }

    int errors = best_score;
    a = alignment(pos1, pos2, errors, gaps, lengths);

    return best_score;
}






