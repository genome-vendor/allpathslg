///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "paths/Sepdev.h"

template digraphE<sepdev>::digraphE(digraphE<sepdev> const&, vec<int> const&);
template void digraphE<sepdev>::ToLeft ( vec<int>& to_left ) const;
template void digraphE<sepdev>::ToRight( vec<int>& to_right ) const;
template void digraphE<sepdev>::ComponentEdges( vec< vec<int> >& edges ) const;
template void digraphE<sepdev>::Initialize( const vec< vec<int> >& from, const vec< vec<int> >& to, const vec<sepdev>& edges, const vec< vec<int> >& to_edge_obj, const vec< vec<int> >& from_edge_obj );
template void digraphE<sepdev>::writeBinary( BinaryWriter& ) const;
template void digraphE<sepdev>::readBinary( BinaryReader& );

template digraphE<fsepdev>::digraphE(digraphE<fsepdev> const&, vec<int> const&);
template void digraphE<fsepdev>::AddVertices( int n_vertices );
template void digraphE<fsepdev>::AddEdge( const int v, const int w, const fsepdev & E );
template void digraphE<fsepdev>::ToLeft ( vec<int>& to_left ) const;
template void digraphE<fsepdev>::ToRight( vec<int>& to_right ) const;
template void digraphE<fsepdev>::ComponentEdges( vec< vec<int> >& edges ) const;

template const Tsepdev<double>& digraphE<Tsepdev<double> >::EdgeObject(int) const;
template void digraphE<Tsepdev<double> >::writeBinary( BinaryWriter& ) const;
template void digraphE<Tsepdev<double> >::readBinary( BinaryReader& );
