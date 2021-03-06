// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#pragma once


/// options for SVLocusGraph edge iteration and noise edge filtration
struct EdgeOptions
{
    EdgeOptions() :
        binCount(1),
        binIndex(0),
        isLocusIndex(false),
        locusIndex(0),
        graphNodeMaxEdgeCount(10)
    {}

    unsigned binCount; ///< divide all edges in the graph into binCount bins of approx equal complexity
    unsigned binIndex; ///< out of binCount bins, iterate through the edges in this bin only

    bool isLocusIndex; ///< if true, generate candidates for a specific SVgraph locus only, and ignore binCount/binIndex
    unsigned locusIndex; ///< if isLocusIndex, report this locus only

    unsigned graphNodeMaxEdgeCount; ///< if both nodes of an edge have an edge count higher than this, then skip evaluation of this edge, set to 0 to turn this filtration off
};
