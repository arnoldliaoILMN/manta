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
// <https://github.com/downloads/sequencing/licenses/>.
//

#include "GetAlignmentStats.hh"

#include "AlignmentStatsOptions.hh"

#include "blt_util/log.hh"
#include "common/OutStream.hh"
#include "manta/ReadGroupStatsSet.hh"

#include "boost/foreach.hpp"
#include "boost/archive/text_oarchive.hpp"
#include "boost/archive/text_iarchive.hpp"
#include "boost/serialization/shared_ptr.hpp"
#include "boost/serialization/vector.hpp"
#include "boost/serialization/hash_map.hpp"

#include <cstdlib>
#include <fstream>



static
void
runAlignmentStats(const AlignmentStatsOptions& opt)
{
    // calculate fragment size statistics for all read groups in all bams

    // instantiate early to test for filename/permissions problems
    OutStream outs(opt.outputFilename);
    if (opt.alignmentFilename.empty())
    {
        log_os << "ERROR: No input files specified.\n";
        exit(EXIT_FAILURE);
    }

    ReadGroupStatsSet rstats;
    // debug...
    ReadGroupStatsSet rstatsDeserized;
    char serizedStatsFile[] = "/home/xchen/projs/manta/results/stats/temp.stats";

    BOOST_FOREACH(const std::string& file, opt.alignmentFilename)
    {
        ReadGroupStats rgs = ReadGroupStats(file);
    	rstats.setStats(file,rgs);
    	// debug...
    	std::ofstream outSerialized(serizedStatsFile);
    	boost::archive::text_oarchive oa(outSerialized);
    	oa << rgs;
/*
    	ReadGroupStats rgsNew;
    	std::ifstream inSerialized(serizedStatsFile);
    	boost::archive::text_iarchive ia(inSerialized);
    	ia >> rgsNew;

    	rstatsDeserized.setStats(file, rgsNew);

    	std::cerr << "numOfFragSize=" << rgsNew.fragSize.numOfFragSize
    			  << "\tquantileNum=" << rgsNew.fragSize.quantileNum;
    			  */
    }

    rstats.write(outs.getStream());
    // debug...
    //rstatsDeserized.write(outs.getStream());

}


void
GetAlignmentStats::
runInternal(int argc, char* argv[]) const
{
    AlignmentStatsOptions opt;

    parseAlignmentStatsOptions(*this,argc,argv,opt);
    runAlignmentStats(opt);
}
