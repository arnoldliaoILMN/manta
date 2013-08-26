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

#include "ReadGroupStatsSet.hh"

#include "blt_util/log.hh"
#include "blt_util/io_util.hh"
#include "blt_util/parse_util.hh"
#include "blt_util/string_util.hh"

#include <fstream>
#include <iostream>


// serialization
void
ReadGroupStatsSet::
save(const char* filename) const
{
	using namespace boost::archive;

	assert(NULL != filename);
    std::ofstream ofs(filename, std::fstream::binary);
    boost::archive::text_oarchive oa(ofs);

	const unsigned numGroups(_group.size());
	for (unsigned i(0); i<numGroups; ++i)
	{
		oa << _group.get_key(i);
		oa << getStats(i);
	}
}

// restore from serialization
void
ReadGroupStatsSet::
load(const char* filename)
{
	using namespace boost::archive;

	clear();

	assert(NULL != filename);
	std::ifstream ifs(filename, std::fstream::binary);
	boost::archive::text_iarchive ia(ifs);

	while (ifs.peek() != EOF)
	{
		std::string bamFile;
		ReadGroupStats rgs;
		ia >> bamFile;
		ia >> rgs;
		setStats(bamFile, rgs);
	}
}


void
ReadGroupStatsSet::
read(const char* filename)
{

}



void
ReadGroupStatsSet::
write(std::ostream& os) const
{
    const unsigned n_groups(_group.size());
    for (unsigned i(0); i<n_groups; ++i)
    {
        os << "# Bam_Path\t" << i << "\t" << _group.get_key(i) << '\n';
    }
    // write column header for better readability
    os << "# index"
       << "\treadOrientation"
       << "\tsample-count\tnumber-of-fragment-sizes"
       << '\n';

    for (unsigned i(0); i<n_groups; ++i)
    {
        os << i << '\t';
        getStats(i).write(os);
        os << '\n';
    }
}

