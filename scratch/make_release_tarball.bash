#!/usr/bin/env bash

set -o nounset
set -o pipefail
set -o xtrace


package_name=manta

pname_root=""
if [ $# -gt 1 ]; then
    cat <<END
usage: $0 [tarball_rootname]

This script makes the manta source tarball

- script assumes that it is located in the release checkout version
- the tarball is written to the caller's working directory
END
    exit 2
elif [ $# == 1 ]; then
    pname_root=$1
fi

rel2abs() {
    (cd $1; pwd -P)
}

thisdir=$(rel2abs $(dirname $0))
outdir=$(pwd)

cd $thisdir
gitversion=$(git describe | sed "s/^v//")
if [ $? != 0 ]; then
    echo "ERROR: 'git describe' failed" 1>&2
    exit 1
fi

if [ "$pname_root" == "" ]; then
    pname_root=${package_name}-$gitversion
fi

pname=$outdir/$pname_root

if [ -d $pname ]; then
    echo "ERROR: directory already exists: $pname"
    exit 1
fi
mkdir -p $pname

cd ..
git archive --prefix=$pname_root/ HEAD: | tar -x -C $outdir

# make version number substitutions:
source_dir=$(pwd)/starka
tmp_file=$(mktemp)
cml=$pname/src/CMakeLists.txt
awk -v gver=$gitversion '{if($1=="set\(MANTA_VERSION") printf "set(MANTA_VERSION \"%s\")\n",gver; else print;}' $cml >| $tmp_file
mv $tmp_file $cml

# tar it up:
(
cd $outdir
tar -f $pname_root.tar.bz2 -cj $pname_root
)

rm -rf $pname
