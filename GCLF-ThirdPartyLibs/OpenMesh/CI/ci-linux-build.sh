#!/bin/bash
source CI/ci-linux-prepare.sh

echo -e "${OUTPUT}"
echo ""
echo "======================================================================"
echo "Basic configuration details:"
echo "======================================================================"
echo -e "${NC}"

echo "Compiler:     $COMPILER"
echo "Options:      $OPTIONS"
echo "Language:     $LANGUAGE"
echo "Make Options: $OPTIONS"
echo "BuildPath:    $BUILDPATH"
echo "Path:         $PATH"
echo "Language:     $LANGUAGE"

echo -e "${OUTPUT}"
echo ""
echo "======================================================================"
echo "Building $BUILD_TYPE version with vectorchecks enabled"
echo "======================================================================"
echo -e "${NC}"

if [ ! -d build-$BUILDPATH-Vector-Checks ]; then
  mkdir build-$BUILDPATH-Vector-Checks
fi

cd build-$BUILDPATH-Vector-Checks

cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DOPENMESH_BUILD_UNIT_TESTS=TRUE -DSTL_VECTOR_CHECKS=ON $OPTIONS ../

#build it
make $MAKE_OPTIONS

cd ..