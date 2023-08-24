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
echo "Building $BUILD_TYPE version unittests"
echo "======================================================================"
echo -e "${NC}"

if [ ! -d build-$BUILDPATH-Vector-Checks ]; then
  mkdir build-$BUILDPATH-Vector-Checks
fi

cd build-$BUILDPATH-Vector-Checks

#build the unit tests
make  $MAKE_OPTIONS unittests

echo -e "${OUTPUT}"
echo ""
echo "======================================================================"
echo "Running unittests $BUILD_TYPE version with vectorchecks enabled"
echo "======================================================================"
echo -e "${NC}"

cd Unittests

#execute tests
./unittests --gtest_color=yes --gtest_output=xml

echo -e "${OUTPUT}"
echo ""
echo "======================================================================"
echo "Running unittests $BUILD_TYPE version with custom vector type"
echo "======================================================================"
echo -e "${NC}"

./unittests_customvec --gtest_color=yes --gtest_output=xml

echo -e "${OUTPUT}"
echo ""
echo "======================================================================"
echo "Running unittests $BUILD_TYPE version with double vector type"
echo "======================================================================"
echo -e "${NC}"

#execute tests
./unittests_doublevec --gtest_color=yes --gtest_output=xml

cd ..
cd ..