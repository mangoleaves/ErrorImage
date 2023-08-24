#!/bin/bash

COMPILER=$1
LANGUAGE=$2
BUILD_TYPE=$3

# Exit script on any error
set -e 

OPTIONS=""
MAKE_OPTIONS=""
BUILDPATH=""

# set GTEST path
OPTIONS="-DGTEST_ROOT=/usr/src/gtest/"

if [ "$COMPILER" == "gcc" ]; then
  echo "Building with GCC";
  BUILDPATH="gcc"

  # without icecc: no options required
  OPTIONS="$OPTIONS -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DCMAKE_C_COMPILER=/usr/bin/gcc"
  MAKE_OPTIONS="-j16"
  export ICECC_CXX=/usr/bin/g++ ; export ICECC_CC=/usr/bin/gcc

elif [ "$COMPILER" == "clang" ]; then

  OPTIONS="$OPTIONS -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang"
  echo "Building with CLANG";
  BUILDPATH="clang"  
fi  

if [ "$LANGUAGE" == "cpp98" ]; then
  echo "Building with C++98";
  BUILDPATH="$BUILDPATH-cpp98"
elif [ "$LANGUAGE" == "cpp11" ]; then
  echo "Building with C++11";
  OPTIONS="$OPTIONS -DCMAKE_CXX_FLAGS='-std=c++11' "
  BUILDPATH="$BUILDPATH-cpp11"  
elif [ "$LANGUAGE" == "cpp14" ]; then
  echo "Building with C++14";
  OPTIONS="$OPTIONS -DCMAKE_CXX_FLAGS='-std=c++14' "
  BUILDPATH="$BUILDPATH-cpp14"  
fi  

#=====================================
# Color Settings:
#=====================================
NC='\033[0m'
OUTPUT='\033[0;32m'
WARNING='\033[0;93m'

if [ "$BUILD_TYPE" == "release" ]; then
    export BUILD_TYPE=release
    BUILDPATH="$BUILDPATH-release"  
else
    export BUILD_TYPE=debug
    BUILDPATH="$BUILDPATH-debug"  
fi