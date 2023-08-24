#!/bin/bash

#Exit on any error
set -e 

LANGUAGE=$1
BUILD_TYPE=$2

PATH=$PATH:/opt/local/bin
export PATH

OPTIONS=""

# set GTEST path
OPTIONS="$OPTIONS -DGTEST_ROOT=~/sw/gtest-1.7.0/"

if [ "$LANGUAGE" == "C++98" ]; then
  echo "Building with C++98";
  BUILDPATH="cpp98"
elif [ "$LANGUAGE" == "C++11" ]; then
  echo "Building with C++11";
  OPTIONS="$OPTIONS -DCMAKE_CXX_FLAGS='-std=c++11' "
  BUILDPATH="cpp11"  
elif [ "$LANGUAGE" == "C++14" ]; then
  echo "Building with C++14";
  OPTIONS="$OPTIONS -DCMAKE_CXX_FLAGS='-std=c++14' "
  BUILDPATH="cpp14"
fi  

#=====================================
# Color Settings:
#=====================================
NC='\033[0m'
OUTPUT='\033[0;32m'
WARNING='\033[0;93m'

if [ "$BUILD_TYPE" == "release" ]; then
    export BUILD_TYPE=Release
    export BUILD_TYPE_L=release
else
    export BUILD_TYPE=Debug
    export BUILD_TYPE_L=debug
fi
