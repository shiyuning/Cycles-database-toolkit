#!/bin/sh

TEMP_DIR=$1
ARCHIVE=$2
SEVEN_ZIP=$4
PREFIX=$5

cd $TEMP_DIR
../$SEVEN_ZIP a -sdel ../$ARCHIVE $PREFIX_0*.soil
../$SEVEN_ZIP a -sdel ../$ARCHIVE $PREFIX_1*.soil
../$SEVEN_ZIP a -sdel ../$ARCHIVE $PREFIX_2*.soil
../$SEVEN_ZIP a -sdel ../$ARCHIVE $PREFIX_3*.soil
