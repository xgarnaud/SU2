#!/bin/sh

find . -name "*.dat" -exec rm {} \;
find . -name "*.vtk" -exec rm {} \;
find . -name "*.csv" -exec rm {} \;
find . -name "*~" -exec rm {} \;
find . -name "*.pyc" -exec rm {} \;
