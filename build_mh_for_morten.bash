#!/bin/bash

meson --buildtype=debugoptimized build
cd build
ninja
cd ..
cp mh_script_for_morten/run.bash


echo " "
echo "==============================================================="
echo " "
echo "Running self tests"
echo " "
echo "==============================================================="
echo " "
ninja test




echo " "
echo "==============================================================="
echo " "
echo "If the above worked (more or less) do the following:
echo " "
echo "    cd build/Tests/EVP"
echo "    ./run_morten.bash"
echo " "
echo "Results are in directory NEW_RUNS"
echo " "
echo "==============================================================="
echo " "



