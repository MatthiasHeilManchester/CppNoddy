#!/bin/bash


echo " "
echo "==============================================================="
echo " "
echo "Before you start building this make sure that you have lapack"
echo "installed in a publically visible place. On ubuntu do this with"
echo " " 
echo "    sudo  apt-get install liblapack-dev"
echo " "
echo "==============================================================="
echo " "
read -p "Hit return to continue"

echo " "
echo "Righto, let's build this thing!"
echo " "

meson --buildtype=debugoptimized build
cd build
ninja
cd ..
cp mh_script_for_morten/run_morten.bash build/Tests/EVP


echo " "
echo "==============================================================="
echo " "
echo "Running self tests"
echo " "
echo "==============================================================="
echo " "
cd build
ninja test




echo " "
echo "==============================================================="
echo " "
echo "If the above worked (more or less) do the following:"
echo " "
echo "    cd build/Tests/EVP"
echo "    ./run_morten.bash"
echo " "
echo "Results are in directory NEW_RUNS (but note that this script"
echo "needs some paraview scripts in ~/bin)"
echo " "
echo "==============================================================="
echo " "



