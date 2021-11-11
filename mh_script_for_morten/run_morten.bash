#! /bin/bash

# Stuff to move
important_files=" "

cd ../..
ninja
cd -

# Setup directories
main_dir=NEW_RUNS
if [ -e $main_dir ]; then
    echo "Directory " $main_dir "already exists -- please move it out of the way
."
    exit
fi
mkdir $main_dir

cp $important_files $main_dir
cd $main_dir



half_width_list="10 20 30"
n_waves_list="1 2 3"
for half_width in `echo $half_width_list`; do
    
    for n_waves in `echo $n_waves_list`; do
        ../EVPOrrSommerfeld_lapack_morten --n_waves_in_box $n_waves --y_min -$half_width --y_max $half_width
        oomph-convert DATA/two_d_field.dat
        mv DATA DATA_n_waves`echo $n_waves`_half_width`echo $half_width`
    done
    
    ~/bin/make_combined_pvd.bash DATA_n_waves*_half_width`echo $half_width`/two*vtu
    mv combined.pvd combined_half_width`echo $half_width`.pvd
    
done

for n_waves in `echo $n_waves_list`; do
    ~/bin/make_combined_pvd.bash DATA_n_waves`echo $n_waves`_half_width*/two*vtu
    mv combined.pvd combined_n_waves`echo $n_waves`.pvd
done

