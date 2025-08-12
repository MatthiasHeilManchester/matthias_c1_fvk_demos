#! /bin/bash


executable=matthias_disk_explorer_fvk
make $executable

main_dir=RESLT
if [ -e $main_dir ]; then
    echo " "
    echo "WARNING: Directory " $main_dir " already exists!"
    read -p "         remove it and continue? [Y/n] " yn
    case $yn in
        ''|[Yy]* ) rm -rf $main_dir;;
        [Nn]* ) echo "Can't continue until you move $main_dir"; exit;;
    esac
fi
mkdir $main_dir

# Do it
./$executable --test_damped_solve > OUTPUT 
                            
cd $main_dir
oomph-convert soln*.dat; makePvd soln soln.pvd 
oomph-convert full_soln*.dat; makePvd full_soln full_soln.pvd

echo " "
echo " "
echo "Done!"
echo " "
echo " "

exit
