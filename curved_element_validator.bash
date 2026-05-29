#! /bin/bash


# dummy comment

# Setup directories YOU MUST PICK A NAME FOR YOUR OURPUT DIRECTORY.
main_dir=CURVED_ELEMENT_VALIDATOR
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


stem=curved_element_validator
executable=$stem

echo "Executable : "$executable
make $executable

# Stuff to move
important_files="$executable  $stem.cc
                 validate_dpsidn_bubble.gp
                 validate_dpsidn_nodal.gp
                 validate_psi_nodal.gp  
                 validate_psi_bubble.gp"

# Move across
cp $important_files $main_dir

# Go there and create dirs
cd $main_dir

dir_list="
RESLT_boundary_order3_phi0.0_unrotated_coordinates 
RESLT_boundary_order3_phi0.3_unrotated_coordinates 
RESLT_boundary_order3_phi0.0_rotated_coordinates 
RESLT_boundary_order3_phi0.3_rotated_coordinates 
RESLT_boundary_order5_phi0.0_unrotated_coordinates 
RESLT_boundary_order5_phi0.3_unrotated_coordinates 
RESLT_boundary_order5_phi0.0_rotated_coordinates 
RESLT_boundary_order5_phi0.3_rotated_coordinates"

for dir in $dir_list; do
    mkdir $dir
done


# Do it
./$executable $args > OUTPUT 

home_dir=`pwd`
for dir in $dir_list; do
    cd $dir
    oomph-convert soln0.dat
    oomph-convert test_curved_element.dat
    oomph-convert -z test_curved_bell_nodal_basis*.dat
    makePvd test_curved_bell_nodal_basis test_curved_bell_nodal_basis.pvd
    oomph-convert -z test_curved_bell_bubble_basis*.dat
    makePvd test_curved_bell_nodal_basis test_curved_bell_bubble_basis.pvd
    gnuplot -c ../validate_dpsidn_bubble.gp
    gnuplot -c ../validate_dpsidn_nodal.gp
    gnuplot -c ../validate_psi_nodal.gp
    gnuplot -c ../validate_psi_bubble.gp
    cd $home_dir
done

echo " "
echo " "
echo "Done!"
echo " "
echo " "


exit
