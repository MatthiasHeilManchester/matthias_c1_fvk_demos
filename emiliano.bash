#! /bin/bash


# Setup directories YOU MUST PICK A NAME FOR YOUR OURPUT DIRECTORY.
main_dir=NEW_RUNS
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


stem=matthias_disk_explorer
postfix_list="_fvk" #  _ks"
for postfix in `echo $postfix_list`; do
    
    executable=$stem$postfix

    echo "Executable : "$executable
    make $executable

    # Stuff to move
    important_files="$executable  $stem.cc"
    
    case_list="_balance_on_edge " # " _free " # _clamped _pinned  _free _balance_on_edge"
    echo $case_list
    for the_case in `echo $case_list`; do

        inner_boundary_list="straight" # " curved"
        for inner_boundary in `echo $inner_boundary_list`; do

            boundary_postfix="_straight"
            boundary_command_line_arg=" --use_polyline_for_internal_boundaries "
            if [ "$inner_boundary" == "curved" ]; then
                boundary_postfix="_curved"
                boundary_command_line_arg=" "
            fi
            
            t_shape_list="not_t_shaped" # "t_shaped not_t_shaped"
            for t_shape in `echo $t_shape_list`; do
                
                t_shape_postfix="_use_t_shape"
                t_shape_command_line_arg=" --use_t_shape_internal_boundaries "
                if [ "$t_shape" == "not_t_shaped" ]; then
                    t_shape_postfix="_no_t_shape"
                    t_shape_command_line_arg=" "
                fi              
                
                rotate_list="rotate not_rotate"
                for rotate in `echo $rotate_list`; do
                    
                    rotate_postfix="_rotated_coords"
                    rotate_command_line_arg=" "
                    if [ "$rotate" == "not_rotate" ]; then
                        rotate_postfix="_unrotated_coords"
                        rotate_command_line_arg=" --do_not_rotate_coords_on_curved_boundaries "
                    fi              
                    
                    the_dir=$main_dir/Case$boundary_postfix$rotate_postfix$t_shape_postfix$postfix$the_case
                    
                    mkdir $the_dir
                    home_dir=`pwd`
                    
                    # Transfer the important files to the main directory and go there
                    cp $important_files $the_dir
                    cd $the_dir
                    
                    reslt_dir=RESLT
                    mkdir $reslt_dir
                    
                    
                    
                    args="--nplot 50 "$boundary_command_line_arg" "$t_shape_command_line_arg" "$rotate_command_line_arg
                    
                    if [ "$the_case" == "_clamped" ]; then
                        args=" --use_clamped_bc "
                    elif [ "$the_case" == "_pinned" ]; then
                        args=" --use_pinned_bc "
                    elif [ "$the_case" == "_balance_on_edge" ]; then
                        args=" --use_balance_on_edge_bc "
                    fi
                    
                    echo " "
                    echo " " 
                    echo "Running: "
                    ls -l $executable
                    echo $args
                    echo " in "`pwd`
                    echo " "
                    echo " "
                    
                    # Do it
                    ./$executable $args > OUTPUT 
                    
                    
                    oomph-convert mesh_black_box_upgrade.dat
                    oomph-convert mesh_before_black_box_upgrade.dat
                    
                    oomph-convert -p2 duplicated_nodes.dat rotated_nodes.dat 
                    oomph-convert elements_upgraded_to_curved.dat rotated_elements.dat split_elements.dat
                    
                    oomph-convert -p2 boundary_coordinate*.dat
                    oomph-convert test_bulk_elements_on_boundary*.dat
                    oomph-convert -p2 test_face_elements_on_boundary*.dat
                    
                    cd RESLT
                    oomph-convert soln*.dat; makePvd soln soln.pvd 
                    oomph-convert full_soln*.dat; makePvd full_soln full_soln.pvd
                    
                    # Next
                    cd $home_dir

                done
            done 
        done
    done
done

echo " "
echo " "
echo "Done!"
echo " "
echo " "


exit
