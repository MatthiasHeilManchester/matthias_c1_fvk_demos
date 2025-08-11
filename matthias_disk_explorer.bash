#! /bin/bash


# dummy comment

# Setup directories YOU MUST PICK A NAME FOR YOUR OURPUT DIRECTORY.
main_dir=MATTHIAS_RUNS
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
    
    outer_boundary_list="curved_curved"  # "curved_curved straight_curved straight_poly"
    for outer_boundary in `echo $outer_boundary_list`; do
        
        outer_boundary_postfix="_outer_curved_curved"
        outer_boundary_command_line_arg=""
        if [ "$outer_boundary" == "straight_curved" ]; then
            outer_boundary_postfix="_outer_straight_curved"
            outer_boundary_command_line_arg=" --outer_boundary_straight_curved"
        fi
        if [ "$outer_boundary" == "straight_poly" ]; then
            outer_boundary_postfix="_outer_straight_poly"
            outer_boundary_command_line_arg=" --outer_boundary_straight_poly"
        fi
        
        el_area_list="0.5 0.1 0.05 0.01" # 0.05"
        for el_area in `echo $el_area_list`; do
            
            case_list=" _free " # _clamped _pinned  _free _balance_on_edge"
            for the_case in `echo $case_list`; do
                
                inner_boundary_list="curved" # "straight" # " curved"
                for inner_boundary in `echo $inner_boundary_list`; do
                    
                    inner_boundary_postfix="_straight"
                    inner_boundary_command_line_arg=" --use_polyline_for_internal_boundaries "
                    if [ "$inner_boundary" == "curved" ]; then
                        inner_boundary_postfix="_curved"
                        inner_boundary_command_line_arg=" "
                    fi
                    
                    t_shape_list="t_shaped" # "not_t_shaped" # "t_shaped not_t_shaped"
                    for t_shape in `echo $t_shape_list`; do
                        
                        t_shape_postfix="_use_t_shape"
                        t_shape_command_line_arg=" --use_t_shape_internal_boundaries "
                        if [ "$t_shape" == "not_t_shaped" ]; then
                            t_shape_postfix="_no_t_shape"
                            t_shape_command_line_arg=" "
                        fi              
                        
                        rotate_list="rotate" # "rotate not_rotate"
                        for rotate in `echo $rotate_list`; do
                            
                            rotate_postfix="_rotated_coords"
                            rotate_command_line_arg=" "
                            if [ "$rotate" == "not_rotate" ]; then
                                rotate_postfix="_unrotated_coords"
                                rotate_command_line_arg=" --do_not_rotate_coords_on_curved_boundaries "
                            fi              
                            
                            the_dir=$main_dir/Case$outer_boundary_postfix$inner_boundary_postfix"_el_area"$el_area$rotate_postfix$t_shape_postfix$postfix$the_case
                            
                            mkdir $the_dir
                            home_dir=`pwd`
                            
                            # Transfer the important files to the main directory and go there
                            cp $important_files $the_dir
                            cd $the_dir
                            
                            reslt_dir=RESLT
                            mkdir $reslt_dir
                            
                            args="--nplot 50 --el_area "$el_area" "$inner_boundary_command_line_arg" "$outer_boundary_command_line_arg" "$t_shape_command_line_arg" "$rotate_command_line_arg
                            
                            if [ "$the_case" == "_clamped" ]; then
                                args=$args" --use_clamped_bc "
                            elif [ "$the_case" == "_pinned" ]; then
                                args=$args" --use_pinned_bc "
                            elif [ "$the_case" == "_balance_on_edge" ]; then
                                args=$args" --use_balance_on_edge_bc "
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
    done
done

echo " "
echo " "
echo "Done!"
echo " "
echo " "


exit
