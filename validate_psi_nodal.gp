reset
set terminal pngcairo size 2400,1200   # or 'wxt' for interactive window
set output 'validate_psi_nodal.png'           # optional


set termoption noenhanced
set multiplot layout 3,6 title GPVAL_PWD . "\n test_curved_bell_curved_edge_nodal_basis (col 8 vs col 7)"
set grid
set xlabel "s"
set ylabel "psi"
#set label GPVAL_PWD at graph 0.01, 0.01
do for [i=0:17] {
    filename = sprintf("test_curved_bell_curved_edge_nodal_basis%d.dat", i)
    set title sprintf("basis%d", i)
    plot filename using 7:8 with lines lw 2 title sprintf("basis %d", i)
}

unset multiplot
