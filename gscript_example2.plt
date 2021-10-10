# Gnuplot script file for emf_multilayer example2.out
set terminal postscript eps color enhanced "Arial" 26 size 7in,14in
set output "RTI_example2.eps"

file1="power_coefficients.txt"
file2="intensity_distributions.txt"

set multiplot layout 3,1

set key outside
set key spacing 1.5

#set xrange [-2 : 2]
set yrange [-0.1 : 1.1]
#set xtics -2, 2.0, 2
set ytics -0, 0.25, 1

set xlabel "{incident angle} ({/Symbol \260})"
set ylabel "power coefficient"

plot file1 u 1:2 t "{/Arial-Italic R}_{/Arial-Italic p}" w l lw 4, \
     file1 u 1:3 t "{/Arial-Italic R}_{/Arial-Italic s}" w l lw 4, \
     file1 u 1:4 t "{/Arial-Italic T}_{/Arial-Italic p}" w l lw 4, \
     file1 u 1:5 t "{/Arial-Italic T}_{/Arial-Italic s}" w l lw 4

set pm3d map
set nokey
set yrange [-5 : 5]
set ytics -5, 1, 5
set xlabel "{incident angle} ({/Symbol \260})"
set ylabel "{/Arial-Italic z}"

set title "electric field intensity on the z-axis"
splot file2 u 1:2:3

set title "magnetic field intensity on the z-axis"
splot file2 u 1:2:4

unset multiplot
set terminal x11
