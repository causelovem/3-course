set title "L2"
set term postscript eps enhanced color
set style line 1 lt 1 lw 50
set xrange [0:4]
set xtics ("32x32 ikj" 1, "32x32 ijk" 2, "opbl ijk" 3)
set term "postscript" eps
set output "plots/l2.eps"
plot "plot_data/l2_data.dat" using 1:2 with imp ls 1

