set title "Diagram"
set xrange [1:6]
set xtics ("ijk" 1, "ikj" 2, "jki" 3, "jik" 4, "kij" 5, "kji" 6)
set term "postscript" eps
set output "plot.eps"
plot "plot.dat" u 1:2 w l