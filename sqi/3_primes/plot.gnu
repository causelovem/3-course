set title "Angel time in sec"
set xrange [1:3]
set xtics ("2 cores" 1, "4 cores" 2, "8 cores" 3)
set term "postscript" eps
set output "plot.eps"
plot "plot.dat" u 1:2 w l
