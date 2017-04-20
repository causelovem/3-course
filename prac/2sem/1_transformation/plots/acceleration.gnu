set title "Acceleration"
set term postscript eps enhanced color
set style line 1 lt 1 lw 50
set xrange [1:4]
set yrange [0:3]
set xtics ("1" 1, "2" 2, "4" 3, "8" 4)
set term "postscript" eps
set output "acceleration.eps"
plot "acceleration.dat" u 1:2 w l title "acceleration N = 20, K = 20"
