set title "Efficiency"
set term postscript eps enhanced color
set style line 1 lt 1 lw 50
set xrange [1:4]
set yrange [0:5]
set xtics ("1" 1, "8" 2, "64" 3, "125" 4)
set term "postscript" eps
set output "efficiency.eps"
plot "efficiency1.dat" u 1:2 w l title "1204", \
"efficiency2.dat" u 1:2 w l title "2048"
