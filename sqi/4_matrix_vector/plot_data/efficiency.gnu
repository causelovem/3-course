set title "Efficiency"
set term postscript eps enhanced color
set style line 1 lt 1 lw 50
set xrange [0:6]
set yrange [0:1.1]
set xtics ("1" 1, "32" 2, "64" 3, "128" 4, "256" 5, "512" 6)
set term "postscript" eps
set output "efficiency.eps"
plot "efficiency1.dat" u 1:2 w l title "1 cores", \
"efficiency2.dat" u 1:2 w l title "32 cores", \
"efficiency3.dat" u 1:2 w l title "64 cores", \
"efficiency4.dat" u 1:2 w l title "128 cores", \
"efficiency5.dat" u 1:2 w l title "256 cores", \
"efficiency6.dat" u 1:2 w l title "512 cores"