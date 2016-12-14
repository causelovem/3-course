set title "Time of computation"
set term postscript eps enhanced color
set style line 1 lt 1 lw 50
set xrange [1:4]
set yrange [0:1160]
set xtics ("1" 1, "8" 2, "64" 3, "125" 4)
set term "postscript" eps
set output "time.eps"
plot "time1.dat" u 1:2 w l title "1024", \
"time2.dat" u 1:2 w l title "2048", \
"time3.dat" u 1:2 w l title "4096", \
"time4.dat" u 1:2 w l title "8192"
