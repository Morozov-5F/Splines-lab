set grid

set key outside
set xlabel "x"
set ylabel "y(x)"

plot "output.txt" title "Interpolated function",\
"output.txt" using 1:3 title "Exact function" w l
