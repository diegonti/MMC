reset

set title "Time evolution of A, B and C"
set xl "time (day)"
set yl "%"
set xr [0:28]
set yr [0:100]
set key center bottom bmargin outside horizontal
set border 1+2
set xtics 3 nomirror
set ytics 25 nomirror
set size ratio 0.5

set linetype 1 lc rgb "red" lw 1 pt 5 ps 2
set linetype 2 lc rgb "green" lw 1 pt 7 ps 2
set linetype 3 lc rgb "blue" lw 1 pt 9 ps 2
set linetype 4 lc rgb "blue" lw 3


cfit(x) = a*x**2+b*x
fit cfit(x) "Exercise_3.txt" u 1:($4)*100/($2+$3+$4) via a,b	

set terminal postscript enhanced color eps
set output "ex3.eps"
plot "Exercise_3.txt" u 1:($2)*100/($2+$3+$4) t "A" w lp lt 1,\
"" u 1:($3)*100/($2+$3+$4) t "B" w lp lt 2,\
"" u 1:($4)*100/($2+$3+$4) t "C" w lp lt 3,\
 cfit(x) t "Cfit" w l lt 4

