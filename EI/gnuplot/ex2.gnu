reset
set sample 100
set xr [0:10]
set yr [-0.6:1]
set title "Figure 2"
set xl "x"
set yl "y"
set key top right reverse box title "Functions"
set border lw 2
set label "Arrow to (4,0.5)" at 6,0.7
set arrow from 6,0.7 to 4,0.5

set terminal postscript enhanced color eps
set output "ex2.eps"
plot  sin(x)*cos(x) w lp lt 1 pt 5, 2/(x+1) w l lw 2 lc "sea-green"
