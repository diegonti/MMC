reset
set xr [1:5]
set yr [-1.5:4]
set title "Figure 1"
set xl "x"
set yl "y"
set key bmargin horizontal center box
set border lw 2

set linetype 1 lc rgb "dark-violet" lw 1
set linetype 2 lc rgb "sea-green" lw 3
set linetype 3 lc rgb "cyan" lw 3

set terminal postscript enhanced color eps
set output "ex1.eps"
plot cos(x**2) lt 1, atan(x)/x**2 lt 2, sin(x)*exp(1/x) lt 3
