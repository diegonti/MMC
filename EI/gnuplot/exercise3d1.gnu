reset
set isosamples 100
set samples 100
set view 47,345
set border 4095
set key left outside

set xr [-20:20];set yr [0:8];set zr [-2:2]
set ztics 1

set xl "x"; set yl "y"; set zl "z"
set title "erf(x)+cos(y), View: 47,345"

set hidden3d
set pm3d
set cntrparam levels incremental -2,0.5,2
set contour base

set terminal postscript enhanced color eps
set output "ex3d1.eps
splot erf(x)+cos(y)
