reset
set terminal qt font ",13"
set size ratio 1
set size 0.8, 0.8
set title "Evolution of the distance with respect to the initial position"
set xl "time"
set yl "distance"
set yr [0:25]
set key bottom right box inside 

set linetype 1 lc rgb "black" lw 5
set linetype 2 lc rgb "orange" lw 5
set linetype 3 lc rgb "green" pt 5 ps 1.5

fit1(x) = a1*x+b1
fit2(x) = a2*x+b2
fit [:20] fit1(x) "Exercise_4.dades.txt" u ($1):((($3)-50)**2+(($2)-50)**2)**0.5  via a1,b1 
fit [60:120] fit2(x) "Exercise_4.dades.txt" u ($1):((($3)-50)**2+(($2)-50)**2)**0.5   via a2,b2 

set terminal postscript enhanced color eps
set output "ex4.eps"
plot "Exercise_4.dades.txt" u ($1):((($3)-50)**2+(($2)-50)**2)**0.5 t "dist" w p lt 3,\
fit1(x) t "fitting 0 to 20" w l lt 1, fit2(x) t "fitting 60 to 120" w l lt 2
