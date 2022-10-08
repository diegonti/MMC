###Exercise 2.3
import numpy as np
import matplotlib.pyplot as plt 

def f(x): return np.sin(x)/x

range = 100
points1 = 1000
x1 = np.linspace(-range,range, num=points1)  # x array with higher number of points
x2 = np.linspace(-range,range)               # default is 50 points

#Plotting Settings
fig,plots = plt.subplots(2,1)

plots[0].plot(x1,f(x1), color ="teal",label = fr"f($x_1$), points = {points1}")
plots[1].plot(x2,f(x2), color = "orange", label = r"f($x_2$), points = 50")

for i,plot in enumerate(plots):
    plot.set_ylabel(f"y{i+1}")
    plot.legend()
plt.xlabel("x")

plt.show()

#When plotting, the program cannot show the infiniy of rational numbers for the function
#it uses a sample of a certain number of points that "do the job" and the function appears smooth.
#But we have to be carefoul with this fact because as in the example before, using the default number of points with large ranges,
#may not be enough in some cases and important information of the function may be missing.