###Exercise 2.4
import numpy as np
import matplotlib.pyplot as plt

#Arrays and Mesh creation
range = 2                             
x = np.linspace(-range,range)             
y = np.linspace(-range,range)            #np.linspace produces linear arrays.
X,Y = np.meshgrid(x,y)                   #To plot a 3d function, a mesh (grid/matrix) has to be created
Z = np.exp(-(X**2+Y**2))

#Plot Settings
ax=plt.axes(projection = "3d")
ax.plot_surface(X,Y,Z, cmap="plasma")
plt.xlabel("x");plt.ylabel("y")
ax.set_zlabel("z")
plt.title(r"f(x,y) = exp(-($x^2$+$y^2$))")
plt.show()

#If you plot with linear arrays (x,y --> np.linspace()) the function would only be evaluated in those (x,y) linear pairs.
#To do a surface plot, all the combination between the two linear arrays have to be created (X,Y --> matrix), to "see the full picture".