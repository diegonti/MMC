###Exercise 2.6
import numpy as np
import matplotlib.pyplot as plt

#Random Numbers Section
n = 10000                                   #Number of random numbers generated
n_bins = 20                                 #Number of bins (boxes) in the histogram
random_numbers = np.random.normal(size=n)   #Array with the random numbers

mean = np.mean(random_numbers)              #Average of the random numbers
std = np.std(random_numbers)                #Standard Deviation of the random numbers
print(f"Average of the {n} random numbers: {mean:.3f}")
text = fr"$\mu$ = {mean:.3f}"+"\n"+fr"$\sigma$ = {std:.3f}"          #Text to plot

#Normal distribution Section
x = np.linspace(-4,4,100)
gaussian = 1/(std*(2*np.pi)**0.5)*np.exp(-0.5*((x-mean)/std)**2)

#Plotting and Plot settings
plt.hist(random_numbers,rwidth=0.95,density=True,bins=n_bins,color="teal")     #Histogram of the random numbers
plt.plot(x,gaussian, color="red", alpha=0.5)                                   #Plot of Gaussian function with the random numbers mean and std
plt.axvline(mean, ls="--",lw=1,color="red")                                    #Vertical line at sample mean

plt.xlim(-4,4)                                       
plt.xlabel("Random value");plt.ylabel("Probabiity density")

props = dict(boxstyle='round', facecolor='wheat', alpha=0.85)
plt.text(0.05, 0.85, text, bbox=props,transform=plt.gca().transAxes)
plt.title(f"{n} random generated numbers")
plt.show()

#When a high number of sample points and bins are taken, one can see that the nummpy.random.normal() function
#does generate a gaussian distribution centered in 0.