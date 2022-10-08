import numpy as np
import matplotlib.pyplot as plt

n = 10000                                   #Number of random numbers generated
n_bins = 20                                 #Number of bins (boxes) in the histogram
random_numbers = np.random.uniform(size=n)  #Array with the random numbers

mean = np.mean(random_numbers)              #Average of the random numbers
print(f"Average of the {n} random numbers: {mean:.3f}")
text = f"Mean = {mean:.3f}"

#Plotting and Plot settings
plt.hist(random_numbers,rwidth=0.95,density=True,bins=n_bins,color="teal")
plt.axhline(1, color = "red")                       #density=True returns the probability density 
plt.xlim(0,1)                                       #calculated with counts/(sum(counts)*bin.width)
plt.xlabel("Random value")                          #such that the area under the histogram is 1.
plt.ylabel("Probability density")
props = dict(boxstyle='round', facecolor='wheat', alpha=0.85)
plt.text(0.05, 0.15, text,verticalalignment='top', bbox=props)
plt.title(f"{n} random generated numbers")
plt.show()