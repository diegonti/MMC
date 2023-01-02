"""
Reads output from MC program and performs block average calculations (bining)
to estimate the statistical error. Saves the plots of the sigma vs. m.
Plots also timeseries of E and M for the first 1e3 and 1e5 MCS.

Diego Ontiveros
"""
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import time as cpu_time
import warnings
warnings.filterwarnings("ignore") # To avoid scipy optimize warnings
to = cpu_time.time()

###################### FUNCTION DEFINITION ######################

def blockAverage(data, maxBlockSize=None):
	"""Computes the block average of a timeseries "x", and 
	provides error bounds for the estimated mean <x>. 

	Parameters
	--------------------
	`data`  : Time series array of an observable X
	`maxBlocksize` : Maximum number of observations/block

	Returns
	--------------------
	`m_points` : Points used (observations/block) for computing averages
	`blockVar` : Variances for each blocksize array 
	`blockMean` : Variances for each blocksize array 	
	"""
 
	Nobs = len(data)           # total number of observations in data
	minBlockSize = 1           # min 1 observation/block
 
	if maxBlockSize is None: maxBlockSize = int(Nobs/4)   # max: 4 blocks (otherwise can't calc variance)
  
	# m_points = 2**n until being less of the inputed maxblocksize
	power = np.arange(int(np.log(maxBlockSize)/np.log(2)))
	m_points = 2**power

	NumBlocks = len(m_points)   				# total number of block sizes

	blockMean = np.zeros(NumBlocks)             # mean for each blocksize
	blockVar  = np.zeros(NumBlocks)             # variance associated with each blockSize


	# Loop for all considered blocksizes (m)
	for k,m in enumerate(m_points):

		Nblock    = int(Nobs/m)               # Number of blocks
		obsProp   = np.zeros(Nblock)          # Container for parcelling block 

		# Loop to chop datastream into blocks and take average
		for i in range(1,Nblock+1):
			
			i1 = (i-1) * m
			i2 =  i1 + m
			obsProp[i-1] = np.mean(data[i1:i2])

		blockMean[k] = np.mean(obsProp)
		blockVar[k]  = np.var(obsProp)/(Nblock - 1)
	
	return m_points, blockVar, blockMean


def plotBlockAverage(m_points,blockVar,blockMean,fit_params=None,save = True,save_name="BlockAverage.jpg",label="x"):
	"""
	Returns a plot of the square root of the block variance and the block means as a function of block size m.

	Parameters
	---------------
	`m_points` : array with the number of samples/block used.
	`blockVar` : array with the block variances gathered.
	`blockMean` : array with the block means gathered.
	`fit_params`: (optional) optimized parameters of the fitting function a-b*exp(-m/t).
	`save`: (optional) Boolean to save or not the image. Defauls to True.
	`save_name`: (optional) File name of the saved image.

	"""

	fig,(ax1,ax2) = plt.subplots(1,2, figsize = (10,5))

	ax1.scatter(m_points, np.sqrt(blockVar),marker="x",c="k",lw=1)
	if fit_params is not None:
		x = np.linspace(m_points[0],m_points[-1])
		ax1.plot(x,fit_function(x,*fit_params),"k:",label="fit")

	ax1.set_xlabel('m')
	ax1.set_ylabel('$\sigma_m$')
	ax1.legend()

	ax2.errorbar(m_points, blockMean, yerr=np.sqrt(blockVar),fmt="k.-",capsize=5)
	ax2.set_ylabel(f'<{label}>')
	ax2.set_xlabel('m')

	fig.tight_layout()

	if save:
		fig.savefig(save_name,dpi=600)


def fit(fit_function,m_points,blockSTD):
	"""Fitting the computed block statistichal errors to the fit_function given"""
	params,cov = sp.optimize.curve_fit(fit_function,m_points,blockSTD)
	return params,cov

def fit_function(x,a,b,tau):
	"""Fit function for the statistichal error points"""
	return a-b*np.exp(-x/tau)


###################### MAIN PROGRAM ######################

# Loading data from test files
file = "out.dat"
L = int(open(file).readline().split("=")[-1])       # Reads L from the first line
dataT = np.loadtxt(file,skiprows=2)                 # Reads the rest of the file
data = dataT.T      # Each parameter in a column
labels = ["E","|M|"]

N = L*L				# Number of spins
t,E,M = data		# data (time,energy,magnetization)
E = E/N				# normalizing E by N
M = M/N 			# normalizing M by N
nMCS = len(t)		# number of montecarlo steps
nmeas = t[1]-t[0]	# MCS between measures

outFile = open("stats_out.dat","a")
# 
start = 1000                                    # If we want to start from a specific datapoint (discardin 1000 first, non equilibrium)
m_points,blockVars,blockMeans = [],[],[]        # Initial empty arrays for each parameter
params = []                                     # Initial empty array for optimized fitting parameters
for i,xi in enumerate(data[1:]/N):              # For each data set (avoiding time)
    n_tot = len(xi)
    # Calculate Block Averages (Bining)
    if labels[i] == "|M|": xi = np.abs(xi)
    m_points_i,blockVar_i,blockMean_i = blockAverage(xi[start:],maxBlockSize=int(n_tot/100))
    m_points.append(m_points_i)
    blockVars.append(blockVar_i)
    blockMeans.append(blockMean_i)

    # Calulate fitting to the block std computed
    params_i, cov_i = fit(fit_function,m_points_i,np.sqrt(blockVar_i))      # Fitting --> Getting optimized params
    x_m = np.linspace(m_points_i[0],m_points_i[-1])                         # linspace for the fitted function (smoother)
    fitted_sigma = fit_function(x_m,*params_i)                              # values for the fitted function
    mean = blockMean_i.mean()                                               # Calculation of the total mean (average of BockMeans)
    std = fitted_sigma[-1]                                                  # Correlated STD will be the one at large BlockSizes (m) (last values of fitted function, plateau)
    params.append(params_i)
    
    # Printing Results
    print(f"\nAverage: <{labels[i]}>/N = {mean:10.5f} +/- {std:.5f}")
    print(f"Optimized fit parameters for {labels[i]}: {params_i}")
    print(f"Statistical Error:    {std:10.8f}")
    print(f"Autocorrelation time: {params_i[-1]:10.8f}")

    outFile.write(f"\nAverage: <{labels[i]}>/N = {mean:10} +/- {std:}\n")
    outFile.write(f"Optimized fit parameters for {labels[i]}: {params_i}\n")
    outFile.write(f"Statistical Error:    {std:10}\n")
    outFile.write(f"Autocorrelation time: {params_i[-1]:10}\n")

# Print first few values of Block Variances
print()
outFile.write("\n")
few = 5
for i,var_i in enumerate(blockVars):
    std_i = np.sqrt(var_i)
    print(f"<{labels[i]}>: First {few} values of Block STD: {std_i[:few]}")
    outFile.write(f"<{labels[i]}>: First {few} values of Block STD: {std_i[:few]}\n")

# Saves Plots
for i in range(len(labels)):
    plotBlockAverage(m_points[i],blockVars[i],blockMeans[i],fit_params=params[i],label = labels[i], save_name = f"BlockAverage{i+1}.jpg")
outFile.close()

# Plotting data time series
fig, ax = plt.subplots(2,2, figsize = (9,8.5))
# Choosing starting and finishing points
start = 0
finish = (np.array((1000,100000))/nmeas).astype(int)    # 1e3 and 1e5 MCS
for i in range(2):
    if i==1: step = 100 
    else: step = 1
    time = t[start:finish[i]:step]
    energy = E[start:finish[i]:step]
    magnetization = M[start:finish[i]:step]

    ax[0,i].plot(time,energy,c="red", lw=1)
    ax[1,i].plot(time,magnetization,c="orange", lw=1)

    ax[0,i].set_xlim(0,time[-1]+step*nmeas);ax[1,i].set_xlim(0,time[-1]+step*nmeas)

    ax[0,i].set_ylabel("Energy");ax[1,i].set_ylabel("Magnetization")
    ax[0,i].set_xlabel("Time (MCS)");ax[1,i].set_xlabel("Time (MCS)")
    ax[0,i].legend(["E"]);ax[1,i].legend(["M"])


fig.tight_layout()
fig.savefig("timeSeries.jpg", dpi=600) # Saving the resulting plot


tf = cpu_time.time()

print(f"\nProcess finished in {tf-to:.2f}s.")
