###Exercise 2.2

import sys

#Epsilon calculation
print("Calculated epsilon and baundaries:\n")
eps = 1
while eps+1.0>1.0:
    eps = eps*0.5
print(f"Calculated epsilon = {eps}")

#In Python, int numbers have no boundaries, they are limited by the memory of the computer.
#Nevertheless, floating numbers do have an upper and lower bounds
n = 0    #exponent
a = 10   #base
while True:
    try: x_positive = float(a**n)
    except OverflowError: 
        print(f"Upper bound of float numbers: ~ +10^{n}")
        break
    n +=1

n = 0
while True:
    x_negative = float(a**-n)
    if x_negative == 0.0:
        print(f"Lower bound of float numbers: ~ 10^-{n}")
        break
    n += 1
   

#In fact the exact boundaries and epsilon can be called with sys.float_ino:
print("\nInformation for float numbers (from sys module):\n")
epsilon = sys.float_info.epsilon/2 #
print(f"Epsilon: {epsilon}")

boundaries = (-sys.float_info.max,sys.float_info.max)
print(f"Exact boundaries: {boundaries}")

#And the exact minimum representable positive non-normalized float:
print(f"Minimum representable positive float: {sys.float_info.min*sys.float_info.epsilon}")
#(sys.float_info.min is normalized with epsilon, multiply to get the "real" min value)
