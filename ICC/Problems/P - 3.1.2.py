###Exercici 3.1.2
from numpy import e

n=20

#Recursion, stable loop, unstable loop:
E_unstable = 1/e
for i in range(1,n+1):
    E_unstable = (1-i*E_unstable)

#Inverse recursion, stable loop:
E_stable = 0
for i in range(n,1,-1):
    E_stable = (1-E_stable)/(i)
    
#Printing results
print(f"With {n} iterations:")
print(f"Recursive iteration leads to: {E_unstable} (unstable)")
print(f"Inverse recursion leads to: {E_stable} (stable)")
print(f"1/e = {1/e} (from numpy)")

#Relative error between the iterated e and the "real" one
e_iter = 1/E_stable
e_rel = abs(e_iter-e)/(e)
print(f"\nRelative error from stable inverse recursion: {e_rel}")
#With higher n, e_rel will be smaller. From n=18, e_rel=0

