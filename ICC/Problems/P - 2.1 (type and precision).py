###Exercise 2.1

def case1(b,c):
    """Case 1: a,b,c=int"""
    b,c= int(b),int(c)
    a = int(b/c)
    print(f"Case 1: When a,b,c = int, a = {a}")

def case2(b,c):
    """Case 2: a,b,c=float"""
    b,c = float(b),float(c)
    a = float(b/c)
    print(f"Case 2: When a,b,c = float, a = {a}")

def case3(b,c):
    """Case 3: a=float,b,c=int"""
    b,c = int(b),int(c)
    a = float(b/c)
    print(f"Case 3: When a = float and b,c = int, a = {a}")

#Main Program
b,c = 1,2
case1(b,c) #Same as b//c
case2(b,c) #Same as b/c
case3(b,c)

#Se puede ver como a diferencia de Fortran, en Python, aunque dos numeros sean enteros
#al dividirse entre ellos, el resultado pasa a ser racional (a no ser que sea forzado con la función int()).
#Es importante tener eso en cuenta ya que si no, puede haber truncamientos donde no queremos, o al reves.
#En este sentido, Python es más "humano" ya que lo hace más como alguien con una calculadora, en cambio 
#Fortran, al declarar el tipo de variable, es mas fiel a esa declaración, por lo tanto int/int da int, no real.

#En ambos lenguajes, poner un punto despues del numero hace que sea float (o real), dandole más precisión y evitando 
#errores de precisión.