###Exercici 3.1.5

def changeBase(x,n):
    """Takes number (x) and changes it to base n"""
    new_x = ""                      #Changed base number
    quotient = x
    while quotient != 0:            #Loop that gets the quiotient and reminder
        remainder = quotient%n      #for deviding the number by the new base.
        quotient = quotient//n      #The reminder sets the digits of the new number.
        
        new_x += str(remainder)     #Adding the digits to the new number

    return new_x[::-1]              #Returns the reversed string (numbers in the right order)

number = 10
for base in [2,4,8,16,60]:
    new_number = changeBase(number,base)
    print(f"{number} in base {base} is {new_number}")