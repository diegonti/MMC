#/bin/bash
read -p "Input a number: " number1

if [ $number1 -gt 0 ]; then
    echo "Number is positive"
elif [ $number1 -lt 0 ]; then
    echo "Number is negative"
elif [ $number1 -eq 0 ]; then
    echo "Number is 0."
else
    echo "Number not definedi"
    exit 1
fi

read -p "Input a second number:" number2

if [ $number2 -gt $number1 ]; then
    echo "$number2 greater than $number1."
elif [ $number2 -lt $number1 ]; then
    echo "$number2 less than $number1."
elif [ $number2 -eq $number1 ]; then
    echo "$number2 equal to $number1."
else
    echo "Number not defined. Enter integer number."
    exit 1
fi

exit 0
