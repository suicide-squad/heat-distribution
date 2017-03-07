#!/bin/bash 
#указываем где у нас хранится bash-интерпретатор 
parametr1=$1 #присваиваем переменной parametr1 значение первого параметра скрипта
script_name=$0 #присваиваем переменной script_name значение имени скрипта

data=20
data2=7

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
do
echo 20^$i
./JacobiMethod $i
echo
done #цикл окончен
exit 0



