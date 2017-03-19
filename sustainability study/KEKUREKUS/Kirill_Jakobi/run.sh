#!/bin/bash 
PATH_PROGRAM=$(find ../../../modules -name 'implicit_euler_SpMatrix')
echo "$PATH_PROGRAM"
cp "$PATH_PROGRAM" .
for i in 3 4 5 6 7 8 9 10
do
	echo 20^$i
	./implicit_euler_SpMatrix $i
	echo
done
