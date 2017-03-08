#!/bin/bash 
PATH_PROGRAM=$(find ../../../modules -name 'implicit_euler_SpMatrix')
echo "$PATH_PROGRAM"
cp "$PATH_PROGRAM" .
for i in 11 12 13 14 15 16 17 18 19 20 21 22
do
	echo 20^$i
	./implicit_euler_SpMatrix $i
	echo
done
