#!/bin/sh

Build () {
	Clear

	local dir_cmake=$1
	local target=$2

	echo $dir_cmake
	echo $target

	mkdir build
	cd build

	mkdir modules
	cd modules

	cmake -DCMAKE_BUILD_TYPE=RELEASE ../../$dir_cmake

	cp ./../initial/INPUT.txt ./
	make -j2 $target

	cd -- "$(find -name "$target" -type f -printf '%h' -quit)"
	./$target
}

Clear() {
	if [ -d build ]
		then
		rm -rf ./build
	fi
}

OPT=$1
MODUL=$2

case $OPT in
	-K|-k)
		case $MODUL in
			1)
				Build modules/Kirill/Euler Euler
				;;
			2)
				Build modules/Kirill/RungeKutta RungeKutta
				;;
			3)
				Build modules/Kirill/SpMatrix euler_SpMatrix
				;;
			4)
				Build modules/Kirill/SpMatrix runge_SpMatrix
				;;
			*)
				echo "'1' - Euler; '2' - RungeKutta; '3' - SpMatrix Euler; '4' - SpMatrix RungeKutta" 

		esac
		;;
	*)
		echo "'-K' - Kirill" 
    	exit 1
    	;; 
esac
