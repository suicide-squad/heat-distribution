#!/bin/bash

sh ./clear.sh

echo "build..."

mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=RELEASE ../modules/Kirill/SpMatrix

cp ./../initial/INPUT.txt ./app/

echo "done."

echo "comliler..."

make

echo "done."

cd app/euler
./euler_SpMatrix