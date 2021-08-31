rm ./include/common.c
rm *.dat
cp ./common.c ./include/common.c
cp ./common.h ./include/common.h

cd ./include
make
cd ../

./include/longwaves.exe


mv ./include/*.dat ./
