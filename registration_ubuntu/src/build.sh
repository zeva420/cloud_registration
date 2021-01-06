cd Thirdparty/g2o

echo "Configuring and building Thirdparty/g2o ..."

mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j

cd ../../../


echo "Configuring and building cloud_rigestration ..."

mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j2
