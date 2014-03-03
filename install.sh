g++ ./src/ModelDist2.cpp -o ./bin/ModelDist2 -std=gnu++0x -fopenmp
g++ ./src/Overlap19.cpp -o ./bin/Overlap19 -std=gnu++0x -fopenmp
g++ ./src/OverlapRegion2.cpp -o ./bin/OverlapRegion2 -std=gnu++0x -fopenmp
g++ ./src/RUFUSv5.Filter.cpp -o ./bin/RUFUSv5.Filter -std=gnu++0x -fopenmp
g++ ./src/RUFUSv6.BuildHash.cpp -o ./bin/RUFUSv6.BuildHash -std=gnu++0x -fopenmp
cd ./bin/jellyfish-2.0.0/

path=$( pwd )
echo "path=$path"
./configure --prefix=$path
make 
make install
cd ../

