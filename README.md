#This is the dark matter repo for work done @ CUHK under supervision of Professor M.C. Chu and Dr. P.K. Leung

#Compilation:
g++ -std=c++11 -O3 -fopenmp main.cpp rk.cpp -lgsl -lgslcblas -o a.out

#Running:
export OMP_NUM_THREADS=1 # Use only one therad for testing
./a.out

#Notes:
# Testing:
# Create the file directory ./checks in this folder to run the ./a.out checks
mkdir checks
./a.out checks

# Scattered files:
# The only c++ files are main.cpp and the rk.h and rk.cpp files; the .nb files are a little unorganized currently but they involve various derivations for e.g. 4-velocity,  density integrals, limits,  and they are mainly documented by their names.



