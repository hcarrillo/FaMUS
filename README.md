Introduction
---------------
This repository contains the accompanying code of our IROS'12 paper Fast Minimum Uncertainty Search on a Graph Map Representation.

In the paper we takle the problem of how to plan the minimum uncertainty path in a roadmap like structure?.

We have proposed a fast path planning algorithm capable of obtaining the minimum uncertainty path according to a reduced representation of the environment using a determinant-based criterion (e.g. D-opt).

The code in this repository is based on the code used in the experiments report in the paper. Moreover, It can be used as a starting point in an minimum uncertainty navigation framework based on SLAM algorithms.

If you use this work, please cite our corresponding paper :

@INPROCEEDINGS{Carrillo2012,
author = {H. Carrillo and Y. Latif and J. Neira and J. A. Castellanos},
title = {{Fast Minimum Uncertainty Search on a Graph Map Representation}},
booktitle = {IEEE / RSJ International Conference on Intelligent Robots and Systems, IROS’12},
year = {2012},
address = {Vilamoura, Algarve, Portugal},
month = {October}
}

Requirements
---------------
Suitesparse ( Avalible in http://www.cise.ufl.edu/research/sparse/SuiteSparse/ Version 3.4.0 )
Also with sudo apt-get install libsuitesparse-dev
g2o (avalible in http://openslam.org/g2o Rev. 30)
If you get "No rule to make target `/usr/lib/x86_64-linux-gnu/libGL.so'" (probably becasue you are using 64 bits linux and a NVIDIA card). Check the symbolic link of libGL.so and point it to the correct place.
Eigen3 (avalible in http://eigen.tuxfamily.org/index.php?title=Main_Page Version 3.1.1)
Nanoflann (avalible in http://code.google.com/p/nanoflann/ Version 1.1.3)
Boost graph 1.42 (avalible in http://www.boost.org/)
Version 1.42 is strongly recommended
Although g2o was used for the experiment, as stated in the paper others graph-SLAM based implementations could be used (e.g. iSAM), as long as the correct interface is coded.

Installation
---------------
Install all the requirements.
Go to the build folder
cd build
compile tesp.cpp
make
Usage
Usage: ./FaMuS g2o_filename startVertexID endVertexID
-./FaMuS ../Data/intel_opt.g2o 10 500
The outputs in the command line are:
The nodes of the shortest path and the minimum uncertainty path.
The percentage of overlap between the shortest path and the minimum uncertainty path.
Overall time of the reduction process.
Final uncertainty value of the shortest path.
Also if WRITE_MATLAB is defined, a MATLAB script with the resulting paths will be generated. In the folder named "Results" is available a MATLAB script (i.e. plottingFunctions.m) to plot figures presented in the paper.

Contact
---------------
If you have any queries please contact:

Henry Carrillo L.
http://webdiis.unizar.es/~hcarri/
