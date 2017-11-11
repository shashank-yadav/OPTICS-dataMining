# OPTICS-dataMining #
Implementation of OPTICS algorithm for density based clustering

Requirements:
Boost 1.5.9 or above

##How to run:##
compilation: sh 2013CS50799_compile.sh
running: 2013CS50799_run.sh ( Output in result.txt )
plot: 2013CS50799_plot.sh ( Output as plot.png )

#Files included:
nanoflann.hpp: Header only library, used for index tree structure
KDTreeVectorOfVectorsAdaptor.h: Contains definition for making kd tree on a n dimensional vector

#Source files:
optics.h: Contains definitions of optics class
optics.f: Contains all function initializations
main.cpp: Used to call optics
plot.py: For creating the reachability graph


