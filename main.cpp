#include <iostream>
#include <sstream>
#include <fstream>
#include "optics.h"

using namespace std;

int main(int argc, char *argv[])
{
    
	if (argc != 3){
		cout<<"Usage: ./optics minPts epsilon"<<endl;
		return 1;
	}

	int minPts = stoi(argv[1]);
	float eps = strtof( argv[2] , 0 );

	cout<<minPts<<' '<<eps<<endl;

    optics O;
    
    O.readFile("data.tsv");
    O.runAlgorithm(eps, minPts);
    O.print();

    return 0;
}
