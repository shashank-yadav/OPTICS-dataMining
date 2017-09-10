#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <list>
#include <math.h>
#include <algorithm>
#include <time.h>
#include <boost/heap/fibonacci_heap.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include "nanoflann.hpp"

using namespace nanoflann;
using namespace std;

#include "KDTreeVectorOfVectorsAdaptor.h"

using namespace std;
using namespace boost::heap;

typedef std::vector<std::vector<float> > my_vector_of_vectors_t;
typedef KDTreeVectorOfVectorsAdaptor< my_vector_of_vectors_t, float >  my_kd_tree_t;

struct heap_data
{
    int id;
    float reachability_distance;

    heap_data(int id, float reachability_distance) : id(id), reachability_distance(reachability_distance) {}
};

struct compare_heap_data
{
    bool operator()(const heap_data& n1, const heap_data& n2) const
    {
        return n1.reachability_distance > n2.reachability_distance;
    }
};

using Heap = boost::heap::fibonacci_heap< heap_data, boost::heap::compare< compare_heap_data> >;

struct optics_descriptor
{
	std::size_t index = -1;
	float coreDistance = -1;
	float reachability_distance = -1;
	bool isProcessed = false;
	Heap::handle_type handle;
};

class optics
{
	int dim = 0;
	my_vector_of_vectors_t dataPts;
	vector<optics_descriptor> descriptors;
	vector<int> ordered_list;
	my_kd_tree_t  *mat_index;

public:

	optics();
	~optics();

	void readFile(string filename );
	void createKDTree();
	float getDist(int id1, int id2);
	vector<int> getNeighbors( int id , float eps);
	vector<int> getNeighbors_indexing( int id , float eps, int minPts );
	bool isCore( int id , float eps , int minPts);
	float coreDistance( int id, vector<int> &N , int minPts);
	void runAlgorithm( float eps, int minPts);
	void update( vector<int> &N, int &id, Heap &seeds, float &eps, int &minPts );
	void print();
};