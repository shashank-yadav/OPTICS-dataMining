#include "optics.h"

using namespace std;
using namespace boost::heap;

optics::optics(){}
optics::~optics(){}

void optics::readFile(string filename ){

	ifstream f(filename.c_str());
	string line;
	while (getline(f, line))
	{
	    std::istringstream iss(line);
	    float temp;
	    std::vector<float> temp_v;

	    while( iss>>temp ){
	    	temp_v.push_back(temp);
	    }
		
		dataPts.push_back(temp_v);
	}

	dim = dataPts[0].size();
	descriptors.resize( dataPts.size() );
	ordered_list.reserve( dataPts.size() );
	createKDTree();
}

void optics::createKDTree(){

	const size_t nSamples = dataPts.size();
	
	// construct a kd-tree index:
	mat_index = new my_kd_tree_t(dim /*dim*/, dataPts, 10 /* max leaf */ );
	mat_index->index->buildIndex();
}


float optics::getDist(int id1, int id2){
	
	float dist = 0;
	for( int i=0; i<dim; i++){
		float temp = (dataPts[id1][i] - dataPts[id2][i]); 
		dist += temp*temp;
	}
	return dist;
}


vector<int> optics::getNeighbors( int id , float eps ){
	vector<int> nbrs;
	float squared_eps = eps*eps;
	for( int i=0; i<dataPts.size(); i++){
		if( i == id )
			continue;

		if( getDist(id, i) <= squared_eps ){
			nbrs.push_back( i );
		}
	}

	return nbrs;
}

vector<int> optics::getNeighbors_indexing( int id , float eps, int minPts ){
	vector<int> nbrs;
	float squared_eps = eps*eps;
	
	int prod = 4;
	
	while( 1 ){

		const size_t num_results = minPts*prod;
		std::vector<size_t> ret_indexes(num_results);
		std::vector<float> out_dists_sqr(num_results);

		nanoflann::KNNResultSet<float> resultSet(num_results);

		resultSet.init(&ret_indexes[0], &out_dists_sqr[0] );
		mat_index->query( &dataPts[id][0] , num_results, &ret_indexes[0], &out_dists_sqr[0]);
		
		auto it = lower_bound( out_dists_sqr.begin(), out_dists_sqr.end() , squared_eps );

		if( it != out_dists_sqr.end() ){
			int size = std::distance( out_dists_sqr.begin(), it );
			nbrs.resize(size);
			for (int i = 1; i < size; ++i){
				nbrs[i] = ret_indexes[i];
			}
			return nbrs;	
		
		}else{
			
			prod = prod<<1;
		}
	}
}

float optics::coreDistance( int id, vector<int> &N , int minPts){
	
	float dist = getDist( id, N[minPts-2] );
	return dist;
}

void optics::runAlgorithm( float eps, int minPts){

	for( int i=0; i<descriptors.size(); i++){
		
		if( descriptors[i].isProcessed )
			continue;

		vector<int> N = getNeighbors_indexing( i, eps, minPts );

		descriptors[i].isProcessed = true;
		ordered_list.push_back(i);

		if( N.size() >= (minPts - 1) ){
			Heap seeds;
			update( N, i, seeds, eps, minPts );

			while( seeds.size() ){
				// cout<<seeds.size()<<endl;
				heap_data top = seeds.top();
				seeds.pop();						
				vector<int> N_p = getNeighbors_indexing( top.id, eps, minPts );
				descriptors[top.id].isProcessed = true;
				ordered_list.push_back( top.id );

				if( N_p.size() >= minPts - 1 ){
					update( N_p, top.id, seeds, eps, minPts );
				}
			}
		}
	}
}


void optics::update( vector<int> &N, int &id, Heap &seeds, float &eps, int &minPts ){

	float coreDist = coreDistance( id, N, minPts );
	for( int i=0; i<N.size(); i++){
		
		int Ni = N[i]; 
		
		if( descriptors[Ni].isProcessed )
			continue;

		float newReachDistance = sqrt( max( coreDist , getDist(id, Ni) ) );

		if(  descriptors[Ni].reachability_distance == -1 ){
			
			descriptors[Ni].reachability_distance = newReachDistance;
			descriptors[Ni].handle = seeds.push( heap_data(Ni, descriptors[Ni].reachability_distance ) );
		}else{

			if( newReachDistance < descriptors[Ni].reachability_distance ){
				descriptors[Ni].reachability_distance = newReachDistance;
				(*descriptors[Ni].handle).reachability_distance = newReachDistance;
				seeds.update( descriptors[Ni].handle );
			}
		}
	}
}

void optics::print(){
	ofstream f("result.txt");
	
	for (int i = 0; i < ordered_list.size(); ++i){
		f<<descriptors[ordered_list[i]].reachability_distance<<'\n';
	}

	f.close();
}

