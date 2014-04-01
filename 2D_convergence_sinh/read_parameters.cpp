   // C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <sstream>

// Basic include file needed for the mesh functionality.
#include "libmesh.h"

#include "assemble.h"

using namespace std;


void read_parameters(EquationSystems& es,  int& argc, char**& argv){

if (argc >2){
		es.parameters.set<Real> ("n_timesteps") = atoi( argv[1] );
		es.parameters.set<Real> ("N_eles") = atoi( argv[2] );
		es.parameters.set<std::string> ("output_file_name") = argv[3] ;
		es.parameters.set<std::string> ("result_file_name") =  argv[4] ;

	if (argc >5){
   es.parameters.set<Real> ("DELTA") = atof( argv[5] );
	}else{
	es.parameters.set<Real> ("DELTA") = 1;
	}
	
	if (argc >6){
   es.parameters.set<Real> ("DELTA_BC") = atof( argv[6] );
   
	}else{
	es.parameters.set<Real> ("DELTA_BC") = 100;
	}

}else{

	es.parameters.set<Real> ("n_timesteps") = 4;
	es.parameters.set<Real> ("N_eles") = 12;
	es.parameters.set<std::string> ("output_file_name") = "data/test_2D.mat";
	es.parameters.set<std::string> ("result_file_name") = "data/Test3D_neumann";
	es.parameters.set<Real> ("DELTA") = 0.001;

	
  es.parameters.set<Real> ("DELTA_BC") = 100000000;

}


  std::cout<<"n_timesteps "<< es.parameters.get<Real>("n_timesteps") <<" \n";
  std::cout<<"N_eles "<< es.parameters.get<Real>("N_eles") <<" \n";
  std::cout<<"output_file_name "<< es.parameters.get<std::string>("output_file_name") <<" \n";
  std::cout<<"result_file_name "<< es.parameters.get<std::string>("result_file_name") <<" \n";
  std::cout<<"DELTA "<< es.parameters.get<Real>("DELTA") <<" \n";
  std::cout<<"DELTA_BC "<< es.parameters.get<Real>("DELTA_BC") <<" \n";


}
