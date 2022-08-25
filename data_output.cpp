#include "data_output.h"

#include <iostream>

#include "TString.h"
#include "TVector3.h"


// constructor
data_output::data_output(const char* output_file_name, const bool include_timings, const bool include_reflected, const char* original_file_name) : include_reflected{include_reflected} {

	// create file

	output_file = new TFile(original_file_name, "UPDATE", "Output File");
	if (output_file->IsOpen()) {
		std::cout << "Output file created successfully." << std::endl << std::endl;
	}
	else {
		std::cout << "Output file could not be opened, check file name is valid." << std::endl << std::endl; exit(1);
	}

	// create trees
	test_tree = new TTree("ScintSim_tree", "Scintillation Light Simulation");
	test_tree->Branch("total_time_vuv", &total_time_vuv); // A vector of three-vectors of hit locations
	//test_tree->Branch("LightYield", &LightYield);

}

// destructor
data_output::~data_output(){
	// deleting output_file also deletes all trees properly
	delete output_file;
}


void data_output::add_data_till(const std::vector<std::vector<double>> &times_vuv) {
    total_time_vuv = times_vuv;
    test_tree->Fill();
}
