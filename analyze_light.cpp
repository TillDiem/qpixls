// main driver for scintillation toy mc simulation code
#include<string>
#include<iostream>
#include<fstream>
#include<chrono>
#include <sstream>
#include <vector>
#include <algorithm>
#include "TH1.h"
#include "TRandom.h"
#include "TVector3.h"

#include "data_output.h"
#include "semi_analytic_hits.h"
#include "time_parameterisation.h"
#include "utility_functions.h"
#include "radiological_parameters.h"

// include parameter file
#include "simulation_parameters.h"
using namespace std;

int main() {

	gRandom->SetSeed(0);

	// -------- Initialise semi-analytic hits class ---------
	semi_analytic_hits hits_model;

	// -------- Initialise timing parametrisation class ---------
	time_parameterisation times_model(parameters::timing_discretisation_step_size);

	// -------- Initialise utility/energy spectrum class ---------
	utility_functions utility;
	utility.initalise_scintillation_functions_argon(parameters::t_singlet, parameters::t_triplet, parameters::singlet_fraction_electron, parameters::triplet_fraction_electron,
        parameters::singlet_fraction_alpha, parameters::triplet_fraction_alpha, parameters::scint_time_window);
        utility.initalise_scintillation_functions_xenon(parameters::t_singlet_Xe, parameters::t_triplet_Xe, parameters::singlet_fraction_Xe, parameters::triplet_fraction_Xe, parameters::scint_time_window);

	// ------- Read photon detector positions and types --------
	std::vector<std::vector<int>> opdet_type;
	std::vector<int> opdet_direction;
	std::vector<std::vector<double>> opdet_position;

       	std::cout << "Loading Photon Detector positions..." << std::endl;
        std::ifstream detector_positions_file;
	//detector_positions_file.open("optical_detectors_dune1x2x6.txt");
        //detector_positions_file.open("optical_detectors_QPIX4mm.txt");
        detector_positions_file.open("out");
        if(detector_positions_file.is_open()) std::cout << "File opened successfully" << std::endl;
        else {std::cout << "File not found." << std::endl; exit(1);}
        while(!detector_positions_file.eof()) {
		int num_opdet, type_opdet, direction; double x_opdet, y_opdet, z_opdet;
		if(detector_positions_file >> num_opdet >> x_opdet >> y_opdet >> z_opdet >> type_opdet >> direction) {
		    std::vector<int> type({num_opdet, type_opdet});
		    std::vector<double> position({x_opdet, y_opdet, z_opdet});
		    opdet_type.push_back(type);
		    opdet_position.push_back(position);
		    opdet_direction.push_back(direction);
		}
		else{ break; }
        }
	detector_positions_file.close();
	int number_opdets = opdet_type.size();
	std::cout << "Positions Loaded: " << number_opdets << " optical detectors." << std::endl << std::endl;

	// ------- Read G4 simulation data --------
	// READ IN THE G4 SIMULATION  - CURRENTLY HARDCODED
	char* G4OutputFileName = "single_electron.root";
	TFile * G4OutputFile = new TFile(G4OutputFileName);

	TTree *G4OutputTree = (TTree*)G4OutputFile->Get("event_tree");

	double energy_deposit;
	vector <double> *hit_start_x = nullptr;
	vector <double> *hit_start_y = nullptr;
	vector <double> *hit_start_z = nullptr;
	vector <double> *hit_start_t = nullptr;
	vector <double> *hit_energy_deposit = nullptr;
	vector <double> *hit_length = nullptr;
	double pixel_size;

	G4OutputTree->SetBranchAddress("energy_deposit", &energy_deposit);
	G4OutputTree->SetBranchAddress("hit_start_x", &hit_start_x);
	G4OutputTree->SetBranchAddress("hit_start_y", &hit_start_y);
	G4OutputTree->SetBranchAddress("hit_start_z", &hit_start_z);
	G4OutputTree->SetBranchAddress("hit_length", &hit_length);
	G4OutputTree->SetBranchAddress("hit_start_t", &hit_start_t);
	G4OutputTree->SetBranchAddress("hit_energy_deposit", &hit_energy_deposit);
	G4OutputTree->SetBranchAddress("hit_length", &hit_length);

	int NEventsToLoopOver =  100; //G4OutputTree->GetEntries(); // 100000

	data_output output_file(parameters::output_file_name, parameters::include_timings, parameters::include_reflected, G4OutputFileName );

	for (int EventIt=0; EventIt < NEventsToLoopOver; EventIt++)
	{
		//if(EventIt!=1894) continue;
		std::cout << "Event: " << EventIt << std::endl;
		G4OutputTree->GetEntry(EventIt);

		int max_events = hit_start_x->size();

		// Vector inlcuding all the hit positions in (x,y,z)
		std::vector<std::vector<double>> position_list(max_events, std::vector<double>(3,-999.9));

		// Vector including the number of photons per SiPM
		vector<int> num_VUV_array;

		// Vector including the current hits x,y,z positoin
		vector<TVector3> ScintPoint_array;

		// Vector of vectors.
		// [ SiPM1[t0, t1, t2,t3,...], SiPM2[t0, t1, t2,t3,...], SiPM3[t0, t1, t2,t3,...], ...]
		// For each SiPM there should be a vector. This vector will contain the time of each photon.
		vector<vector<double> > total_time_vuv_array ;//(number_opdets, vector<double>);
		vector<vector<double>> op_channel_pos(number_opdets, vector<double>(3,0.0));
		//total_time_vuv_array[i] = vector_ofHits;

		// Go through every SiPM
		for(int op_channel = 0; op_channel < number_opdets; op_channel++) {

			// get optical detector type - rectangular or disk aperture
			int op_channel_type = opdet_type[op_channel][1];
			// get optical detector direction - (x, y, z)
			int op_direction = opdet_direction[op_channel];

			// get detection channel coordinates (in cm)
			TVector3 OpDetPoint(opdet_position[op_channel][0],opdet_position[op_channel][1],opdet_position[op_channel][2]);
			vector<double> op_det_pos = {OpDetPoint.X(), OpDetPoint.Y(), OpDetPoint.Z()};

			std::vector<double> total_time_vuv;


			// Go through every hit in the event
			for (int i=0; i < hit_start_x->size(); i++)
			{
		   		position_list[i][0] = hit_start_x->at(i);
		    		position_list[i][1] = hit_start_y->at(i);
		    		position_list[i][2] = hit_start_z->at(i);



		    		double time_hit = hit_start_t->at(i); // in ns
				time_hit*=0.001; // in us
		    		double energy_deposit = hit_energy_deposit->at(i);
				double light_yield = hits_model.LArQL(energy_deposit, hit_length->at(i), parameters::electric_field);

		    		int number_photons;
		    		number_photons = utility.poisson(light_yield * energy_deposit, gRandom->Uniform(1.), energy_deposit);

		    		double singlet_fraction = parameters::singlet_fraction_electron;
		    		double triplet_fraction = parameters::triplet_fraction_electron;

		    		// Setting the point of energy deposition
		    		TVector3 ScintPoint(position_list[i][0],position_list[i][1],position_list[i][2]);

				// TODO: Make sure that normalisation fits!
				int num_VUV = hits_model.VUVHits(number_photons, ScintPoint, OpDetPoint, op_channel_type, 0, op_direction);


				if(num_VUV== 0) { continue; } // forces the next iteration

					// Calculate the angle between the scinitllation point and the optical detector
				double distance_to_pmt = (OpDetPoint-ScintPoint).Mag();
				double cosine = -999.0;

				// TODO: Correct and check
				if(op_direction == 1)  cosine = sqrt(pow(ScintPoint[0] - OpDetPoint[0],2)) / distance_to_pmt;
				else if(op_direction == 2)  cosine = sqrt(pow(ScintPoint[1] - OpDetPoint[1],2)) / distance_to_pmt;
				else if(op_direction == 3)  cosine = sqrt(pow(ScintPoint[2] - OpDetPoint[2],2)) / distance_to_pmt;
				else { std::cout << "Error: Optical detector direction not defined." << std::endl; }

				double theta = acos(cosine)*180./3.14159;

				// Due to rounding errors, this can return a value like 90.000x01 which will throw an error in the angle_bin histogram.
				// Cap the angle to 90 degrees.
				if (theta >= 90){ theta=89.99999; }
				int angle_bin = theta/45;       // 45 deg bins


				// For each photon we will store its arival time in the SiPM currently being looped over.
				// Returns the transport time that the photon takes from the sicntillation point to the detector
				std::vector<double> transport_time_vuv = times_model.getVUVTime(distance_to_pmt, angle_bin, num_VUV);

				// For each photon get the time info
				for(auto& x: transport_time_vuv) {
					// emission time
					double emission_time;
					emission_time = utility.get_scintillation_time_electron()*1000000.0; // in us
					double total_time = time_hit+(x*0.001 + emission_time + 2.5*0.001); // in microseconds
					total_time_vuv.push_back(total_time);
				}// End for transporttime
			// combine timings into single vectors for Direct and Reflected light
		    } // end of hit
		total_time_vuv_array.push_back(total_time_vuv);

		} // End of SiPM Loop

	output_file.add_data_till(total_time_vuv_array);
	} // end event loop

     //Write to OUTPUT FILE
     output_file.write_output_file();
     std::cout << "Program finished." << std::endl;
}
