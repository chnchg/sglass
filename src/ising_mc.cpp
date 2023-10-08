#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <H5Cpp.h>
#include <algorithm>
#include <string>
#include <exception>
#include <memory>
#include <sstream>

// using namespace std;
// Calculation parameters
const int warm_up_time = 10000;
const int calculation_time = 10000;
const int chosen_num = 64;
const int pair_num = chosen_num*(chosen_num-1)/2;
// System parameters
const int linear_size = 32; // Linear size
double beta; // Inverse temperature
double field; // Magnetic field

// State
const int num_spins = linear_size*linear_size;
int spins[num_spins];

struct Selection {
	std::vector<int> chosen_indices;
	std::vector<double> M;
	std::vector<double> aC;	
};
std::vector<Selection> selection_list;

struct Statics{
	std::vector<int> config;
	std::vector<double> E;
	std::vector<double> E2;
	std::vector<double> ME;
	std::vector<double> Mtotal;
};
Statics Statics_list;

std::mt19937 rng(1234);

// Auxiliary
std::uniform_real_distribution<> uniform;
int neighbors[num_spins][4];

std::string format(float f, int digits) {
    std::ostringstream ss;
    ss.precision(digits);
    ss << f;
    return ss.str();
}


void show_sys() {
	for (int i = 0; i < linear_size; i ++) {
		if (i) std::cout << '\n';
		for (int j = 0; j < linear_size; j ++) std::cout << (spins[i * linear_size + j] > 0 ? '+' : '-');
	}
	std::cout << '\n';
}

void initialize() {
	for (int i = 0; i < num_spins; i ++) spins[i] = uniform(rng) < 0.5 ? 1 : -1;

	int lz = linear_size - 1;
	int lnz = linear_size * lz;
	for (int i = 0; i <  num_spins; i ++) {
		// Left
		neighbors[i][0] = i % linear_size ? i - 1 : i + lz;
		// Right
		neighbors[i][1] = (i + 1) % linear_size ? i + 1 : i - lz;
		// Up
		neighbors[i][2] = i < linear_size ? i + lnz : i - linear_size;
		// Down
		neighbors[i][3] = i < lnz ? i + linear_size : i - lnz;
	}
}

void update_spin() {
	int ran_site = num_spins * uniform(rng); // Random site
	int total_neighbor_spins = 0;
	for (int n = 0; n < 4; n++) total_neighbor_spins += spins[neighbors[ran_site][n]];
	double delta_energy = 2 * spins[ran_site] * (total_neighbor_spins + field);
	// Metropolis algorithm
	if (uniform(rng) < exp(-beta * delta_energy)) spins[ran_site] = - spins[ran_site];
}

void mc_step() {
	for (int i = 0; i < num_spins; i++) update_spin();
}

void warm_up(int time)
{
	for (int t = 0; t < time; t++) mc_step();
}
void setup(){
	int i,j;
	Selection square_selection;
	//square choose
	for(i=0;i<sqrt(chosen_num);i++)	for(j=0;j<sqrt(chosen_num);j++) {
		square_selection.chosen_indices.push_back(i*linear_size+j);
		// std::cout<<"square_selected position is "<<i*linear_size+j<<std::endl;
	}
	square_selection.M.resize(chosen_num);

	square_selection.aC.resize(pair_num);
	selection_list.push_back(square_selection);

	Selection uniform_selection;
	//uniformly chosen
	for (i = 0; i<linear_size; i += 4) for (j = 0; j<linear_size; j += 4) {
		uniform_selection.chosen_indices.push_back(i*linear_size+j);
	}
	uniform_selection.M.resize(chosen_num);
	uniform_selection.aC.resize(pair_num);
	selection_list.push_back(uniform_selection);

	for (int seed = 0; seed<16; seed ++) {
		Selection random_selection;
		std::mt19937 rng(seed);
		bool selected[num_spins];
		for (auto &s: selected) s = false;
		while (random_selection.chosen_indices.size()<chosen_num) {
			int random_site = uniform(rng)*num_spins;
			if (not selected[random_site]) {
				selected[random_site] = true;
				random_selection.chosen_indices.push_back(random_site);
			}
		}
		random_selection.M.resize(chosen_num);
		random_selection.aC.resize(pair_num);
		selection_list.push_back(random_selection);
	}
	Statics_list.Mtotal.resize(1);
	Statics_list.ME.resize(1);
	Statics_list.E.resize(1);
	Statics_list.E2.resize(1);

}

//choose random spins and calculate cij
void calc_cij_mag(){
	int i,j,k=0,m=0;
	double energy;
	double total_spin = 0;
	double total_pair_product = 0;
	warm_up(warm_up_time);
	std::cout<<"end warm up"<<std::endl;
	std::cout<<"mc step"<<std::endl;

	for (int m = 0; m< calculation_time; m++) {
		mc_step();
		// if (m%1000==0){
		// 	std::cout<<"calculating correlation for m="<<m<<std::endl;
		// }
		Statics_list.config.insert(std::end(Statics_list.config), std::begin(spins), std::end(spins));
		total_spin = 0;
		total_pair_product = 0;
		for(i=0;i<num_spins;i++){
			total_spin += spins[i];
			total_pair_product += spins[i] * (spins[neighbors[i][0]] + spins[neighbors[i][2]]); // Only consider left and up
		}
		energy = -(field * total_spin + total_pair_product);
		Statics_list.Mtotal[0] += total_spin;
		Statics_list.E[0] += energy;
		Statics_list.E2[0] += energy*energy;
		Statics_list.ME[0] += total_spin*energy;

		for (auto &selection: selection_list) {
			// selection.config.insert(std::end(selection.config), std::begin(spins), std::end(spins));
			int aC_index = 0;
			int aC_index2 = 0;

			for (i = 0; i<chosen_num; i++) {
				selection.M[i] += spins[selection.chosen_indices[i]];
				for (j = 0; j<i; j++) {
					selection.aC[aC_index++] += spins[selection.chosen_indices[i]]*spins[selection.chosen_indices[j]];
				}

			}
		// if (m%1000==0){
		// 	std::cout<<"calculating end, m="<< m <<std::endl;
		// }
		}
	}
}
void normalize(){
	for (auto &selection: selection_list) {
		for (auto & m: selection.M) m /= calculation_time;
		for (auto & aC: selection.aC) aC /= calculation_time;
	}
	Statics_list.Mtotal[0] /=calculation_time;
	Statics_list.E[0] /=calculation_time;
	Statics_list.E2[0]/=calculation_time;
	Statics_list.ME[0]/=calculation_time;
}
void clear_struct(){
	for (int i = 0; i<selection_list.size(); i++) {
		for (auto &selection: selection_list) {
			selection.aC.clear();
			selection.M.clear();
			selection.M.resize(chosen_num);
			selection.aC.resize(pair_num);
		}
	}
	Statics_list.config.clear();
	Statics_list.E.clear();
	Statics_list.E2.clear();
	Statics_list.Mtotal.clear();
	Statics_list.ME.clear();

	Statics_list.Mtotal.resize(1);
	Statics_list.ME.resize(1);
	Statics_list.E.resize(1);
	Statics_list.E2.resize(1);


}

int main()
{
	int i,j;
	field = 0.1;
	double temp=0.05;
	beta = 0.1;

	initialize();
	setup();

	for(beta=0.1;beta<=2;beta+=temp){
		std::cout<<"beta="<< beta <<std::endl;
		clear_struct();
		calc_cij_mag();
		normalize();
		std::string Beta=format(beta,4);
		std::string h=format(field,4);
		for (int i = 0; i<selection_list.size(); i++) {

			auto &selection = selection_list[i];
			std::string fn = "selection_mc/l" +std::to_string(linear_size)+"Beta"+(Beta)+"h"+(h)+"_"+ std::to_string(i) + ".h5";
			H5::H5File file(fn,H5F_ACC_TRUNC);

			hsize_t l = pair_num;
			H5::DataSet ds = H5::DataSet(file.createDataSet("c",H5::PredType::NATIVE_DOUBLE,H5::DataSpace(1,&l)));
			ds.write(selection.aC.data(),H5::PredType::NATIVE_DOUBLE);
			
			l = chosen_num;
			ds = H5::DataSet(file.createDataSet("m",H5::PredType::NATIVE_DOUBLE,H5::DataSpace(1,&l)));
			ds.write(selection.M.data(),H5::PredType::NATIVE_DOUBLE);


			ds = H5::DataSet(file.createDataSet("chosen_indices",H5::PredType::NATIVE_INT,H5::DataSpace(1,&l)));
			ds.write(selection.chosen_indices.data(), H5::PredType::NATIVE_INT);

		}
		std::string fn = "selection_mc/config_l" +std::to_string(linear_size)+"Beta"+Beta+"h"+h+"_"+ std::to_string(i) + ".h5";
		H5::H5File file2(fn,H5F_ACC_TRUNC);

		hsize_t l = Statics_list.config.size();
		H5::DataSet ds = H5::DataSet(file2.createDataSet("config",H5::PredType::NATIVE_INT,H5::DataSpace(1,&l)));
		ds.write(Statics_list.config.data(),H5::PredType::NATIVE_INT);

		l= Statics_list.Mtotal.size();
		ds = H5::DataSet(file2.createDataSet("Mtotal",H5::PredType::NATIVE_DOUBLE,H5::DataSpace(1,&l)));
		ds.write(Statics_list.Mtotal.data(),H5::PredType::NATIVE_DOUBLE);

		l= Statics_list.ME.size();
		ds = H5::DataSet(file2.createDataSet("ME",H5::PredType::NATIVE_DOUBLE,H5::DataSpace(1,&l)));
		ds.write(Statics_list.ME.data(),H5::PredType::NATIVE_DOUBLE);
		
		l= Statics_list.E.size();
		ds = H5::DataSet(file2.createDataSet("E",H5::PredType::NATIVE_DOUBLE,H5::DataSpace(1,&l)));
		ds.write(Statics_list.E.data(),H5::PredType::NATIVE_DOUBLE);
		
		l= Statics_list.E2.size();
		ds = H5::DataSet(file2.createDataSet("E2",H5::PredType::NATIVE_DOUBLE,H5::DataSpace(1,&l)));
		ds.write(Statics_list.E2.data(),H5::PredType::NATIVE_DOUBLE);
}
	return 0;
}
