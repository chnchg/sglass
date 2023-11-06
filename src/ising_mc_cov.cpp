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
#include "arg.hh"

// using namespace std;
// Calculation parameters
const int warm_up_time = 10000;
const long long calculation_time = 100;
const long long chosen_num = 64;
const int pair_num = chosen_num*(chosen_num-1)/2;
// System parameters
long long linear_size ; // Linear size
const int draw_list_size = 2000; // Linear size
double beta; // Inverse temperature
double field; // Magnetic field

// State
long long num_spins = linear_size*linear_size;
// int *spinspointer;
std::vector<int> spins;

struct Selection {
	std::vector<int> chosen_indices;
	std::vector<double> M;
	std::vector<double> aC;	
	std::vector<double> aC3;	
	std::vector<double> aC4;
	std::vector<double> pk;
    std::vector<double> sub_cov;
	
};
std::vector<Selection> selection_list;

struct Statics{
	// std::vector<int> config;
	double E;
	double E2;
	double ME;
	double Mtotal;
	// std::vector<double> aC2;
	std::vector<int> draw_list3;
	std::vector<int> draw_list4;
	std::vector<double> pktotal;
    std::vector<double> total_cov;

};
Statics Statics_list;

std::mt19937 rng(1234);

// Auxiliary
std::uniform_real_distribution<> uniform;

std::vector<std::vector<int>> neighbors;

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
	num_spins = linear_size*linear_size;
	spins.resize(num_spins);
	neighbors.resize(num_spins);
	Statics_list.pktotal.resize(num_spins);

	for (int i = 0; i < num_spins; i ++) spins[i] = uniform(rng) < 0.5 ? 1 : -1;

	int lz = linear_size - 1;
	int lnz = linear_size * lz;
	for (int i = 0; i <  num_spins; i ++) {
		// periodic boundary
		// Left
		i % linear_size ? neighbors[i].push_back(i - 1) : neighbors[i].push_back(i + lz);
		// Right
		(i + 1) % linear_size ? neighbors[i].push_back(i + 1) : neighbors[i].push_back(i - lz);
		// Up
		i < linear_size ? neighbors[i].push_back(i + lnz) : neighbors[i].push_back(i - linear_size);
		// Down
		i < lnz ?  neighbors[i].push_back(i + linear_size) : neighbors[i].push_back(i - lnz);
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
	for(i=0;i<draw_list_size;i++){
		bool selected3[chosen_num],selected4[chosen_num];
		int draw_list3_tem=0;
		int draw_list4_tem=0;

		for (auto &s: selected3) s = false;
		for (auto &s: selected4) s = false;
		while (draw_list3_tem<3) {
			int random_site = uniform(rng)*chosen_num;
			if (not selected3[random_site]) {
				selected3[random_site] = true;
				Statics_list.draw_list3.push_back(random_site);
				draw_list3_tem+=1;
			}
		}
		while (draw_list4_tem<4) {
			int random_site = uniform(rng)*chosen_num;
			if (not selected4[random_site]) {
				selected4[random_site] = true;
				Statics_list.draw_list4.push_back(random_site);
				draw_list4_tem+=1;
			}
		}
	}
	// std::cout<<"Statics_list.draw_list3 size is "<<Statics_list.draw_list3.size();

	//square choose
	for(i=0;i<sqrt(chosen_num);i++)	for(j=0;j<sqrt(chosen_num);j++) {
		square_selection.chosen_indices.push_back(i*linear_size+j);
		// std::cout<<"square_selected position is "<<i*linear_size+j<<std::endl;
	}


	// for (i=0;i<30;i++){
	// 	std::cout<<draw3[0][i]<<std::endl;
	// }
	square_selection.M.resize(chosen_num);
	square_selection.aC.resize(pair_num);
	square_selection.aC3.resize(draw_list_size);
	square_selection.aC4.resize(draw_list_size);
	square_selection.pk.resize(chosen_num);
	selection_list.push_back(square_selection);

	Selection uniform_selection;
	//uniformly chosen
	for (i = 0; i<linear_size; i += linear_size/8) for (j = 0; j<linear_size; j += linear_size/8) {
		uniform_selection.chosen_indices.push_back(i*linear_size+j);
	}



	uniform_selection.M.resize(chosen_num);
	uniform_selection.aC.resize(pair_num);
	uniform_selection.aC3.resize(draw_list_size);
	uniform_selection.aC4.resize(draw_list_size);
	uniform_selection.pk.resize(chosen_num);
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
		random_selection.aC3.resize(draw_list_size);
		random_selection.aC4.resize(draw_list_size);
		random_selection.pk.resize(chosen_num);
		selection_list.push_back(random_selection);

	}
	// Statics_list.Mtotal.resize(1);
	// Statics_list.ME.resize(1);
	// Statics_list.E.resize(1);
	// Statics_list.E2.resize(1);
	// Statics_list.aC2.resize(pair_num);
    Statics_list.total_cov.resize(num_spins*num_spins);
    for (auto &selection: selection_list)selection.sub_cov.resize(chosen_num*chosen_num);

	
}

//choose random spins and calculate cij
void calc_cij_mag(){
	int i,j,k=0,m=0;
    std::vector<double> total_cov_config;
    std::vector<double> sub_cov_config;
    std::vector<double> total_cov_aC;
    std::vector<double> sub_cov_aC;
    total_cov_aC.resize(num_spins*num_spins);
    sub_cov_aC.resize(chosen_num*chosen_num*18);

    long long total_config_len=num_spins*num_spins*calculation_time;
    long long sub_config_len=chosen_num*chosen_num*calculation_time*18;
    sub_cov_config.resize(sub_config_len);
    total_cov_config.resize(total_config_len);

    std::cout<<"mc step start"<<std::endl;
	for (m=0; m< calculation_time; m++) {
		mc_step();
        
		for(i=0;i<num_spins;i++){
            for (j=0;j<num_spins;j++){
                total_cov_aC[i*num_spins+j]+=spins[i]*spins[j];
                total_cov_config.push_back(spins[i]*spins[j]);
            }
		}


        int listnum=0;
		for (auto &selection: selection_list) {
			for (i = 0; i<chosen_num; i++) {
				selection.M[i] += spins[selection.chosen_indices[i]];
				for (j = 0; j<chosen_num; j++) {
                    sub_cov_aC[chosen_num*chosen_num*listnum+i*chosen_num+j]+=spins[selection.chosen_indices[i]]*spins[selection.chosen_indices[j]];
                    sub_cov_config.push_back(spins[selection.chosen_indices[i]]*spins[selection.chosen_indices[j]]);
				}
			}

        listnum+=1;
		}
	}
    std::cout<<"mc step end"<<std::endl;

    for(i=0;i<num_spins;i++){
        for (j=0;j<num_spins;j++){
            for (m=0; m< calculation_time; m++){
                Statics_list.total_cov[i*num_spins+j]+=total_cov_config[m*num_spins*num_spins+i*num_spins+j]-total_cov_aC[i*num_spins+j]/calculation_time;


            }                
        }
    }
    std::cout<<"total covariance end"<<std::endl;
    int listnum=0;
    for (auto &selection: selection_list){
        for(i=0;i<chosen_num;i++){
            for (j=0;j<chosen_num;j++){
                for (m=0; m< calculation_time; m++){
                    selection.sub_cov[i*chosen_num+j]+=sub_cov_config[m*listnum*chosen_num*chosen_num+i*chosen_num+j]-sub_cov_aC[i*chosen_num+j]/calculation_time;
                }                
            }
        }
        listnum+=1;
    }
    std::cout<<"sub covariance end"<<std::endl;
}

void normalize(){
	for (auto &selection: selection_list) {
		for (auto & m: selection.sub_cov) m /= calculation_time;
	}
	for(auto & cov: Statics_list.total_cov) cov/= calculation_time;
}

void save_data()
{
    std::cout<<"save sub data"<<std::endl;
	std::string Beta=format(beta,4);
	std::string h=format(field,4);
	for (int i = 0; i<selection_list.size(); i++) {
		auto &selection = selection_list[i];

		std::string fn = "selection_mc/cov/subsample_l" +std::to_string(linear_size)+"Beta"+(Beta)+"h"+(h)+"_"+ std::to_string(i)+ + ".h5";
		H5::H5File file(fn,H5F_ACC_TRUNC);
		
		hsize_t l = selection.sub_cov.size();
		H5::DataSet ds = H5::DataSet(file.createDataSet("covariance",H5::PredType::NATIVE_DOUBLE,H5::DataSpace(1,&l)));
		ds.write(selection.sub_cov.data(),H5::PredType::NATIVE_DOUBLE);

        l=selection.chosen_indices.size();
        ds = H5::DataSet(file.createDataSet("chosen_indices",H5::PredType::NATIVE_INT,H5::DataSpace(1,&l)));
		ds.write(selection.chosen_indices.data(), H5::PredType::NATIVE_INT);

	}
    std::cout<<"save total data"<<std::endl;
	std::string fn = "selection_mc/cov/statics_l" +std::to_string(linear_size)+"Beta"+(Beta)+"h"+(h)+ ".h5";
	
    H5::H5File file2(fn,H5F_ACC_TRUNC);
    hsize_t l= Statics_list.total_cov.size();
	H5::DataSet ds = H5::DataSet(file2.createDataSet("total_covariance",H5::PredType::NATIVE_DOUBLE,H5::DataSpace(1,&l)));
	ds.write(Statics_list.total_cov.data(),H5::PredType::NATIVE_DOUBLE);
}



int main(int argc, char ** argv)
{
	using namespace arg;
	field = -0.01;
	beta = 0.1;
	linear_size=32;
	Parser psr;
	psr.add_opt("linear_size",linear_size).desc("linear size of ising model");
	psr.add_opt("beta",beta).desc("Inverse temperature");
	psr.add_opt("field",field).desc("field");
	try {
		psr.parse(argc,argv);
	}
	catch (std::logic_error &e) {
		std::cerr << e.what() << "\n\n";
		std::cerr << psr.get_help();
		std::cerr << std::endl;
		return -1;
	}

	std::cout<<"linear size = "<<linear_size<<" beta="<< beta <<" field = "<<field<<std::endl;


	initialize();
	setup();
    std::cout<<"set up system" <<std::endl;
	// clear_struct();
	calc_cij_mag();
	normalize();
	save_data();

	return 0;
}
