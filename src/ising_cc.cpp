#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
// Calculation parameters
const int warm_up_time = 10000;
const int calculation_time = 10000;

// System parameters
const int linear_size = 32; // Linear size
double beta; // Inverse temperature
double field; // Magnetic field

// State
const int num_spins = linear_size*linear_size;
int spins[num_spins];
double correlation[num_spins*num_spins];
bool selected[num_spins];

//store parameters
std::vector<double> E;
std::vector<double> M;
std::vector<double> C;
std::vector<double> T;

std::mt19937 rng(1234);

// Auxiliary
std::uniform_real_distribution<> uniform;
int neighbors[num_spins][4];

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

void calc_energy_magnetization(double &energy, double &magnetization) {
	double total_spin = 0;
	double total_pair_product = 0;
	for (int i = 0; i < num_spins; i++) {
		total_spin += spins[i];
		total_pair_product += spins[i] * (spins[neighbors[i][0]] + spins[neighbors[i][2]]); // Only consider left and up
	}
	energy = -(field * total_spin + total_pair_product);
	magnetization = total_spin;
}

void warm_up(int time)
{
	for (int t = 0; t < time; t++) mc_step();
}


void calc_mean_statistics(
	int time,
	double &energy,
	double &energy_sqr,
	double &magnetization,
	double &magnetization_sqr,
	double &energy_magnetization,
	double &Cv
) {
	double acc_energy = 0;
	double acc_energy_sqr = 0;
	double acc_magnetization = 0;
	double acc_magnetization_sqr = 0;
	double acc_energy_magnetization = 0;
	for (int i = 0; i< time; i++) {
		mc_step();
		double energy, magnetization;
		calc_energy_magnetization(energy, magnetization);
		acc_energy += energy;
		acc_energy_sqr += energy * energy;
		acc_magnetization += magnetization;
		acc_magnetization_sqr += magnetization * magnetization;
		acc_energy_magnetization += energy * magnetization;
	}
	energy = acc_energy / time;
	energy_sqr = acc_energy_sqr / time;
	magnetization = acc_magnetization / time;
	magnetization_sqr = acc_magnetization_sqr / time;
	energy_magnetization = acc_energy_magnetization / time;
	Cv=energy_sqr-energy*energy;
}
//choose random spins and calculate cij
void calc_cij(){
	int i,j,k=0,l=0,m=0;

	// int random_site=uniform(rng)*1023;
	// bool selected[1024];
	int num_selected =0;
	//choose random spins
	std::cout<<"Calculating correlation"<<std::endl;
	while(num_selected<8){
		int random_site=uniform(rng)*(num_spins-1);
		if(selected[random_site]== false){
			selected[random_site]= true;
			num_selected=num_selected+1;
			// std::cout<<"selected position is "<<random_site<<std::endl;
			std::cout<<"Number selected "<<num_selected<<std::endl;
		}
	}
	std::cout<<"warm up"<<std::endl;
	warm_up(warm_up_time);
	std::cout<<"end warm up"<<std::endl;
	std::cout<<"mc step"<<std::endl;
	for (int m = 0; m< 1000; m++) {
		mc_step();
		if (m%100==0){
			std::cout<<"calculating correlation for m="<<m<<std::endl;
		}
		
		for (i=0;i<num_spins;i++){
			// std::cout<<"spin "<<i<<std::endl;
			if(selected[i]== true){
				// std::cout<<"spin i is true "<<i<<std::endl;
				k=k+1;
				for(j=0;j<num_spins;j++){
					// std::cout<<"spin "<<j<<std::endl;
					if(selected[j]== true){
						// std::cout<<"spin j is true "<<j<<std::endl;
						l=l+1;
						correlation[i*num_spins+j]=correlation[i*num_spins+j]+spins[i]*spins[j];
					}
				}
			}
		}
		if (m%100==0){
			std::cout<<"calculating end, m="<< m <<std::endl;
		}
	}
}

int main()
{
	int i,j;
	field = 0.;
	initialize();
	for (double temp = 2.5; temp > 1.5; temp -= 0.01) {
		beta = 1./temp;
		warm_up(warm_up_time);
		double energy, energy_sqr, magnetization, magnetization_sqr, energy_magnetization,Cv;
		calc_mean_statistics(
			calculation_time, energy, energy_sqr, magnetization,
			magnetization_sqr, energy_magnetization,Cv
		);
		E.push_back(energy);
        M.push_back(magnetization);
        C.push_back(energy_sqr-energy*energy);
        T.push_back(temp);
		std::cout << temp << ' ' << energy << ' ' << Cv << ' '<< std::endl;
		// std::cout << magnetization << ' ' << magnetization_sqr << ' ';
		// std::cout << energy_magnetization <<' '<< Cv << std::endl;
	}
	// save data
	std::fstream newFile;
	newFile.open("M.txt", std::ios::out);
    for(auto &v : M){
        newFile << v << " ";
    }
    newFile.close();

    newFile.open("E.txt", std::ios::out);
    for(auto &v : E){
        newFile << v << " ";
    }
    newFile.close();

	newFile.open("C.txt", std::ios::out);
    for(auto &v : C){
        newFile << v << " ";
    }
    newFile.close();

    newFile.open("T.txt", std::ios::out);
    for(auto &v : T){
        newFile << v << " ";
    }
    newFile.close();

	std::vector<double>::iterator max_c;
    
    max_c = std::max_element(C.begin(), C.end());
    std::cout << "Max capacity: " << *max_c << "\n";
    


	// beta=1.0/2.35;
	// calc_cij();

	// for (int i = 0; i < num_spins*num_spins; i++) {
    //     std::cout << correlation[i] << " ";
	// 	if(i==num_spins){
	// 		std::cout << std::endl;
	// 	}
    // }

	return 0;
}
