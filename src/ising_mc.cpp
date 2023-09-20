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

// Calculation parameters
const int warm_up_time = 10000;
const int calculation_time = 100000;
const int choosen_num=64;
// System parameters
const int linear_size = 32; // Linear size
double beta; // Inverse temperature
double field; // Magnetic field

// State
const int num_spins = linear_size*linear_size;
int spins[num_spins];
// double correlation[num_spins*num_spins];

bool selected[num_spins];


//store parameters
std::vector<double> E;
std::vector<double> E_squr;
std::vector<double> M;
std::vector<double> C;
std::vector<double> T;

std::vector<double> aC_tem(choosen_num*choosen_num);  //temporary correlation 
std::vector<double> aC; //correlation
std::vector<int> choosen_index;
std::vector<double> mag(choosen_num); //magnitization for indiviual spin



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

void warm_up(int time)
{
	for (int t = 0; t < time; t++) mc_step();
}

//choose random spins and calculate cij
void calc_cij_mag(){
	int i,j,k=0,l=0,m=0;

	// int random_site=uniform(rng)*1023;
	// bool selected[1024];
	int num_selected =0;
	//choose random spins
	while(num_selected<choosen_num){
		int random_site=uniform(rng)*(num_spins);
		if(selected[random_site]== false){
			selected[random_site]= true;
			num_selected=num_selected+1;
			choosen_index.push_back(random_site);
			// std::cout<<"selected position is "<<random_site<<std::endl;
		}
	}
	warm_up(warm_up_time);
	std::cout<<"end warm up"<<std::endl;
	std::cout<<"mc step"<<std::endl;
	for (int m = 0; m< calculation_time; m++) {
		mc_step();
		if (m%1000==0){
			std::cout<<"calculating correlation for m="<<m<<std::endl;
		}
		for(i=0;i<choosen_num;i++){
			mag[i]=mag[i]+spins[choosen_index[i]];
			for(j=0;j<choosen_num;j++){
				aC_tem[i*choosen_num+j]=aC_tem[i*choosen_num+j]+spins[choosen_index[i]]*spins[choosen_index[j]];
			}
		}
		// k=0;
		// for (i=0;i<num_spins;i++){
		// 	// std::cout<<"spin "<<i<<std::endl;
		// 	if(selected[i]== true){
		// 		// std::cout<<"spin i is true "<<i<<std::endl;
		// 		mag[k]=mag[k]+spins[i]; //choosen mag
		// 		// std::cout<<"i = "<< i <<std::endl;
		// 		l=0;
		// 		for(j=0;j<num_spins;j++){
		// 			// std::cout<<"spin "<<j<<std::endl;
		// 			if(selected[j]== true){
		// 				// std::cout<<"spin j is true "<<j<<std::endl;
						
		// 				// correlation[k*num_spins+l]=correlation[k*num_spins+l]+spins[i]*spins[j];
		// 				aC_tem[k*choosen_num+l]=aC_tem[k*choosen_num+l]+spins[i]*spins[j];
		// 				l=l+1;
		// 				// std::cout<<"l = "<< l <<std::endl;
		// 			}
		// 		}
		// 		k=k+1;
		// 	}
		// }
		if (m%1000==0){
			std::cout<<"calculating end, m="<< m <<std::endl;
		}
	}
}

int main()
{
	int i,j;
	field = 0.;
	double temp=2.34;
	beta = 1./temp;

	initialize();


	warm_up(warm_up_time);
	calc_cij_mag();

	for (j=0;j<choosen_num;j++){
		for(i=0;i<j;i++){
			aC.push_back(aC_tem[j*choosen_num+i]/calculation_time);
		}
	}
	for (i=0;i<choosen_num;i++){
		mag[i]=mag[i]/calculation_time;
	}
	std::cout<<std::endl<<"vector aC"<<std::endl;
	for (auto &v : aC) {
        std::cout << v << " ";
    }

	// std::cout<<std::endl<<"selected spins"<<std::endl;
	// for (auto &v : selected) {
    //     std::cout << v << " ";
    // }

	std::cout<<std::endl<<"vector aC_tem"<<std::endl;
	for (int i = 0; i < choosen_num*choosen_num; i++) {
        std::cout << aC_tem[i] << " ";
		if(i%(choosen_num)==choosen_num-1){
			std::cout << std::endl;
		}
    }

	for(int i = 0; i < choosen_num; i++){
		std::cout << mag[i] << " ";
	}

	// save data
	std::string fn = "mag_corre_file.h5" ;
	
	H5::H5File file(fn,H5F_ACC_TRUNC);
	// std::unique_ptr<H5::H5File> f;
	// H5::Group g(f->createGroup("result1"));
	// H5::Group g(file.createGroup("result1"));

	hsize_t l = aC.size();
    H5::DataSet ds= H5::DataSet(file.createDataSet("c",H5::PredType::NATIVE_DOUBLE,H5::DataSpace(1,&l)));
	ds.write(aC.data(),H5::PredType::NATIVE_DOUBLE);


	l = mag.size();
    ds= H5::DataSet(file.createDataSet("m",H5::PredType::NATIVE_DOUBLE,H5::DataSpace(1,&l)));
	ds.write(mag.data(),H5::PredType::NATIVE_DOUBLE);


	return 0;
}
