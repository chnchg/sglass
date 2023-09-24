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
const int calculation_time = 10000;
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

std::vector<double> random_aC_tem(choosen_num*choosen_num);  //temporary correlation 
std::vector<double> random_aC; //correlation
std::vector<double> uniform_aC_tem(choosen_num*choosen_num);
std::vector<double> uniform_aC;
std::vector<double> square_aC_tem(choosen_num*choosen_num);
std::vector<double> square_aC;

std::vector<int> random_choosen_index;
std::vector<int> uniform_choosen_index;
std::vector<int> square_choosen_index;

std::vector<double> random_mag(choosen_num); //magnitization for random indiviual spin
std::vector<double> uniform_mag(choosen_num); 
std::vector<double> square_mag(choosen_num); 


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
			random_choosen_index.push_back(random_site);
			// std::cout<<"selected position is "<<random_site<<std::endl;
		}
	}
	//uniformly chhose
	num_selected=0;
	uniform_choosen_index.push_back(0);
	while(num_selected<choosen_num){
		num_selected=num_selected+1;
		uniform_choosen_index.push_back(num_selected*4-1);
		
	}
	//square choose
	square_choosen_index.push_back(0);
	for(i=0;i<sqrt(choosen_num);i++){
		for(j=0;j<sqrt(choosen_num);j++){
			square_choosen_index.push_back(i*linear_size+j);
			std::cout<<"square_selected position is "<<i*linear_size+j<<std::endl;
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
			random_mag[i]=random_mag[i]+spins[random_choosen_index[i]];
			uniform_mag[i]=uniform_mag[i]+spins[uniform_choosen_index[i]];
			square_mag[i]=square_mag[i]+spins[square_choosen_index[i]];
			for(j=0;j<choosen_num;j++){
				random_aC_tem[i*choosen_num+j]=random_aC_tem[i*choosen_num+j]+spins[random_choosen_index[i]]*spins[random_choosen_index[j]];
				uniform_aC_tem[i*choosen_num+j]=uniform_aC_tem[i*choosen_num+j]+spins[uniform_choosen_index[i]]*spins[uniform_choosen_index[j]];
				square_aC_tem[i*choosen_num+j]=square_aC_tem[i*choosen_num+j]+spins[square_choosen_index[i]]*spins[square_choosen_index[j]];
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
void normalize(){
	int i,j;
	for (j=0;j<choosen_num;j++){
		random_mag[j]=random_mag[j]/calculation_time;
		uniform_mag[j]=uniform_mag[j]/calculation_time;
		square_mag[j]=square_mag[j]/calculation_time;
		for(i=0;i<j;i++){
			random_aC.push_back(random_aC_tem[j*choosen_num+i]/calculation_time);
			uniform_aC.push_back(uniform_aC_tem[j*choosen_num+i]/calculation_time);
			square_aC.push_back(square_aC_tem[j*choosen_num+i]/calculation_time);
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
	calc_cij_mag();
	normalize();

	// std::cout<<std::endl<<"vector aC"<<std::endl;
	// for (auto &v : random_aC) {
    //     std::cout << v << " ";
    // }

	// std::cout<<std::endl<<"vector aC_tem"<<std::endl;
	// for (int i = 0; i < choosen_num*choosen_num; i++) {
    //     std::cout << random_aC_tem[i] << " ";
	// 	if(i%(choosen_num)==choosen_num-1){
	// 		std::cout << std::endl;
	// 	}
    // }

	// for(int i = 0; i < choosen_num; i++){
	// 	std::cout << random_mag[i] << " ";
	// }

	// save data
	std::string fn = "random_mag_corre_file.h5" ;
	
	H5::H5File file1(fn,H5F_ACC_TRUNC);

	hsize_t l = random_aC.size();
    H5::DataSet ds= H5::DataSet(file1.createDataSet("c",H5::PredType::NATIVE_DOUBLE,H5::DataSpace(1,&l)));
	ds.write(random_aC.data(),H5::PredType::NATIVE_DOUBLE);


	l = random_mag.size();
    ds= H5::DataSet(file1.createDataSet("m",H5::PredType::NATIVE_DOUBLE,H5::DataSpace(1,&l)));
	ds.write(random_mag.data(),H5::PredType::NATIVE_DOUBLE);

	fn = "uniform_mag_corre_file.h5" ;
	
	H5::H5File file2(fn,H5F_ACC_TRUNC);

	l = uniform_aC.size();
    ds= H5::DataSet(file2.createDataSet("c",H5::PredType::NATIVE_DOUBLE,H5::DataSpace(1,&l)));
	ds.write(uniform_aC.data(),H5::PredType::NATIVE_DOUBLE);


	l = uniform_mag.size();
    ds= H5::DataSet(file2.createDataSet("m",H5::PredType::NATIVE_DOUBLE,H5::DataSpace(1,&l)));
	ds.write(uniform_mag.data(),H5::PredType::NATIVE_DOUBLE);

	fn = "square_mag_corre_file.h5" ;
	
	H5::H5File file3(fn,H5F_ACC_TRUNC);

	l = square_aC.size();
    ds= H5::DataSet(file3.createDataSet("c",H5::PredType::NATIVE_DOUBLE,H5::DataSpace(1,&l)));
	ds.write(square_aC.data(),H5::PredType::NATIVE_DOUBLE);


	l = square_mag.size();
    ds= H5::DataSet(file3.createDataSet("m",H5::PredType::NATIVE_DOUBLE,H5::DataSpace(1,&l)));
	ds.write(square_mag.data(),H5::PredType::NATIVE_DOUBLE);



	return 0;
}
