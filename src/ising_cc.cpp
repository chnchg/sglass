#include <random>
#include <cmath>
#include <iostream>

// Calculation parameters
const int warm_up_time = 10000;
const int calculation_time = 100000;

// System parameters
const int linear_size = 10; // Linear size
double beta; // Inverse temperature
double field; // Magnetic field

// State
const int num_spins = linear_size*linear_size;
int spins[num_spins];
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
	double &energy_magnetization
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
}

int main()
{
	field = 0.;
	initialize();
	for (double temp = 5; temp > 1; temp *= .9) {
		beta = 1./temp;
		warm_up(warm_up_time);
		double energy, energy_sqr, magnetization, magnetization_sqr, energy_magnetization;
		calc_mean_statistics(
			calculation_time, energy, energy_sqr, magnetization,
			magnetization_sqr, energy_magnetization
		);
		std::cout << temp << ' ' << energy << ' ' << energy_sqr << ' ';
		std::cout << magnetization << ' ' << magnetization_sqr << ' ';
		std::cout << energy_magnetization << std::endl;
	}
}
