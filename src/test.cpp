#include <vector>
#include <random>
#include <iostream>
#include <fstream>
#include <memory>
#include <cassert>
#include <H5Cpp.h>
#include <cstdlib>
#include <time.h>   /* 時間相關函數 */


using namespace std;
#define L 10
#define SIZE L *L

void initialize(double spins[L][L]);
double calcenergy(double spins[L][L]);
double calcmag(double spins[L][L]);
void montecarlo(double spins[L][L],double T);

vector<double> E;
vector<double> M;
vector<double> C;
vector<double> X;

std::uniform_real_distribution<double> prob;
std::mt19937 rng;
mt19937 gen(958431198);                             // Mersenne Twister RNG
uniform_int_distribution<int> brandom(10);        // Get any random integer
uniform_int_distribution<int> ran_pos(0, SIZE - 1); // Get any random integer
uniform_real_distribution<double> ran_u(0.0, 1.0);  // Our uniform variable generator

int main()
{
    // random_device rd;
    srand( time(NULL) );


    double spins[L][L];
    int i,j;
    int neighs[L][L][4];                     // To store nearest neighbours
    double tstar;                            // Control parameter
    double deltat, deltat_crit;              // Change in control parameter by iteration
    double tmax, tmin, tcrit_up, tcrit_down; // Max and min temperature, and interval where we apply Wolff
    double energy;                           // Value of the energy of the system
    double mag;  
    double e_tem;
    double m_tem;
    double m2_tem;
    double e2_tem;
    int N; // number of iteration

    tmax = 5.0;
    tmin = 0.1;
    tcrit_up = 2.4;
    tcrit_down = 2.2;

    deltat = 0.1;
    deltat_crit = 0.01;

    initialize(spins);
    // get_neighbors(neighs);              // Get neighbour table
    energy = calcenergy(spins); // Compute initial energy
    
    N = 1e4;
    // cout<< energy << endl;


    for (tstar = tmax; tstar > tcrit_up; tstar -= deltat)
    {
        // warmup
        e_tem = 0;
        e2_tem = 0;
        m_tem = 0;
        m2_tem = 0;
        for (i = 0; i < N; i++)
        {
            montecarlo(spins,tstar);
        }
        // frame
        for (i = 0; i < N; i++)
        {
            montecarlo(spins,tstar);
            energy = calcenergy(spins);
            mag = calcmag(spins);

            e_tem += energy;
            e2_tem += energy*energy;
            m_tem += mag;
            m2_tem += mag*mag;
        }

        E.push_back(e_tem/N);
        M.push_back(m_tem/N);
        C.push_back(((e2_tem-e_tem*e_tem)/N)/tstar/tstar);
        X.push_back((m2_tem-m_tem*m_tem/N)/N/tstar);

        cout << tstar << endl;
    }

    // save_data();
    fstream newFile;
    newFile.open("M.txt", ios::out);

    // Write to the file
    for(auto &v : M){
        newFile << v << " ";
    }
    newFile.close();

    newFile.open("E.txt", ios::out);

    // Write to the file
    for(auto &v : E){
        newFile << v << " ";
    }
    newFile.close();
    for(auto &v : C){
        newFile << v << " ";
    }
    newFile.close();

    newFile.open("C.txt", ios::out);

    for(auto &v : X){
        newFile << v << " ";
    }
    newFile.close();

    newFile.open("X.txt", ios::out);

    return 0;
}

void initialize(double spins[L][L])
{
    int i, j;

    // Init spins with a random distribution
    for (i = 0; i < L * L; i++)
    {
        for (j=0;j<L;j++){
            spins[i][j] = 2*(rand()%2)-1; // Generate numbers

        }
    }

    return;
}

double calcenergy(double spins[L][L])
{
    int i,j; // Counters
    int neighboresum;

    int en = 0; // Sum

    // For every spin,
    for (i = 0; i < L; i++){
        for(j=0;j<L;j++){
            neighboresum = spins[(i+1)%L][j]+spins[(i-1)%L][j]+spins[i][(j+1)%L]+spins[i][(j-1)%L];
        // And compute the energy change
          en = en + spins[i][j] +spins[i][j]*neighboresum;

        }
        // Get sum of the neighbours
 
    }

    return en ; // Return the energy
}

double calcmag(double spins[L][L])
{
    int i,j;

    double sum = 0.0;
    // Sum all the values of the spins
    for (i = 0; i < L; i++){
        for (j = 0; j < L; j++){
            sum += spins[i][j];
        }
        
    }
    // And then return them
    return abs(sum)/L/L;
}

void montecarlo(double spins[L][L],double T){
	// random spin
    int i,j;
    int ran_posi;
    int ran_posj;
    int neighborsum;
    double cost=0 ;

	for(i=0;i<L;i++){
        for(j=0;j<L;j++){
            ran_posi=rand()%L;
            ran_posj=rand()%L;
            // cout<<"i posotion = "<<ran_posi;
            // cout<<" j posotion = "<<ran_posj;
            // cout<<endl;
            neighborsum = spins[(ran_posi+1)%L][ran_posj]+spins[(ran_posi-1)%L][ran_posj]+spins[ran_posi][(ran_posj+1)%L]+spins[ran_posi][(ran_posj-1)%L];
            cost= 2*spins[ran_posi][ran_posj]*neighborsum;
	    // dE = spin[i]*heff*2
            if(cost<0){
                spins[ran_posi][ran_posj] = -spins[ran_posi][ran_posj];
            }
	        else if ((prob(rng))<exp(-1/T*cost)) {
        
                spins[ran_posi][ran_posj] = -spins[ran_posi][ran_posj];
        
            }

        }


    }

    return;
}
// void mc_step()
// {
//     int i,j;

// 	for(i=0;i<L;i++){
//         for(j=0;j<L;j++){
//             metropolis();
//         }


//     }
// }

void flip_spin(double spins[SIZE], int neighs[SIZE][4], double h[5], double& energy, mt19937& gen, uniform_real_distribution<double>& ran_u, uniform_int_distribution<int>& ran_pos,double T){
    int index = ran_pos(gen); //Get a random position to flip
    int i;
    //Compute the sum of neighbours
    int sum_neigh = spins[neighs[index][0]] + spins[neighs[index][1]] + spins[neighs[index][2]] + spins[neighs[index][3]];
    //Use this to get the energy change (depending on the value of my spin)
    int change = spins[index] ? 2.0 * (sum_neigh) - 4.0 : 4.0 - 2.0 * (sum_neigh);

    for (i=-4; i <= 4; i += 2)
	{
		h[(i+4)/2] =  min(1.0, exp(- 2.0 * i / T));
	}
    //Apply Metropolis skim
    if (ran_u(gen) < h[(change+4)/2])
    {
        spins[index] = !spins[index];
        energy += (2.0*change)/(1.0*SIZE);
        //cout << change << "  " << (2.0*change)/(1.0*SIZE) << "  " << energy << endl;
    }

    return;
}

// void flip_spin(bool spins[SIZE], int neighs[SIZE][4], double h[5], double &energy, mt19937 &gen, uniform_real_distribution<double> &ran_u, uniform_int_distribution<int> &ran_pos, double &T)
// {
//     int index = ran_pos(gen); // Get a random position to flip
//     double heff = h[index];
//     // Compute the sum of neighbours
//     int sum_neigh = spins[neighs[index][0]] + spins[neighs[index][1]] + spins[neighs[index][2]] + spins[neighs[index][3]];
//     // Use this to get the energy change (depending on the value of my spin)
//     int change = spins[index] ? 2.0 * (sum_neigh)-4.0 : 4.0 - 2.0 * (sum_neigh);

//     // Apply Metropolis skim
//     if (ran_u(gen) < exp(-(1 / T) * spins[index] * heff * 2))
//     {
//         spins[index] = !spins[index];
//         energy += (2.0 * change) / (1.0 * SIZE);
//         // cout << change << "  " << (2.0*change)/(1.0*SIZE) << "  " << energy << endl;
//     }
// }

void mcmove(bool spins[SIZE], int neighs[SIZE][4], double tstar, int N, double &energy, mt19937 &gen, uniform_real_distribution<double> &ran_u, uniform_int_distribution<int> &ran_pos, double &T)
{
    int i;
    int index = ran_pos(gen);
    int sum_neigh = spins[neighs[index][0]] + spins[neighs[index][1]] + spins[neighs[index][2]] + spins[neighs[index][3]];
    int difference = 2 * spins[index] * sum_neigh;

    if (difference < 0)
    {
        spins[index] = -spins[index];
    }
    else if (ran_u(gen) < exp(-difference / T))
    {
        spins[index] = -spins[index];
    }
}