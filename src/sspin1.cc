/**
   @file sspin1.cc
   @brief Spin-glass model in symmetric-spin representation.
*/
#include "arg.hh"
#include <H5Cpp.h>
#include <exception>
#include <vector>
#include <random>
#include <iostream>
#include <memory>
#include <cassert>

class Logger
{
	bool on;
public:
	Logger() : on(false) {}
	template<class T>
    Logger& operator<<(const T& out)
    {
		if (on) std::cerr<<out;
		return *this;
    }
} info;

///\name System variables
/// Parameters and state of the system
///@{
std::vector<double> h; ///<Field
std::vector<double> j; ///<Pair coupling, condensed format
// runtime thermo
double beta; ///<Inverse temperature
// state
std::vector<int> spin; ///<Spin states 
// RNG
std::mt19937 rng; ///<Random number generator
///@}
///\name Auxiliary variables
/// Variables that can be recalculated
///@{
int sz; ///<System size
int link_num; ///<Number of links
/// Connection between nodes
struct Link {
	int n0; ///<Node0
	int n1; ///<Node1
};
std::vector<Link> links; ///<All connections of the system
/// Connection of a spin
struct SpinBond {
	int l; ///<Index of the connection
	int n; ///<Node it's connecting to
};
std::vector<std::vector<SpinBond>> slnk; ///<Connections of all spins
std::uniform_int_distribution<int> rnd_site; ///<Random picking of a site
std::uniform_real_distribution<double> prob; ///<Uniform random value
///@}
///\name Ensemble averages
/// Results of calculation
///@{
double cnt; ///<Count for measurements
double aE;  ///<Energy averager/accumulator
double aE2; ///<Accumulator for energy square
double aM; ///<Magnetization averager
double aM2; ///<Accumulator for magnetization square
double aME; ///<Accumulator for energy-magnetization product
std::vector<double> aS; ///<Averager for spins
std::vector<double> aC; ///<Averager for correlations
///@}
/// Load a 1D array of double
/// @param ds Dataset to load from
/// @param v Storage of data
void load_1d_array(H5::DataSet const & ds, std::vector<double> & v)
{
	H5::DataSpace s = ds.getSpace();
	int r = s.getSimpleExtentNdims();
	if (r!=1) throw std::exception();
	hsize_t n;
	s.getSimpleExtentDims(&n);
	v.resize(n);
	ds.read(v.data(),H5::PredType::NATIVE_DOUBLE);
}
/// Load a 1D array of int
/// @param ds Dataset to load from
/// @param v Storage of data
void load_1d_array(H5::DataSet const & ds, std::vector<int> & v)
{
	H5::DataSpace s = ds.getSpace();
	int r = s.getSimpleExtentNdims();
	if (r!=1) throw std::exception();
	hsize_t n;
	s.getSimpleExtentDims(&n);
	v.resize(n);
	ds.read(v.data(),H5::PredType::NATIVE_INT);
}
/// Load parameters from a file
/// @param fn File name
void load_params(std::string fn)
{
	// read parameters
	H5::H5File f(fn,H5F_ACC_RDONLY);
	load_1d_array(f.openDataSet("h"),h);
	load_1d_array(f.openDataSet("j"),j);
	beta = 1.0;
}
/// Setup system from parameters
void setup_system()
{
	sz = h.size();   // number of nodes
	link_num = j.size(); // number of pairs
	spin.resize(sz);
	slnk.resize(sz); // links of each node
	links.clear();    // links in system
	for (int & s:spin) s = -1;
	for (auto & i:slnk) i.clear();
	int l = 0;
	// setup fully connected network: i<j, little endian
	for (int j = 0; j<sz; j++) for (int i = 0; i<j; i++) {
		slnk[i].push_back(SpinBond{l,j});
		slnk[j].push_back(SpinBond{l,i});
		links.push_back(Link{i,j});
		l += 1;
	}
	if (l!=link_num) throw std::exception();
	rnd_site.param(std::uniform_int_distribution<int>::param_type(0,sz-1));
	// measurements
	aE = 0;
	aE2 = 0;
	aM = 0;
	aM2 = 0;
	aME = 0;
	aS.resize(sz);
	std::fill(aS.begin(),aS.end(),0);
	aC.resize(link_num);
	std::fill(aC.begin(),aC.end(),0);
}
/// Perform one Metropolis updating
void metropolis()
{
	// random spin
	int i = rnd_site(rng);
	double heff = h[i];
	for (auto b:slnk[i]) heff += j[b.l]*spin[b.n];
	// dE = spin[i]*heff*2
	if (prob(rng)<exp(-beta*spin[i]*heff*2)) spin[i] = -spin[i];
}
/// Perform one MC step
void mc_step()
{
	for (int i = 0;i<sz;i++) metropolis();
}
/// Perform measurement
void measure()
{
	double E = 0; // total energy
	double M = 0; // total magnetization
	// field
	for (int i = 0;i<sz;i++) {
		E += h[i]*spin[i];
		M += spin[i];
		aS[i] += spin[i];
	}
	// pair coupling
	for (int l = 0;l<link_num;l++) {
		int ss = spin[links[l].n0]*spin[links[l].n1];
		E += j[l]*ss;
		aC[l] += ss;
	}
	E = -E;
	aE += E;
	aE2 += E*E;
	aM += M;
	aM2 += M*M;
	aME += M*E;
	cnt ++;
}
/// Save a scalar variable
/// @param loc Where to save
/// @param name Name of the scalar
/// @param val The scalar data
void save_scalar(H5::H5Location & loc,std::string const & name,double val)
{
    H5::DataSet ds(loc.createDataSet(name,H5::PredType::NATIVE_DOUBLE,H5::DataSpace()));
    ds.write(&val,H5::PredType::NATIVE_DOUBLE);
}
/// Save a 1D array of double
/// @param loc Where to save
/// @param name Name of the data array
/// @param val Storage of the array of data
void save_1d_array(H5::H5Location & loc,std::string const & name,std::vector<double> const & val)
{
    hsize_t l = val.size();
    H5::DataSet ds(loc.createDataSet(name,H5::PredType::NATIVE_DOUBLE,H5::DataSpace(1,&l)));
    ds.write(val.data(),H5::PredType::NATIVE_DOUBLE);
}
/// Save system data
/// @param fn Name of a HDF5 file to save
void save_data(std::string fn = "result_file.h5")
{
	std::unique_ptr<H5::H5File> f;
	H5::Exception::dontPrint();
	try {
		f.reset(new H5::H5File(fn,H5F_ACC_RDWR));
		info << "opened " << fn << '\n';
	}
	catch (H5::FileIException const & e) {
		f.reset(new H5::H5File(fn,H5F_ACC_TRUNC));
		info << "created " << fn << '\n';
	}
    if (f->exists("result1")) f->unlink("result1");
    H5::Group g(f->createGroup("result1"));
    save_scalar(g,"cnt",cnt);
    save_scalar(g,"aE",aE);
    save_scalar(g,"aE2",aE2);
    save_scalar(g,"aM",aM);
    save_scalar(g,"aM2",aM2);
    save_scalar(g,"aME",aME);
    save_1d_array(g,"aS",aS);
    save_1d_array(g,"aC",aC);
}
/// Main function for performing simulation
int main(int argc, char ** argv)
{
	using std::string, std::cerr, std::endl;
	using namespace arg;

	string infn;
	string oufn;
	double count;
	double warmup = 1024;
	double pm_beta;
	unsigned pm_seed;

	Parser psr;
	psr.add_prm(infn).name("infn").desc("Input parameter file");
	psr.add_prm(oufn).name("oufn").desc("Output result file");
	psr.add_prm(count).name("count").desc("Measurement count");
	psr.required();
	psr.add_prm(warmup).name("warmup").desc("Warmup-step count");
	psr.add_opt("beta",pm_beta).desc("Inverse temperature");
	psr.add_opt("seed",pm_seed).desc("Seed for random number generator");
	try {
		psr.parse(argc,argv);
	}
	catch (std::logic_error &e) {
		cerr << e.what() << "\n\n";
		cerr << psr.get_help();
		cerr << endl;
		return -1;
	}

	load_params(infn);
	setup_system();
	// parameter override
	if (psr.opt_is_present("beta"))	beta = pm_beta;
	rng.seed(psr.opt_is_present("seed")?pm_seed:0);
	for (auto i = 0;i<warmup;i++) mc_step();
	while (cnt<count) {
		mc_step();
		measure();
	}
	save_data(oufn);
	return 0;
}
