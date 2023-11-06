all: sspin1
sspin1.o: sspin1.cc arg.hh
	c++ -I/usr/include/hdf5/serial -c sspin1.cc -g -O3 -Wall -Wextra -Werror -static
sspin1: sspin1.o
	c++ -I/usr/include/hdf5/serial sspin1.o -o sspin1 -g -O3 -Wall -Wextra -Werror -static \
	-L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 \
	-lpthread -lsz -lz -ldl -lm -laec -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial

all: ising_mc
ising_mc.o: ising_mc.cpp 
	c++ -I/usr/include/hdf5/serial -c ising_mc.cpp -g -O3  -static
ising_mc: ising_mc.o
	c++ -I/usr/include/hdf5/serial ising_mc.o -o ising_mc -g -O3 -static \
	-L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 \
	-lpthread -lsz -lz -ldl -lm -laec 

all: sspin2
sspin2.o: sspin2.cc arg.hh
	c++ -I/usr/include/hdf5/serial -c sspin2.cc -g -O3 -Wall -Wextra -Werror -static
sspin2: sspin2.o
	c++ -I/usr/include/hdf5/serial sspin2.o -o sspin2 -g -O3 -Wall -Wextra -Werror -static \
	-L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 \
	-lpthread -lsz -lz -ldl -lm -laec -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial

all: ising_mc_cb
ising_mc_cb.o: ising_mc_cb.cpp 
	c++ -I/usr/include/hdf5/serial -c ising_mc_cb.cpp -g -O3  -static
ising_mc_cb: ising_mc_cb.o
	c++ -I/usr/include/hdf5/serial ising_mc_cb.o -o ising_mc_cb -g -O3 -static \
	-L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 \
	-lpthread -lsz -lz -ldl -lm -laec 

all: ising_mc_cov
ising_mc_cov.o: ising_mc_cov.cpp 
	c++ -I/usr/include/hdf5/serial -c ising_mc_cov.cpp -g -O3  -static
ising_mc_cov: ising_mc_cov.o
	c++ -I/usr/include/hdf5/serial ising_mc_cov.o -o ising_mc_cov -g -O3 -static \
	-L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 \
	-lpthread -lsz -lz -ldl -lm -laec 

