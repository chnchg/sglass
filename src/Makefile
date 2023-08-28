all: sspin1
sspin1.o: sspin1.cc arg.hh
	c++ -I/usr/include/hdf5/serial -c sspin1.cc -g -O3 -Wall -Wextra -Werror -static
sspin1: sspin1.o
	c++ -I/usr/include/hdf5/serial sspin1.o -o sspin1 -g -O3 -Wall -Wextra -Werror -static \
	-L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 \
	-lpthread -lsz -lz -ldl -lm -laec -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial
