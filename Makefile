#CCFLAGS=-O3 -march=native -msse2
#TODO: find out correct optimization flags

# with OpenMP acceleration
CCFLAGS=-O3 -fopenmp
LDFLAGS=-lgomp
# without OpenMP acceleration
#CCFLAGS=-O3

all: libhremd.so

hremd-boost.o: hremd-boost.c boost.h hremd-boost.h Makefile
	gcc -Wall $(CCFLAGS) -fPIC -c -o hremd-boost.o hremd-boost.c

boost.o: boost.c boost.h hremd-boost.h Makefile
	gcc -Wall $(CCFLAGS) -fPIC -c -o boost.o boost.c

libhremd.so: hremd-boost.o boost.o Makefile
	gcc --shared $(LDFLAGS) -o libhremd.so hremd-boost.o boost.o

clean:
	rm hremd-boost.o boost.o libhremd.so
