CC = g++
WLECEL = -Wall
CFLAGS = -c $(WLEVEL) -O3 -I../include -I/usr/local/include/g2o/ -I/usr/local/include/eigen3/ -I/usr/include/suitesparse/ -L/usr/local/lib/ -lg2o_core -lg2o_types_slam2d -lg2o_solver_csparse

all: test.o

test.o:
	$(CC) test.cpp -I../include -I/usr/local/include/g2o/ -I/usr/local/include/eigen3/ -I/usr/include/suitesparse/ -L/usr/local/lib/ -lg2o_core -lg2o_types_slam2d -lg2o_solver_csparse -o FaMuS
	
clean:
	rm -f *.o FaMuS

#g++ ../test.cpp -I/usr/local/include/eigen3/ -I/usr/include/suitesparse/ -L/usr/local/lib/ -lg2o_core -lg2o_types_slam2d -lg2o_solver_csparse
#$(CC) $(CFLAGS) test.cpp -o test
