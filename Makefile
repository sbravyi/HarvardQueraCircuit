all : ex ex_parallel

CC = g++ -std=c++11
CFLAGS = -O3

ex : sim.cpp $(DEPS)
	${CC} ${CFLAGS} sim.cpp -o $@

ex_parallel : sim.cpp $(DEPS)
	${CC} ${CFLAGS} sim_parallel.cpp -o $@ -lpthread

clean :
	rm ex ex_parallel
