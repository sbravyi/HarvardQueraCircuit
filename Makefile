all : ex_k4 ex_k5 ex ex_parallel

CC = g++ -std=c++11
CFLAGS = -O3

ex : sim.cpp $(DEPS)
	${CC} ${CFLAGS} sim.cpp -o $@

ex_parallel : sim_parallel.cpp $(DEPS)
	${CC} ${CFLAGS} sim_parallel.cpp -o $@ -lpthread

ex_k4 : sim_k4.cpp $(DEPS)
	${CC} ${CFLAGS} sim_k4.cpp -o $@ -lpthread

ex_k5 : sim_k5.cpp $(DEPS)
	${CC} ${CFLAGS} sim_k5.cpp -o $@ -lpthread

clean :
	rm ex ex_parallel ex_k4 ex_k5
