all : ex_k4_worst_case ex_k5_worst_case ex_k4_random ex_k5_random ex ex_parallel ex_k5_extra_cnot_layers

CC = g++ -std=c++11
CFLAGS = -O3

ex : sim.cpp $(DEPS)
	${CC} ${CFLAGS} sim.cpp -o $@

ex_parallel : sim_parallel.cpp $(DEPS)
	${CC} ${CFLAGS} sim_parallel.cpp -o $@ -lpthread

ex_k4_worst_case : sim_k4_worst_case.cpp $(DEPS)
	${CC} ${CFLAGS} sim_k4_worst_case.cpp -o $@ -lpthread

ex_k5_worst_case : sim_k5_worst_case.cpp $(DEPS)
	${CC} ${CFLAGS} sim_k5_worst_case.cpp -o $@ -lpthread

ex_k4_random : sim_k4_random.cpp $(DEPS)
	${CC} ${CFLAGS} sim_k4_random.cpp -o $@ -lpthread

ex_k5_random : sim_k5_random.cpp $(DEPS)
	${CC} ${CFLAGS} sim_k5_random.cpp -o $@ -lpthread

ex_k5_extra_cnot_layers : sim_k5_extra_cnot_layers.cpp $(DEPS)
	${CC} ${CFLAGS} sim_k5_random.cpp -o $@ -lpthread

clean :
	rm ex ex_parallel ex_k4_worst_case ex_k5_worst_case ex_k4_random ex_k5_random ex_k5_extra_cnot_layers
