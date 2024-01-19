all : ex ex_parallel ex_symm

CC = g++ -std=c++11
CFLAGS = -O3

ex : sim.cpp $(DEPS)
	${CC} ${CFLAGS} sim.cpp -o $@

ex_parallel : sim_parallel.cpp $(DEPS)
	${CC} ${CFLAGS} sim_parallel.cpp -o $@ -lpthread

ex_symm : sim_symmetries.cpp $(DEPS)
	${CC} ${CFLAGS} sim_symmetries.cpp iqp_swap_symmetries.cpp -o $@

test_symmetries : iqp_swap_symmetries.cpp test_symmetries.cpp $(DEPS)
	${CC} -g test_symmetries.cpp -o $@

clean :
	rm ex ex_parallel ex_symm test_symmetries
