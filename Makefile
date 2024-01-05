all : ex

CC = g++ 
#CC = g++ -std=c++11
CFLAGS = -O3



ex : ExponentialSumReal.cpp $(DEPS)
	${CC} ${CFLAGS} ExponentialSumReal.cpp -o $@


clean :
	rm *.o
