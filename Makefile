CC = mpicc
CFLAGS = -g -Wall -Wextra 
OBJS = lcs-mpi.o lcsAlgorithm.o readFileUtils.o utils.o

lcs-mpi: $(OBJS)
	$(CC) $(CFLAGS) -o lcs-mpi $(OBJS)

lcs-mpi.o: lcs-mpi.c structs.h readFileUtils.h readFileUtils.c lcsAlgorithm.h lcsAlgorithm.c utils.h utils.c
	$(CC) -c lcs-mpi.c

lcsAlgorithm.o: structs.h utils.h utils.c
	$(CC) -c lcsAlgorithm.c

readFileUtils.o: structs.h utils.h utils.c
	$(CC) -c readFileUtils.c

utils.o: structs.h
	$(CC) -c utils.c

clean:
	rm -f *.o *~ lcs-mpi lcsAlgorithm readFileUtils utils
