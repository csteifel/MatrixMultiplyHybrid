all: main.c
	mpicc -Wall -lpthread main.c
