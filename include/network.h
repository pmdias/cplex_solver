#ifndef NETWORK_H
#define NETWORK_H 1

#include <stdlib.h>
#include <stdio.h>

typedef struct network_basis {
	int * cstat;
	int * rstat;
} BASIS;

BASIS * create_network_basis(int, int);
void free_and_null_basis(BASIS **);

typedef struct network_solution {
	int solstat;
	double objval;
	double * x;
	double * dj;
	double * pi;
	double * slack;
	BASIS * netbasis;
} SOLUTION;

SOLUTION * create_network_solution(int, int);
void free_and_null_solution(SOLUTION **);


#endif /* NETWORK_H */
