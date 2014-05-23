#include "network.h"

BASIS * create_network_basis(int narcs, int nnodes)
{
	if(!narcs || !nnodes) {
		return NULL;
	}
	
	BASIS * base = malloc(sizeof(BASIS));
	if(!base) {
		fprintf(stderr, "Unable to allocate memory for basis.\n");
		return NULL;
	}
	
	base->cstat = malloc(narcs * sizeof(double));
	base->rstat = malloc(nnodes * sizeof(double));
	if(base->cstat == NULL || base->rstat == NULL) {
		fprintf(stderr, "unable to allocate memory for basis arrays.\n");
		free(base);
		return NULL;
	}
	
	return base;
}

void free_and_null_basis(BASIS ** base)
{
	if(base) {
		free((*base)->cstat);
		free((*base)->rstat);
		free(*base);
		base = NULL;
	}
}

SOLUTION * create_network_solution(int narcs, int nnodes)
{
	if(!narcs || !nnodes) {
		return NULL;
	}

	SOLUTION * sol = malloc(sizeof(SOLUTION));
	if(!sol) {
		fprintf(stderr, "Unable to allocate memory for solution.\n");
		return NULL;
	}
	
	sol->x = malloc(narcs * sizeof(double));
	sol->dj = malloc(narcs * sizeof(double));
	sol->pi = malloc(nnodes * sizeof(double));
	sol->slack = malloc(nnodes * sizeof(double));
	if(!(sol->x) || !(sol->dj) || !(sol->pi) || !(sol->slack)) {
		fprintf(stderr, "Unable to allocate memory for solution arrays.\n");
		free(sol);
		return NULL;
	}
	
	sol->solstat = 0;
	sol->objval = 0.0;
	
	sol->netbasis = create_network_basis(narcs, nnodes);
	if(sol->netbasis == NULL) {
		fprintf(stderr, "Error allocating memory for basis on solution object.\n");
		free(sol);
		return NULL;
	}
	
	return sol;
}

void free_and_null_solution(SOLUTION ** sol)
{
	if(sol) {
		free_and_null_basis((BASIS **) &((*sol)->netbasis));
		free(*sol);
		sol = NULL;
	}
}
