/*** CPLEX Interface ***/
#include <ilcplex/cplex.h>

/*** System Interfaces ***/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*** Public Interfaces ***/
#include "network.h"


/*** Constants Definitions ***/
static SOLUTION * solution1 = NULL;
static SOLUTION * solution2 = NULL;


/*** Forward Declarations ***/
static void free_and_null(void **);
static double objective_value(double *, double *, int);



/*** Main Routine ***/
int
main(int argc, char **argv)
{
	/* Main Variables */
	CPXENVptr	env1 = NULL;
	CPXENVptr	env2 = NULL;
	CPXNETptr	net1 = NULL;
	CPXNETptr	net2 = NULL;
	CPXLPptr	lp1 = NULL;
	CPXLPptr	lp2 = NULL;
	int status = 0;
	int i, j;

	/* Problem Info Variables */
	int narcs, nnodes;
	double * costs1 = NULL;
	double * costs2 = NULL;

	/* Sanity Check to Command Line Args */
	if(argc != 3) {
		fprintf(stderr, "Usage: ./solver [NETWORK1] [NETWORK2]\n");
		goto TERMINATE;
	}

	/*** Initialize CPLEX Environment ***/
	env1 = CPXopenCPLEX(&status);
	if(!env1) {
		char errmsg[CPXMESSAGEBUFSIZE];
		CPXgeterrorstring(env1, status, errmsg);
		fprintf(stderr, "Unable to start CPLEX, %d, %s\n", status, errmsg);
		goto TERMINATE;
	}

	env2 = CPXopenCPLEX(&status);
	if(!env2) {
		char errmsg[CPXMESSAGEBUFSIZE];
		CPXgeterrorstring(env2, status, errmsg);
		fprintf(stderr, "Unable to start CPLEX, %d, %s\n", status, errmsg);
		goto TERMINATE;
	}




	/*** Create NET problem objects ***/
	net1 = CPXNETcreateprob(env1, &status, "network1");
	if(!net1) {
		fprintf(stderr, "Unable to create NET problem object 1.\n");
		goto TERMINATE;
	}
	
	net2 = CPXNETcreateprob(env2, &status, "network2");
	if(!net2) {
		fprintf(stderr, "Unable to create NET problem object 2.\n");
		goto TERMINATE;
	}




	/*** Fill the problem data ***/
	status = CPXNETreadcopyprob(env1, net1, argv[1]);
	if(status) {
		fprintf(stderr, "Unable to copy problem to NET object 1.\n");
		goto TERMINATE;
	}

	status = CPXNETreadcopyprob(env2, net2, argv[2]);
	if(status) {
		fprintf(stderr, "Unable to copy problem to NET object.\n");
		goto TERMINATE;
	}





	/* Copy to LP objects */
	lp1 = CPXcreateprob(env1, &status, "lp1");
	if(!lp1) {
		fprintf(stderr, "Unable to create LP problem object 1.\n");
		goto TERMINATE;
	}

	status = CPXcopynettolp(env1, lp1, net1);
	if(status) {
		fprintf(stderr, "Unable to copy problem 1.\n");
		goto TERMINATE;
	}

	lp2 = CPXcreateprob(env2, &status, "lp2");
	if(!lp2) {
		fprintf(stderr, "Unable to create LP problem object 2.\n");
		goto TERMINATE;
	}

	status = CPXcopynettolp(env2, lp2, net2);
	if(status) {
		fprintf(stderr, "Unable to copy problem 2.\n");
		goto TERMINATE;
	}




	/* Set CPLEX Parameters */
	status = CPXsetintparam(env1, CPX_PARAM_SCRIND, CPX_OFF);
	if(status) {
		fprintf(stderr, "Unable to set screen output.\n");
		goto TERMINATE;
	}

	status = CPXsetintparam(env1, CPX_PARAM_ITLIM, 1);
	if(status) {
		fprintf(stderr, "Unable to set iteration limit to 1.\n");
		goto TERMINATE;
	}

	status = CPXsetintparam(env2, CPX_PARAM_SCRIND, CPX_OFF);
	if(status) {
		fprintf(stderr, "Unable to set screen output.\n");
		goto TERMINATE;
	}

	status = CPXsetintparam(env2, CPX_PARAM_ITLIM, 1);
	if(status) {
		fprintf(stderr, "Unable to set iteration limit to 1.\n");
		goto TERMINATE;
	}




/*** ********************** MEMORY ALLOCATION STAGE ************************ ***/
	narcs = CPXgetnumcols(env1, lp1);
	nnodes = CPXgetnumrows(env1, lp1);

	solution1 = create_network_solution(narcs, nnodes);
	if(!solution1) {
		fprintf(stderr, "Error on solution 1 alloc.\n");
		goto TERMINATE;
	}
	
	solution2 = create_network_solution(narcs, nnodes);
	if(!solution2) {
		fprintf(stderr, "Error on solution 2 alloc.\n");
		goto TERMINATE;
	}

/*** *********************** OPTIMIZATION STAGE **************************** ***/
	int itcnt = 0;
	double obj2_val = 0.0;
	narcs = CPXgetnumcols(env1, lp1);
	
	while(solution1->solstat != 1) {
		/* Read basis file if not in initial iteration */
		if(itcnt != 0) {
			status = CPXreadcopybase(env1, lp1, "basis");
			if(status) {
				fprintf(stderr, "Error reading basis.\n");
				goto TERMINATE;
			}
		}

		/* Optimize objective function 1 */
		status = CPXprimopt(env1, lp1);
		if(status) {
			fprintf(stderr, "Error during optimization.\n");
			goto TERMINATE;
		}

		/* Get objective function 1 solution */
		status = CPXsolution(env1, lp1, &(solution1->solstat), &(solution1->objval), 
		                     solution1->x, solution1->pi, solution1->slack, solution1->dj);
		if(status) {
			fprintf(stderr, "Unable to get solution.\n");
			goto TERMINATE;
		}

		fprintf(stdout, "\n\n------------------------------------------------------\n");
		fprintf(stdout, "Iteration: %d\n", itcnt);
		fprintf(stdout, "Objective Function 1: %lf\n", solution1->objval);
		int i, j;
		for(i = 0; i < narcs; i++) {
			fprintf(stdout, "X: %lf\tDJ: %lf\n", solution1->x[i], solution1->dj[i]);
		}

		for(j = 0; j < nnodes; j++) {
			fprintf(stdout, "PI: %lf\tSLACK: %lf\n", solution1->pi[j], solution1->slack[j]);
		}

		/* Get objective value of function 2 */
		costs2 = malloc(narcs * sizeof(double));
		if(!costs2) {
			fprintf(stderr, "Unable to alloc costs2 array.\n");
			goto TERMINATE;
		}

		status = CPXgetobj(env2, lp2, costs2, 0, narcs - 1);
		if(status) {
			fprintf(stderr, "Unable to get objective function 2 costs.\n");
			goto TERMINATE;
		}

		obj2_val = objective_value(costs2, solution1->x, narcs);

		fprintf(stdout, "Function 2: %lf\n", obj2_val);

		/* Print new basis to basis file */
		status = CPXmbasewrite(env1, lp1, "basis");
		if(status) {
			fprintf(stderr, "Error printing new basis to file.\n");
			goto TERMINATE;
		}

		itcnt++;
	}





//------------------------------------------------------------------------------	
	//while(/* Control variable to assert solution status */) {
		/* 1) Check if a base exists and if yes, pass it to CPLEX */
		//if(/* Base existance boolean */) {
			/* Pass it to CPLEX */
			
		//}
		
		/* 2) CPLEX performs an iteration only */
		//status = CPXNETprimopt(env, net);
		//if(status) {
			//fprintf(stderr, "Unable to optimize problem.\n");
			//goto TERMINATE;
		//}
		
		/* 3) Get the new base from CPLEX NET object and store it in memory */
		
		

	//}
//------------------------------------------------------------------------------



	/* Get Solution Data:
	 * Alloc memory for the solution arrays
	 * Get the solution data using CPXsolution routine
	 */


TERMINATE:

	/* Free alloc'd memory */
	free_and_null_solution(&solution1);
	free_and_null_solution(&solution2);

	free_and_null((void *) &costs1);
	free_and_null((void *) &costs2);


	/* Free problem structure */
	status = CPXNETfreeprob(env1, &net1);
	if(net1) {
		fprintf(stderr, "Unable to free NET1 problem object.\n");
	}
	
	status = CPXNETfreeprob(env2, &net2);
	if(net2) {
		fprintf(stderr, "Unable to free NET2 problem object.\n");
	}




	/* Close CPLEX */
	status = CPXcloseCPLEX(&env1);
	if(env1) {
		fprintf(stderr, "Unable to close CPLEX.\n");
	}
	
	status = CPXcloseCPLEX(&env2);
	if(env2) {
		fprintf(stderr, "Unable to close CPLEX.\n");
	}

	return status;
} /* END MAIN */






/*** Functions Definitions ***/
void free_and_null(void **ptr)
{
	if(ptr) {
		free(*ptr);
		ptr = NULL;
	}
} /* END FREE_AND_NULL */


/* The function objective_value receives two arrays with the objective values of
 * the objective function and the flow for each variable and the number of
 * variables (arcs) and returns the objective value.
 */
double objective_value(double * objs, double * flow, int narcs)
{
	if(!objs || !flow) {
		fprintf(stderr, "Unable to calculate objective flow due to NULL array.\n");
		return 0.0;
	}

	double value = 0.0;
	int i = 0;

	while(i < narcs) {
		value += objs[i] * flow[i];
		i++;
	}

	return value;
}
