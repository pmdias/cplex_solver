#include <ilcplex/cplex.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*** Constants Definitions ***/
#define SOLVER_NETWORK_SIMPLEX	1
#define SOLVER_LP_SIMPLEX		2
#define SOLVER_HYBRID_SIMPLEX	3

/*** Forward Declarations ***/
static void free_and_null(void **);

/*** Main Routine ***/
int
main(int argc, char **argv)
{
	/* Main Variables */
	CPXENVptr	env = NULL;
	CPXNETptr	net = NULL;
	CPXLPptr	lp  = NULL;
	int status = 0;
	int i, j;

	/* Problem Info Variables */
	int narcs, nnodes;

	/* Solution Variables */
	int solstat;
	double objval;
	double * x     = NULL;
	double * pi    = NULL;
	double * dj    = NULL;
	double * slack = NULL;

	/* Basis Variables */
	int * cstat = NULL;
	int * rstat = NULL;

	/* Sanity Check to Command Line Args */
	if(argc != 2) {
		fprintf(stderr, "Usage: ./solver [FILE]\n");
		goto TERMINATE;
	}

	/* Initialize CPLEX Environment */
	env = CPXopenCPLEX(&status);
	if(!env) {
		char errmsg[CPXMESSAGEBUFSIZE];
		CPXgeterrorstring(env, status, errmsg);
		fprintf(stderr, "Unable to start CPLEX, %d, %s\n", status, errmsg);
		goto TERMINATE;
	}

	/* Create NET problem object */
	net = CPXNETcreateprob(env, &status, "network");
	if(!net) {
		fprintf(stderr, "Unable to create NET problem object.\n");
		goto TERMINATE;
	}

	/* Fill the problem data */
	status = CPXNETreadcopyprob(env, net, argv[1]);
	if(status) {
		fprintf(stderr, "Unable to copy problem to NET object.\n");
		goto TERMINATE;
	}



	/* Create the LP objects that are going to be used in the optimization stage
	 * and copy their data from the NET object.
	 */
	lp = CPXcreateprob(env, &status, "network_lp");
	if(!lp) {
		fprintf(stderr, "Unable to create LP problem object.\n");
		goto TERMINATE;
	}

	if(net) {
		status = CPXcopynettolp(env, lp, net);
		if(status) {
			fprintf(stderr, "Unable to copy NET object to LP object.\n");
			goto TERMINATE;
		}
	}



	/* Set the CPLEX parameters, by the following order:
	 * - Set CPLEX screen output ON
	 * - Turn off the CPLEX aggregator
	 */
	status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
	if(status) {
		fprintf(stderr, "Unable to set screen output.\n");
		goto TERMINATE;
	}

	status = CPXsetintparam(env, CPX_PARAM_AGGIND, 0);
	if(status) {
		fprintf(stderr, "Unable to set aggregator fill to 0.\n");
		goto TERMINATE;
	}



	/* Optimization Stage */
	status = CPXprimopt(env, lp);
	if(status) {
		fprintf(stderr, "Unable to optimize problem.\n");
		goto TERMINATE;
	}


	/* Get Solution Data:
	 * Alloc memory for the solution arrays
	 * Get the solution data using CPXsolution routine
	 */
	narcs = CPXgetnumcols(env, lp);
	nnodes = CPXgetnumrows(env, lp);

	x     = (double *) malloc(narcs * sizeof(double));
	dj    = (double *) malloc(narcs * sizeof(double));
	pi    = (double *) malloc(nnodes * sizeof(double));
	slack = (double *) malloc(nnodes * sizeof(double));

	if(!x || !dj || !pi || !slack) {
		fprintf(stderr, "Unable to alloc solution memory.\n");
		goto TERMINATE;
	}

	rstat = (int *) malloc(nnodes * sizeof(int));
	cstat = (int *) malloc(narcs * sizeof(int));

	if(!rstat || !cstat) {
		fprintf(stderr, "Unable to alloc basis memor.\n");
		goto TERMINATE;
	}

	status = CPXsolution(env, lp, &solstat, &objval, x, pi, slack, dj);
	if(status) {
		fprintf(stderr, "Unable to get solution.\n");
		goto TERMINATE;
	}

	status = CPXgetbase(env, lp , cstat, rstat);
	if(status) {
		fprintf(stderr, "Unable to get solution basis.\n");
		goto TERMINATE;
	}

	/* Print Solution Data */
	fprintf(stdout, "Solution Status: %d\n", solstat);
	fprintf(stdout, "Solution Value: %f\n", objval);

	fprintf(stdout, "\nARCS:\n");
	for(i = 0; i < narcs; i++) {
		fprintf(stdout, "Arc %d: ", i);
		fprintf(stdout, "Value: %10f ", x[i]);
		fprintf(stdout, "Reduced Costs: %10f ", dj[i]);
		/* Print the basis status */
		if(cstat[i] == CPX_BASIC) {
			fprintf(stdout, " AT BASIS\n");
		} else if(cstat[i] == CPX_AT_LOWER) {
			fprintf(stdout, " AT LOWER\n");
		} else if(cstat[i] == CPX_AT_UPPER) {
			fprintf(stdout, " AT UPPER\n");
		}
	}

	fprintf(stdout, "\nNODES:\n");
	for(j = 0; j < nnodes; j++) {
		fprintf(stdout, "Node %d: ", j);
		fprintf(stdout, "PI: %10f ", pi[j]);
		fprintf(stdout, "Slack: %10f", slack[j]);
		/* Print the basis status */
		if(rstat[j] == CPX_BASIC) {
			fprintf(stdout, " AT BASIS\n");
		} else if(rstat[j] == CPX_AT_LOWER) {
			fprintf(stdout, " AT LOWER\n");
		} else if(rstat[j] == CPX_AT_UPPER) {
			fprintf(stdout, " AT UPPER\n");
		}
	}

TERMINATE:

	/* Free alloc'd memory */
	free_and_null((void *) &x);
	free_and_null((void *) &pi);
	free_and_null((void *) &dj);
	free_and_null((void *) &slack);

	free_and_null((void *) &rstat);
	free_and_null((void *) &cstat);

	/* Free problem structure */
	status = CPXNETfreeprob(env, &net);
	if(net) {
		fprintf(stderr, "Unable to free NET problem object.\n");
	}

	status = CPXfreeprob(env, &lp);
	if(lp) {
		fprintf(stderr, "Unable to free LP problem object.\n");
	}

	/* Close CPLEX */
	status = CPXcloseCPLEX(&env);
	if(env) {
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
