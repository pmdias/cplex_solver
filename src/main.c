#include <ilcplex/cplex.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*** Constants Definitions ***/




/*** Forward Declarations ***/
static void free_and_null(void **);



/*** Data Structures ***/
struct network_basis {
	int * astat; /* Basis status for problem arcs */
	int * nstat; /* Basis status for problem nodes */
};

struct network_solution {
	int netstat;	/* Solution Status */
	double objval;	/* Objective value of solution */
	double * x;		/* Flow values */
	double * pi;	/* Pi values for nodes */
	double * dj;	/* Reduced costs for arcs */
	double * slack;	/* Slack values for nodes */
};





/*** Main Routine ***/
int
main(int argc, char **argv)
{
	/* Main Variables */
	CPXENVptr	env = NULL;
	CPXNETptr	net = NULL;
	CPXLPptr	lp = NULL;
	int status = 0;
	int i, j;

	/* Problem Info Variables */
	int narcs, nnodes;

	/* Solution Variables */
	int itcnt;
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


	/* Copy to LP */
	lp = CPXcreateprob(env, &status, argv[1]);
	if(!lp) {
		fprintf(stderr, "Unable to create LP problem object.\n");
		goto TERMINATE;
	}

	status = CPXcopynettolp(env, lp, net);
	if(status) {
		fprintf(stderr, "Unable to copy problem.\n");
		goto TERMINATE;
	}

	/* Set the CPLEX parameters, by the following order:
	 * - Set CPLEX screen output ON
	 * - Turn off the CPLEX aggregator
	 * - Set iteration limit to 1
	 */
	status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
	if(status) {
		fprintf(stderr, "Unable to set screen output.\n");
		goto TERMINATE;
	}
/*
	status = CPXsetintparam(env, CPX_PARAM_AGGIND, 0);
	if(status) {
		fprintf(stderr, "Unable to set aggregator fill to 0.\n");
		goto TERMINATE;
	}
*/
	status = CPXsetintparam(env, CPX_PARAM_ITLIM, 1);
	if(status) {
		fprintf(stderr, "Unable to set iteration limit to 1.\n");
		goto TERMINATE;
	}


	/* Read basis file */
	if(argc == 3) {
		status = CPXreadcopybase(env, lp, argv[2]);
		if(status) {
			fprintf(stderr, "Error reading basis.\n");
			goto TERMINATE;
		}
	}

	/* Optimization Stage */
	status = CPXprimopt(env, lp);
	if(status) {
		fprintf(stderr, "Unable to optimize problem.\n");
		goto TERMINATE;
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
	itcnt = CPXgetitcnt(env, lp);
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

	memset(rstat, 0, nnodes * sizeof(int));
	memset(cstat, 0, narcs * sizeof(int));

	status = CPXsolution(env, lp, &solstat, &objval, x, pi, slack, dj);
	if(status) {
		fprintf(stderr, "Unable to get solution.\n");
		goto TERMINATE;
	}

	status = CPXmbasewrite(env, lp, "basis");
	if(status) {
		fprintf(stderr, "Unable to print current basis to file\n");
		goto TERMINATE;
	}

	status = CPXgetbase(env, lp, cstat, rstat);
	if(status) {
		fprintf(stderr, "Unable to get solution basis.\n");
		goto TERMINATE;
	}

	/* Print Solution Data */
	fprintf(stdout, "Solution Status: %d\n", solstat);
	fprintf(stdout, "Solution Value: %f\n", objval);

	fprintf(stdout, "Number of arcs: %d\n", narcs);
	fprintf(stdout, "Number of nodes: %d\n", nnodes);
	fprintf(stdout, "Iteration Count: %d\n", itcnt);

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
