/***********************
 *** CPLEX Interface ***
 ***********************/
#include <ilcplex/cplex.h>

/*************************
 *** System Interfaces ***
 *************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>




/************************ 
 *** Type Definitions ***
 ************************/

/** The NET_BASIS struct is made with two arrays of specific size and holds the
 ** information of a basis associated to a given problem. The first array must
 ** have a length equal to the number of arcs and the second array must have a
 ** length equal to the number of nodes in the network problem.
 **/
typedef struct net_basis_struct {
	int * arc_basis;
	int * node_basis;
} NET_BASIS;

/*** Methods to handle the NET_BASIS struct ***/
static NET_BASIS * create_basis(int, int);
static void free_basis(NET_BASIS **);


/** The NET_SOLUTION struct holds all the information of a specific solution to
 ** a network problem, including the basis associated. The arrays x and dj must
 ** have a length equal to the number of arcs in the network and the arrays pi
 ** and slack must have a length equal to the number of nodes. The integer value
 ** solstat is an indicator of the status of the solution, with 1 being an
 ** optimal solution and 10 indicating the optimization was stopped due to the
 ** limit of iterations being reached.
 **/
typedef struct net_sol_struct {
	double * x;
	double * dj;
	double * pi;
	double * slack;
	double objval;
	int solstat;
	NET_BASIS * basis;
} NET_SOLUTION;

/*** Methods to handle the NET_SOLUTION struct ***/
static NET_SOLUTION * create_solution(int, int);
static void free_solution(NET_SOLUTION **);
static void print_solution(NET_SOLUTION *, int, int);

/*****************************
 *** Constants Definitions ***
 *****************************/
static NET_SOLUTION * solution1 = NULL;
static NET_SOLUTION * solution2 = NULL;
static NET_SOLUTION * perturbsol = NULL;



/****************************
 *** Forward Declarations ***
 ****************************/
static int usage(int);
static void free_and_null(void **);
static double objective_value(double *, double *, int);
static int copy_cplex_problem(CPXENVptr, CPXNETptr, CPXLPptr, char *);
static int update_solution(CPXENVptr, CPXLPptr, char *, NET_SOLUTION *, char *);
static int entering_arc(double *, double *, int *, int);

static NET_SOLUTION * get_initial_objective(char * net_file);
static NET_SOLUTION * get_perturbation_solution(CPXENVptr, CPXENVptr, CPXLPptr, CPXLPptr);



/********************
 *** Main Routine ***
 ********************/
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

	/* Timestamps */
	double start, end;
	int start_st, end_st;

	/* Sanity Check to Command Line Args */
	if(usage(argc)) {
		goto TERMINATE;
	}


	


	/*** CPLEX INITIALIZATION:
	 *** The CPLEX environments are initialized here, with the LP objects also
	 *** set from the NET objects generated from the input files passed to the
	 *** solver at startup.
	 ***
	 *** During the initialization steps, the NET objects are also set free as
	 *** they aren't needed during the rest o the optimization. At the end, only
	 *** the environment and the LP objects need to be set free.
	 ***/
	
	/* env1 is the CPLEX environment of the first objective function */
	env1 = CPXopenCPLEX(&status);
	if(!env1) {
		char errmsg[CPXMESSAGEBUFSIZE];
		CPXgeterrorstring(env1, status, errmsg);
		fprintf(stderr, "Unable to start CPLEX environment 1, %d, %s\n", status, errmsg);
		goto TERMINATE;
	}

	/* the NET and LP problem objects are also created now ***/
	net1 = CPXNETcreateprob(env1, &status, "network1");
	if(!net1) {
		fprintf(stderr, "Unable to create NET problem object 1.\n");
		goto TERMINATE;
	}

	lp1 = CPXcreateprob(env1, &status, "lp1");
	if(!lp1) {
		fprintf(stderr, "Unable to create LP problem object 1.\n");
		goto TERMINATE;
	}

	/* finally the objective function problem data is copied from the file
	 * provided as argument 1 to the solver.
	 */
	status = copy_cplex_problem(env1, net1, lp1, argv[1]);
	if(status) {
		fprintf(stderr, "An error ocurred copying problem 1 data.\n");
		goto TERMINATE;
	}
	
	/***** PROCESS REPEATED FOR OBJECTIVE FUNCTION 2 ******/
	
	/* env2 is the CPLEX environment of the second objective function */
	env2 = CPXopenCPLEX(&status);
	if(!env2) {
		char errmsg[CPXMESSAGEBUFSIZE];
		CPXgeterrorstring(env2, status, errmsg);
		fprintf(stderr, "Unable to start CPLEX environment 2, %d, %s\n", status, errmsg);
		goto TERMINATE;
	}

	/* the NET and LP problem objects are also created now ***/
	net2 = CPXNETcreateprob(env2, &status, "network2");
	if(!net2) {
		fprintf(stderr, "Unable to create NET problem object 2.\n");
		goto TERMINATE;
	}

	lp2 = CPXcreateprob(env2, &status, "lp2");
	if(!lp2) {
		fprintf(stderr, "Unable to create LP problem object 2.\n");
		goto TERMINATE;
	}

	/* finally the objective function problem data is copied from the file
	 * provided as argument 2 to the solver.
	 */
	status = copy_cplex_problem(env2, net2, lp2, argv[2]);
	if(status) {
		fprintf(stderr, "An error ocurred copying problem 2 data.\n");
		goto TERMINATE;
	}






	/*** CPLEX PARAMETERS SETTINGS:
	 *** parameters are set for both environments in this stage. Below is a list
	 *** of the settings apllied to both the CPLEX environments. If there is any
	 *** change to the code, please post it here also:
	 *** 	- CPLEX SCREEN OUTPUT   = OFF
	 ***	- CPLEX ITERATION LIMIT = 1
	 ***/
	
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






	/*** MEMORY ALLOCATION STAGE:
	 *** memory is alloc'd for the various number of objects that are used during
	 *** the optimization, like structs to save information about basis and
	 *** solution. If this stage fails, the solver MUST be immediatly terminated
	 *** and if not the user takes the risk of facing potentially serious
	 *** SEGFAULTS furing the optimization stage as CPLEX will try to access
	 *** memory that wasn't alloc'd to the solver.
	 ***
	 *** We make the assumption that the number of rows and colums in both
	 *** objective functions is the same. If not, the risk of a SEGFAULT is also
	 *** very probable due to the same reasons as stated above.
	 ***/

	/* Geting the number of arcs and number of nodes with:
	 *		n. of arcs  = n. of cols
	 *		n. of nodes = n. of rows
	 */
	narcs  = CPXgetnumcols(env1, lp1);
	nnodes = CPXgetnumrows(env1, lp1);


	/* solution1 object holds the status of the solution relative to the first
	 * objective function.
	 */
	solution1 = create_solution(narcs, nnodes);
	if(!solution1) {
		fprintf(stderr, "Error on solution 1 alloc.\n");
		goto TERMINATE;
	}
	
	/* solution2 object holds the information of the solution relative to the
	 * second objective function. Notice that the data stored in solution2->basis
	 * must be the same ALWAYS as the data stored in solution1->basis. If this
	 * does not hold, something is wrong with the optimization.
	 */
	solution2 = create_solution(narcs, nnodes);
	if(!solution2) {
		fprintf(stderr, "Error on solution 2 alloc.\n");
		goto TERMINATE;
	}

	/* the perturbsol object holds the initial basis for the bi-objective problem.
	 * this object is used only two times: when finding the initial basis and
	 * when passing it to the bi-objective solving loop. As its information is
	 * necessary only during a brief period, with may be used for other purposes
	 * if needed, as its only set free at termination like the rest of the
	 * objects.
	 */
	perturbsol = create_solution(narcs, nnodes);
	if(!perturbsol) {
		fprintf(stderr, "Error on perturbation solution alloc.\n");
		goto TERMINATE;
	}









	/*** OPTIMIZATION STAGE:
	 *** the optimization stage is made in several steps and must be followed
	 *** with precision or the solver fails to generate the correct results (at
	 *** least in a computational sense). The steps are the following:
	 *** 	- Get the global minimum for objective function 1 and 2
	 ***	- Get the initial basis using the perturbation of both objective functions.
	 ***	  The method of the perturbation is Z(x) = 0.99*z1(x) + 0.01*z2(x).
	 *** 	- Start the while loop with the stop condition of solution2->value === temp_sol2->value
	 ***	- Get the reduced costs of both networks and identify the arcs that may
	 ***	  enter the basis next.
	 ***    - Perform ratio testing to finding the entering arc.
	 ***    - Build the new basis (not a valid basis but may be used to tell CPLEX
	 ***      which arc is entering).
	 ***	- Print information to screen and/or file.
	 ***    - Check stop criteria.
	 ***/

	NET_SOLUTION * initial_sol2 = get_initial_objective(argv[2]);
	if(!initial_sol2) {
		fprintf(stderr, "Failed to get global objective 2 minimum.\n");
		goto TERMINATE;
	}

	printf("Objective 2 Objective Value Min: %lf\n\n\n", initial_sol2->objval);


	perturbsol = get_perturbation_solution(env1, env2, lp1, lp2);
	if(!perturbsol) {
		fprintf(stderr, "Error on perturbation method..\n");
		goto TERMINATE;
	}
	


	/* Copy the initial basis found by the perturbation and stored in perturbsol
	 * to both solution1 and solution2.
	 */
	memcpy(perturbsol, solution1, sizeof(NET_SOLUTION));
	memcpy(perturbsol, solution2, sizeof(NET_SOLUTION));
	
	
	
	status = update_solution(env1, lp1, "pbasis", solution1, "basis1");
	if(status) {
		goto TERMINATE;
	}
	
	print_solution(solution1, narcs, nnodes);
	
	status = update_solution(env2, lp2, "pbasis", solution2, "basis2");
	if(status) {
		goto TERMINATE;
	}
	
	print_solution(solution2, narcs, nnodes);
	
	
	
	/* Start the optimization loop for the bi-objective problem now if everything
	 * went right and no failures where detected.
	 */

	int it_test_cnt = 0;
	while(solution2->objval > initial_sol2->objval) {
	
		status = CPXgettime(env2, &start);
		if(status) {
			fprintf(stderr, "Unable to get time.\n");
			goto TERMINATE;
		}
	
		/* Find the entering arc */
		int arc = entering_arc(solution1->dj, solution2->dj, solution2->basis->arc_basis, narcs);
		if(arc == -1) {
			fprintf(stderr, "Error calculating entering arc. Break.\n");
			goto TERMINATE;
			//break;
		}
		
		printf("Entering arc: %d\n", arc);
		
		/* Turn off presolve and set parameters for CPLEX accept advanced basis */
		status = CPXsetintparam(env2, CPX_PARAM_ADVIND, 2);
		if(status) {
			fprintf(stderr, "Unable to set advanced start switch to 1.\n");
			goto TERMINATE;
		}
		
		status = CPXsetintparam(env2, CPX_PARAM_PREIND, CPX_OFF);
		if(status) {
			fprintf(stderr, "Error turning off presolve. ERROR %d\n", status);
			goto TERMINATE;
		}
		
		status = CPXsetintparam(env2, CPX_PARAM_AGGIND, 0);
		if(status) {
			fprintf(stderr, "Error turning off aggregator. ERROR: %d\n", status);
			goto TERMINATE;
		}
		
		status = CPXsetintparam(env2, CPX_PARAM_DEPIND, 0);
		if(status) {
			fprintf(stderr, "Error turning off DEPIND. ERROR: %d\n", status);
			goto TERMINATE;
		}
				
		status = CPXsetintparam(env2, CPX_PARAM_PREDUAL, -1);
		if(status) {
			fprintf(stderr, "Error turning off PREDUAL. ERROR: %d\n", status);
			goto TERMINATE;
		}
		
		status = CPXsetintparam(env2, CPX_PARAM_PREPASS, 0);
		if(status) {
			fprintf(stderr, "Error turning off PREPASS. ERROR: %d\n", status);
			goto TERMINATE;
		}
		
		status = CPXsetintparam(env2, CPX_PARAM_SCAIND, -1);
		if(status) {
			fprintf(stderr, "Error turning off SCAIND. ERROR: %d\n", status);
			goto TERMINATE;
		}
		
		status = CPXsetintparam(env2, CPX_PARAM_SIMDISPLAY, 2);
		if(status) {
			fprintf(stderr, "Error turning off SCAIND. ERROR: %d\n", status);
			goto TERMINATE;
		}
		
		
		
		/* Enter the arc using CPXpivot */
		status = CPXpivot(env2, lp2, arc, CPX_NO_VARIABLE, CPX_AT_LOWER);
		if(status) {
			fprintf(stderr, "CPXpivot failed.\n");
			goto TERMINATE;
		}
		
		
		status = CPXgettime(env2, &end);
		if(status) {
			fprintf(stderr, "Unable to get time.\n");
			goto TERMINATE;
		}
		fprintf(stdout, "Time Elapsed: %lf s\n", end - start);
	
	
		
		/* Get the solution */
		status = CPXsolution(env2, lp2, &solution2->solstat, &solution2->objval,
							 solution2->x, solution2->pi, solution2->slack, solution2->dj);
		if(status) {
			fprintf(stderr, "Error getting solution at end of loop.\n");
			goto TERMINATE;
		}

		status = CPXgetbase(env2, lp2, solution2->basis->arc_basis, solution2->basis->node_basis);
		if(status) {
			fprintf(stderr, "Error getting base at end of loop.\n");
			goto TERMINATE;
		}
		
		
		/* Update all the stuff and print the solutions */
		status = CPXmbasewrite(env2, lp2, "basis2");
		if(status) {
			fprintf(stderr, "Error %d\n", status);
			goto TERMINATE;
		}
		
		status = update_solution(env2, lp2, "basis2", solution2, "basis2");
		if(status) {
			goto TERMINATE;
		}
		
		status = update_solution(env1, lp1, "basis2", solution1, "basis1");
		if(status) {
			goto TERMINATE;
		}
	
		print_solution(solution1, narcs, nnodes);
		print_solution(solution2, narcs, nnodes);
	
		it_test_cnt++;
	}


	/* Get Solution Data:
	 * Alloc memory for the solution arrays
	 * Get the solution data using CPXsolution routine
	 */





TERMINATE:

	/* Free alloc'd memory */
	free_solution(&solution1);
	free_solution(&solution2);
	free_solution(&perturbsol);

	free_and_null((void *) &costs1);
	free_and_null((void *) &costs2);





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

/** Function: usage
 ** A function that receives the number of inputs given to the main function
 ** and if an error is detected, it prints the necessary usage information to
 ** the standard error output stream.
 **/
int usage(int argc) {
	if(argc != 3) {
		fprintf(stderr, "Usage: ./solver [NETWORK1] [NETWORK2]\n");
		return 1;
	}
	
	return 0;
}


/** Function: free_and_null
 ** The function receives a double pointer to some previously allocated memory
 ** and frees that memory also setting all the necessary pointers to NULL
 **/
void free_and_null(void **ptr)
{
	if(ptr) {
		free(*ptr);
		ptr = NULL;
	}
}


/** The function objective_value receives two arrays with the objective values of
 ** the objective function and the flow for each variable and the number of
 ** variables (arcs) and returns the objective value.
 **/
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


/** Function: copy_cplex_problem
 ** The function receives three CPLEX objects (environment, net and lp) and a
 ** string with the filename of the network to be used and copies the network to
 ** the objects performing all the initialization necessary. If it fails, the
 ** function returns a value of -1. Otherwise, it returns a value of 0 or the
 ** error code of the error that happened. The CPLEX objetcs must be passed
 ** to the function already initialized by the respective CPLEX methods.
 **/
int copy_cplex_problem(CPXENVptr env, CPXNETptr net, CPXLPptr lp, char * filename)
{
	if(!env) {
		fprintf(stderr, "Unable to copy problem due to NULL environment.\n");
		return -1;
	}

	if(!net || !lp) {
		fprintf(stderr, "Unable to copy problem due to NULL Net or Lp object.\n");
		return -1;
	}

	if(!filename) {
		fprintf(stderr, "Unable to copy problem due to NULL file.\n");
		return -1;
	}

	int status = 0;
	
	status = CPXNETreadcopyprob(env, net, filename);
	if(status) {
		fprintf(stderr, "Unable to copy problem to NET object.\n");
		return status;
	}
	
	status = CPXcopynettolp(env, lp, net);
	if(status) {
		fprintf(stderr, "Unable to copy problem to LP object.\n");
		return status;
	}
	
	status = CPXNETfreeprob(env, &net);
	if(net) {
		fprintf(stderr, "Unable to free NET problem object.\n");
	}
	
	return status;
}


/** The function update_solution performs a single iteration that updates the
 ** information relative to the current basis stored in the LP problem object
 ** that is passed as an argument to the function.
 ** If any error is detected, the function returns a non-zero value and an error
 ** message is logged into standard error output stream.
 **/
int update_solution(CPXENVptr env, CPXLPptr lp, char * basis_file, NET_SOLUTION * sol, char * output_basis)
{
	int status = 0;
	
	/* Check if the input arguments are valid */
	if(!env || !lp || !basis_file || !sol || !output_basis) {
		fprintf(stderr, "Error: NULL pointer\n");
		status = 1;
		goto TERMINATE;
	}
	
	/* Set the iteration limit to zero so that CPLEX doesn't change the current
	 * basis
	 */
	status = CPXsetintparam(env, CPX_PARAM_ITLIM, 0);
	if(status) {
		fprintf(stderr, "Unable to set iteration limit to 0.\n");
		goto TERMINATE;
	}
	
	
	status = CPXreadcopybase(env, lp, basis_file);
	if(status) {
		fprintf(stderr, "Error reading basis from file %s\n", basis_file),
		status = 2;
		goto TERMINATE;
	}
	
	status = CPXprimopt(env, lp);
	if(status) {
		fprintf(stderr, "Error during optimization. ERROR %d\n", status);
		goto TERMINATE;
	}
	
	status = CPXsolution(env, lp, &sol->solstat, &sol->objval, sol->x, sol->pi, sol->slack, sol->dj);
	if(status) {
		fprintf(stderr, "Error getting solution. ERROR %d\n", status);
		goto TERMINATE;
	}
	
	status = CPXgetbase(env, lp, sol->basis->arc_basis, sol->basis->node_basis);
	if(status) {
		fprintf(stderr, "Error getting basis. ERROR %d\n", status);
		goto TERMINATE;
	}

	status = CPXmbasewrite(env, lp, output_basis);
	if(status) {
		fprintf(stderr, "Error printing basis. ERROR %d\n", status);
		goto TERMINATE;
	}
	
	status = CPXsetintparam(env, CPX_PARAM_ITLIM, 1);
	if(status) {
		fprintf(stderr, "Unable to set iteration limit to 1.\n");
		goto TERMINATE;
	}
	
TERMINATE:

	return status;
}


/** Function: entering_arc
 ** This method receives three arrays with the reduced costs of both objective
 ** functions and a valid basis and also receives an integer equal to the number
 ** of arcs in the problem (size of the arrays). It returns an integer equal to
 ** the index of the arc that will enter the basis.
 **/
int entering_arc(double * dj1, double * dj2, int * basis, int size)
{
	int arc = 0;
	int i;
	
	/* Sanity check of input data arrays --- if one of them its NULL the function
	 * exists with a value of -1
	 */
	if(!dj1 || !dj2 || !basis) {
		fprintf(stderr, "Error due to NULL pointer when calculating entering arc.\n");
		return -1;
	}
	
	double * ratios = malloc(size * sizeof(double));
	if(!ratios) {
		fprintf(stderr, "Unable to alloc ratios array.\n");
		return -1;
	}
	
	/* Calculate the ratio between the reduced costs of objective function 1 and
	 * objective function 2.
	 */
	for(i = 0; i < size; i++) {
		/* The ratio test must follow:
		 * if basis[i] != 1 AND basis[i] == 0 (arc at lower bound) ----> dj2[i] must be < 0
		 * if basis[i] != 1 AND basis[i] == 2 (arc at upper bound) ----> dj2[i] must be > 0
		 */
		if(basis[i] != 1 && basis[i] == 0 && dj2[i] < 0) {
			ratios[i] = dj2[i] / dj1[i];
			printf("DJ1: %lf DJ: %lf\n", dj1[i], dj2[i]);
			printf("Ratio %d: %lf\n", i, ratios[i]);
		} else if(basis[i] != 1 && basis[i] == 2 && dj2[i] > 0) {
			ratios[i] = dj2[i] / dj1[i];
			printf("DJ1: %lf DJ: %lf\n", dj1[i], dj2[i]);
			printf("Ratio %d: %lf\n", i, ratios[i]);
		} else {
			ratios[i] = 0;
		}
	}

	/* Find the best ratio, i.e., the slope that gives the highest rate of
	 * decrease in the objective function 2.
	 */
	
	for(i = 0; i < size; i++) {
		if(ratios[i] != 0) {
			if(ratios[i] < ratios[arc]) {
				arc = i;
			}
		}
	}
	
	free(ratios);
	
	return arc;
}

/** Function: create_basis
 ** the function receives two variables of integer type equal to the number of
 ** arcs and nodes of the problem and returns a pointer to a NET_BASIS object.
 ** If any error is detected, an error is printed and NULL is returned.
 **/
NET_BASIS * create_basis(int narcs, int nnodes)
{
	NET_BASIS * basis = malloc(sizeof(NET_BASIS));
	if(!basis) {
		fprintf(stderr, "Unable to alloc basis object.\n");
		return NULL;
	}

	basis->arc_basis = malloc(narcs * sizeof(int));
	basis->node_basis = malloc(nnodes * sizeof(int));
	if(!(basis->arc_basis) || !(basis->node_basis)) {
		fprintf(stderr, "Unable to alloc basis arrays.\n");
		return NULL;
	}

	return basis;
}


/** Function: free_basis
 ** The function receives a double pointer to a NET_BASIS object that was
 ** previously allocated into dynamic memory and frees that memory.
 **/
void free_basis(NET_BASIS ** basis)
{
	if(basis) {
		free(*basis);
		basis = NULL;
	}
}


/** Function: create_solution
 ** This function receives two integers equal to the number of arcs and nodes in
 ** the network and creates an object to hold information of solutions of that
 ** problem, carrying out all the initialization procedurs necessary. If any
 ** error is found, the function returns NULL, otherwise it returns a pointer to
 ** the object.
 **/
NET_SOLUTION * create_solution(int narcs, int nnodes)
{
	NET_SOLUTION * solution = (NET_SOLUTION *) malloc(sizeof(NET_SOLUTION));
	if(!solution) {
		fprintf(stderr, "Unable to alloc solution.\n");
		return NULL;
	}

	solution->x = malloc(narcs * sizeof(double));
	solution->dj = malloc(narcs * sizeof(double));
	solution->pi = malloc(nnodes * sizeof(double));
	solution->slack = malloc(nnodes * sizeof(double));
	solution->objval = 0.0;
	solution->solstat = 0;
	solution->basis = create_basis(narcs, nnodes);

	if(!(solution->x) || !(solution->dj) || !(solution->pi) || !(solution->slack) || !(solution->basis)) {
		fprintf(stderr, "Error on alloc of solution arrays.\n");
		return NULL;
	}

	return solution;
}


/** Function: free_solution
 ** the function receives a double pointer to a NET_SOLUTION object that was
 ** created using the create_solution method and frees the memory of that object
 ** also making all pointers NULL
 **/
void free_solution(NET_SOLUTION ** solution)
{
	if(solution) {
		//free_basis((NET_BASIS **) &(*solution)->basis);
		free(*solution);
		solution = NULL;
	}
}


/** Function: print_solution
 ** The function receives a NET_SOLUTION object and prints the solution stored
 ** in that object to the standard output stream.
 **/
void print_solution(NET_SOLUTION * solution, int narcs, int nnodes)
{
	int i;
	
	if(!solution) {
		return;
	}
	
	fprintf(stdout, "********************************************************************\n");
	fprintf(stdout, "Printing Solution Data:\n\n");
	
	fprintf(stdout, "Objective Value:\t\t%lf\n", solution->objval);
	fprintf(stdout, "Objective Status:\t\t%d\n\n", solution->solstat);
	
	fprintf(stdout, "Objective Arc Data:\n");
	for(i = 0; i < narcs; i++) {
		fprintf(stdout, "Arc %d\tx: %lf\t reduced cost: %lf\t\tbasis: %d\n", i, solution->x[i], solution->dj[i], solution->basis->arc_basis[i]);
	}
	
	fprintf(stdout, "Objective Node Data:\n");
	for(i = 0; i < nnodes; i++) {
		fprintf(stdout, "Node %d\tpi: %lf\t slack: %lf\t\tbasis: %d\n", i, solution->pi[i], solution->slack[i], solution->basis->node_basis[i]);
	}
	
	fprintf(stdout, "\n");
}


/** Function: get_initial_objective
 ** The function receives a string with a filename to a file that holds information
 ** of a network problem and than it solves that problem until an optimal
 ** solution is reached and returns a pointer to a NET_SOLUTION object with all
 ** the information of that solution. If any error is detected, the function 
 ** returns NULL.
 **/
NET_SOLUTION * get_initial_objective(char * net_file)
{
	if(!net_file) {
		return NULL;
	}
	
	int status = 0;
	CPXENVptr	env = NULL;
	CPXNETptr	net = NULL;
	CPXLPptr	lp = NULL;
	
	env = CPXopenCPLEX(&status);
	if(!env) {
		char errmsg[CPXMESSAGEBUFSIZE];
		CPXgeterrorstring(env, status, errmsg);
		fprintf(stderr, "Unable to start CPLEX free env, %d, %s\n", status, errmsg);
		goto TERMINATE;
	}

	net = CPXNETcreateprob(env, &status, "network_free");
	if(!net) {
		fprintf(stderr, "Unable to create NET free problem object.\n");
		goto TERMINATE;
	}

	lp = CPXcreateprob(env, &status, "lp_free");
	if(!lp) {
		fprintf(stderr, "Unable to create LP free problem object.\n");
		goto TERMINATE;
	}

	status = CPXNETreadcopyprob(env, net, net_file);
	if(status) {
		fprintf(stderr, "Unable to copy problem to NET free object.\n");
		goto TERMINATE;
	}
	
	status = CPXcopynettolp(env, lp, net);
	if(status) {
		fprintf(stderr, "Unable to copy free problem.\n");
		goto TERMINATE;
	}

	status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	if(status) {
		fprintf(stderr, "Unable to set screen output.\n");
		goto TERMINATE;
	}

	



	int narcs = CPXgetnumcols(env, lp);
	int nnodes = CPXgetnumrows(env, lp);

	NET_SOLUTION * solution = create_solution(narcs, nnodes);
	if(!solution) {
		fprintf(stderr, "Error on free solution alloc.\n");
		goto TERMINATE;
	}

	status = CPXprimopt(env, lp);
	if(status) {
		fprintf(stderr, "Error during optimization of free problem.\n");
		goto TERMINATE;
	}

	status = CPXsolution(env, lp, &(solution->solstat), &(solution->objval), 
		                     solution->x, solution->pi, solution->slack, solution->dj);
	if(status) {
		fprintf(stderr, "Unable to get solution.\n");
		goto TERMINATE;
	}

	status = CPXgetbase(env, lp, solution->basis->arc_basis, solution->basis->node_basis);
	if(status) {
		fprintf(stderr, "Unable to get basis.\n");
		goto TERMINATE;
	}

TERMINATE:
	
	status = CPXNETfreeprob(env, &net);
	if(net) {
		fprintf(stderr, "Unable to free NET free problem object.\n");
	}
	
	status = CPXfreeprob(env, &lp);
	if(lp) {
		fprintf(stderr, "Unable to free LP free problem object.\n");
	}

	status = CPXcloseCPLEX(&env);
	if(env) {
		fprintf(stderr, "Unable to close CPLEX.\n");
	}

	return solution;
}



NET_SOLUTION * get_perturbation_solution(CPXENVptr env1, CPXENVptr env2, CPXLPptr lp1, CPXLPptr lp2)
{
	int status = 0;
	CPXENVptr penv = NULL;
	CPXLPptr plp = NULL;
	NET_SOLUTION * solution = NULL;
	
	double * costs1 = NULL;
	double * costs2 = NULL;
	double * costs3 = NULL;
	int narcs, nnodes;

	if(!env1 || !env2 || !lp1 || !lp2) {
		fprintf(stderr, "Unable to perform perturbation.\n");
		return NULL;
	}

	penv = CPXopenCPLEX(&status);
	if(!penv) {
		char errmsg[CPXMESSAGEBUFSIZE];
		CPXgeterrorstring(penv, status, errmsg);
		fprintf(stderr, "Unable to start CPLEX free env, %d, %s\n", status, errmsg);
		goto TERMINATE;
	}

	plp = CPXcloneprob(penv, lp1, &status);
	if(!plp) {
		fprintf(stderr, "Unable to create LP free problem object.\n");
		goto TERMINATE;
	}	


	/* Get the number of arcs and nodes ***/
	narcs = CPXgetnumcols(env1, lp1);
	nnodes = CPXgetnumrows(env1, lp1);


	/* Get the cost arrays for each objective function ***/
	costs1 = malloc(narcs * sizeof(double));
	if(!costs1) {
		fprintf(stderr, "Unable to alloc costs 1 array.\n");
		goto TERMINATE;
	}

	status = CPXgetobj(env1, lp1, costs1, 0, narcs - 1);
	if(status) {
		fprintf(stderr, "Unable to get objective function 1 costs.\n");
		goto TERMINATE;
	}

	costs2 = malloc(narcs * sizeof(double));
	if(!costs2) {
		fprintf(stderr, "Unable to alloc costs 2 array.\n");
		goto TERMINATE;
	}

	status = CPXgetobj(env2, lp2, costs2, 0, narcs - 1);
	if(status) {
		fprintf(stderr, "Unable to get objective function 2 costs.\n");
		goto TERMINATE;
	}
	
	/* Build the new cost array and pass it to the new LP object */
	costs3 = malloc(narcs * sizeof(double));
	if(!costs3) {
		fprintf(stderr, "Unable to alloc costs 2 array.\n");
		goto TERMINATE;
	}

	int * index_list = malloc(narcs * sizeof(int));
	if(!index_list) {
		fprintf(stderr, "Error on index_list alloc.\n");
		goto TERMINATE;
	}
	
	int i;
	for(i = 0; i < narcs; i++) {
		costs3[i] = 0.999 * costs1[i] + 0.001 * costs2[i];
		index_list[i] = i;
	}
	
	status = CPXchgobj(penv, plp, narcs, index_list, costs3);
	if(status) {
		fprintf(stderr, "Unable to change objective values.\n");
		goto TERMINATE;
	}
	
	/* Alloc the solution object */
	solution = create_solution(narcs, nnodes);
	if(!solution) {
		fprintf(stderr, "Unable to alloc solution.\n");
		goto TERMINATE;
	}
	
	/* Optimize */
	status = CPXprimopt(penv, plp);
	if(status) {
		fprintf(stderr, "Unable to optimize problem.\n");
		goto TERMINATE;
	}
	
	/* Get solution */
	status = CPXsolution(penv, plp, &(solution->solstat), &(solution->objval), 
		                     solution->x, solution->pi, solution->slack, solution->dj);
	if(status) {
		fprintf(stderr, "Unable to get solution.\n");
		goto TERMINATE;
	}
	
	/* Get the basis */
	status = CPXmbasewrite(penv, plp, "pbasis");
	if(status) {
		fprintf(stderr, "Unable to print basis.\n");
		goto TERMINATE;
	}
	
	status = CPXgetbase(penv, plp, solution->basis->arc_basis, solution->basis->node_basis);
	if(status) {
		fprintf(stderr, "Error getting basis from perturbation into struct.\n");
		goto TERMINATE;
	}
	

TERMINATE:

	free(costs1);
	free(costs2);
	
	status = CPXfreeprob(penv, &plp);
	if(status) {
		fprintf(stderr, "Unable to free LP object.\n");
	}
	
	status = CPXcloseCPLEX(&penv);
	if(penv) {
		fprintf(stderr, "Unable to close CPLEX.\n");
	}	
	
	return solution;
}
