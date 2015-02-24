#include"petscInterface.h"

void petscSolve(double** matrix,int N,double* rhs, double* x1){
  Vec		 x, b, u;      /* approx solution, RHS, exact solution */
  Mat            A;            /* linear system matrix */
  KSP            ksp;         /* linear solver context */
  PC             pc;           /* preconditioner context */
  PetscReal      norm,tol=1.e-14;  /* norm of solution error */
  PetscErrorCode ierr;
  PetscInt       i,n = N;
  PetscMPIInt    size;
  //PetscScalar    one = 1.0, two=0, zero=0;
  PetscBool      nonzeroguess = PETSC_FALSE;

/*  PetscInitialize(&argc,&args,(char*)0,help);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");
  ierr = PetscOptionsGetInt(NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,"-nonzero_guess",&nonzeroguess,NULL);CHKERRQ(ierr);
*/

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Compute the matrix and right-hand-side vector that define
         the linear system, Ax = b.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
i
  /*
     Create vectors.  Note that we form 1 vector from scratch and
     then duplicate as needed.
  */
  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
  ierr = VecSetSizes(x,PETSC_DECIDE,n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&b);CHKERRQ(ierr);


// index where value will be added

int* lines;
lines=malloc(N*sizeof(int));
for (i=0;i<N;i++){
	lines[i]=i;
}

ierr = VecSetValues(b,n,lines,rhs, INSERT_VALUES);//	Set vector b from the input function(RHS)  
//ierr = VecDuplicate(x,&u);CHKERRQ(ierr);

  /*
     Create matrix.  When using MatCreate(), the matrix format can
     be specified at runtime.

     Performance tuning note:  For problems of substantial size,
     preallocation of matrix memory is crucial for attaining good
     performance. See the matrix chapter of the users manual for details.
  */


  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  /*
     Assemble matrix
  */

double* valueM=malloc((N*N)*sizeof(double));
double* linesm=malloc((N*N)*sizeof(double));
double* colsm=malloc((N*N)*sizeof(double));

// Set the entry of matrix become array 1D
int i,j;
int k=0;
for(i=0;i<N;i++){
	for(j=0;j<N;j++){
		valueM[k]=matrix[i][j];
		k++;
	}
}

//Build the array to configure lines and coloums that the value of matrix will be added

k=0;
for(i=0;i<N;i++){
        for(j=0;j<N;j++){
                linesm[k]=i;
                colsm[k]=j;
		k++;
        }
}

CHKERRQ(ierr);
  
  ierr = MatSetValues(A,N,linesm,N,colsm,valueM,INSERT_VALUES);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /*
     Set exact solution; then compute right-hand-side vector.// Set the entry of matrix become array 1D
  */
  //ierr = VecSet(u,one);CHKERRQ(ierr);
 // ierr = MatMult(A,u,b);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the linear solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create linear solver context
  */
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

  /*
     Set operators. Here the matrix that defines the linear system
     also serves as the preconditioning matrix.
  */
  ierr = KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);

  /*
     Set linear solver defaults for this problem (optional).
     - By extracting the KSP and PC contexts from the KSP context,
       we can then directly call any KSP and PC routines to set
       various options.
     - The following four statements are optional; all of these
       parameters could alternatively be specified at runtime via
       KSPSetFromOptions();
  */
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,1.e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);

  /*
    Set runtime options, e.g.,
        -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    These options will override those specified above as long as
    KSPSetFromOptions() is called _after_ any other customization
    routines.
  */
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

  if (nonzeroguess) {
    PetscScalar p = .5;
    ierr = VecSet(x,p);CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Solve linear system
  */
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);

  /*
     View solver info; we could instead use the option -ksp_view to
     print this info to the screen at the conclusion of KSPSolve().
  */
  ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Check solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Check the error
  */
 
/*
 ierr = VecAXPY(x,neg_one,u);CHKERRQ(ierr);
  ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
  if (norm > tol) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %G, Iterations %D\n",norm,its);CHKERRQ(ierr);
  }
*/
  /*
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
  */
/*  ierr = VecDestroy(&x);CHKERRQ(ierr); ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr); ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
*/
  /*
     Always call PetscFinalize() before exiting a program.  This routine
       - finalizes the PETSc libraries as well as MPI
       - provides summary and diagnostic information if certain runtime
         options are chosen (e.g., -log_summary).
  */
  ierr = PetscFinalize();
}
