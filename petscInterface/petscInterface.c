#include"petscInterface.h"

void petscSolve(double* matrix,int N,double* rhs, double* x1){
	Vec		 x, b;      /* approx solution, RHS, exact solution */
	Mat            A;            /* linear system matrix */
	KSP            ksp;         /* linear solver context */
	PC             pc;           /* preconditioner context */
	PetscErrorCode ierr;
	PetscInt       n = N;

	int i;

	PetscInitialize(NULL,NULL,(char*)0,NULL);

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Compute the matrix and right-hand-side vector that define
	   the linear system, Ax = b.
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	ierr = VecCreate(PETSC_COMM_WORLD,&x);
	ierr = PetscObjectSetName((PetscObject) x, "Solution");
	ierr = VecSetSizes(x,PETSC_DECIDE,n);
	ierr = VecSetFromOptions(x);
	ierr = VecDuplicate(x,&b);


	PetscInt* linesv=malloc(N*sizeof(PetscInt));

	for(i=0;i<N;i++){
		linesv[i]=i;
	}
	ierr = VecSetValues(b,n,linesv,rhs, INSERT_VALUES);//	Set vector b from the input function(RHS)  

	/*
	   Create matrix.  When using MatCreate(), the matrix format can
	   be specified at runtime.

	   Performance tuning note:  For problems of substantial size,
	   preallocation of matrix memory is crucial for attaining good
	   performance. See the matrix chapter of the users manual for details.
	   */

	ierr = MatCreate(PETSC_COMM_WORLD,&A);
	ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);
	ierr = MatSetFromOptions(A);
	ierr = MatSetUp(A);

	/*
	   Assemble matrix
	   */

	int* idx=malloc(N*sizeof(int));

	//Build the array to configure lines and coloums that the value of matrix will be added

	for(i=0;i<N;i++){
		idx[i]=i;
	}


	ierr = MatSetValues(A,N,idx,N,idx,matrix,INSERT_VALUES);
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);



	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Create the linear solver and set various options
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/*
	   Create linear solver context
	   */
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);

	/*
	   Set operators. Here the matrix that defines the linear system
	   also serves as the preconditioning matrix.
	   */
	ierr = KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);

	/*
	   Set linear solver defaults for this problem (optional).
	   - By extracting the KSP and PC contexts from the KSP context,
	   we can then directly call any KSP and PC routines to set
	   various options.
	   - The following four statements are optional; all of these
	   parameters could alternatively be specified at runtime via
	   KSPSetFromOptions();
	   */
	ierr = KSPGetPC(ksp,&pc);
	ierr = PCSetType(pc,PCJACOBI);
	ierr = KSPSetTolerances(ksp,1.e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

	ierr = KSPSetFromOptions(ksp);

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Solve the linear system
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	ierr = KSPSolve(ksp,b,x);

	PetscScalar *temp;

	ierr = VecGetArray(x,&temp);

	for(i=0;i<N;i++)
		x1[i]=temp[i];


	ierr = PetscFinalize();
}
