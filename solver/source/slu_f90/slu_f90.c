#include "slu_ddefs.h"

void slu_solve(int n, int nnz, int isrow, int shiftindex, double val[], int colind[], int rowptr[], double rhs[], int *ierr) {
  SuperMatrix A, b, L, U;
  superlu_options_t options;
  SuperLUStat_t stat;
  int *perm_c = NULL, *perm_r = NULL;
  double *valcol = NULL;
  int *colptr = NULL, *rowind = NULL;

  if (shiftindex) {
    int i;
    for (i = 0; i <= n; i++) rowptr[i] -= 1;
    for (i = 0; i < nnz; i++) colind[i] -= 1;
  }

  if (isrow) {
    dCompRow_to_CompCol(n, n, nnz, val, colind, rowptr, &valcol, &rowind, &colptr);
    dCreate_CompCol_Matrix(&A, n, n, nnz, valcol, rowind, colptr, SLU_NC, SLU_D, SLU_GE);
  }
  else {
    dCreate_CompCol_Matrix(&A, n, n, nnz, val, colind, rowptr, SLU_NC, SLU_D, SLU_GE);
  }
  dCreate_Dense_Matrix(&b, n, 1, rhs, n, SLU_DN, SLU_D, SLU_GE);

  set_default_options( &options );
  StatInit( &stat );

  if (!(perm_c = intMalloc(n))) {
    *ierr = -1;
    goto error;
  }
  if (!(perm_r = intMalloc(n))) {
    *ierr = -1;
    goto error;
  }


  dgssv(&options, &A, perm_c, perm_r, &L, &U, &b, &stat, ierr);

error:
  StatFree( &stat );

  SUPERLU_FREE(perm_c);
  SUPERLU_FREE(perm_r);
  SUPERLU_FREE(valcol);
  SUPERLU_FREE(colptr);
  SUPERLU_FREE(rowind);

  Destroy_SuperNode_Matrix( &L );
  Destroy_CompCol_Matrix( &U );
  Destroy_SuperMatrix_Store(&A);
  Destroy_SuperMatrix_Store(&b);
}
