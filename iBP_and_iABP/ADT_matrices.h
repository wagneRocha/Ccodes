/* *********************************************************************************** */
/* ---------------------------------- ADT double matrices ---------------------------- */
/* *********************************************************************************** */
int **alloc_mat_d(int m, int n);
/* *********************************************************************************** */
void dealloc_mat_d(int **M, int m);
/* *********************************************************************************** */
void print_mat_d(int **M, int m, int n);
/* *********************************************************************************** */
double **alloc_mat_lf(int m, int n);
/* *********************************************************************************** */
void dealloc_mat_lf(double **M, int m);
/* *********************************************************************************** */
void dealloc_mat_lf_2(double **M, int m);
/* *********************************************************************************** */
void print_mat_lf(double **M, int m, int n);
/* *********************************************************************************** */
double **zeros_mat_lf(int m, int n);
/* *********************************************************************************** */
double **ones_mat_lf(int m, int n);
/* *********************************************************************************** */
void col_of_mat_lf(double **M, int m, int col, double *v);
/* *********************************************************************************** */
double ***alloc_3Darray_lf(int n, int *vec);
/* *********************************************************************************** */
void dealloc_3Darray_lf(double ***T, int *vec_kv, int n);
/* *********************************************************************************** */
double ***alloc_fix3Darray_lf(int l, int m, int n);
/* *********************************************************************************** */
void dealloc_fix3Darray_lf(double ***T, int l, int m);
/* *********************************************************************************** */
double *get_column_mat_lf(double **M, int m, int c);
/* *********************************************************************************** */
double **get_columns_mat_lf(double **M, int m, int *cols, int ncols);
/* *********************************************************************************** */
int **mat_lf_2_mat_d(double **M, int m, int n);
/* *********************************************************************************** */
void l_t_mat_lf(double l, double **M, int m, int n);
/* *********************************************************************************** */
void mat_p_mat_lf(double **M, double **M1, double **M0, int m, int n);
/* *********************************************************************************** */
void mat_m_mat_lf(double **M, double **M1, double **M0, int m, int n);
/* *********************************************************************************** */
void mat_e_mat_lf(double **M1, double **M0, int m, int n);
/* *********************************************************************************** */
double **trans_mat_lf(double **M, int m, int n);
/* *********************************************************************************** */
double **mat_times_mat_lf(double **M1, double **M2, int m, int p, int n);
/* *********************************************************************************** */
double frobenius_norm(double **A, int m, int n);
/* *********************************************************************************** */
void mat_2_vec_lf(double **M, int m, int n, double *v);
/* *********************************************************************************** */
void vec_2_mat_lf(double **M, int m, int n, double *v);
/* *********************************************************************************** */
void vec_2_diagmat_lf(double **M, double *v, int p);
/* *********************************************************************************** */
void svd_row_major_matrix(double **M, int m0, int n0, double **U, double **S, double **Vt);
/* *********************************************************************************** */
/* ---------------------------------- ADT string matrices ---------------------------- */
/* *********************************************************************************** */
char **alloc_mat_s(int m, int n);
/* *********************************************************************************** */
char ***alloc_3Darray_s(int m, int n, int l);
/* *********************************************************************************** */
void dealloc_mat_s(char **T, int m);
/* *********************************************************************************** */
void dealloc_3Darray_s(char ***T, int m, int n);
/* *********************************************************************************** */
void print_mat_s(char ***T, int m, int n);
