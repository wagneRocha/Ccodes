/* *********************************************************************************** */
/* ---------------------------------- ADT files -------------------------------------- */
/* *********************************************************************************** */
int num_lines_file(char filename[]);
/* *********************************************************************************** */
int num_columns_file_lines(char filename[], char sep_char);
/* *********************************************************************************** */
void mat_lf_2_file(double **M, int m, int n, char filename[], char sep_char);
/* *********************************************************************************** */
double **file_2_mat_lf(char filename[], char sep_char);
/* *********************************************************************************** */
int **file_2_mat_d(char filename[], char sep_char);
/* *********************************************************************************** */
void vec_lf_2_file(double *v, int m, char filename[]);
/* *********************************************************************************** */
char ***txtfile_2_mat_s(char filename[], char sep_char);
/* *********************************************************************************** */
dcinputfile *dcInputFile_2_dcDataVector(char *fname);
/* *********************************************************************************** */
void read_input_file(char *fname, char **structure_id, char **structure_chain, char **method, char **fname_dc, char **fname_T0, char **fname_X0, double *timeLimit, double *tolerance, double *angularResolution, int *numOfSols, int *sampleSize);
/* *********************************************************************************** */
void save_protein_PDBformat(char method[], proteinstructure *protein, int n, char structure_id[], char model[], char structure_chain[], char fname[]);
