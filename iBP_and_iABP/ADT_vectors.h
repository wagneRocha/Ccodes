/* *********************************************************************************** */
/* ---------------------------------- ADT int vectors -------------------------------- */
/* *********************************************************************************** */
int *alloc_vec_d(int n);
/* *********************************************************************************** */
int *zeros_vec_d(int n);
/* *********************************************************************************** */
int *ones_vec_d(int n);
/* *********************************************************************************** */
void print_vec_d(int *v, int n);
/* *********************************************************************************** */
void print_vec_d_t(int *v, int n);
/* *********************************************************************************** */
void vec_e_vec_d(int *v, int *u, int n);
/* *********************************************************************************** */
void vec_p_vec_d(int *w, int *u, int *v, int n);
/* *********************************************************************************** */
int sum_vec_d(int *v, int n);
/* *********************************************************************************** */
int *find_val_vec_d(int *v, int n, int val);
/* *********************************************************************************** */
int max_vec_d(int *v, int n);
/* *********************************************************************************** */
int *find_ones_position_vec_d(int *v, int n, int num1);
/* *********************************************************************************** */
int minAB_d(int a, int b);
/* *********************************************************************************** */
int maxAB_d(int a, int b);
/* *********************************************************************************** */
/* ---------------------------------- ADT double vectors ----------------------------- */
/* *********************************************************************************** */
double *alloc_vec_lf(int n);
/* *********************************************************************************** */
double *zeros_vec_lf(int n);
/* *********************************************************************************** */
double *ones_vec_lf(int n);
/* *********************************************************************************** */
void vec_e_vec_lf(double *v, double *u, int n);
/* *********************************************************************************** */
void vec_p_vec_lf(double *w, double *u, double *v, int n);
/* *********************************************************************************** */
void vec_m_vec_lf(double *w, double *u, double *v, int n);
/* *********************************************************************************** */
void l_t_vec_lf(double *v, double l, double *u, int n);
/* *********************************************************************************** */
void print_vec_lf(double *v, int n);
/* *********************************************************************************** */
double dot_product_lf(double *u, double *v, int n);
/* *********************************************************************************** */
double norm2_lf(double *v, int n);
/* *********************************************************************************** */
void vec_lf_zeros(double *v, int n);
/* *********************************************************************************** */
int *find_val_vec_lf(double *v, int n, double val);
/* *********************************************************************************** */
double sum_vec_lf(double *v, int n);
/* *********************************************************************************** */
double min_val_vec_lf(double *v, int n);
/* *********************************************************************************** */
double max_val_vec_lf(double *v, int n);
/* *********************************************************************************** */
double minAB_lf(double a, double b);
/* *********************************************************************************** */
double maxAB_lf(double a, double b);
/* *********************************************************************************** */
double dot_product_lf(double *u, double *v, int n);
/* *********************************************************************************** */
int *vec_lf_2_vec_d(double *vec_lf, int n);
/* *********************************************************************************** */
void swap_lf(double *a, double *b);
/* *********************************************************************************** */
void swap_d(int *a, int *b);
/* *********************************************************************************** */
int partition(double array[], int array_index[], int low, int high);
/* *********************************************************************************** */
void quicksort(double array[], int array_index[], int low, int high);
/* *********************************************************************************** */
/* ---------------------------------- ADT 3D double vectors -------------------------- */
/* *********************************************************************************** */
double *alloc_vec_lf3();
/* *********************************************************************************** */
double *zeros_vec_lf3();
/* *********************************************************************************** */
void vec_e_vec_lf3(double *v, double *u);
/* *********************************************************************************** */
void vec_p_vec_lf3(double *w, double *u, double *v);
/* *********************************************************************************** */
void vec_m_vec_lf3(double *w, double *u, double *v);
/* *********************************************************************************** */
void l_t_vec_lf3(double *v, double l, double *u);
/* *********************************************************************************** */
void print_vec_lf3(double *v);
/* *********************************************************************************** */
double dot_product_lf3(double *u, double *v);
/* *********************************************************************************** */
double norm2_lf3(double *v);
/* *********************************************************************************** */
void cross_product0_lf3(double *w, double *u, double *v);
/* *********************************************************************************** */
double *cross_product_lf3(double *u, double *v);
/* *********************************************************************************** */
void vec_lf3_zeros(double *v);
/* *********************************************************************************** */
/* ---------------------------------- ADT string vectors ----------------------------- */
/* *********************************************************************************** */
char *alloc_vec_s(int n);
/* *********************************************************************************** */
void blank_string(char *str, int L);
/* *********************************************************************************** */
void remove_spaces(char *input, char *output);
/* *********************************************************************************** */
/* -------------------------------- ADT dcfile vectors --------------------------- */
/* *********************************************************************************** */
typedef struct typedcInputFile{
	int i;
	int j;
	int ri;
	int rj;
	double dl;
	double du;
	char ai[5];
	char aj[5];
	char rti[4];
	char rtj[4];
} dcinputfile;
/* *********************************************************************************** */
dcinputfile *alloc_vec_dcInputFile(int n);
/* *********************************************************************************** */
/* ---------------------------------- ADT protein strucute --------------------------- */
/* *********************************************************************************** */
typedef struct typeProteinStructure{
	int v_k;
	int res_seq;
	char atom_name[5];
	char res_name[4];
	double x;
	double y;
	double z;
} proteinstructure;
/* *********************************************************************************** */
proteinstructure *alloc_vec_proteinstructure(int n);
/* *********************************************************************************** */
/* ---------------------------------- ADT edges set ---------------------------------- */
/* *********************************************************************************** */
typedef struct typePruneEdgesSet{
	int cardUkP;
	double **precise;
	int cardUkI;
	double **interval;
} prune_edges_set;
/* *********************************************************************************** */
prune_edges_set *alloc_vec_pruneedgesset(int n);
/* *********************************************************************************** */
/* ---------------------------------- ADT DDGP parameters ---------------------------- */
/* *********************************************************************************** */
typedef struct typeMatrixInterval{
	int is_point;
	int num_branches;
	int num_positive_intervals;
	int num_negative_intervals;
	double **positive_intervals_matrix;
	double **negative_intervals_matrix;
} matrix_intervals;
/* *********************************************************************************** */
typedef struct typeDDGPinterval{
	int num_pos_inter;
	int num_neg_inter;
	double **pos_inter;
	double **neg_inter;
} ddgp_interval;
/* *********************************************************************************** */
typedef struct typeDdgpParameters{
	double kappa;
	double kappabar;
	double rho2;
	double Rho;
	matrix_intervals tau;
} ddgp_parameters;
/* *********************************************************************************** */
ddgp_parameters *alloc_vec_ddgparameters(int n);
