/* *********************************************************************************** */
int *adjacent_predecessors_cardinality(int *vertices, int m, int n);
/* *********************************************************************************** */
void dcinputfile_2_instance(dcinputfile *dc_vec, double ****instance_adjlist, proteinstructure **protein, int **kv, int m, int *n);
/* *********************************************************************************** */
double dij_from_discretizable_edges_set(double ***discretizableEdges, int i, int j, char l_or_u[]);
/* *********************************************************************************** */
double dij_from_adjlist(double ***adjlist, int *cardUi, int i, int j, char l_or_u[]);
/* *********************************************************************************** */
void adjlist_2_discretation_edges(double ***instance_adjlist, int *kv, int n, int **cliques, double ****discretationEdges);
/* *********************************************************************************** */
void adjlist_2_prune_edges(double ***instance_adjlist, int *kv, int n, int **cliques, prune_edges_set **pruneEdges);
/* *********************************************************************************** */
void adjlist_2_graph_parameters(double ***instance_adjlist, int *kv, int n, int **cliques, double ****discretationEdges, prune_edges_set **pruneEdges);
/* *********************************************************************************** */
void referential_x1_x2_x3_0(double ***instance_adjlist, int *kv, double ***X);
/* *********************************************************************************** */
void referential_x1_x2_x3(double ***X, double ***discretationEdges);
/* *********************************************************************************** */
void discretation_edges_2_ddgp_parameters(double ***discretationEdges, int **cliques, int n, ddgp_parameters **ddgpparameters);
/* *********************************************************************************** */
double dxixj_lf3(double *x1, double *x2);
/* *********************************************************************************** */
double d2xixj_lf3(double *x1, double *x2);
/* *********************************************************************************** */
double kappai_d2(double dii2_2, double di1i2, double dii1_2);
/* *********************************************************************************** */
double rhoi2_d2(double dii2_2, double ki);
/* *********************************************************************************** */
double round_double(double v, int n);
/* *********************************************************************************** */
double abs_torsion_angle_with_constants_d2(double p, double q, double dii3_2);
/* *********************************************************************************** */
double abs_torsion_angle_with_distances_d2(double d12_2, double d13_2, double d14_2, double d23, double d24_2, double d34_2);
/* *********************************************************************************** */
double torsion_angle_with_points_d2(double *x1, double *x2, double *x3, double *x4);
/* *********************************************************************************** */
int sign_torsion_angle(double *x1, double *x2, double *x3, double *x4);
/* *********************************************************************************** */
int sign_double(double v);
/* *********************************************************************************** */
void circunference_parameters(int i, double ***A, double ***B, double ***C, double *xi3, double *xi2, double *xi1, double kappa_210, double kappa_213, double rho2_210, double q_30, double di1i2);
/* *********************************************************************************** */
double changeReferential(double phase, double tau_ref1);
/* *********************************************************************************** */
void interIntervals(ddgp_interval *tau, ddgp_interval tau_k, double myZero);
/* *********************************************************************************** */
void cap2intervals(double *interA, double *interB, int *pos_interC, double ***interC, double myZero);
/* *********************************************************************************** */
void tree_backtracking(int *i, int **exploredVertex, double ***branches, int **branchNum, int sampleSize);
/* *********************************************************************************** */
void vertex_embedding(int i, double ***X, double *A, double *B, double *C, double tauAngle);
/* *********************************************************************************** */
void change_referential_precise_case(ddgp_interval *tau_k, double tauk_m, double phase);
/* *********************************************************************************** */
void change_referential_interval_case(ddgp_interval *tau_k, double tauk_l, double tauk_u, double phase, double myZero);
/* *********************************************************************************** */
void sampleInterval(int i, double ***branches, int **branchNum, ddgp_interval tau, int sampleSize, double tolerance);
/* *********************************************************************************** */
void sample_interval_0(double **J, int num_inter, int *vec_pos, double **sample, int sample_size, double tolerance);
/* *********************************************************************************** */
void sampleSimpleInterval(int i, double ***branches, int **branchNum, ddgp_interval tau, int sampleSize, double myZero);
/* *********************************************************************************** */
void tauii3_with_parameters_d2(ddgp_interval *tau, double p30, double q30, double dii3_2_l, double dii3_2_u, double myZero);
/* *********************************************************************************** */
void compute_tauii3(int i, double **X, double ***discretationEdges, double myZero, double **xi1, double **xi2, double **xi3, double *di1i2, double *kappa_210, double *rho2_210, double *kappa_213, double *rho2_213, double *q30, ddgp_interval *tau);
/* *********************************************************************************** */
void compute_tauiik_precise(int i, int k, double **X, prune_edges_set *pruneEdges, double kappa_210, double rho2_210, double kappa_213, double rho2_213, double di1i2, double *xi3, double *xi2, double *xi1, ddgp_interval *tau_k);
/* *********************************************************************************** */
void compute_tauiik_interval(int i, int k, double **X, prune_edges_set *pruneEdges, double kappa_210, double rho2_210, double kappa_213, double rho2_213, double di1i2, double *xi3, double *xi2, double *xi1, double myZero, ddgp_interval *tau_k, int *isTauIntervalEmpty);
/* *********************************************************************************** */
void dealloc_tau(ddgp_interval tau);
/* *********************************************************************************** */
void iABP(double ***X, int n, int *lev, double *cpu_time_used, double *numIt, double ***discretationEdges, prune_edges_set *pruneEdges, int sampleSize, double myZero, double timeMAX, double tolerance);
/* *********************************************************************************** */
void iBP(double ***X, int n, int *lev, double *cpu_time_used, double *numIt, double ***discretationEdges, prune_edges_set *pruneEdges, int sampleSize, double myZero, double timeMAX, double tolerance);
/* *********************************************************************************** */
void satisfiesInstanceQM(double **X, dcinputfile *dc_vec, int m, double myZero);
/* *********************************************************************************** */
double RMSD_calculation(double **X0, double **Xr, int n);

