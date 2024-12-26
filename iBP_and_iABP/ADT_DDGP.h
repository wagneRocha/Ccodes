typedef struct ddgp_solution_metrics{
        double npp;
        int lev;
        long nos;
        double cpu_time_used;
        stack *MDE;
        stack *LDE;
        stack *RMSD;
} solution_metrics;
/* *********************************************************************************** */
int *adjacent_predecessors_cardinality(int *vertices, int m, int n);
/* *********************************************************************************** */
void dcinputfile_2_instance(dcinputfile *dc_vec, double ****instance_adjlist, proteinstructure **protein, int **kv, int m, int *n);
/* *********************************************************************************** */
double dij_from_discretizable_edges_set(double ***discretizableEdges, int i, int j, char l_or_u[]);
/* *********************************************************************************** */
double dij_from_adjlist(double ***adjlist, int *cardUi, int i, int j, char l_or_u[]);
/* *********************************************************************************** */
void adjlist_2_discretization_edges_2(double ***instance_adjlist, int *kv, int n, int **cliques, double ****discretizationEdges_2);
/* *********************************************************************************** */
void adjlist_2_prune_edges_2(double ***instance_adjlist, int *kv, int n, int **cliques, prune_edges_set **pruneEdges_2);
/* *********************************************************************************** */
void adjlist_2_graph_parameters(double ***instance_adjlist, int *kv, int n, int **cliques, double ****discretizationEdges, prune_edges_set **pruneEdges);
/* *********************************************************************************** */
void referential_x1_x2_x3_0(double ***instance_adjlist, int *kv, double ***X);
/* *********************************************************************************** */
void referential_x1_x2_x3(double ***X, double ***discretizationEdges);
/* *********************************************************************************** */
void discretization_edges_2_ddgp_parameters(double ***discretizationEdges, int **cliques, int n, ddgp_parameters **ddgpparameters);
/* *********************************************************************************** */
double dxixj_lf3(double *x1, double *x2);
/* *********************************************************************************** */
double d2xixj_lf3(double *x1, double *x2);
/* *********************************************************************************** */
double kappai_d2(double dii2_2, double di1i2, double dii1_2);
/* *********************************************************************************** */
double rhoi2_d2(double dii2_2, double ki);
/* *********************************************************************************** */
double round_lf(double val, int n);
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
double change_referential_1_to_0(double phase, double tau_ref1);
/* *********************************************************************************** */
void cap_interAk(ddgp_interval *tau, ddgp_interval tau_k);
/* *********************************************************************************** */
void interA_cap_interB(double *interA, double *interB, double **interC, int *isCapEmpty);
/* *********************************************************************************** */
void interA_cap_interB_matrix(double *interA, double *interB, double ***interC_matrix, int *matrix_row);
/* *********************************************************************************** */
void tree_backtracking(int *i, int **exploredVertex, double ***branches, int **branchNum, int sampleSize);
/* *********************************************************************************** */
void vertex_embedding(int i, double ***X, double *A, double *B, double *C, double tauAngle);
/* *********************************************************************************** */
void change_referential_precise_case(ddgp_interval *tau_k, double tauk_m, double phase, double resolution);
/* *********************************************************************************** */
void change_referential_interval_case(ddgp_interval *tau_k, double tauk_l, double tauk_u, double phase);
/* *********************************************************************************** */
void sample_interval_union(double **J, int numInter, double **sample, int *sampleSize, int maxSampleSize, double resolution);
/* *********************************************************************************** */
void sample_interval_union_DDGP(int i, double ***branches, int **branchNum, ddgp_interval tau, int sampleSize, double resolution);
/* *********************************************************************************** */
void sample_simple_DDGPinterval(int i, double ***branches, int **branchNum, ddgp_interval tau, int sampleSize, double resolution);
/* *********************************************************************************** */
void sample_simple_interval(double **sample, double *interval, int sampleSize, double tolerance, int *numPoints);
/* *********************************************************************************** */
void sort_sample_indices(double **sortedSample, double *sample, int sampleSize);
/* *********************************************************************************** */
void tauii3_with_parameters_d2(ddgp_interval *tau, double p30, double q30, double dii3_2_l, double dii3_2_u);
/* *********************************************************************************** */
void compute_tauii3_parameters(int i, double **X, double ***discretizationEdges, int *i1, int *i2, int *i3, double *di1i2, double *kappa_210, double *rho2_210, double *kappa_213, double *rho2_213, double *q30);
/* *********************************************************************************** */
void compute_tauii3(int i, double **X, double ***discretizationEdges, int *i1, int *i2, int *i3, double *di1i2, double *kappa_210, double *rho2_210, double *kappa_213, double *rho2_213, double *q30, ddgp_interval *tau);
/* *********************************************************************************** */
void compute_tauiik_precise(int i, int k, double **X, prune_edges_set *pruneEdges, double kappa_210, double rho2_210, double kappa_213, double rho2_213, double di1i2, int i3, int i2, int i1, ddgp_interval *tau_k, double resolution);
/* *********************************************************************************** */
void compute_tauiik_interval(int i, int k, double **X, prune_edges_set *pruneEdges, double kappa_210, double rho2_210, double kappa_213, double rho2_213, double di1i2, int i3, int i2, int i1, ddgp_interval *tau_k, int *isTauIntervalEmpty);
/* *********************************************************************************** */
void dealloc_tau(ddgp_interval tau);
/* *********************************************************************************** */
void get_torsion_angle(ddgp_interval *tau, double tauAst, double deltaTau);
/* *********************************************************************************** */
void satisfiesInstanceQM(double **X, dcinputfile *dc_vec, int m, double resolution);
/* *********************************************************************************** */
void compute_mde_and_lde(double **X, dcinputfile *dc_vec, int m, double *mde, double *lde);
/* *********************************************************************************** */
double compute_rmsd(double **X, double **Y, int n);
/* *********************************************************************************** */
void iABP(int n, double ***discretizationEdges_2, int *tauSign, double *givenTau, double *givenTauDeviation, prune_edges_set *pruneEdges_2, int sampleSize, double resolution, double tolerance, double timeLimit, proteinstructure **protein, int GivenNumOfSols, solution_metrics *solutionMetrics, dcinputfile *dc_vec, int num_dc, int is_X0_known, double **X0);
/* *********************************************************************************** */
void iBP(int n, double ***discretizationEdges_2, prune_edges_set *pruneEdges_2, int sampleSize, double tolerance, double timeLimit, proteinstructure **protein, int GivenNumOfSols, solution_metrics *solutionMetrics, dcinputfile *dc_vec, int num_dc, int is_X0_known, double **X0);
/* *********************************************************************************** */
void copy_protein_structure(proteinstructure **p1, proteinstructure *p0, int n);
