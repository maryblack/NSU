void psi_newlay_1(int N, int M, double **ksi, double **ksi_1_2, double tau, double alpha);

void psi_newlay_2(int N, int M, double **ksi_1, double **ksi_1_2, double tau, double alpha);

void psi_final(int N, int M, double **psi_1, double **psi, double **ksi_1, double tau);

void initial_ksi(int N, int M, vector<vector<double>>& v1, vector<vector<double>>& v2, double **ksi, double **psi, double tau, double alpha);

double lambda11(int N, int M, double **v, double h, int i, int j);

double lambda22(int N, int M, double **v, double h, int i, int j);