#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <fstream>
using namespace std;

/*double x0 = 0., xN = 1.;
double y10 = 0., yM = 1.;
const int N = 10;
const int M = 10;
double T0 = 0., Tn = 1.;
const int n = 1;//разбиение по времени
double tau = (Tn - T0) / n;
const double p0 = 1.;
double z11 = 1., z12 = 1., z21 = 1., z22 = 1.;//мы изначально рассматриваем единичный квадрат, поэтому параметризаци€ тождественна€
double Re = 100;//число –ейнольдса
double mu = 1 / Re;
int temp = 1;//температура
double g = 9.81;
double eps = 0.01;
double alpha = 0.6;*/


double lambda11(int N, int M, double **v, double h, int i, int j)
{
	double sum = 0;
	sum = ((*v)[i + 1 + j*(N + 2)] - 2 * (*v)[i + j*(N + 2)] + (*v)[i - 1 + j*(N + 2)]) / (h*h);
	return sum;
}

double lambda22(int N, int M, double **v, double h, int i, int j)
{
	double sum = 0;
	sum = ((*v)[i + (j + 1)*(N + 2)] - 2 * (*v)[i + j*(N + 2)] + (*v)[i + (j - 1)*(N + 2)]) / (h*h);
	return sum;
}

void initial_ksi(int N, int M, vector<vector<double>>& v1, vector<vector<double>>& v2, double **ksi, double **psi, double tau, double alpha)
{
	*ksi = (double*)malloc(sizeof(double)* (M + 2) * (N + 2));
	*psi = (double*)malloc(sizeof(double)* (M + 2) * (N + 2));
	double   h1 = 2.0 / N;
	double   h2 = 2.0 / M;
	double v_i, v_j;

	for (int j = 0; j <= M + 1; ++j)
	{
		for (int i = 0; i <= N + 1; ++i){
			(*psi)[i + j*(N + 2)] = 0;
		}
	}


	for (int j = 1; j <= M; ++j)
	{
		for (int i = 1; i <= N; ++i)
		{
			v_i = (v2[i + 1][j] - v2[i - 1][j]) / (2 * h1);
			v_j = (v1[i][j + 1] - v1[i][j - 1]) / (2 * h2);
			(*ksi)[i + j*(N + 2)] = lambda11(N, M, psi, h1, i, j) + lambda22(N, M, psi, h2, i, j) + v_i + v_j;
		}
	}

	for (int j = 0; j <= M + 1; j++)
	{
		int i = 0; //лева€ граница
		(*ksi)[i + j* (N + 2)] = (*ksi)[i + 1 + j* (N + 2)];
		i = N + 1;//права€ граница
		(*ksi)[i + j* (N + 2)] = (*ksi)[i - 1 + j* (N + 2)];
	}

	for (int i = 0; i <= N + 1; i++)
	{
		int j = 0;//нижн€€ граница
		(*ksi)[i + j*(N + 2)] = (*ksi)[i + (j + 1) * (N + 2)];

		j = M + 1;//верхн€€ граница
		(*ksi)[i + j* (N + 2)] = (*ksi)[i + (j - 1) * (N + 2)];
	}
}

void psi_newlay_1(int N, int M, double **ksi, double **ksi_1_2, double tau, double alpha)
{
	*ksi_1_2 = (double*)malloc(sizeof(double)* (M + 2) * (N + 2));
	double   h1 = double(2.0 / N);
	double   h2 = double(2.0 / M);
	double *c = (double*)malloc(sizeof(double)* (N + 1));
	double *b = (double*)malloc(sizeof(double)* (N + 1));
	double *a = (double*)malloc(sizeof(double)* (N + 1));
	double *beta = (double*)malloc(sizeof(double)* (N + 1));
	double *alph = (double*)malloc(sizeof(double)* (N + 1));
	double *func = (double*)malloc(sizeof(double)* (N + 1));


	for (int j = 1; j <= M; ++j) {
		for (int i = 1; i <= N; ++i){

			alph[0] = -1;
			beta[0] = 0;
			for (int i = 0; i <N; i++) {

				c[i] = tau*alpha / (h1*h1);
				//std::cout << tau << " " << alpha << " " << h2 << "\n";
				b[i] = -1 - 2 * tau*alpha / (h1*h1);
				a[i] = tau*alpha / (h1*h1);
				func[i] = (*ksi)[i + j*(N + 2)];
				alph[i + 1] = -c[i] / (a[i] * alph[i] + b[i]);
				beta[i + 1] = (func[i] - a[i] * beta[i]) / (a[i] * alph[i] + b[i]);
				//std::cout << (*ksi)[i + j*(N + 2)] << " " << a[i] << " " << b[i] << " " << c[i] << " " << func[i] << " " << alph[i + 1] << " " << beta[i + 1] << "\n";
			}

			(*ksi_1_2)[N + 1 + j*(N + 2)] = beta[N] / (alph[N] + 1);
			//std::cout << "\n" << alph[N] << " " << beta[N] << "\n";
			for (int i = N; i >= 0; i--) {
				(*ksi_1_2)[i + j*(N + 2)] = alph[i] * (*ksi_1_2)[i + 1 + j*(N + 2)] + beta[i];
				//std::cout << (*ksi_1_2)[i + (j)* (N + 2)] << "  " << i <<"  " << j << "\n\n";
			}
		}
	}

	for (int i = 0; i <= N + 1; i++){
		int j = 0;//нижн€€ граница

		(*ksi_1_2)[i + j* (N + 2)] = -(*ksi_1_2)[i + (j + 1) * (N + 2)];

		j = M + 1;//верхн€€ граница
		(*ksi_1_2)[i + j* (N + 2)] = -(*ksi_1_2)[i + (j - 1) * (N + 2)];
	}
}

void psi_newlay_2(int N, int M, double **ksi_1, double **ksi_1_2, double tau, double alpha){
	*ksi_1 = (double*)malloc(sizeof(double)* (M + 2) * (N + 2));
	double   h1 = double(2.0 / N);
	double   h2 = double(2.0 / M);

	double *c = (double*)malloc(sizeof(double)* (M + 1));
	double *b = (double*)malloc(sizeof(double)* (M + 1));
	double *a = (double*)malloc(sizeof(double)* (M + 1));
	double *beta = (double*)malloc(sizeof(double)* (M + 1));
	double *alph = (double*)malloc(sizeof(double)* (M + 1));
	double *func = (double*)malloc(sizeof(double)* (M + 1));



	for (int i = 1; i <= N; ++i){
		for (int j = 1; j <= M; ++j) {

			alph[0] = -1;
			beta[0] = 0;

			for (int i = 0; i <N; i++) {


				c[i] = tau*alpha / (h2*h2);
				//std::cout << tau << " " << alpha << " " << h2 << "\n";
				b[i] = -1 - 2 * tau*alpha / (h2*h2);
				a[i] = tau*alpha / (h2*h2);
				func[i] = (*ksi_1_2)[i + j*(N + 2)];
				alph[i + 1] = -c[i] / (a[i] * alph[i] + b[i]);
				beta[i + 1] = (func[i] - a[i] * beta[i]) / (a[i] * alph[i] + b[i]);
				//cout << (*ksi_1_2)[i + j*(N + 2)] << " " << a[i] << " " << b[i] << " " << c[i] << " " << func[i] << " " << alph[i + 1] << " " << beta[i + 1] << "\n";
			}

			(*ksi_1)[N + 1 + j*(N + 2)] = beta[N] / (alph[N] + 1);
			//std::cout << "\n" << alph[N] << " " << beta[N] << "\n";
			for (int i = N; i >= 0; i--) {
				(*ksi_1)[i + j*(N + 2)] = alph[i] * (*ksi_1)[i + 1 + j*(N + 2)] + beta[i];
				//std::cout << (*ksi_1_2)[i + (j)* (N + 2)] << "  " << i <<"  " << j << "\n\n";
			}
		}
	}

	for (int j = 0; j <= M + 1; j++){
		int i = 0;
		(*ksi_1)[i + j* (N + 2)] = (*ksi_1)[i + 1 + (j)* (N + 2)];

		i = N + 1;
		(*ksi_1)[i + j* (N + 2)] = (*ksi_1)[i - 1 + (j)* (N + 2)];
	}
}


void psi_final(int N, int M, double **psi_1, double **psi, double **ksi_1, double tau){
	*psi_1 = (double*)malloc(sizeof(double)* (M + 2) * (N + 2));
	for (int i = 0; i <= N + 1; ++i){
		for (int j = 0; j <= M + 1; ++j) {
			(*psi_1)[i + j* (N + 2)] = (*psi)[i + j* (N + 2)] + tau*(*ksi_1)[i + j* (N + 2)];
		}
	}
}