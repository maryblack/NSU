#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <fstream>
using namespace std;

#include "isoline.h"





double x0 = 0., xN = 2.;
double y10 = 0., yM = 2.;
const int N = 20;
const int M = 20;
double T0 = 0., Tn = 1.;
const int n = 10000;//разбиение по времени
double tau = (Tn - T0) / n;
const double p0 = 1.;
//double z11 = N, z12 = 0., z21 = 0., z22 = M;//мы изначально рассматриваем единичный квадрат, поэтому параметризаци€ тождественна€
double z11 = N/2., z12 = 0., z21 = 0., z22 = M;
double Re = 100;//число –ейнольдса
double mu = 1 / Re;
int temp = 1;//температура
double g = 9.81;
double eps = 0.01;
double alpha = 1.;





typedef struct box{//информаци€ о €чейке
	double w;//площадь €чейки
}box;

typedef struct el{
	double x;//координата по х
	double y;//координата по у
}el;

double dist(double x1, double y1, double x2, double y2){
	return sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2));
}

double square(int N, int i, int j, el **V){//площадь i,j-ой €чейки
	double a = 0, b = 0, c = 0, d = 0, p = 0;
	a = dist((*V)[i + j*(N + 1)].x, (*V)[i + j*(N + 1)].y, (*V)[i + (j + 1)*(N + 1)].x, (*V)[i + (j + 1)*(N + 1)].y);
	b = dist((*V)[i + (j + 1)*(N + 1)].x, (*V)[i + (j + 1)*(N + 1)].y, (*V)[i + 1 + (j + 1)*(N + 1)].x, (*V)[i + 1 + (j + 1)*(N + 1)].y);
	c = dist((*V)[i + 1 + (j + 1)*(N + 1)].x, (*V)[i + 1 + (j + 1)*(N + 1)].y, (*V)[i + 1 + j*(N + 1)].x, (*V)[i + 1 + j*(N + 1)].y);
	d = dist((*V)[i + 1 + j*(N + 1)].x, (*V)[i + 1 + j*(N + 1)].y, (*V)[i + j*(N + 1)].x, (*V)[i + j*(N + 1)].y);
	p = (double)(a + b + c + d) / 2;
	return sqrt((p - a)*(p - b)*(p - c)*(p - d));
}

void fragmentation(int N, int M, el **V, box **B){
	*V = (el*)malloc(sizeof(el)* (M + 1) * (N + 1));// в V хран€тс€ координаты вершин
	*B = (box*)malloc(sizeof(box)* (M)* (N));
	double h1 = (xN - x0) / N;
	double h2 = (yM - y10) / M;
	el k;
	double xh = 0, yh = 0;
	for (int j = 0; j <= M; j++){//массив с координатами вершин
		yh = y10 + h2*j;
		k.y = yh;//дл€ равномерной сетки
		for (int i = 0; i <= N; i++){

			xh = x0 + h1*i;
			k.x = xh;//дл€ равномерной сетки
			(*V)[i + j*(N + 1)] = k;
		}
	}

	for (int j = 0; j <M; j++){//массив с информацией о €чейках
		for (int i = 0; i <N; i++){
			(*B)[i + j*N].w = square(N, i, j, V);//площадь 
		}
	}

}


void projection(int N, int M, double **S1i, double **S2i, double **S1j, double **S2j, el **V){
	*S1i = (double*)malloc(sizeof(double)* (M)* (N + 1));//
	*S2i = (double*)malloc(sizeof(double)* (M)* (N + 1));//
	*S1j = (double*)malloc(sizeof(double)* (M + 1) * (N));//проекции горизонтальной грани по y
	*S2j = (double*)malloc(sizeof(double)* (M + 1) * (N));//проекции вертикальной грани по у

	for (int j = 0; j <= M; j++){//горизонтальные
		for (int i = 0; i < N; i++){
			(*S1j)[i + j*(N)] = abs((*V)[i + j*(N + 1)].y - (*V)[i + 1 + (j)*(N + 1)].y);
			(*S2j)[i + j*(N)] = abs((*V)[i + j*(N + 1)].x - (*V)[i + 1 + (j)*(N + 1)].x);

		}
	}


	for (int j = 0; j <M; j++){//вертикальные
		for (int i = 0; i <= N; i++){
			(*S1i)[i + j*(N + 1)] = abs((*V)[i + j*(N + 1)].y - (*V)[i + (j + 1)*(N + 1)].y);
			(*S2i)[i + (j)*(N + 1)] = abs((*V)[i + j*(N + 1)].x - (*V)[i + (j + 1)*(N + 1)].x);
		}
	}
}


void initial(int N, int M, vector<vector<double>>& p, vector<vector<double>>& v1, vector<vector<double>>& v2, el **V){//начальные услови€

	double y;
	double x;
	for (int j = 1; j <= M; j++){
		for (int i = 0; i <= N; i++){
			v2[i][j] = 0;
			//v1[i][j] = 0;
			//p[i][j] = 1;
			x = (1 / 10.)*(i + 1 / 2.);
			p[i][j] = 1-2*mu*x;
			y = (1 / 10.)*(j - 1 / 2.);
			v1[i][j] = y*(2 - y);

		}
	}


	for (int j = 0; j <= M + 1; j++){
		int i = 0;
		//p[i][j] = p[i + 1][j];
		p[i][j] = 1 - 2 * mu*(1 / 20.);

		v2[i][j] = v2[i + 1][j];
		y = (1 / 10.)*(j - 1 / 2.);
		v1[i][j] = y*(2-y);

		i = N + 1;
		p[i][j] = 1-2*mu*(2+2/N);

		v2[i][j] = v2[i - 1][j];
		v1[i][j] = v1[i - 1][j];

	}



	for (int i = 0; i <= N + 1; i++){
		int j = 0;

		p[i][j] = p[i][j + 1];
		v2[i][j] = -v2[i][j + 1];
		v1[i][j] = -v1[i][j + 1];

		j = M + 1;
		v1[i][j] = - v1[i][j - 1];
		p[i][j] = p[i][j - 1];
		v2[i][j] = -v2[i][j - 1];

	}
	v1[0][0] = -v1[0][1];//(0,0)
	v1[0][M + 1] = -v1[0][M];
	/*v2[0][0] = 2 * v2[1][0] - v2[2][0];
	v1[N + 1][0] = 2 * v1[N][0] - v1[N - 1][0];//(N+1,0)
	v2[N + 1][0] = 2 * v2[N][0] - v2[N - 1][0];
	v1[0][M + 1] =0;//(0,M+1)
	v2[0][M + 1] = 2 * v2[1][M + 1] - v2[2][M + 1];
	v1[N + 1][M + 1] = 2 * v1[N][M + 1] - v1[N - 1][M + 1];//(N+1,M+1)
	v2[N + 1][M + 1] = 2 * v2[N][M + 1] - v2[N - 1][M + 1];*/

}

void predict(int N, int M, vector<vector<double>>& p0, vector<vector<double>>& v10, vector<vector<double>>& v20, vector<vector<double>>& p1, vector<vector<double>>& v11, vector<vector<double>>& v21, box **B, double **S1i, double **S2i, double **S1j, double **S2j, el **V){
	vector<vector<double>> T1(N + 1, vector<double>(M + 1, 0));
	vector<vector<double>> T2(N + 1, vector<double>(M + 1, 0));

	/*vector<vector<double>> sigma11(N + 2, vector<double>(M + 2, 0));
	vector<vector<double>> sigma12(N + 2, vector<double>(M + 2, 0));
	vector<vector<double>> sigma21(N + 2, vector<double>(M + 2, 0));
	vector<vector<double>> sigma22(N + 2, vector<double>(M + 2, 0));

	int i = 0, j = 0;
	double w = 0, sum1 = 0, sum2 = 0, sum3 = 0;
	double d = 0;
	for (j = 1; j <= M; j++){
	for (i = 1; i <= N; i++){

	sigma11[i][j] = mu*z11*(v10[i + 1][j] - v10[i - 1][j]);
	//cout << sigma11[i][j] << "\n";
	sigma22[i][j] = mu*z22*(v20[i][j+1] - v20[i][j-1]);
	//cout << v10[i + 1][j] << " " << v10[i - 1][j] << "\n";
	//cout << mu*z11*(v10[i + 1][j] - v10[i - 1][j]) << "*\n\n";
	//std::cout << "*" << (*v20)[i + (j + 1)*(N + 2)] << " " << (*v20)[i + (j - 1)*(N + 2)] << "**\n";

	sigma21[i][j] = mu*(z22*(v10[i][j + 1] - v10[i][j - 1]) + z11*(v20[i + 1][j] - v20[i - 1][j])) / 2;
	sigma12[i][j] = sigma21[i][j];
	//cout <<"*"<<i<<" "<<j<<" "<< sigma11[i][j]<< "  " << sigma21[i][j] << "  " << sigma22[i][j] << "*\n";

	}
	}
	for (int j = 1; j <= M; j++){
	i = 0;
	sigma11[i][j] = sigma11[i + 1][j];
	sigma22[i][j] = sigma22[i + 1][j];
	sigma21[i][j] = sigma21[i + 1][j];
	sigma12[i][j] = sigma21[i][j];
	//cout << "*" << i << " " << j << " " << sigma11[i][j] << "  " << sigma21[i][j] << "  " << sigma22[i][j]<<"\n";
	i = N + 1;
	sigma11[i][j] = sigma11[i - 1][j];
	sigma22[i][j] = sigma22[i - 1][j];
	sigma21[i][j] = sigma21[i - 1][j];
	sigma12[i][j] = sigma21[i][j];
	//cout << "*" << i << " " << j << " " << sigma11[i][j] << "  " << sigma21[i][j] << "  " << sigma22[i][j] << "\n";
	}
	for (int i = 1; i <= N; i++){
	j = 0;
	sigma11[i][j] = sigma11[i][j+1];
	sigma22[i][j] = sigma22[i][j+1];
	sigma21[i][j] = sigma21[i][j+1];
	sigma12[i][j] = sigma21[i][j];
	//cout << "*" << i << " " << j << " " << sigma11[i][j] << "  " << sigma21[i][j] << "  " << sigma22[i][j] << "\n";
	j = M + 1;
	sigma11[i][j] = -sigma11[i][j - 1];
	sigma22[i][j] = sigma22[i][j - 1];
	sigma21[i][j] = sigma21[i][j - 1];
	sigma12[i][j] = sigma21[i][j];
	//cout << "*" << i << " " << j << " " << sigma11[i][j] << "  " << sigma21[i][j] << "  " << sigma22[i][j] << "\n";
	}
	*/
	double v1ij = 0, v1i1j = 0, v1i_1j = 0, v2ij = 0, v2i1j = 0, v2i_1j = 0, v1ij1 = 0, v1ij_1 = 0, v2ij1 = 0, v2ij_1 = 0;
	double v1i1j1 = 0, v1i1j_1 = 0, v1i_1j1 = 0, v1i_1j_1 = 0, v2i1j1 = 0, v2i1j_1 = 0, v2i_1j1 = 0, v2i_1j_1 = 0;
	int i = 0, j = 0;
	double w = 0, sum1 = 0, sum2 = 0, sum3 = 0;
	double d = 0;


	for (j = 1; j <= M; j++){
		for (i = 1; i <= N; i++){
			w = (*B)[i - 1 + (j - 1)*N].w;
			//	w = 1 / (N*M);
			d = tau*alpha / w;

			//cout << i << " " << j << "\n";
			v1ij = v10[i][j];
			v1i1j = v10[i + 1][j];
			v1i_1j = v10[i - 1][j];

			v2ij = v20[i][j];
			v2i1j = v20[i + 1][j];
			v2i_1j = v20[i - 1][j];

			v1ij1 = v10[i][j + 1];
			v1ij_1 = v10[i][j - 1];

			v2ij1 = v20[i][j + 1];
			v2ij_1 = v20[i][j - 1];

			v1i1j1 = v10[i + 1][j + 1];
			v1i1j_1 = v10[i + 1][j - 1];
			v1i_1j1 = v10[i - 1][j + 1];
			v1i_1j_1 = v10[i - 1][j - 1];


			v2i1j1 = v20[i + 1][j + 1];
			v2i1j_1 = v20[i + 1][j - 1];
			v2i_1j1 = v20[i - 1][j + 1];
			v2i_1j_1 = v20[i - 1][j - 1];

			/*std::cout << v1ij << " " << v1i1j << " " << v1i_1j << "\n";
			std::cout << v2ij << " " << v2i1j << " " << v2i_1j << "\n";
			std::cout << v1ij1 << " " << v1ij_1 << "\n";
			std::cout << v2ij1 << " " << v2ij_1 << "\n";
			std::cout << v1i1j1 << " " << v1i1j_1 << " " << v1i_1j1 << " " << v1i_1j_1 << "\n";
			std::cout << v2i1j1 << " " << v2i1j_1 << " " << v2i_1j1 << " " << v2i_1j_1 << "\n\n";*/

			//std::cout << "d=" << d << "\n";
			/*T1[i][j] = 2 * mu*z11*(*S1i)[(i)+(j - 1)*(N + 1)] * (v1i1j - 2 * v1ij + v1i_1j) +
				+mu*z22*(*S2j)[i - 1 + j*(N)] * (v1ij1 - 2 * v1ij + v1ij_1 + (1 / 4.)*(v2i1j1 - v2i1j_1 - v2i_1j1 + v2i_1j_1));
			T2[i][j] = mu*z11*(*S1i)[(i)+(j - 1)*(N + 1)] * ((1 / 4.)*(v1i1j1 - v1i1j_1 - v1i_1j1 + v1i_1j_1) + v2i1j - 2 * v2ij + v2i_1j) +
				+2 * mu*z22*(*S2j)[i - 1 + j*(N)] * (v2ij1 - 2 * v2ij + v2ij_1);*/
			T1[i][j] = mu*z11*(*S1i)[(i)+(j - 1)*(N + 1)] * (v1i1j - 2 * v1ij + v1i_1j) + mu*z22*(*S2j)[i - 1 + j*(N)] * (v1ij1 - 2 * v1ij + v1ij_1);
			T2[i][j] = mu*z11*(*S1i)[(i)+(j - 1)*(N + 1)] * (v2i1j - 2 * v2ij + v2i_1j) + mu*z22*(*S2j)[i - 1 + j*(N)] * (v2ij1 - 2 * v2ij + v2ij_1);//cout << "*" << mu*z11*(*S1i)[(i)+(j - 1)*(N + 1)] << "  " << (1 / 4.)*(v1i1j1 - v1i1j_1 - v1i_1j1 + v1i_1j_1) << " " << i << " " << j << "\n";
			//cout << mu*z22*(*S2j)[i - 1 + j*(N)] << " " << v2ij << " " << v2ij1 << " " << v2ij_1 << " " << v2ij1 - 2 * v2ij + v2ij_1 << "\n";
			//std::cout << v1i1j1 << " " << v1i1j_1<< " " << v1i_1j1 << " " << v1i_1j_1 <<  "\n\n";



			/*for (j = 1; j <= M; j++){
			for (i = 1; i <= N; i++){
			w = (*B)[i - 1 + (j - 1)*N].w;
			//	w = 1 / (N*M);
			d = tau*alpha / w;
			//std::cout << "d=" << d << "\n";
			T1[i][j] = ((*S1i)[(i)+(j - 1)*(N + 1)] * sigma11[i + 1][j] - (*S1i)[i - 1 + (j - 1)*(N + 1)] * sigma11[i-1][j] +
			+(*S2j)[i - 1 + j*(N)] * sigma21[i][j + 1] - (*S2j)[i - 1 + (j - 1)*(N)] * sigma21[i][j-1]) / 2;
			T2[i][j] = ((*S1i)[(i)+(j - 1)*(N + 1)] * sigma12[i + 1][j] - (*S1i)[i - 1 + (j - 1)*(N + 1)] * sigma12[i-1][j] +
			+(*S2j)[i - 1 + j*(N)] * sigma22[i][j + 1] - (*S2j)[i - 1 + (j - 1)*(N)] * sigma22[i][j-1]) / 2;*/
			//cout << "*" << T1[i][j] << "  " << T2[i][j] << " " << i << " " << j << "\n\n";

			sum1 = (*S1i)[(i)+(j - 1)*(N + 1)] * v10[i + 1][j] - (*S1i)[i - 1 + (j - 1)*(N + 1)] * v10[i - 1][j] +
				+(*S2j)[i - 1 + j*(N)] * v20[i][j + 1] - (*S2j)[i - 1 + (j - 1)*(N)] * v20[i][j - 1];

			p1[i][j] = p0[i][j] - ((d*p0[i][j]) / (eps))*sum1 / 2;

			sum2 = (*S1i)[(i)+(j - 1)*(N + 1)] * (v10[i + 1][j] * v10[i + 1][j] + p0[i + 1][j]) - (*S1i)[i - 1 + (j - 1)*(N + 1)] * (v10[i - 1][j] * v10[i - 1][j] + p0[i - 1][j]) +
				+(*S2j)[i - 1 + j*(N)] * (v10[i][j + 1] * v20[i][j + 1]) - (*S2j)[i - 1 + (j - 1)*(N)] * (v10[i][j - 1] * v20[i][j - 1]);

			//std::cout << (*S1i)[(i)+(j - 1)*(N + 1)] << " " << (*p0)[i + 1 + j*(N + 2)] << " " << (*v10)[i + 1 + j*(N + 2)]<< "\n";
			//std::cout << (*S1i)[i - 1 + (j - 1)*(N + 1)] << " " << (*p0)[i - 1 + j*(N + 2)] << " " << (*v10)[i - 1 + j*(N + 2)] << "\n\n";
			//std::cout << (*S1i)[(i)+(j - 1)*(N + 1)] * ((*v10)[i + 1 + j*(N + 2)] * (*v10)[i + 1 + j*(N + 2)] + (*p0)[i + 1 + j*(N + 2)]) -(*S1i)[i - 1 + (j - 1)*(N + 1)] * ((*v10)[i - 1 + j*(N + 2)] * (*v10)[i - 1 + j*(N + 2)] + (*p0)[i - 1 + j*(N + 2)]) << " " << (*p0)[i + 1 + j*(N + 2)] << " " << (*v10)[i + 1 + j*(N + 2)] << "\n";
			//std::cout << (*S1i)[i - 1 + (j - 1)*(N + 1)] << " " << (*p0)[i - 1 + j*(N + 2)] << " " << (*v10)[i - 1 + j*(N + 2)] << "\n\n";

			v11[i][j] = v10[i][j] - (d / 2)*(-2 * T1[i][j] + sum2);

			//std::cout << (*v11)[i + j*(N + 2)] << " " << (*v10)[i + j*(N + 2)] << " " << T1[i + j*(N + 1)] << " " << sum2 << "\n\n";
			//std::cout <<i<<" "<<j<<" "<< (*v11)[i + j*(N + 2)] << " " << (*v10)[i + j*(N + 2)] << " " << T1[i + j*(N + 1)] << " " << sum2<<"\n\n" ;

			sum3 = (*S2j)[i - 1 + j*(N)] * (v20[i][j + 1] * v20[i][j + 1] + p0[i][j + 1]) - (*S2j)[i - 1 + (j - 1)*(N)] * (v20[i][j - 1] * v20[i][j - 1] + p0[i][j - 1]) +
				+(*S1i)[(i)+(j - 1)*(N + 1)] * (v10[i + 1][j] * v20[i + 1][j]) - (*S1i)[i - 1 + (j - 1)*(N + 1)] * (v10[i - 1][j] * v20[i - 1][j]);

			v21[i][j] = v20[i][j] - (d / 2)*(-2 * T2[i][j] + sum3);
			//cout << v21[i][j] << " " << p0[i][j + 1] << " " << p0[i][j - 1] <<  "\n";
			//cout << v1i1j1 << " " << v1i1j_1 << " " << v1i_1j1 << " " << v1i_1j_1 << " " << v1i1j1 - v1i1j_1 - v1i_1j1 + v1i_1j_1  <<"\n\n";
			sum1 = 0;
			sum2 = 0;
			sum3 = 0;
		}
		//printf("\n\n");
	}

	for (int i = 0; i <= N + 1; i++){
		int j = 0;

		p1[i][j] = p1[i][j + 1];
		v21[i][j] = -v21[i][j + 1];
		v11[i][j] = -v11[i][j + 1];

		j = M + 1;
		v11[i][j] = -v11[i][j - 1];
		p1[i][j] = p1[i][j - 1];
		v21[i][j] = -v21[i][j - 1];

	}
	double y,x;
	for (int j = 0; j <= M + 1; j++){
		int i = 0;
		x = (1 / 20.);
		p1[i][j] = 1-2*mu*x;

		v21[i][j] = 0;

		y = (1 / 10.)*(j - 1 / 2.);
		v11[i][j] = y*(2 - y);

		i = N + 1;
		p1[i][j] = 1 - 2 * mu*(2 + 2 / N);

		v21[i][j] = v21[i - 1][j];
		v11[i][j] = v11[i - 1][j];

	}
	v11[0][0] = -v11[0][1];//(0,0)
	v11[0][M + 1] = -v11[0][M];
	/*v11[0][0] = 2 * v11[1][0] - v11[2][0];//(0,0)
	v21[0][0] = 2 * v21[1][0] - v21[2][0];
	v11[N + 1][0] = 2 * v11[N][0] - v11[N - 1][0];//(N+1,0)
	v21[N + 1][0] = 2 * v21[N][0] - v21[N - 1][0];
	v11[0][M + 1] = 2 * v11[1][M + 1] - v11[2][M + 1];//(0,M+1)
	v21[0][M + 1] = 2 * v21[1][M + 1] - v21[2][M + 1];
	v11[N + 1][M + 1] = 2 * v11[N][M + 1] - v11[N - 1][M + 1];//(N+1,M+1)
	v21[N + 1][M + 1] = 2 * v21[N][M + 1] - v21[N - 1][M + 1];*/
	


}


int main(void){
	el *V = NULL;
	box *B = NULL;
	double *S1i = NULL;
	double *S2i = NULL;
	double *S1j = NULL;
	double *S2j = NULL;
	vector<vector<double>> p0(N + 2, vector<double>(M + 2, 0));
	vector<vector<double>> v10(N + 2, vector<double>(M + 2, 0));
	vector<vector<double>> v20(N + 2, vector<double>(M + 2, 0));
	vector<vector<double>> p1(N + 2, vector<double>(M + 2, 0));
	vector<vector<double>> v11(N + 2, vector<double>(M + 2, 0));
	vector<vector<double>> v21(N + 2, vector<double>(M + 2, 0));
	/*vector<vector<double>> v11_4;
	vector<vector<double>> v21_4;
	vector<vector<double>> p11_4;
	vector<vector<double>> v11_2;
	vector<vector<double>> p11_2;
	vector<vector<double>> v21_2;
	vector<vector<double>> ksi;
	vector<vector<double>> ksi_1;
	vector<vector<double>> psi;
	vector<vector<double>> psi_1;
	vector<vector<double>> ksi_1_2;*/
	double *ksi=NULL;
	double *ksi_1 = NULL;
	double *psi = NULL;
	double *psi_1 = NULL;
	double *ksi_1_2 = NULL;
	double maxv1, maxv2, maxp;
	FILE  *out, *lines, *component, *v1, *v2, *p;

	fopen_s(&component, "component.txt", "w+");

	fopen_s(&v1, "v1.xlsх", "w+");
	fopen_s(&v2, "v2.xlsх", "w+");
	fopen_s(&p, "p.xlsх", "w+");

	fopen_s(&out, "out.txt", "w+");
	fopen_s(&lines, "isolines.txt", "w+");

	fragmentation(N, M, &V, &B);
	projection(N, M, &S1i, &S2i, &S1j, &S2j, &V);
	initial(N, M, p0, v10, v20, &V);
	for (int l = 0; l <= 1000; l++){
		fprintf(out, "new lay :%d\n", l);
		fprintf(v1, "new lay :%d\n", l);
		fprintf(v2, "new lay :%d\n", l);
		fprintf(p, "new lay :%d\n", l);
		fprintf(component, "new lay :%d\n", l);


		for (int j = M + 1; j >= 0; --j){
			for (int i = 0; i <= N + 1; ++i){
				fprintf(v1, "%.5f   ", v10[i][j]);
				fprintf(v2, "%.5f   ", v20[i][j]);
				fprintf(p, "%.5f   ", p0[i][j]);
			}
			fprintf(v1, "\n");
			fprintf(v2, "\n");
			fprintf(p, "\n");
		}

		predict(N, M, p0, v10, v20, p1, v11, v21, &B, &S1i, &S2i, &S1j, &S2j, &V);



		for (int j = 0; j <= M + 1; j++){
			for (int i = 0; i <= N + 1; i++){
				//fprintf(component, "i=%d      j=%d\n", i, j);
				fprintf(component, "%.5f       %.5f      %.5f     |", p1[i][j], v11[i][j], v21[i][j]);
				fprintf(component, "%.5f       %.5f      %.5f     \n ", p0[i][j], v10[i][j], v20[i][j]);

			}
			fprintf(component, "\n\n");
			/*fprintf(v1, "\n");
			fprintf(v2, "\n");
			fprintf(p, "\n");*/
		}
		maxv1 = abs(v10[0][0] - v11[0][0]);
		maxv2 = abs(v20[0][0] - v21[0][0]);
		//maxv2 = abs(v21[0][0]);
		maxp = abs(p0[0][0] - p1[0][0]);
		for (int i = 1; i <= N; i++){
			for (int j = 1; j <= M; j++){
				if (abs(v10[i][j] - v11[i][j]) > maxv1){
					maxv1 = abs(v10[i][j] - v11[i][j]);
				}
				if (abs(v20[i][j] - v21[i][j]) > maxv2){
					maxv2 = abs(v20[i][j] - v21[i][j]);
				}
				if (abs(p0[i][j] - p1[i][j]) > maxp){
					maxp = abs(p0[i][j] - p1[i][j]);
				}
			}
		}

		fprintf(out, "%.15f        %.15f        %.15f\n", maxv1, maxv2, maxp);
		for (int j = 0; j <= M + 1; j++){
			for (int i = 0; i <= N + 1; i++){
				p0[i][j] = p1[i][j];
				v10[i][j] = v11[i][j];
				v20[i][j] = v21[i][j];//присваивание
				//cout << p1[i][j] << " " << v11[i][j] << " " << v21[i][j] << "\n";


			}
		}
		fprintf(v1, "\n\n");
		fprintf(v2, "\n\n");
		fprintf(p, "\n\n");
		fprintf(p, "\n\n");
		initial_ksi(N, M, v11, v21, &ksi, &psi, tau, alpha);
		psi_newlay_1(N, M, &ksi, &ksi_1_2, tau, alpha);
		psi_newlay_2(N, M,  &ksi_1, &ksi_1_2, tau, alpha);
		psi_final(N, M, &psi_1, &psi, &ksi_1, tau);

		int j = 0;
		for (int j = 1; j <= M; j++){
			for (int i = 1; i <= N; i++){
				//fprintf(lines, "%.5f          %.5f\n", ksi_1_2[i + j*(N + 2)], ksi_1[i + j*(N + 2)]);
				fprintf(lines, "%.5f        %.5f        %.5f\n", V[i - 1 + (j - 1)*(N + 1)].x + 1 * (2.0 / (2 * N)), V[i - 1 + (j - 1)*(N + 1)].y + 1 * (2.0 / (2 * M)), psi_1[i + (j)*(N + 2)]);
				//fprintf(lines, "%.5f", psi_1[i + j*(N + 2)]);
				//fprintf(lines, "%.5f            %.5f\n", y0 + (i-1)*1.0 / M, psi_1[i + j*(N + 2)]);
			}
			fprintf(lines, "\n\n");

		}
	}

	fclose(out);
	fclose(component);
	fclose(lines);
	fclose(v1);
	fclose(v2);
	fclose(p);

	return 0;
}
