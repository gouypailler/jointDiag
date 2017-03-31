#include <R.h>
#include <math.h>
#include <stdio.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define size 3				/* dimension of matrix */

extern void dgeev_(char *JOBVL,char *JOBVR,int *N,double *A,int *LDA,
	double *WR,double *WI,double *VL,int *LDVL,double *VR,int *LDVR,
	double *WORK,int *LWORK,int *INFO);

// be careful because we do not check anything
void eigen3(double *x, double *V, double *D) {
	int n=size;
	int i,j; 
	int	lwork, info;
	double *work, wI[size], *left;
	double xvals[size*size], tmp;
	char jobVL[1], jobVR[1];
	n = size;
	/* work on a copy of x */
	memcpy(xvals, x, n * n*sizeof(double));
	jobVL[0] = jobVR[0] = 'N';
	left = (double *) 0;
	jobVR[0] = 'V';
	/* ask for optimal size of work array */
	lwork = -1;
	dgeev_(jobVL, jobVR, &n, xvals, &n, D, wI,
		    left, &n, V, &n, &tmp, &lwork, &info);
	if (info != 0) fprintf(stderr, "failure with error %d\n", info);
	lwork = (int) tmp;
	work = (double*)malloc(lwork*sizeof(double));
	dgeev_(jobVL, jobVR, &n, xvals, &n, D, wI,
		    left, &n, V, &n, work, &lwork, &info);
	for (i=0; i<size; i++) { 
		if (fabs(wI[i])>1e-12) 
			fprintf(stderr, "failure: imaginary parts in evd\n");
	}
	if (info != 0) fprintf(stderr, "failure with error %d\n", info);
	double auxv[size],aux;
	if (fabs(D[1])<fabs(D[2])) {
		aux = D[1];
		D[1] = D[2];
		D[2] = aux;
		memcpy(auxv,&V[3],3*sizeof(double));
		memcpy(&V[3],&V[6],3*sizeof(double));
		memcpy(&V[6],auxv,3*sizeof(double));
	}
	if (fabs(D[0])<fabs(D[1])) {
		aux = D[0];
		D[0] = D[1];
		D[1] = aux;
		memcpy(auxv,V,3*sizeof(double));
		memcpy(V,&V[3],3*sizeof(double));
		memcpy(&V[3],auxv,3*sizeof(double));
	}
	if (fabs(D[1])<fabs(D[2])) {
		aux = D[1];
		D[1] = D[2];
		D[2] = aux;
		memcpy(auxv,&V[3],3*sizeof(double));
		memcpy(&V[3],&V[6],3*sizeof(double));
		memcpy(&V[6],auxv,3*sizeof(double));
	}
	free(work);
}	

double m_max(double a, double b) {
	if (a>b) return a;
	return b;
}

void iterJDI(double *C, int *pMatSize, int *pMatNumber, int *ptn,int *ptm,
				double *s_max, double *sh_max, double *A) {
	int n=*ptn-1;
	int m=*ptm-1;
	int MN=*pMatNumber;
	int MS=*pMatSize;
	int i,j,k;
	int d3,d2;

	double tmm[MN];
	double tnn[MN];
	double tmn[MN];
	for (i=0,d3=0;i<MN;i++,d3+=MS*MS) {
		tmm[i] = C[d3+m*MS+m];
		tnn[i] = C[d3+n*MS+n];
		tmn[i] = C[d3+n*MS+m];
	}

	// here we evd
	double G[MN][3];
	double evectors[9], evalues[3];
	for (i=0;i<MN;i++) {
		G[i][0] = 0.5*(tmm[i]+tnn[i]);	
		G[i][1] = 0.5*(tmm[i]-tnn[i]);
		G[i][2] = tmn[i];	
	}
	double GG[9];
	for (i=0;i<3;i++) {
		for (j=0;j<=i;j++) {
			GG[3*j+i] = 0;
			for (k=0;k<MN;k++) GG[3*j+i] += G[k][i]*G[k][j];
			GG[3*i+j] = GG[3*j+i];
		}
	}
	for (i=0;i<3;i++) GG[3*i] *= -1;
	eigen3(GG,evectors,evalues);

	double aux[9];
	double diagaux[3];
	double v[3];
	for (i=0;i<9;i++) aux[i] = evectors[i];
	for (i=0;i<3;i++) aux[3*i] *= -1;
	for (i=0;i<3;i++) {
		diagaux[i]	= 0;
		for (j=0;j<3;j++) diagaux[i] += evectors[3*i+j] * aux[3*i+j];
		diagaux[i] = 1/sqrt(fabs(diagaux[i]));
	}
	int ind_min=2;
	for (i=0;i<3;i++) v[i] = evectors[3*ind_min+i]*diagaux[ind_min];
	if (v[2]<0) for (i=0;i<3;i++) v[i]*=-1;
	double ch = sqrt( (1+sqrt(1+v[0]*v[0]))/2 );
	double sh = v[0]/(2*ch);
	double c = sqrt( ( 1 + v[2]/sqrt(1+v[0]*v[0]) )/2 );
	double s = -v[1]/( 2*c*sqrt( (1+v[0]*v[0]) ) );
	*s_max = m_max(*s_max,fabs(s));
	*sh_max = m_max(*sh_max,fabs(sh));
	
	double rm11=c*ch - s*sh;
	double rm12=c*sh - s*ch;
	double rm21=c*sh + s*ch;
	double rm22=c*ch + s*sh;

	double h_slice1,h_slice2;
	double buf[MS][MN];
	for (i=0;i<MS;i++) {
		for (j=0,d3=0;j<MN;j++,d3+=MS*MS) {
			h_slice1 = C[d3+i*MS+m];
			h_slice2 = C[d3+i*MS+n];
			buf[i][j] = rm11*h_slice1 + rm21*h_slice2;
			C[d3+i*MS+n] = rm12*h_slice1 + rm22*h_slice2;
			C[d3+i*MS+m] = buf[i][j];
		}
	}

	for (i=0;i<MS;i++) {
		for (j=0,d3=0;j<MN;j++,d3+=MS*MS) {
			C[d3+m*MS+i] = buf[i][j];
		}
	}	
	for (i=0;i<MS;i++) {
		for (j=0,d3=0;j<MN;j++,d3+=MS*MS) {
			C[d3+n*MS+i] = C[d3+i*MS+n];
		}
	}
	for (i=0;i<MN;i++) {
		C[i*MS*MS+m*MS+m] = (rm11*rm11)*tmm[i]+(2*rm11*rm21)*tmn[i]+(rm21*rm21)*tnn[i];
		C[i*MS*MS+n*MS+m] = (rm11*rm12)*tmm[i]+(rm11*rm22+rm12*rm21)*tmn[i]+(rm21*rm22)*tnn[i];
		C[i*MS*MS+m*MS+n] = C[i*MS*MS+n*MS+m];
		C[i*MS*MS+n*MS+n] = (rm12*rm12)*tmm[i]+(2*rm12*rm22)*tmn[i]+(rm22*rm22)*tnn[i];
	}

	double col1;
	double col2;
	for (i=0;i<MS;i++) {
		col1 = A[m*MS+i];
		col2 = A[n*MS+i];
		A[m*MS+i] = rm22*col1 - rm12*col2;
		A[n*MS+i] = -rm21*col1 + rm11*col2;
	}
}

void sweepjedi(double *C, int *pMatSize, int *pMatNumber,
				double *s_max, double *sh_max, double *A) {
	for (int n=2;n<=*pMatSize;n++) 
		for (int m=1;m<n;m++) 
			iterJDI(C, pMatSize, pMatNumber, &n, &m, s_max, sh_max, A);
	int MS=*pMatSize;
	int MN=*pMatNumber;
	double col_norm[MS];
	for (int i=0;i<MS;i++) {
		col_norm[i] = 0;
		for (int j=0;j<MS;j++) col_norm[i] += A[i*MS+j]*A[i*MS+j];
		col_norm[i] = sqrt(col_norm[i]);
	}
	double daux=1;
	double d[MS];
	for (int i=0;i<MS;i++) daux*=col_norm[i];
	daux=pow((double)daux,(double) 1/MS);
	for (int i=0;i<MS;i++) d[i]=daux/col_norm[i];
	for (int i=0;i<MS;i++) 
		for (int j=0;j<MS;j++)
			A[j*MS+i] *= d[j];
	for (int k=0;k<MN;k++) {
		for (int i=0;i<MS;i++) {
			for (int j=0;j<MS;j++) {
				C[k*MS*MS+i*MS+j] *= (1/d[i])*(1/d[j]);
			}
		}
	}
}	

















