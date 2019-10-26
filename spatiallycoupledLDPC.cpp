//****************************************************************************************
//*																						 *
//*								EE 388 - Spatial Coupling LDPC		     				 *
//*																						 *
//****************************************************************************************

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <malloc.h>

#define m_PROTOGRAPH				100
#define L_SPATIALCOUPLING			500
#define k_LDPC_reg_ensemble			10
#define l_LDPC_reg_ensemble			5		// k/l = q, an integer

#define MAX_BP_ITER					200
#define MAX_BP_UNCHANGED_ITER		10
#define SIMULATION_NUM				20
double eps_vals[] = { 0, 0.2, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.475, 0.5, 0.6, 0.7, 0.8, 0.9, 1 };


// Create Fisher Yates random permutation of 0,1,...,n-1
// Output in vector[]
//
void randperm(int n, int* vector) {
	int k, j, temp;
	for (k = 0; k < n; k++) {
		vector[k] = k;
	}
	for (k = 0; k < n; k++) {
		j = rand() % (n - k) + k;
		temp = vector[j];
		vector[j] = vector[k];
		vector[k] = temp;
	}
}


// Create regular (l,k) LDPC code
//_ #var node = n
//_ m = n*l/k integer
//_ varConnect[i][.] = a if (a,i) is an edge, and -1 otherwise (assumed at initialization), always ends with -1!!!
//
void regularLDPC(int l, int k, int n, int **varConnect) {
	int i, a, *vector;	
	int m = n*l / k, numEdge = n*l;
	vector = (int*)calloc(numEdge, sizeof(int)); if (vector == NULL) exit(1);
	randperm(numEdge, vector);
	for (i = 0; i < n; i++) {
		for (a = 0; a < l; a++) {
			varConnect[i][a] = vector[a + i*l]/k;
		}
	}
	free(vector);
}


// Create L-coupling for regular(l,k), q = k/l integer
// varConnect[i][.] = a if (a,i) is an edge, and -1 otherwise (assumed at initialization), always ends with -1!!!
//
void coupling_regular(int l, int q, int L, int **varConnect) {
	int i, a, cnt, r, l1, l2;
	if ((l % 2) == 1) {
		l1 = (l - 1) / 2;
		l2 = (l - 1) / 2;
	}
	else {
		l1 = (l - 2) / 2;
		l2 = l / 2;
	}
	for (i = 0; i < L; i++) {
		cnt = 0;
		for (a = i; a < i + l; a++) {
			for (r = 0; r < q; r++) {
				varConnect[r*L+i][cnt] = a;
			}
			cnt++;
		}
	}
}


// Create m-fold protograph code, from a protograph
// Protograph: #variable nodes = n, #check nodes = p, allow double edges
// varConnect[i][.] = a if (a,i) is an edge, and -1 otherwise (assumed at initialization), always ends with -1!!!
//
void protograph(int m, int p, int n, int **varConnectProtograph, int **varConnect) {
	int a, i, cnt, b;	
	int *vector, *posVar;
	vector = (int*)calloc(m, sizeof(int)); if (vector == NULL) exit(1);
	posVar = (int*)calloc(m*n, sizeof(int)); if (posVar == NULL) exit(1);
	for (i = 0; i < n; i++) {
		cnt = 0;
		while (1) {
			if (varConnectProtograph[i][cnt] != -1) {
				a = varConnectProtograph[i][cnt];
				randperm(m, vector);
				for (b = 0; b < m; b++) {
					varConnect[b + m*i][posVar[b + m*i]] = vector[b] + a*m;
					posVar[b + m*i]++;
				}
				cnt++;
			}
			else {
				break;
			}
		}
	}
	free(vector);
	free(posVar);
}


// Get rid of double edges in varConnect & complete other information about the code based on varConnect
// n = # variable nodes
// chkConnect[a][.] = i if (a,i) is an edge, and -1 otherwise (assumed at initialization), always ends with -1!!!
// degVar[i] = degree of variable node i, must be initialized to 0
// degChk[a] = degree of check node a, must be initialized to 0
// varConnectPosition[i][j] = k, for chkConnect[varConnect[i][j]][k] = i (i.e. trace the corresponding index in chkConnect), always ends with -1!!!
// chkConnectPosition[a][j] = k, for varConnect[chkConnect[a][j]][k] = a (i.e. trace the corresponding index in varConnect), always ends with -1!!!
//
void complete_graph(int n, int **varConnect, int **chkConnect, int *degVar, int *degChk, int **varConnectPosition, int **chkConnectPosition) {
	int i, cnt, j, newlength;

	// Get rid of double edges
	for (i = 0; i < n; i++) {
		newlength = 1;
		cnt = 0;
		while (1) {
			if (varConnect[i][cnt]==-1) {
				break;
			}
			for (j = 0; j < newlength; j++) {
				if (varConnect[i][cnt] == varConnect[i][j]) {
					break;
				}
			}
			if (j == newlength) {
				varConnect[i][newlength] = varConnect[i][cnt];
				newlength++;
			}
			cnt++;
		}
		cnt = newlength;
		while (1) {
			if (varConnect[i][cnt] != -1) {
				varConnect[i][cnt] = -1;
			}
			else {
				break;
			}
		}
	}

	// Complete other info
	for (i = 0; i < n; i++) {
		cnt = 0;
		while (1) {
			if (varConnect[i][cnt] != -1) {
				chkConnect[varConnect[i][cnt]][degChk[varConnect[i][cnt]]] = i;
				varConnectPosition[i][cnt] = degChk[varConnect[i][cnt]];
				chkConnectPosition[varConnect[i][cnt]][degChk[varConnect[i][cnt]]] = cnt;
				degChk[varConnect[i][cnt]]++;
				degVar[i]++;
				cnt++;
			}
			else {
				break;
			}
		}
	}
}


// Return #check nodes in m-folded L-Spatial Coupling code from regular LDPC(k,l)
//
int chkSize_SC_regLDPC(int l, int m, int L) {
	int l1, l2;
	if ((l % 2) == 1) {
		l1 = (l - 1) / 2;
		l2 = (l - 1) / 2;
	}
	else {
		l1 = (l - 2) / 2;
		l2 = l / 2;
	}
	return m*(L + l1 + l2);
}


// BEC(eps) channel simulation
// Binary input {0, 1}
// Output: {0, 1, -1}, where erasure = -1
// If input=NULL, treat as all-0
//
void BEC(double eps, int n, int *input, int *output) {
	int i;	
	if (input == NULL) {
		for (i = 0; i < n; i++) {
			if ((rand()*1.0 / RAND_MAX) < 1 - eps) {
				output[i] = 0;
			}
			else {
				output[i] = -1;
			}
		}
	}
	else {
		for (i = 0; i < n; i++) {
			if ((rand()*1.0 / RAND_MAX) < 1 - eps) {
				output[i] = input[i];
			}
			else {
				output[i] = -1;
			}
		}
	}
}


// Return i if vector[i] = k (for first time), and -1 if none found.
// vector[last entry] = endMarker
//
int findPosition(int k, int *vector, int endMarker) {
	int cnt = 0;
	while (1) {
		if (vector[cnt] == k) {
			return cnt;
		}
		else if (vector[cnt] == endMarker) {
			return -1;
		}
		else {
			cnt++;
		}
	}
}


// Belief Propagation (or rather MP) decoding on BEC
// This assumes all-0 codeword is transmitted
// Output: decodedword
// Code parameters: varConnect, chkConnect, varConnectPosition, chkConnectPosition, degVar, degChk, varNum (#var nodes), 
// chkNum (#chk nodes), memorySizeVar (memory size of array varConnect[.]), memorySizeChk (memory size of array chkConnect[.])
// BP Parameters: MAX_BPIter (max # BP iterations), MAX_UnchangedIter (max # iterations that decoded word remains unchanged)
//
// First do Message Passing version, then convert to whatever that is asked from BP:
// decodedword_without_yi: decoded word without using y[i].
//
void BPDecoder_BEC(int *y, int **varConnect, int **chkConnect, int **varConnectPosition, int **chkConnectPosition, int *degVar, int *degChk, int varNum, int chkNum, 
	int memorySizeVar, int memorySizeChk, int MAX_BPIter, int MAX_UnchangedIter, int *decodedword, int *decodedword_without_yi) {
	int i, a, j, r, iter, unchangedCnt, changed, sumMessage;	
	int **messageVar2Chk, **messageChk2Var, *prevword;
	messageVar2Chk = (int**)calloc(varNum, sizeof(int)); if (messageVar2Chk == NULL) exit(1);
	for (i = 0; i<varNum; i++) {
		messageVar2Chk[i] = (int*)calloc(memorySizeVar, sizeof(int)); if (messageVar2Chk[i] == NULL) exit(1);
	}
	messageChk2Var = (int**)calloc(chkNum, sizeof(int)); if (messageChk2Var == NULL) exit(1);
	for (i = 0; i<chkNum; i++) {
		messageChk2Var[i] = (int*)calloc(memorySizeChk, sizeof(int)); if (messageChk2Var[i] == NULL) exit(1);
		for (j = 0; j < memorySizeChk; j++) {
			messageChk2Var[i][j] = -1;
		}
	}
	prevword = (int *)calloc(varNum, sizeof(int)); if (prevword == NULL) exit(1);
	for (i = 0; i < varNum; i++) {
		prevword[i] = 2;
	}

	iter = 0;
	unchangedCnt = 0;
	while (1) {
		iter++;
		for (i = 0; i < varNum; i++) {
			sumMessage = 0;
			for (j = 0; j < degVar[i]; j++) {
				sumMessage += messageChk2Var[varConnect[i][j]][varConnectPosition[i][j]];
			}
			for (j = 0; j < degVar[i]; j++) {
				if ((y[i] == 0) || ((sumMessage - messageChk2Var[varConnect[i][j]][varConnectPosition[i][j]]) != (1 - degVar[i]))) {
					messageVar2Chk[i][j] = 0;
				} else {
					messageVar2Chk[i][j] = -1;
				}
			}
		}
		for (a = 0; a < chkNum; a++) {
			sumMessage = 0;
			for (j = 0; j < degChk[a]; j++) {
				sumMessage += messageVar2Chk[chkConnect[a][j]][chkConnectPosition[a][j]];
			}
			for (j = 0; j < degChk[a]; j++) {
				if ((sumMessage - messageVar2Chk[chkConnect[a][j]][chkConnectPosition[a][j]]) != 0) {
					messageChk2Var[a][j] = -1;
				}
				else {
					messageChk2Var[a][j] = 0;
				}
			}
		}
		for (i = 0; i < varNum; i++) {
			sumMessage = 0;
			for (j = 0; j < degVar[i]; j++) {
				sumMessage += messageChk2Var[varConnect[i][j]][varConnectPosition[i][j]];
			}
			if (y[i] == 0 || sumMessage != (-degVar[i])) {
				decodedword[i] = 0;
			}
			else {
				decodedword[i] = -1;
			}
			if (sumMessage != (-degVar[i])) {
				decodedword_without_yi[i] = 0;
			}
			else {
				decodedword_without_yi[i] = -1;
			}
		}

		changed = 0;		
		for (i = 0; i < varNum; i++) {
			if (decodedword[i] != prevword[i]) {
				changed = 1;
				break;
			}
		}
		if (changed == 0) {
			unchangedCnt++;
		}
		else {
			unchangedCnt = 0;
		}

		if (iter >= MAX_BPIter || unchangedCnt >= MAX_UnchangedIter) {
			break;
		}
	}

	free(prevword);
	for (i = 0; i < chkNum; i++) free(messageChk2Var[i]);
	free(messageChk2Var);
	for (i = 0; i < varNum; i++) free(messageVar2Chk[i]); 
	free(messageVar2Chk);
}



int _tmain(int argc, _TCHAR* argv[])
{
	srand(time(NULL));

	// ------ Code construction ------
	int m, l, k, L, q;
	m = m_PROTOGRAPH;
	l = l_LDPC_reg_ensemble;
	k = k_LDPC_reg_ensemble;
	L = L_SPATIALCOUPLING;
	q = k / l;

	int **varConnectOriginal, **varConnect, **chkConnect, **varConnectPosition, **chkConnectPosition, *degVar, *degChk, varNum, chkNum;
	int i, j;

	varNum = q*L*m;
	chkNum = chkSize_SC_regLDPC(l, m, L);

	varConnectOriginal = (int**)calloc(q*L, sizeof(int)); if (varConnectOriginal == NULL) exit(1);
	for (i = 0; i<q*L; i++) { 
		varConnectOriginal[i] = (int*)calloc(l+1, sizeof(int)); if (varConnectOriginal[i] == NULL) exit(1);
		for (j = 0; j < l + 1; j++) {
			varConnectOriginal[i][j] = -1;
		}
	}
	varConnect = (int**)calloc(varNum, sizeof(int)); if (varConnect == NULL) exit(1);
	for (i = 0; i<varNum; i++) {
		varConnect[i] = (int*)calloc(l + 1, sizeof(int)); if (varConnect[i] == NULL) exit(1);
		for (j = 0; j < l + 1; j++) {
			varConnect[i][j] = -1;
		}
	}
	chkConnect = (int**)calloc(chkNum, sizeof(int)); if (chkConnect == NULL) exit(1);
	for (i = 0; i<chkNum; i++) {
		chkConnect[i] = (int*)calloc(k + 1, sizeof(int)); if (chkConnect[i] == NULL) exit(1);
		for (j = 0; j < k + 1; j++) {
			chkConnect[i][j] = -1;
		}
	}
	varConnectPosition = (int**)calloc(varNum, sizeof(int)); if (varConnectPosition == NULL) exit(1);
	for (i = 0; i<varNum; i++) {
		varConnectPosition[i] = (int*)calloc(l + 1, sizeof(int)); if (varConnectPosition[i] == NULL) exit(1);
		for (j = 0; j < l + 1; j++) {
			varConnectPosition[i][j] = -1;
		}
	}
	chkConnectPosition = (int**)calloc(chkNum, sizeof(int)); if (chkConnectPosition == NULL) exit(1);
	for (i = 0; i<chkNum; i++) {
		chkConnectPosition[i] = (int*)calloc(k + 1, sizeof(int)); if (chkConnectPosition[i] == NULL) exit(1);
		for (j = 0; j < k + 1; j++) {
			chkConnectPosition[i][j] = -1;
		}
	}
	degVar = (int*)calloc(varNum, sizeof(int)); if (degVar == NULL) exit(1);
	degChk = (int*)calloc(chkNum, sizeof(int)); if (degChk == NULL) exit(1);
	
	coupling_regular(l, q, L, varConnectOriginal);	
	protograph(m, chkNum/m, q*L, varConnectOriginal, varConnect);	
	complete_graph(varNum, varConnect, chkConnect, degVar, degChk, varConnectPosition, chkConnectPosition);
	
	// ------ End of code construction ------	

	int *y, *decodedword, *decodedword_without_yi; 
	int cnt_eps, cnt_simul, error, progress;
	double eps, *BER, *EXIT;
	y = (int*)calloc(varNum, sizeof(int)); if (y == NULL) exit(1);	
	decodedword = (int*)calloc(varNum, sizeof(int)); if (decodedword == NULL) exit(1);
	decodedword_without_yi = (int*)calloc(varNum, sizeof(int)); if (decodedword_without_yi == NULL) exit(1);
	BER = (double*)calloc((sizeof(eps_vals) / sizeof(double)), sizeof(double)); if (BER == NULL) exit(1);
	EXIT = (double*)calloc((sizeof(eps_vals) / sizeof(double)), sizeof(double)); if (EXIT == NULL) exit(1);

	progress = 0;
	for (cnt_eps = 0; cnt_eps < (sizeof(eps_vals) / sizeof(double)); cnt_eps++) {
		for (cnt_simul = 0; cnt_simul < SIMULATION_NUM; cnt_simul++) {
			// Channel simulation - Assume all-0 codeword is transmitted
			eps = eps_vals[cnt_eps];
			BEC(eps, varNum, NULL, y);

			// BP decoding
			BPDecoder_BEC(y, varConnect, chkConnect, varConnectPosition, chkConnectPosition, 
				degVar, degChk, varNum, chkNum, l + 1, k + 1, MAX_BP_ITER, MAX_BP_UNCHANGED_ITER, decodedword, decodedword_without_yi);

			// BER
			error = 0;
			for (i = 0; i < varNum; i++) {
				if (decodedword[i] != 0) {
					error++;
				}
			}
			BER[cnt_eps] += error*1.0 / varNum;

			// EXIT
			error = 0;
			for (i = 0; i < varNum; i++) {
				if (decodedword_without_yi[i] != 0) {
					error++;
				}
			}
			EXIT[cnt_eps] += error*1.0 / varNum;

			// Print progress
			progress++;
			printf("Complete %.3f percents\n", progress*100.0 / SIMULATION_NUM / (sizeof(eps_vals) / sizeof(double)));
		}
		BER[cnt_eps] /= SIMULATION_NUM;
		EXIT[cnt_eps] /= SIMULATION_NUM;
		printf("--- Eps = %.5f,   BER = %.10f,   EXIT = %.10f\n", eps, BER[cnt_eps], EXIT[cnt_eps]);
	}

	printf("\n");
	for (cnt_eps = 0; cnt_eps < (sizeof(eps_vals) / sizeof(double)); cnt_eps++) {	
		eps = eps_vals[cnt_eps];
		printf("--- Eps = %.5f,   BER = %.10f,   EXIT = %.10f\n", eps, BER[cnt_eps], EXIT[cnt_eps]);
	}
	printf("\n");
	

	// Print to file
	FILE* textfile;
	char str[100];
	sprintf_s(str, "Result_%d.txt", L_SPATIALCOUPLING);
	fopen_s(&textfile, str, "w");
	fprintf(textfile, "Eps:\n");
	for (cnt_eps = 0; cnt_eps < (sizeof(eps_vals) / sizeof(double)); cnt_eps++) {
		eps = eps_vals[cnt_eps];
		fprintf(textfile, "%.10f\n", eps);
	}
	fprintf(textfile, "\n");
	fprintf(textfile, "EXIT:\n");
	for (cnt_eps = 0; cnt_eps < (sizeof(eps_vals) / sizeof(double)); cnt_eps++) {
		eps = eps_vals[cnt_eps];
		fprintf(textfile, "%.10f\n", EXIT[cnt_eps]);
	}
	fprintf(textfile, "\n");
	fprintf(textfile, "BER:\n");
	for (cnt_eps = 0; cnt_eps < (sizeof(eps_vals) / sizeof(double)); cnt_eps++) {
		eps = eps_vals[cnt_eps];
		fprintf(textfile, "%.10f\n", BER[cnt_eps]);
	}
	fclose(textfile);


	// Release memory
	free(degVar); free(degChk);
	for (i = 0; i<q*L; i++) {free(varConnectOriginal[i]);} free(varConnectOriginal);
	for (i = 0; i<varNum; i++) {free(varConnect[i]);} free(varConnect);
	for (i = 0; i<chkNum; i++) {free(chkConnect[i]);} free(chkConnect);
	for (i = 0; i<varNum; i++) { free(varConnectPosition[i]); } free(varConnectPosition);
	for (i = 0; i<chkNum; i++) { free(chkConnectPosition[i]); } free(chkConnectPosition);
	free(y); free(decodedword);

	printf("\n\nHit Enter to exit...\n");
	getchar(); return 0;
}

