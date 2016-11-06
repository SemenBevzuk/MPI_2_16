#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include <stddef.h>
#include <cstddef>
#include <random>

#define EPS 0.0001
//#define MSMPI_NO_DEPRECATE_20


int* CreateAndFillMatrixVector(int n, int m) {
	int i;
	int* matrix = (int*)malloc(n * m * sizeof(int));
	for (i = 0; i < n*m; i++) {
		matrix[i] = 0 + rand() % 10;
	}
	return matrix;
}

int FindMaxInMatrix(int* matrix, int n, int m) {
	int i, j;
	int max = matrix[0];
	int *temp = matrix;
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			if (*temp > max) {
				max = matrix[i*n + j];
			}
			temp++;
		}
	}
	return max;
}

void DeleteMatrix(int* matrix) {
	free(matrix);
}

void PrintMatrixVector(int* matrix, int n, int m) {
	int i, j;
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			printf_s("%5.1d ", matrix[i*m + j]);
		}
		printf_s("\n");
	}
}

void PrintVector(int *data, int dataSize) {
	int i;
	for (i = 0; i < dataSize; i++) {
		printf_s("%5d ", data[i]);
		if (i % 10 == 9) {
			printf_s("\n");
		}
	}
	printf_s("\n");
}

int* ConsistentMatrixMultiplication(int* Amatrix, int* Bmatrix, int m, int n, int l)//A*B=C
{
	int* Cmatrix = (int*)malloc(m * l * sizeof(int));
	int i = 0;
	int j = 0;
	int k = 0;
	int temp;
	for (i = 0; i<m; i++) {
		for (j = 0; j<l; j++) {
			temp = 0;
			for (k = 0; k<n; k++) {
				temp += Amatrix[i*n + k] * Bmatrix[k*l + j];
				Cmatrix[i*l + j] = temp;
			}
		}
	}
	return Cmatrix;
}
double StartConsistentMatrixMultiplication(int* Amatrix, int* Bmatrix, int m, int n, int l)
{
	double time1, time2;
	int* Cmatrix = (int*)malloc(m * l * sizeof(int));
	time1 = MPI_Wtime();
	Cmatrix = ConsistentMatrixMultiplication(Amatrix, Bmatrix, m,n,l);
	time2 = MPI_Wtime();
	return time2 - time1;
}

int ScalarMultiplication(int* A, int* B, int n)
{
	int i=0, res=0;
	for (i = 0; i < n; i++)
	{
		res += A[i] * B[i];
	}
	return res;
}
bool CheckResult(int* A, int* B, int m, int l)
{
	for (int i = 0; i < m*l; i++)
	{
		if (A[i]!=B[i])
		{
			return false;
		}
	}
	return true;
}
using namespace std;

int main(int argc, char* argv[]) {
	double time1, time2, delta_time_parallel, delta_time_consistent;
	MPI_Status status;
	FILE *f = NULL;
	int i,j,k;
	int n, m, l;

	int *Amatrix = NULL;
	int *Bmatrix = NULL;
	int *Cmatrix = NULL;
	int *Abuff = NULL;
	int *Bbuff = NULL;
	int thread_count, rank;
	int dataSize, bufferSize;
	int remainingData;

	if (argc >= 4) {
		m = atoi(argv[1]);
		n = atoi(argv[2]);
		l = atoi(argv[3]);
	}
	else {
		printf_s("Error with argv: argc!=4\n");
	}
	int* result_Multiplication = (int*)malloc(l * sizeof(int));
	int* result_matrix = (int*)malloc(m*l* sizeof(int));

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &thread_count);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Datatype StructTypeColomn;
	int count = n;
	int *blocklens = new int[n];
	MPI_Aint *indices = new MPI_Aint[n];
	MPI_Datatype* oldtypes = new MPI_Datatype[n];
	for (i = 0; i<n; i++)
	{
		blocklens[i] = 1;
		indices[i] = i*l*sizeof(int);
		oldtypes[i] = MPI_INT;
	}
	MPI_Type_create_struct(count, blocklens, indices, oldtypes, &StructTypeColomn);
	MPI_Type_commit(&StructTypeColomn);
 
	dataSize = n;
	bufferSize = dataSize;
	int data_slice = m / (thread_count - 1);
	printf_s("data_slice=%d\n", data_slice);
	int delta_data = m - (thread_count - 1)*data_slice;
	printf_s("delta_data=%d\n", delta_data);

	Abuff = new int[bufferSize];
	Bbuff = new int[bufferSize];

	if (rank == 0) {
		Amatrix = CreateAndFillMatrixVector(m, n);
		Bmatrix = CreateAndFillMatrixVector(n, l);
		if (m*n <= 100 && n*l<=100) {
			printf_s("Matrix A:\n");
			PrintMatrixVector(Amatrix, m, n);
			printf_s("\n");
			printf_s("Matrix B:\n");
			PrintMatrixVector(Bmatrix, n, l);
			printf_s("\n");
		}
		Cmatrix = ConsistentMatrixMultiplication(Amatrix, Bmatrix, m, n, l);
		time1 = MPI_Wtime();
		int *temp_start_Amatrix = Amatrix; //рассылка строк
		int count = 0;
		for (count = 0; count < data_slice; count++) {
			for (i = 1; i < thread_count; i++) {
				MPI_Send(temp_start_Amatrix, n, MPI_INT, i, 0, MPI_COMM_WORLD);
				temp_start_Amatrix = temp_start_Amatrix + n;
				MPI_Send(Bmatrix, 1, StructTypeColomn, i, 0, MPI_COMM_WORLD);
			}
			for (i = 1; i < l; i++) {
				for (j = 1; j < thread_count; j++) {
					MPI_Send(Bmatrix + i, 1, StructTypeColomn, j, 0, MPI_COMM_WORLD);
				}
			}
			for (i = 1; i < thread_count; i++) {
				MPI_Recv(result_matrix + count*l*(thread_count - 1) + l*(i - 1), l, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
			}
		}
		if (delta_data!=0)
		{
			for (i = 1; i <= delta_data; i++) {
				MPI_Send(temp_start_Amatrix, n, MPI_INT, i, 0, MPI_COMM_WORLD);
				temp_start_Amatrix = temp_start_Amatrix + n;
				MPI_Send(Bmatrix, 1, StructTypeColomn, i, 0, MPI_COMM_WORLD);
			}
			for (i = 1; i < l; i++) {
				for (j = 1; j <= delta_data; j++) {
					MPI_Send(Bmatrix + i, 1, StructTypeColomn, j, 0, MPI_COMM_WORLD);
				}
			}
			for (i = 1; i <= delta_data; i++) {
				MPI_Recv(result_matrix + count*l*(thread_count - 1) + l*(i - 1), l, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
			}
		}


		PrintMatrixVector(result_matrix, m, l);
		if (CheckResult(result_matrix, Cmatrix, m, l))
		{
			printf_s("Correct.\n");
		}
		else
		{
			printf_s("Error.\n");
		}

	}
	if(rank != 0) {
		for (int count = 0; count < data_slice; count++) {
			MPI_Recv(Abuff, bufferSize, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(Bbuff, n, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			//PrintVector(Bbuff, n);
			result_Multiplication[0] = ScalarMultiplication(Abuff, Bbuff, n);
			for (i = 1; i < l; i++) {
				MPI_Recv(Bbuff, n, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				//PrintVector(Bbuff, n);
				result_Multiplication[i] = ScalarMultiplication(Abuff, Bbuff, n);
			}
			//printf_s("rank = %d, c1 = %d, c2 = %d, c3 = %d, c4 = %d\n", rank, result_Multiplication[0], result_Multiplication[1], result_Multiplication[2], result_Multiplication[3]);
			MPI_Send(result_Multiplication, l, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
		if ((delta_data != 0) && (rank<=delta_data)) {
			MPI_Recv(Abuff, bufferSize, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(Bbuff, n, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			//PrintVector(Bbuff, n);
			result_Multiplication[0] = ScalarMultiplication(Abuff, Bbuff, n);
			for (i = 1; i < l; i++) {
				MPI_Recv(Bbuff, n, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				//PrintVector(Bbuff, n);
				result_Multiplication[i] = ScalarMultiplication(Abuff, Bbuff, n);
			}
			//printf_s("rank = %d, c1 = %d, c2 = %d, c3 = %d, c4 = %d\n", rank, result_Multiplication[0], result_Multiplication[1], result_Multiplication[2], result_Multiplication[3]);
			MPI_Send(result_Multiplication, l, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
	}

	MPI_Type_free(&StructTypeColomn);
	MPI_Finalize();
	return 0;
}
