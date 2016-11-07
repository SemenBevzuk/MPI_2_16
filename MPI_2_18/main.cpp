#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include <cstddef>
#include <random>

#define EPS 0.0001

int* CreateAndFillMatrixVector(int n, int m) {
	int i;
	int* matrix = (int*)malloc(n * m * sizeof(int));
	for (i = 0; i < n*m; i++) {
		matrix[i] = 0 + rand() % 10;
	}
	return matrix;
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
double StartConsistentMatrixMultiplication(int* Amatrix, int* Bmatrix, int **Cmatrix, int m, int n, int l)
{
	double time1, time2;
	time1 = MPI_Wtime();
	*Cmatrix = ConsistentMatrixMultiplication(Amatrix, Bmatrix, m, n, l);
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
	int* result_Multiplication = NULL;
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
	int* result_matrix = (int*)malloc(m*l* sizeof(int)); //матрица для параллельной реализации

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &thread_count);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// структура для столбцов
	int count = n;
	MPI_Datatype StructTypeColomn;
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
 
	//размеры буффера для данных
	int data_slice = m / (thread_count - 1);
	int delta_data = m - (thread_count - 1)*data_slice;
	int* Elements_in_thread = new int[thread_count];
	for (i = 1; i < thread_count - 1; i++)
	{
		Elements_in_thread[i] = data_slice;
	}
	Elements_in_thread[0] = 0;
	Elements_in_thread[thread_count - 1] = data_slice + delta_data;
	dataSize = n;
	bufferSize = dataSize;
	if (rank != 0)
	{
		bufferSize = dataSize*Elements_in_thread[rank];
		result_Multiplication = (int*)malloc(Elements_in_thread[rank] * l * sizeof(int));
	}
	Abuff = new int[bufferSize];
	Bbuff = new int[dataSize];
	

	if (rank == 0) {
		Amatrix = CreateAndFillMatrixVector(m, n);
		Bmatrix = CreateAndFillMatrixVector(n, l);
		Cmatrix = (int*)malloc(m*l* sizeof(int));

		delta_time_consistent = StartConsistentMatrixMultiplication(Amatrix, Bmatrix, &Cmatrix, m, n, l);


		time1 = MPI_Wtime();

		int *temp_start_Amatrix = Amatrix;
		for (i = 1; i < thread_count; i++) {
			MPI_Send(temp_start_Amatrix, n*Elements_in_thread[i], MPI_INT, i, 0, MPI_COMM_WORLD);
			temp_start_Amatrix = temp_start_Amatrix + n*Elements_in_thread[i];
			MPI_Send(Bmatrix, 1, StructTypeColomn, i, 0, MPI_COMM_WORLD);
		}
		for (i = 1; i < l; i++) {
			for (j = 1; j < thread_count; j++) {
				MPI_Send(Bmatrix + i, 1, StructTypeColomn, j, 0, MPI_COMM_WORLD);
			}
		}

		int* temp_resukt_matrix = result_matrix;
		for (i = 1; i < thread_count; i++) {
			MPI_Recv(temp_resukt_matrix, Elements_in_thread[i] * l, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
			temp_resukt_matrix = temp_resukt_matrix + Elements_in_thread[i] * l;
		}

		time2 = MPI_Wtime();
		delta_time_parallel = time2 - time1;
	}

	if(rank != 0) {
			MPI_Recv(Abuff, bufferSize, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(Bbuff, n, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			for (i = 0; i < Elements_in_thread[rank]; i++)
			{
				result_Multiplication[0+l*i] = ScalarMultiplication(Abuff + n*i, Bbuff, n);
			}

			for (i = 1; i < l; i++) {
				MPI_Recv(Bbuff, n, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				for (j = 0; j < Elements_in_thread[rank]; j++) {
					result_Multiplication[i + l*j] = ScalarMultiplication(Abuff + n*j, Bbuff, n);
				}
			}
			MPI_Send(result_Multiplication, Elements_in_thread[rank] * l, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}

	if (rank == 0)
	{
		if (m*n <= 50 && n*l <= 50) {
			printf_s("Matrix A:\n");
			PrintMatrixVector(Amatrix, m, n);
			printf_s("\n");
			printf_s("Matrix B:\n");
			PrintMatrixVector(Bmatrix, n, l);
			printf_s("\n");
			printf_s("Matrix C = A * B:\n");
			PrintMatrixVector(Cmatrix, m, l);
			printf_s("\n");
		}
		if (CheckResult(result_matrix, Cmatrix, m, l)) {
			printf_s("Correct.\n");
		}
		else {
			printf_s("Error.\n");
		}
		printf_s("time_parallel = %f\n", delta_time_parallel);
		printf_s("time_consistent = %f\n", delta_time_consistent);
		printf_s("Acceleration(parallel):  %f\n", (delta_time_consistent / delta_time_parallel));
	}

	MPI_Type_free(&StructTypeColomn);
	MPI_Finalize();
	DeleteMatrix(Amatrix);
	DeleteMatrix(Bmatrix);
	DeleteMatrix(Cmatrix);
	return 0;
}
