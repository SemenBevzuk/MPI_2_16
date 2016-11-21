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
	int data_slice_row = m / (thread_count - 1);
	int delta_data_row = m - (thread_count - 1)*data_slice_row;
	int data_slice_colomn = l / (thread_count - 1);
	int delta_data_colomn = l - (thread_count - 1)*data_slice_colomn;
	
	dataSize = n;
	//bufferSize = dataSize*data_slice_row;
	result_Multiplication = (int*)malloc((thread_count-1)*data_slice_colomn * data_slice_row * sizeof(int));
	Abuff = new int[dataSize*data_slice_row];
	Bbuff = new int[dataSize*data_slice_colomn];
	
	//структура для нескольких столбцов
	int count_colomn = data_slice_colomn;
	MPI_Datatype StructTypeNColomn;
	int *blocklens_colomn = new int[data_slice_colomn];
	MPI_Aint *indices_colomn = new MPI_Aint[data_slice_colomn];
	MPI_Datatype* oldtypes_colomn = new MPI_Datatype[data_slice_colomn];
	for (i = 0; i<data_slice_colomn; i++) {
		blocklens_colomn[i] = 1;
		indices_colomn[i] = i*sizeof(int);
		oldtypes_colomn[i] = StructTypeColomn;
	}
	MPI_Type_create_struct(count_colomn, blocklens_colomn, indices_colomn, oldtypes_colomn, &StructTypeNColomn);
	MPI_Type_commit(&StructTypeNColomn);

	if (rank == 0) {
		Amatrix = CreateAndFillMatrixVector(m, n);
		Bmatrix = CreateAndFillMatrixVector(n, l);
		Cmatrix = (int*)malloc(m*l* sizeof(int));

		delta_time_consistent = StartConsistentMatrixMultiplication(Amatrix, Bmatrix, &Cmatrix, m, n, l);

		time1 = MPI_Wtime();

		int *temp_start_Amatrix = Amatrix;
		for (i = 1; i < thread_count; i++) {
			MPI_Send(temp_start_Amatrix, n*data_slice_row, MPI_INT, i, 0, MPI_COMM_WORLD);
			temp_start_Amatrix = temp_start_Amatrix + n*data_slice_row;
			MPI_Send(Bmatrix, 1, StructTypeNColomn, i, 0, MPI_COMM_WORLD);
		}
		for (i = 1; i < thread_count-1; i++) {
			for (j = 1; j < thread_count; j++) {
				MPI_Send(Bmatrix + i*data_slice_colomn, 1, StructTypeNColomn, j, 0, MPI_COMM_WORLD);
			}
		}

		int* temp_result_matrix = result_matrix;
		for (i = 1; i < thread_count; i++) {
			for (j = 0; j < data_slice_row; j++)
			{
				MPI_Recv(temp_result_matrix, l-delta_data_colomn, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
				//PrintVector(temp_result_matrix, l - delta_data_colomn);
				temp_result_matrix = temp_result_matrix + l;
			}	
		}
		//обрабоать остатки
		int res = 0;
		//PrintMatrixVector(result_matrix, m, l);
		if (delta_data_colomn!=0)
		{
			for (i = 0; i < m; i++) {
				for (k = l - delta_data_colomn; k < l; k++) {
					res = 0;
					for (j = 0; j < n; j++) {
						res += Amatrix[n*i + j] * Bmatrix[k + j*l];
					}

					result_matrix[k + l*i] = res;
				}
			}
		}
		//PrintMatrixVector(result_matrix, m, l);
		if (delta_data_row!=0)
		{
			for (i = m - delta_data_row; i < m; i++) {
				for (j = 0; j < l; j++) {
					res = 0;
					for (k = 0; k < n; k++) {
						res += Amatrix[n*i + k] * Bmatrix[j + k*l];
					}
					//printf_s("res = %d\n", res);
					result_matrix[i*l + j] = res;
				}
			}
		}
		//printf_s("res: \n");
		//PrintMatrixVector(result_matrix, m, l);
		//printf_s("l - delta_data_colomn = %d, l = %d, delta_data_colomn = %d\n", l - delta_data_colomn, l, delta_data_colomn);
		time2 = MPI_Wtime();
		delta_time_parallel = time2 - time1;
	}

	if(rank != 0) {
		MPI_Recv(Abuff, dataSize*data_slice_row, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		//PrintVector(Abuff, dataSize*data_slice_row);
		//PrintVector(result_Multiplication, (thread_count - 1)*data_slice_colomn * data_slice_row);
		for (k = 0; k < thread_count-1; k++)
			{
				MPI_Recv(Bbuff, n*data_slice_colomn, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
					for (i = 0; i < data_slice_row; i++) {
						for (j = 0; j < data_slice_colomn; j++) {
							result_Multiplication[k*data_slice_colomn + i*(l - delta_data_colomn) + j] = ScalarMultiplication(Abuff + n*i, Bbuff + j*n, n);
							/*if (rank == 1)
							{
								printf_s("A: \n");
								PrintVector(Abuff + n*i, n);
								printf_s("B: \n");
								PrintVector(Bbuff + j*n, n);
								printf_s("res = %d: in %d\n", result_Multiplication[k*data_slice_colomn + i*(l - 1) + j], k*data_slice_colomn + i*(l - 1) + j);
							}*/
						}
					}
			}
		//PrintVector(result_Multiplication, (thread_count - 1)*data_slice_colomn * data_slice_row);
		int* temp = result_Multiplication;
		for (i = 0; i < data_slice_row; i++)
		{
			MPI_Send(temp, l-delta_data_colomn, MPI_INT, 0, 0, MPI_COMM_WORLD);
			temp = temp + (l - delta_data_colomn);
		}
		//PrintVector(result_Multiplication, (thread_count - 1)*data_slice_colomn * data_slice_row);
		//MPI_Send(result_Multiplication, (thread_count - 1)*data_slice_colomn * data_slice_row, MPI_INT, 0, 0, MPI_COMM_WORLD);
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
			//printf_s("Matrix C = A * B:\n");
			//PrintMatrixVector(Cmatrix, m, l);
			//printf_s("\n");
		}
		if (CheckResult(result_matrix, Cmatrix, m, l)) {
			printf_s("Correct.\n");
		}
		else {
			printf_s("Error.\n");
		}
		printf_s("\ntime_parallel = %f\n", delta_time_parallel);
		printf_s("time_consistent = %f\n", delta_time_consistent);
		printf_s("Acceleration(parallel):  %f\n", (delta_time_consistent / delta_time_parallel));
	}

	MPI_Type_free(&StructTypeColomn);
	MPI_Type_free(&StructTypeNColomn);
	MPI_Finalize();
	DeleteMatrix(Amatrix);
	DeleteMatrix(Bmatrix);
	DeleteMatrix(Cmatrix);
	return 0;
}
