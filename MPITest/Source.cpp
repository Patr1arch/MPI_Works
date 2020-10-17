#include <cstdio>
#include "mpi.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h> 

using namespace std;

double algebraic_complement(double** data, int i_ex, int j_ex, int dim) {
	if (dim == 3) {
		double sum = 
			(data[0][0] * data[1][1] * data[2][2] + data[0][1] * data[1][2] * data[2][0]
			+ data[1][0] * data[2][1] * data[0][2] - data[0][2] * data[1][1] * data[2][0]
			- data[2][1] * data[1][2] * data[0][0] - data[1][0] * data[0][1] * data[2][2]);
		//cout << sum << " " << endl;
		return sum;
	}
	else if (dim == 2) {
		return data[0][0] * data[1][1] - data[1][0] * data[0][1];
	}
	else if (dim == 1) {
		return data[0][0];
	}

	double** minor = new double* [dim - 1];
	for (int i = 0; i < dim - 1; i++) {
		minor[i] = new double[dim - 1];
	}

	//for (int i = 0; i < dim; i++) {
	//	for (int j = 0; j < dim; j++) {
	//		cout << data[i][j] << " ";
	//	}
	//	cout << endl;
	//}
	//cout << endl;

	for (int i = 0, row_i = 0; i < dim; i++) {
		if (i == i_ex) continue;
		for (int j = 0, col_j = 0; j < dim; j++) {
			if (j == j_ex) continue;
			minor[row_i][col_j++] = data[i][j];
			//cout << data[i][j] << " ";
		}
		row_i++;
		//cout << endl;
	}
	//cout << endl;

	//for (int i = 0; i < dim - 1; i++) {
	//	for (int j = 0; j < dim - 1; j++) {
	//		cout << minor[i][j] << " ";
	//	}
	//	cout << endl;
	//}
	//cout << endl;

	double sum = 0;
	if (dim - 1 == 3)
		sum = (minor[0][0] * minor[1][1] * minor[2][2] + minor[0][1] * minor[1][2] * minor[2][0]
			+ minor[1][0] * minor[2][1] * minor[0][2] - minor[0][2] * minor[1][1] * minor[2][0]
			- minor[2][1] * minor[1][2] * minor[0][0] - minor[1][0] * minor[0][1] * minor[2][2]);
	else {
		for (int j = 0; j < dim - 1; j++) {
			sum += minor[0][j] * algebraic_complement(minor, 0, j, dim - 1);
		}
	}

	//cout << "Sum = " << sum << endl;

	for (int i = 0; i < dim - 1; i++)
		delete[] minor[i];
	delete[] minor;

	//cout << i_ex + j_ex << endl;

	//cout << (int)pow(-1, i_ex + j_ex) << endl;

	//cout << sum * (int)pow(-1, i_ex + j_ex) << endl;

	return sum * pow(-1, i_ex + j_ex);
}

int main(int* argc, char** argv) {
	int numtasks, rank, n;
	//int* n;
	//printf("Enter the order of the matrix\n");
	//scanf("d", &n);
	//printf("Enter %d strings of %d elemetns of matrix", n, n);
	//vector<vector<int>> *matrix;
	//for (int i = 0; i < *n; i++) {

	//}

	//printf(argv[2]);

	MPI_Init(argc, &argv);
	MPI_Status st;
	double start, end;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	double** matrix;
	double rank_res[1];
	rank_res[0] = 0.0;

	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();

	if (rank == 0) {
		//clock_t t;
		//t = clock();
		ifstream matrix_data("matrix_data.txt");
		matrix_data >> n;
		matrix = new double* [n];

		for (int i = 0; i < n; i++) {
			matrix[i] = new double[n];
			for (int j = 0; j < n; j++) {
				double temp;
				matrix_data >> temp;
				matrix[i][j] = temp;
			}
		}

		matrix_data.close();
		// rank * (n / numtasks + n % numtasks)
		//cout << (n / numtasks + n % numtasks) << endl;
		//cout << numtasks << endl;
		//cout << n << endl;
		for (int p = 1; p < numtasks; p++)
		{
			MPI_Send(&n, 1, MPI_INT, p, 99, MPI_COMM_WORLD);
			for (int k = 0; k < n; k++) {
				MPI_Send(matrix[k], n, MPI_DOUBLE, p, 99, MPI_COMM_WORLD);
			}
		}

		double sum = 0;
			//int i = (rank + 0) * (n / numtasks + n % numtasks); i < (rank + 1)* (n / numtasks + n % numtasks) && i < n; i++
		for (int i = 0; i < n; i += numtasks) {
			//cout << "Call " << i << " for rank " << rank << endl;
			sum += matrix[0][i] * algebraic_complement(matrix, 0, i, n);
		}

		for (int p = 1; p < numtasks; p++) {
			MPI_Recv(&rank_res, 1, MPI_DOUBLE, p, 99, MPI_COMM_WORLD, &st);
			sum += rank_res[0];
		}

		//for (int j = 0; j < n; j++) {
		//	auto temp = algebraic_complement(matrix, 0, j, n);
		//	sum += matrix[0][j] * temp;
		//	//cout << "matrix[0][" << j << "] = " << matrix[0][j] << endl;
		//	//cout << "algebraic_complement(matrix, 0," << j << ", " << n << ") = " << temp << endl;;
		//}

		cout << "RESULT: " << sum << endl;

		for (int i = 0; i < n; i++)
			delete [] matrix[i];
		delete[] matrix;
		//t = clock() - t;
		//cout << "It tooks " << (double)t / CLOCKS_PER_SEC << " seconds" << endl;
	}
	else {
		MPI_Recv(&n, 1, MPI_INT, 0, 99, MPI_COMM_WORLD, &st);
		//cout << "Rank " << rank << "has element " << n << endl;
		matrix = new double* [n];
		for (int i = 0; i < n; i++) {
			matrix[i] = new double[n];
			MPI_Recv(matrix[i], n, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD, &st);
		}

		rank_res[0] = 0;
		//cout << "Rank " << rank << "has matrix[rank][rank] element " << matrix[rank][rank] << endl;
		for (int i = rank; i < n; i += numtasks) {
			//cout << "Call " << i << " for rank " << rank << endl;
			rank_res[0] += matrix[0][i] * algebraic_complement(matrix, 0, i, n);
		}

		MPI_Send(&rank_res, 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD);

		for (int i = 0; i < n; i++)
			delete[] matrix[i];
		delete[] matrix;
	}

	MPI_Barrier(MPI_COMM_WORLD); 
	end = MPI_Wtime();

	//printf("MPItest from process - %d, total number of process - %d\n", rank, numtasks);
	MPI_Finalize();
	if (rank == 1) {
		cout << "It took: " << end - start << endl;
	}
}