#include <cstdio>
#include "mpi.h"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

int algebraic_complement(int** data) {
	return 0;
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

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	int** matrix;

	if (rank == 0) {
		ifstream matrix_data("matrix_data.txt");
		matrix_data >> n;
		matrix = new int* [n];

		for (int i = 0; i < n; i++) {
			matrix[i] = new int[n];
			for (int j = 0; j < n; j++) {
				int temp;
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
				MPI_Send(matrix[k], n, MPI_INT, p, 99, MPI_COMM_WORLD);
			}
		}

		for (int i = (rank + 0) * (n / numtasks + n % numtasks); i < (rank + 1) * (n / numtasks + n % numtasks) && i < n; i++) {
			cout << "Call " << i << "For rank " << rank << endl;
		}

		for (int i = 0; i < n; i++)
			delete [] matrix[i];
		delete[] matrix;
	}
	else {
		MPI_Recv(&n, 1, MPI_INT, 0, 99, MPI_COMM_WORLD, &st);
		cout << "Rank " << rank << "has element " << n << endl;
		matrix = new int* [n];
		for (int i = 0; i < n; i++) {
			matrix[i] = new int[n];
			MPI_Recv(matrix[i], n, MPI_INT, 0, 99, MPI_COMM_WORLD, &st);
		}

		cout << "Rank " << rank << "has matrix[rank][rank] element " << matrix[rank][rank] << endl;
		for (int i = (rank + 0) * (n / numtasks + n % numtasks); i < (rank + 1) * (n / numtasks + n % numtasks) && i < n; i++) {
			cout << "Call " << i << "For rank " << rank << endl;
		}

		for (int i = 0; i < n; i++)
			delete[] matrix[i];
		delete[] matrix;
	}

	//printf("MPItest from process - %d, total number of process - %d\n", rank, numtasks);
	MPI_Finalize();
}