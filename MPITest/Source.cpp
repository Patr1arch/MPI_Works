#include <cstdio>
#include "mpi.h"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

int main(int* argc, char** argv) {
	int numtasks, rank;
	//int* n;
	//printf("Enter the order of the matrix\n");
	//scanf("d", &n);
	//printf("Enter %d strings of %d elemetns of matrix", n, n);
	//vector<vector<int>> *matrix;
	//for (int i = 0; i < *n; i++) {

	//}

	//printf(argv[2]);

	MPI_Init(argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	int ** matrix;

	if (rank == 0) {
		ifstream matrix_data("matrix_data.txt");
		int n;
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

		for (int i = 0; i < n; i++)
			delete [] matrix[i];
		delete[] matrix;
	}

	//printf("MPItest from process - %d, total number of process - %d\n", rank, numtasks);
	MPI_Finalize();
}