
#include "lab2_functions.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

// Sequential Experiment
void experiment(int m) {

	int n = scanner("Please enter the power");

	FILE *fp = fopen("lab2.txt", "w");
	clock_t start, end;
	double time_taken;
	time_t start_t, end_t;
	struct timeval start_g, end_g;

	// Row Wise ====================================================
	printf("Starting Row Wise\n");
	fprintf(fp, "Constant Power\nRow Wise openMP: Clock()\n");
	start = clock(); // Start time
	matrixPowerParallelRowWise(m, n);
	end = clock(); // End time
	time_taken = (double)(end - start) / CLOCKS_PER_SEC;
	fprintf(fp, "Data from m, %d, n, %d\n", m, n);
	fprintf(fp, "Time taken, %f, seconds\n", time_taken);

	fprintf(fp, "\n");
	fprintf(fp, "Row Wise openMP: time()\n");
	time(&start_t); // Start time
	matrixPowerParallelRowWise(m, n);
	time(&end_t); // End time
	time_taken = difftime(end_t, start_t);
	fprintf(fp, "Data from m, %d, n, %d\n", m, n);
	fprintf(fp, "Time taken, %f, seconds\n", time_taken);

	fprintf(fp, "\n");
	fprintf(fp, "Row Wise openMP: gettimeofday()\n");
	gettimeofday(&start_g, NULL); // Start time
	matrixPowerParallelRowWise(m, n);
	gettimeofday(&end_g, NULL); // End time
	time_taken = (end_g.tv_sec + end_g.tv_usec / 1000000) - (start_g.tv_sec + start_g.tv_usec / 1000000);
	fprintf(fp, "Data from m, %d, n, %d\n", m, n);
	fprintf(fp, "Time taken, %f, seconds\n", time_taken);

	// Column Wise ====================================================
	printf("Starting Column Wise\n");
	fprintf(fp, "Constant Power\nColumn Wise openMP: Clock()\n");
	start = clock(); // Start time
	matrixPowerParallelColmWise(m, n);
	end = clock(); // End time
	time_taken = (double)(end - start) / CLOCKS_PER_SEC;
	fprintf(fp, "Data from m, %d, n, %d\n", m, n);
	fprintf(fp, "Time taken, %f, seconds\n", time_taken);

	fprintf(fp, "\n");
	fprintf(fp, "Row Wise openMP: time()\n");
	time(&start_t); // Start time
	matrixPowerParallelColmWise(m, n);
	time(&end_t); // End time
	time_taken = difftime(end_t, start_t);
	fprintf(fp, "Data from m, %d, n, %d\n", m, n);
	fprintf(fp, "Time taken, %f, seconds\n", time_taken);
	
	fprintf(fp, "\n");
	fprintf(fp, "Row Wise openMP: gettimeofday()\n");
	gettimeofday(&start_g, NULL); // Start time
	matrixPowerParallelColmWise(m, n);
	gettimeofday(&end_g, NULL); // End time
	time_taken = (end_g.tv_sec + end_g.tv_usec / 1000000) - (start_g.tv_sec + start_g.tv_usec / 1000000);
	fprintf(fp, "Data from m, %d, n, %d\n", m, n);
	fprintf(fp, "Time taken, %f, seconds\n", time_taken);

	// Element Wise ====================================================
	printf("Starting Element Wise\n");
	fprintf(fp, "Constant Power\nElement Wise openMP: Clock()\n");
	start = clock(); // Start time
	matrixPowerParallelElementWise(m, n);
	end = clock(); // End time
	time_taken = (double)(end - start) / CLOCKS_PER_SEC;
	fprintf(fp, "Data from m, %d, n, %d\n", m, n);
	fprintf(fp, "Time taken, %f, seconds\n", time_taken);

	fprintf(fp, "Element Wise openMP: time()\n");
	time(&start_t); // Start time
	matrixPowerParallelElementWise(m, n);
	time(&end_t); // End time
	time_taken = difftime(end_t, start_t);
	fprintf(fp, "Data from m, %d, n, %d\n", m, n);
	fprintf(fp, "Time taken, %f, seconds\n", time_taken);
 
	fprintf(fp, "Element Wise openMP: gettimeofday()\n");
	gettimeofday(&start_g, NULL); // Start time
	matrixPowerParallelElementWise(m, n);
	gettimeofday(&end_g, NULL); // End time
	time_taken = (end_g.tv_sec + end_g.tv_usec / 1000000) - (start_g.tv_sec + start_g.tv_usec / 1000000);
	fprintf(fp, "Data from m, %d, n, %d\n", m, n);
	fprintf(fp, "Time taken, %f, seconds\n", time_taken);
}

// Parallel Experiments
// Performs all of the experiments in order automatically
void openMPExperiment() {

	//int n = scanner("Please enter the power");
	int n, m;
	FILE *fp = fopen("lab2.txt", "w");
	clock_t start, end;
	double time_taken;
	time_t start_t, end_t;
	struct timeval start_g, end_g;

	// ======================================================================================================
	// Row Wise
	// ======================================================================================================
	// -------------------------------------------------------------
	// Constant Power
	// -------------------------------------------------------------
	n = 100;
	printf("Starting Row Wise Constant Power\n");
	fprintf(fp, "Constant Power\nRow Wise openMP: Clock()\n");
	for (m = 100; m <= 1000; m += 100) {
		start = clock(); // Start time
		matrixPowerParallelRowWise(m, n);
		end = clock(); // End time
		time_taken = (double)(end - start) / CLOCKS_PER_SEC;

		fprintf(fp, "Data from m, %d, n, %d\n", m, n);
		fprintf(fp, "Time taken, %f, seconds\n", time_taken);
	} 
	fprintf(fp, "\n");
	fprintf(fp, "Row Wise openMP: time()\n");
	for (m = 100; m <= 1000; m += 100) {
		time(&start_t); // Start time
		matrixPowerParallelRowWise(m, n);
		time(&end_t); // End time
		time_taken = difftime(end_t, start_t);

		fprintf(fp, "Data from m, %d, n, %d\n", m, n);
		fprintf(fp, "Time taken, %f, seconds\n", time_taken);
	} 
	fprintf(fp, "\n");
	fprintf(fp, "Row Wise openMP: gettimeofday()\n");
	for (m = 100; m <= 1000; m += 100) {
		gettimeofday(&start_g, NULL); // Start time
		matrixPowerParallelRowWise(m, n);
		gettimeofday(&end_g, NULL); // End time
		time_taken = (end_g.tv_sec + end_g.tv_usec / 1000000) - (start_g.tv_sec + start_g.tv_usec / 1000000);

		fprintf(fp, "Data from m, %d, n, %d\n", m, n);
		fprintf(fp, "Time taken, %f, seconds\n", time_taken);
	} 
	
	// -------------------------------------------------------------
	// Constant Size
	// -------------------------------------------------------------
	m = 1000;
	printf("Starting Row Wise Constant Size\n");
	fprintf(fp, "Constant Size\nRow Wise openMP: Clock()\n");
	for (n = 100; n <= 500; n += 100) {
		start = clock(); // Start time
		matrixPowerParallelRowWise(m, n);
		end = clock(); // End time
		time_taken = (double)(end - start) / CLOCKS_PER_SEC;

		fprintf(fp, "Data from m, %d, n, %d\n", m, n);
		fprintf(fp, "Time taken, %f, seconds\n", time_taken);
	} 
	fprintf(fp, "\n");
	fprintf(fp, "Row Wise openMP: time()\n");
	for (n = 100; n <= 500; n += 100) {
		time(&start_t); // Start time
		matrixPowerParallelRowWise(m, n);
		time(&end_t); // End time
		time_taken = difftime(end_t, start_t);

		fprintf(fp, "Data from m, %d, n, %d\n", m, n);
		fprintf(fp, "Time taken, %f, seconds\n", time_taken);
	} 
	fprintf(fp, "\n");
	fprintf(fp, "Row Wise openMP: gettimeofday()\n");
	for (n = 100; n <= 500; n += 100) {
		gettimeofday(&start_g, NULL); // Start time
		matrixPowerParallelRowWise(m, n);
		gettimeofday(&end_g, NULL); // End time
		time_taken = (end_g.tv_sec + end_g.tv_usec / 1000000) - (start_g.tv_sec + start_g.tv_usec / 1000000);

		fprintf(fp, "Data from m, %d, n, %d\n", m, n);
		fprintf(fp, "Time taken, %f, seconds\n", time_taken);
	} 

	// ======================================================================================================
	// Column Wise
	// ======================================================================================================
	// -------------------------------------------------------------
	// Constant Power
	// -------------------------------------------------------------
	n = 100;
	printf("Starting Column Wise Constant Power\n");
	fprintf(fp, "Constant Power\nColumn Wise openMP: Clock()\n");
	for (m = 100; m <= 1000; m += 100) {
		start = clock(); // Start time
		matrixPowerParallelColmWise(m, n);
		end = clock(); // End time
		time_taken = (double)(end - start) / CLOCKS_PER_SEC;

		fprintf(fp, "Data from m, %d, n, %d\n", m, n);
		fprintf(fp, "Time taken, %f, seconds\n", time_taken);
	} 

	fprintf(fp, "Column Wise openMP: time()\n");
	for (m = 100; m <= 1000; m += 100) {
		time(&start_t); // Start time
		matrixPowerParallelColmWise(m, n);
		time(&end_t); // End time
		time_taken = difftime(end_t, start_t);

		fprintf(fp, "Data from m, %d, n, %d\n", m, n);
		fprintf(fp, "Time taken, %f, seconds\n", time_taken);
	} 

	fprintf(fp, "Column Wise openMP: gettimeofday()\n");
	for (m = 100; m <= 1000; m += 100) {
		gettimeofday(&start_g, NULL); // Start time
		matrixPowerParallelColmWise(m, n);
		gettimeofday(&end_g, NULL); // End time
		time_taken = (end_g.tv_sec + end_g.tv_usec / 1000000) - (start_g.tv_sec + start_g.tv_usec / 1000000);

		fprintf(fp, "Data from m, %d, n, %d\n", m, n);
		fprintf(fp, "Time taken, %f, seconds\n", time_taken);
	} 
	// -------------------------------------------------------------
	// Constant Size
	// -------------------------------------------------------------
	m = 1000;
	printf("Starting Column Wise Constant Size\n");
	fprintf(fp, "Constant Size\nColumn Wise openMP: Clock()\n");
	for (n = 100; n <= 500; n += 100) {
		start = clock(); // Start time
		matrixPowerParallelColmWise(m, n);
		end = clock(); // End time
		time_taken = (double)(end - start) / CLOCKS_PER_SEC;

		fprintf(fp, "Data from m, %d, n, %d\n", m, n);
		fprintf(fp, "Time taken, %f, seconds\n", time_taken);
	} 

	fprintf(fp, "Column Wise openMP: time()\n");
	for (n = 100; n <= 500; n += 100) {
		time(&start_t); // Start time
		matrixPowerParallelColmWise(m, n);
		time(&end_t); // End time
		time_taken = difftime(end_t, start_t);

		fprintf(fp, "Data from m, %d, n, %d\n", m, n);
		fprintf(fp, "Time taken, %f, seconds\n", time_taken);
	} 

	fprintf(fp, "Column Wise openMP: gettimeofday()\n");
	for (n = 100; n <= 500; n += 100) {
		gettimeofday(&start_g, NULL); // Start time
		matrixPowerParallelColmWise(m, n);
		gettimeofday(&end_g, NULL); // End time
		time_taken = (end_g.tv_sec + end_g.tv_usec / 1000000) - (start_g.tv_sec + start_g.tv_usec / 1000000);

		fprintf(fp, "Data from m, %d, n, %d\n", m, n);
		fprintf(fp, "Time taken, %f, seconds\n", time_taken);
	} 
	// ======================================================================================================
	// Element Wise
	// ======================================================================================================
	// -------------------------------------------------------------
	// Constant Power
	// -------------------------------------------------------------
	n = 100;
	printf("Starting Element Wise Constant Power\n");
	fprintf(fp, "Constant Power\nElement Wise openMP: Clock()\n");
	for (m = 100; m <= 1000; m += 100) {
		start = clock(); // Start time
		matrixPowerParallelElementWise(m, n);
		end = clock(); // End time
		time_taken = (double)(end - start) / CLOCKS_PER_SEC;

		fprintf(fp, "Data from m, %d, n, %d\n", m, n);
		fprintf(fp, "Time taken, %f, seconds\n", time_taken);
	} 

	fprintf(fp, "Element Wise openMP: time()\n");
	for (m = 100; m <= 1000; m += 100) {
		time(&start_t); // Start time
		matrixPowerParallelElementWise(m, n);
		time(&end_t); // End time
		time_taken = difftime(end_t, start_t);

		fprintf(fp, "Data from m, %d, n, %d\n", m, n);
		fprintf(fp, "Time taken, %f, seconds\n", time_taken);
	} 

	fprintf(fp, "Element Wise openMP: gettimeofday()\n");
	for (m = 100; m <= 1000; m += 100) {
		gettimeofday(&start_g, NULL); // Start time
		matrixPowerParallelElementWise(m, n);
		gettimeofday(&end_g, NULL); // End time
		time_taken = (end_g.tv_sec + end_g.tv_usec / 1000000) - (start_g.tv_sec + start_g.tv_usec / 1000000);

		fprintf(fp, "Data from m, %d, n, %d\n", m, n);
		fprintf(fp, "Time taken, %f, seconds\n", time_taken);
	} 
	// -------------------------------------------------------------
	// Constant Size
	// -------------------------------------------------------------
	m = 1000;

	// Experiment with each element as a thread using clock timing function
	printf("Starting Element Wise Constant Size\n");
	fprintf(fp, "Constant Size\nElement Wise openMP: Clock()\n");
	for (n = 100; n <= 500; n += 100) {
		start = clock(); // Start time
		matrixPowerParallelElementWise(m, n);
		end = clock(); // End time
		time_taken = (double)(end - start) / CLOCKS_PER_SEC;

		fprintf(fp, "Data from m, %d, n, %d\n", m, n);
		fprintf(fp, "Time taken, %f, seconds\n", time_taken);
	} 

	// Experiment with each element as a thread using time timing function
	fprintf(fp, "Element Wise openMP: time()\n");
	for (n = 100; n <= 500; n += 100) {
		time(&start_t); // Start time
		matrixPowerParallelElementWise(m, n);
		time(&end_t); // End time
		time_taken = difftime(end_t, start_t);

		fprintf(fp, "Data from m, %d, n, %d\n", m, n);
		fprintf(fp, "Time taken, %f, seconds\n", time_taken);
	} 

	// Experimemnt with each element as a thread using gettimeofday timing function
	fprintf(fp, "Element Wise openMP: gettimeofday()\n");
	for (n = 100; n <= 500; n += 100) {
		gettimeofday(&start_g, NULL); // Start time
		matrixPowerParallelElementWise(m, n);
		gettimeofday(&end_g, NULL); // End time
		time_taken = (end_g.tv_sec + end_g.tv_usec / 1000000) - (start_g.tv_sec + start_g.tv_usec / 1000000);

		fprintf(fp, "Data from m, %d, n, %d\n", m, n);
		fprintf(fp, "Time taken, %f, seconds\n", time_taken);
	} 

	fclose(fp);

}

// Parallel Matrix Multiplication: Each thread is a row
void matrixPowerParallelRowWise(int m, int n) {

	double *matrix1 = (double *)malloc(m * m * sizeof(double));
	double *matrix2 = (double *)malloc(m * m * sizeof(double));
	double *result  = (double *)malloc(m * m * sizeof(double));

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			*(matrix1 + i *m + j) = (double)rand() / (double)RAND_MAX;
		}
	}

	// Deep copying the contents of matrix1 to the result matrix
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			*(result + i *m + j) = *(matrix1 + i *m + j);
		}
	}

	//printMatrix(matrix1, m);
	for (int x = 1; x < n; x++) {
		// Deep copying the contents of result into matrix2
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				*(matrix2 + i *m + j) = *(result + i *m + j);
			}
		}

		// Setting the contents of result to 0
		memset(result, 0.0, sizeof(result));

		int i, j, k;
		// Matrix Multiplication
		#pragma omp parallel for private(j,k) num_threads(m)
		for (i = 0; i < m; i++) {
			for (j = 0; j < m; j++) {
				for (k = 0; k < m; k++) {
					//printf("result[i][j] = %f, matrix2[i][k] = %f, matrix1[k][j] = %f\n", *(result + i *m + j), *(matrix2 + i *m + k), *(matrix1 + k *m + j));
					*(result + i *m + j) += *(matrix2 + i *m + k) * *(matrix1 + k *m + j);
				}
			}
		}
	}

	//printMatrix(result, m);
	free(matrix1);
	free(matrix2);
	free(result);
}

// Parallel Matrix Multiplication: Each thread is a column
void matrixPowerParallelColmWise(int m, int n) {

	double *matrix1 = (double *)malloc(m * m * sizeof(double));
	double *matrix2 = (double *)malloc(m * m * sizeof(double));
	double *result = (double *)malloc(m * m * sizeof(double));

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			*(matrix1 + i *m + j) = (double)rand() / (double)RAND_MAX;
		}
	}

	// Deep copying the contents of matrix1 to the result matrix
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			*(result + i *m + j) = *(matrix1 + i *m + j);
		}
	}
	
	//printMatrix(matrix1, m);
	for (int x = 1; x < n; x++) {
		// Deep copying the contents of result into matrix2
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				*(matrix2 + i *m + j) = *(result + i *m + j);
			}
		}

		// Setting the contents of result to 0
		memset(result, 0.0, sizeof(result));
		
		int i, j, k;
		// Matrix Multiplication
		#pragma omp parallel for private(i,k) num_threads(m)
		for (j = 0; j < m; j++) {
			for (i = 0; i < m; i++) {
				for (k = 0; k < m; k++) {
					//printf("result[i][j] = %f, matrix2[i][k] = %f, matrix1[k][j] = %f\n", *(result + i *m + j), *(matrix2 + i *m + k), *(matrix1 + k *m + j));
					*(result + i *m + j) += *(matrix2 + i *m + k) * *(matrix1 + k *m + j);
				}
			}
		}
	}
	//printMatrix(result, m);
	free(matrix1);
	free(matrix2);
	free(result);
}

void matrixPowerParallelElementWise(int m, int n) {

	double *matrix1 = (double *)malloc(m * m * sizeof(double));
	double *matrix2 = (double *)malloc(m * m * sizeof(double));
	double *result = (double *)malloc(m * m * sizeof(double));

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			*(matrix1 + i *m + j) = (double)rand() / (double)RAND_MAX;
		}
	}

	// Deep copying the contents of matrix1 to the result matrix
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			*(result + i *m + j) = *(matrix1 + i *m + j);
		}
	}

	//printMatrix(matrix1, m);
	for (int x = 1; x < n; x++) {
		// Deep copying the contents of result into matrix2
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				*(matrix2 + i *m + j) = *(result + i *m + j);
			}
		}

		// Setting the contents of result to 0
		memset(result, 0.0, sizeof(result));


		int i, j, k;
		// Matrix Multiplication
		#pragma omp parallel for private(k) num_threads(m)
		for (j = 0; j < m; j++) {
			#pragma omp parallel for private(k) num_threads(m)
			for (i = 0; i < m; i++) {
				for (k = 0; k < m; k++) {
					//printf("result[i][j] = %f, matrix2[i][k] = %f, matrix1[k][j] = %f\n", *(result + i *m + j), *(matrix2 + i *m + k), *(matrix1 + k *m + j));
					*(result + i *m + j) += *(matrix2 + i *m + k) * *(matrix1 + k *m + j);
				}
			}
		}
	}
	//printMatrix(result, m);
	free(matrix1);
	free(matrix2);
	free(result);
}

// Function Parameters
//		str = string prompt
// Function Returns
//		m = returns an integer greater than 0
// Function Description
//		Scans from the terminal and validates that the input is of type integer and is greater than 0.
//		If an input is not validated, the user is prompted to enter a new input.
int scanner(const char *str) {
	char n[1024]; //Maximum size of an input is 1024 characters
	char *ptr; // Stores the string of the next character in n after the numerical value
	double ret; // The numerical value from string n
	int status = -1; // Initialize termination condition of the while loop
	while (status != 1) {
		printf("%s\n", str); // Prints the prompt
		fgets(n, 1024, stdin); // Reads from the terminal
		ret = strtod(n, &ptr); // Seperates the initial numerical value from the string

		// If the numberical value is larger than 0, an integer, and n contains no characters, then the input is validated
		if (ret > 1 && (int)floor(ret) - ret == 0 && strcmp(ptr, "\n") == 0) {
			status = 1;
		}

		// Else the input is not validated 
		else
			status = -1;
		}

	// Converts the input into an integer
	int m = (int)ret;
	printf("The number entered is %d\n", m);
	return m;
}

void printMatrix(double *matrix, int m) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			printf("%f ", *(matrix + i *m + j));
		}
		printf("\n");
	}
	printf("\n");
}

int isNumber(char number[]) {
	char *ptr; // Stores the string of the next character in n after the numerical value
	double ret; // The numerical value from string n
	ret = strtod(number, &ptr); // Seperates the initial numerical value from the string
								// If the numberical value is larger than 0, an integer, and n contains no characters, then the input is validated
	if (ret > 1 && (int)floor(ret) - ret == 0 && strcmp(ptr, "") == 0) {
		return 1;
	}

	// Else the input is not validated 
	else {
		return -1;
	}
}
