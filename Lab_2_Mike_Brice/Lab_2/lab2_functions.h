#pragma once

#ifndef LAB1_FUNCTIONS_H_
	
	void experiment(int m);
	void openMPExperiment();
	int scanner(const char *str);
	void matrixPowerParallelColmWise(int m, int n);
	void matrixPowerParallelRowWise(int m, int n);
	void matrixPowerParallelElementWise(int m, int n);
	void printMatrix(double *matrix, int m);
	int isNumber(char number[]);

#endif // !LAB1_FUNCTIONS_H_
