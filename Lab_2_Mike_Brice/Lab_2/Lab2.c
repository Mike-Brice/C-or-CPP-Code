// Lab1.c : Defines the entry point for the console application.
//

#include "lab2_functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>



int main(int argc, char *argv[])
{
	int m;
	if (isNumber(argv[1]) == 1) {
		m = atoi(argv[1]);
		printf("The size entered is %d\n", m);
	}
	else {
		// Gets the size of the matrix from the user
		m = scanner("Please enter the size of the matrix");
	}
	
	experiment(m);
	return 0;
}





