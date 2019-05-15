// FinalProject.cpp : Defines the entry point for the console application.

#include "stdafx.h"
#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

int main() {

	/* Function Identifiers */
	double calcDistance(double x2, double y2);
	void process(double* houses[2], int start, int end, double& distance, double& index);

	/* Constants */
	const int n = 10000000; // Number of houses
	const int threads = 8; // The number of threads to be used by openMP

	// An array to store the houses coordinates
	double* houses[2];
	houses[0] = new double[n]; // Allocate a array of type "double" on the heap of size [n], assign to houses[0]
	houses[1] = new double[n]; // Allocate another array on the heap, assign to houses[1]

	// Location of the station is (5,10)
	// Generates the random locations of the houses
	for (int j = 0; j < n; j++) {
		houses[0][j] = rand();
		houses[1][j] = rand();
	}

	// An array to store the shortest distances from the threads and their respective indexes 
	double* finalData[2];
	finalData[0] = new double[threads];
	finalData[1] = new double[threads];

	// Arrays to hold the start and end indexes of the subarrays
	int* startArray;
	int* endArray;
	startArray = new int[threads];
	endArray = new int[threads];

	// The initial start and end indexes of the subarrays
	int start = 0;
	int chunk = (int)ceil((n / threads));
	int end = chunk;

	// Calculates the start and end indexes of the subarrays
	for (int i = 0; i < threads; i++) {
		startArray[i] = start;
		endArray[i] = end;

		// Start is 1 after the end of the previous subarray 
		start = end;

		// The end is the end plus the chunk... End is one larger than the subarray size because of how the process function works
		end = chunk + end;
	}

	if (end != n) {
		endArray[threads - 1] = n;
	}

	clock_t t;
	t = clock(); // Start time

// Parallel For loop with the number of active threads set to threads
#pragma omp parallel for num_threads(threads) 
	for (int i = 0; i < omp_get_num_threads(); i++) {
		process(houses, startArray[i], endArray[i], finalData[0][i], finalData[1][i]); 
	}

	/* Finds the shortest distance of the returned distances */
	double closest = finalData[0][0]; // Grabs the first distance and stores it in closest

	int closestIndex = finalData[1][0]; // Grabs the first index and stores it in closestIndex

	for (int i = 1; i < threads; i++) {

		// If the distance is less than closest, replace closest with the distance
		if (closest > finalData[0][i]) {
			closest = finalData[0][i];
			// Update the index
			closestIndex = (int)finalData[1][i];
		}
	}

	t = clock() - t; // end time
	double time_taken = ((double)t) / CLOCKS_PER_SEC; // Time elapsed in seconds
	printf("OpenMP\nNumber of processes = %d\nNumber of houses = %d\nTrain Station located at "
		"(%d, %d)\nClosest House is at x: %f, y: %f, with a distance of %f meters\nTime: %f seconds\n", threads, n, 5, 10, houses[0][closestIndex], houses[1][closestIndex], sqrt(closest), time_taken);

	/* Garbage Collection */
	delete[] houses[0];
	delete[] houses[1];
	delete[] finalData[0];
	delete[] finalData[1];
	delete[] startArray;
	delete[] endArray;

	return 0;
}

// x2 is the x coordinate of the house
// y2 is the y coordinate of the house
// returns the distance squared
double calcDistance(double x2, double y2) {

	// Distance Squared | Station is at (5,10)
	return pow((5 - x2), 2) + pow((10 - y2), 2);
	//return (5 - x2) * (5 - x2) + (10 - y2) * (10 - y2);
}

// houses is the 2D array of house coordinates
// start is the starting index
// end is the ending index + 1
// distance references the finalData[0][i] index of the array
// index references the finalData[1][i] index of the array
void process(double* houses[2], int start, int end, double& distance, double& index) {

	// i is set to the starting index
	int i = start;

	// Stores the closest house's index
	int closestIndex = 0;

	// Gets the distance squared of the first house to the train station
	double closest = calcDistance(houses[0][i], houses[1][i]);

	for (i = start+1; i < end; i++) {

		// The distance squared is used because you do not need to know the actual distance to compare which one is closer 
		double dist2 = calcDistance(houses[0][i], houses[1][i]);

		// If the distance squared is smaller than the current closest distance, replace the distance and the location of the house. 
		if (closest > dist2) {
			closest = dist2;
			closestIndex = i;
		}
	}

	distance = closest;
	index = closestIndex;
}