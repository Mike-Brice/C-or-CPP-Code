// FinalProjectMPI.cpp : Defines the entry point for the console application.

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include <time.h>

// Stores the local distance and index
struct distance {
	double value;
	int index;
};

// Stores the global distance and index
struct closest {
	double value;
	int index;
};

int main(int argc, char *argv[]) {

	int numberHouses = atoi(argv[1]); // Number of houses entered from command line arguments
	int id; // Process ID
	int p; // Number of processes
	double elapsed_time; // Elapsed Time
	struct distance distance; // local distance structure
	struct closest closest; // global distance structure

	// Function Identifiers
	double calcDistance(double x2, double y2);
	void process(double* houses[2], int start, int end, double& distance, int& index);

	// MPI Initializers
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	// Two dimensional array of the x and y coordinates of the houses
	double* houses[2];
	houses[0] = new double[numberHouses];
	houses[1] = new double[numberHouses];

	// Arrays to hold the start and end indexes of the subarrays
	int* startArray;
	int* endArray;
	startArray = new int[p];
	endArray = new int[p];

	// Generates the x and y coordinates in processes 0
	if (id == 0) {

		// Generates the x any y coordinates of the houses
		for (int i = 0; i < numberHouses; i++) {
			houses[0][i] = (double)rand();
			houses[1][i] = (double)rand();
		}
	}

	// Generates the start and end indexes of the subarrays in processes 1 unless there are only 1 processes than it is generated in processes 0
	if (id == 1 || p == 1) {

		// The initial start and end indexes of the subarrays
		int start = 0;
		int chunk = (int)ceil((numberHouses / p));
		int end = chunk;

		// Calculates the start and end indexes of the subarrays
		for (int i = 0; i < p; i++) {
			startArray[i] = start;
			endArray[i] = end;

			// Start is 1 after the end of the previous subarray 
			start = end;

			// The end is the end plus 
			end = chunk + end;
		}

		// If the last subarray ends before the last house add the remainding houses
		if (end != numberHouses) {
			endArray[p - 1] = numberHouses;
		}
	}

	// Wait for all processes to complete the preprocessing
	MPI_Barrier(MPI_COMM_WORLD);

	// Broadcast if number of processes is greater than 1
	if (p > 1) {
		MPI_Bcast(houses[0], numberHouses, MPI_DOUBLE, 0, MPI_COMM_WORLD); // Broadcasts the x coordinates
		MPI_Bcast(houses[1], numberHouses, MPI_DOUBLE, 0, MPI_COMM_WORLD); // Broadcasts the y coordiantes
		MPI_Bcast(startArray, p, MPI_INT, 1, MPI_COMM_WORLD); // Broadcasts the starting indexes
		MPI_Bcast(endArray, p, MPI_INT, 1, MPI_COMM_WORLD); // Broadcasts the ending indexes
	}

	// Start the clock
	elapsed_time = -MPI_Wtime();

	// Calls the processes
	process(houses, startArray[id], endArray[id], distance.value, distance.index);

	// Wait for all processes to complete before performing the reduction
	MPI_Barrier(MPI_COMM_WORLD);

	// Finds the miniumum distance and the index
	MPI_Reduce(&distance, &closest, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);

	// Stop the clock
	elapsed_time += MPI_Wtime();

	// If the processes is processes 0 print the closest house
	if (!id) {
		printf("MPI\nNumber of processes = %d\nNumber of houses = %d\nTrain Station located at "
			"(%d, %d)\nClosest House is at x: %f, y: %f, with a distance of %f meters\nTime: %f seconds\n", 
			p, numberHouses, 5, 10, houses[0][closest.index], houses[1][closest.index], sqrt(closest.value), elapsed_time);
		fflush(stdout);
	}
	MPI_Finalize();

	// Garbage Collection
	delete[] houses[0];
	delete[] houses[1];
	delete[] startArray;
	delete[] endArray;

    return 0;
}

// x2 is the x coordinate of the house
// y2 is the y coordinate of the house
// 5 is the x coordinate of the train station
// 10 is the y coordinate of the train station
// returns the distance squared
double calcDistance(double x2, double y2) {
	return pow((5 - x2), 2) + pow((10 - y2), 2);
}

// houses is the array of houses coordinates
// start is the starting index of the subarray
// end is the ending index + 1 of the subarray
// distance is referencing the structure distance's value
// index is referencing the structure distance's index
void process(double* houses[2], int start, int end, double& distance, int& index) {

	double closest = calcDistance(houses[0][start], houses[1][start]); // Calculates the initial distance
	int closestIndex = start; // Sets the cloest Index to  the starting index

	for (int i = start + 1; i < end; i++) {
		double dist2 = calcDistance(houses[0][i], houses[1][i]); // Calculates the next distance

		// If the current distance is less than the closest distance, replace closest with the current distance. Set closestIndex to the current index
		if (closest > dist2) {
			closest = dist2;
			closestIndex = i;
		}
	}

	distance = closest; // Assign closest to the structure distance's value
	index = closestIndex; // Assign closestIndex to the structure distance's index
}