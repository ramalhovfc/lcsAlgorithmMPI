
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "structs.h"

void checkNullPointer(void *p) {
	
	if (p == NULL) {
		puts("Error alocing memory");
		fflush(stdout);
		MPI_Abort(MPI_COMM_WORLD, -1);
		MPI_Finalize();
		exit(-1);
	}
}

void free_memory(LcsMatrix *lcsMtx, LcsResult result) {

	int i;

	for (i=0; i<lcsMtx->lines+1; i++) {
		free(lcsMtx->mtx[i]);
	}
	free(lcsMtx->mtx);

	free(lcsMtx->seqLine);
	free(lcsMtx->seqColumn);

	free(lcsMtx);

	free(result.sequence);

	return;
}

void printLcsMatrix(LcsMatrix* lcsMtx) {
	
	int i, j;
	
	printf("\n");
	//printf("seqLine: %s, %d\n", lcsMtx->seqLine, lcsMtx->lines);
	//printf("seqColumn: %s, %d\n", lcsMtx->seqColumn, lcsMtx->cols);
	printf("\n\t\tMATRIX\n");
	
	for (i=-2; i<lcsMtx->lines+1; i++) {
		for (j=-2; j<lcsMtx->cols+1; j++) {
			
			// blanks
			if (i==-2 || i==-1) {
				if (j==-2 || j==-1) {
					printf("%3c", ' ');
					continue;
				}
			}
			
			// column index
			if (i==-2) {
				printf("%3d", j);
				continue;
			
			// column sequence
			} else if(i==-1 && j>=0) {
				printf("%3c", lcsMtx->seqColumn[j]);
				continue;
			}
			
			// line index
			if (j==-2) {
				printf("%3d", i);
				continue;
				
			// line sequence
			} else if (j==-1 && i>=0) {
				printf("%3c", lcsMtx->seqLine[i]);
				continue;
			}
			
			// matrix value
			printf("%3d", lcsMtx->mtx[i][j]);
		}
		printf("\n");
	}
	
}
