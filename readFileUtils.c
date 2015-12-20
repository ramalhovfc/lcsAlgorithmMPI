#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "structs.h"
#include "utils.h"

LcsMatrix* readFile(FILE *fp, int nproc) {
	
	int i;
	LcsMatrix *lcs = (LcsMatrix *)calloc(1, sizeof(LcsMatrix));
	checkNullPointer(lcs);
	
	/*
	Ler dimensões do tamanho total da matriz e guardá-los no processo 0
	Nota: lcs->lines e cols no [0] irá sempre conter o tamanho total da matriz. Só os outros processos terão variação no numero de linhas(colunas é sempre constante também)
	*/
	if (fscanf(fp, "%d %d", &(lcs->lines), &(lcs->cols)) != 2) {
		puts("Error reading file");
		fflush(stdout);
		MPI_Abort(MPI_COMM_WORLD, -1);
		MPI_Finalize();
		exit(-1);
	}

	//Linhas extra para incluir linha e coluna de zeros
	lcs->lines=lcs->lines+1;
	lcs->cols=lcs->cols+1;
	
	/* 
	Alocar tamanho da matriz de [0]. 
	Nota: apesar de conter os tamanhos e sequencias totais do ficheiro, mtx de [0] também será uma sub-matriz com tamanho reduzido
	*/
	lcs->mtx = (int **)calloc((lcs->lines/nproc)+1, sizeof(int **));
	checkNullPointer(lcs->mtx);
	for (i=0; i<(lcs->lines/nproc)+1; i++) {
		lcs->mtx[i] = (int *)calloc(lcs->cols, sizeof(int *));
		checkNullPointer(lcs->mtx[i]);
	}
	
	//Alocar espaço em zero para guardar sequencias totais. +1 é para incluir o \0
	lcs->seqLine = (char *)calloc(lcs->lines+1, sizeof(char));
	checkNullPointer(lcs->seqLine);
	lcs->seqColumn = (char *)calloc(lcs->cols+1, sizeof(char));
	checkNullPointer(lcs->seqColumn);
	
	// read sequences
	if (fscanf(fp, "%s", &lcs->seqLine[1]) != 1) {
		puts("Error reading file");
		fflush(stdout);
		MPI_Abort(MPI_COMM_WORLD, -1);
		MPI_Finalize();
		exit(-1);
	}
	lcs->seqLine[0] = ' ';
	lcs->seqLine[lcs->lines]='\0';
	if (fscanf(fp, "%s", &lcs->seqColumn[1]) != 1) {
		puts("Error reading file");
		fflush(stdout);
		MPI_Abort(MPI_COMM_WORLD, -1);
		MPI_Finalize();
		exit(-1);
	}
	lcs->seqColumn[0] = ' ';
	lcs->seqColumn[lcs->cols]='\0';
	
	return lcs;
}

void create_receive(LcsMatrix *lcs) {

	//Aloca a sub-matriz de cada processo(excepto o 0), tendo em conta as dimensões calculadas atrás
	int i;
	
	/*
	printf("[%d], allocing %d lines, %d columns\n", lcs->id, lcs->lines+1, lcs->cols+1);
	fflush(stdout);
	*/

	lcs->mtx = (int **) calloc(lcs->lines+1, sizeof(int *));
	checkNullPointer(lcs->mtx);
	
	for (i=0; i<=lcs->lines; i++) {
		lcs->mtx[i] = (int *) calloc(lcs->cols+1, sizeof(int)); 
		checkNullPointer(lcs->mtx[i]);
	}
}
