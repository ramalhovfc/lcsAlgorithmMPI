#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>

#include "structs.h"
#include "readFileUtils.h"
#include "lcsAlgorithm.h"
#include "utils.h"

// tamanho do araray de comunicaçao para obter resultado final
#define PR_SIZE 3
// index das coordenadas com informação relevante
#define PR_FOUND_INDEX 0
#define PR_NEXTID_INDEX 1
#define PR_NEXTCOLUMN_INDEX 2

int main(int argc, char *argv[]) {

	int n_proc, id;
	FILE *fp;
	LcsMatrix *lcs;	//estrutura principal que contem todos os ficheiros necessários: nºlinhas/colunas, id, n_processos, sequencias, matriz.
	LcsResult result; // estrutura que contem a soluçao
	int pr[3]; 
	int *sendcount;
	int *displs;
	int i;

    double secs;

	char *resultString;
	int counter=0;

	MPI_Status status;

	//Variavies de MPI e inicialização dos processos. Guarda o id de 0->n-1 em id, o numero de processos n em n_proc
	MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
	
	/*
	printf("Hello %d of %d\n", id, n_proc);
	fflush(stdout);
	*/

	if (argc != 2) {
		printf("Wrong number of arguments\n");
		printf("Usage: %s <input_file>\n", argv[0]);
		fflush(stdout);
		MPI_Abort(MPI_COMM_WORLD, -1);
		MPI_Finalize();
		return -1;
	}

	//so depois de todos os processos abertos e confirmação de ficheiro de entrada é que o tempo começa a contar
	MPI_Barrier (MPI_COMM_WORLD);
    secs = - MPI_Wtime();
	
	/*
	Leitura do ficheiro
	Processo 0: reserva a estrutura principal, e lê do ficheiro todas as informações principais: dimensão total(linhas+colunas), sequencias, e guarda-as na sua struct. Também cria
	a sua matriz mtx
	Outros processos: reservam a struct principal
	*/
	if (id == 0) {
		fp = fopen(argv[1], "r");
		if (fp == NULL) {
			printf("Error openning file %s\n", argv[1]);
			fflush(stdout);
			MPI_Abort(MPI_COMM_WORLD, -2);
			MPI_Finalize();
			exit(-2);
		}
			
		lcs = readFile(fp, n_proc);
	} else {
		lcs = (LcsMatrix *)calloc(1, sizeof(LcsMatrix));
		checkNullPointer(lcs);
	}

	lcs->n_proc = n_proc;
	lcs->id = id;
	
	//[0] Envia numero de linhas e colunas para todos os outros processos
	MPI_Bcast(&lcs->lines, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&lcs->cols, 1, MPI_INT, 0, MPI_COMM_WORLD);

	/*
	Neste projecto fez-se uma divisão horizontal em submatrizes. Assim, cada processo vai ter um bloco horizontal, variando o numero de linhas mas mantendo o numero de colunas.
	Este bloco pega nas dimensões recebidas de 0(lembrar, [0], apesar de ser também uma submatriz, contem sempre as dimensões totais da matriz total) e calcula as dimensões 
	das sub-matrizes para todos os processos(excepto o 0).
	*/	
	if (id != 0) {
		if (id == n_proc-1 && (lcs->lines%n_proc > 0)) {
			// se for o ultimo tem de ficar com o resto
			i = (lcs->lines/n_proc)+((lcs->lines-1)%n_proc);
		} else {
			i = (lcs->lines/lcs->n_proc);
		}
		
		lcs->lines = i+1; //linha extra, para incluir uma linha auxiliar em cada submatriz uma linha inicial auxiliar, que vai receber as informações dos outros processos
		
		/*
		[0] contem a sequencia total. Agora é necessário reparti-la pelas sub-matrizes de todos os processos
		Cada um guarda toda a sequencia de colunas, mas apenas parte da das linhas, de dimensão lines+1. +1 é para incluir o \0
		*/
		lcs->seqLine = (char *)calloc(lcs->lines+1, sizeof(char));
		checkNullPointer(lcs->seqLine);

		lcs->seqColumn = (char *)calloc(lcs->cols+1, sizeof(char));
		checkNullPointer(lcs->seqColumn);
		
		/*
		printf("[%d] done alloc: %d lines, %d columns allocLine-size %d, allocCol-size %d \n", 
			id, lcs->lines, lcs->cols, lcs->lines+1, lcs->cols+1);
		fflush(stdout);
		*/
	}

	MPI_Bcast(&lcs->seqColumn[0], lcs->cols+1, MPI_CHAR, 0, MPI_COMM_WORLD);
		
	/*
	if (!id) {
		printf("[%d], now will send line:%s, column:%s\n",id, lcs->seqLine, lcs->seqColumn);
		fflush(stdout);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	*/

	//Com as dimensões definidas e a memória reservada, enviar todas as sequencias para cada um dos processos. Usar Scatterv para distribuir a sequencia de linhas.
	sendcount = calloc(lcs->n_proc,sizeof(int));
	displs = calloc(lcs->n_proc,sizeof(int));
	displs[0] = 0;
	sendcount[0] = 0;
	for(i=1;i<lcs->n_proc;i++) {
	
		if (i != (lcs->n_proc-1)) {
			sendcount[i] = (lcs->lines/lcs->n_proc)+1;
		} else {
			sendcount[i] = (lcs->lines/lcs->n_proc) + (lcs->lines%lcs->n_proc)+1;
		}

		displs[i] = i*(lcs->lines/lcs->n_proc);
	}
	
	/*
	if(id==0) {
		for (i=0;i<lcs->n_proc;i++) {
			printf("[%d], sendount[%d]=%d\n",id,i,sendcount[i]);
		}
		for (i=0;i<lcs->n_proc;i++) {
			printf("[%d], displs[%d]=%d\n",id,i,displs[i]);
		}
		fflush(stdout);
	}
	*/

	// enviar sequencia linha.
	MPI_Scatterv(&lcs->seqLine[0], sendcount, displs, MPI_CHAR, &lcs->seqLine[0], lcs->lines+1, MPI_CHAR, 0, MPI_COMM_WORLD);
	
	/*
	printf("[%d] antes de entrar no create_receive, lines:%d cols:%d seqLine:%s seqColumn:%s\n", id, lcs->lines, lcs->cols, lcs->seqLine, lcs->seqColumn);
	fflush(stdout);
	*/

	/*
	[0] ja contem todos os dados principais e ja definiu a matriz. Os outros processos ja receberam as informações necessárias, portanto só falta
	reservar a memória de cada sub-matriz (excepto a 0) usando os dados calculados antes
	*/
	if (id != 0) {
		create_receive(lcs);
		/*
		printf("[%d], done create_receive\n",id);
		fflush(stdout);
		*/
	}

	//função que preenche cada uma das sub-matrizes. No fim espera que todos os processos terminem, pois [0] será sempre o primeiro a terminar
	fillMatrix(lcs);
	MPI_Barrier(MPI_COMM_WORLD);

	/*
	Com a matriz preenchida e dividida por sub-matrizes espalhadas pelos processos, usa-se o mesmo sistema para encontrar a resposta.
	Encontar a sub-sequencia percorrendo desde o ultimo processo até ao [0] ... até ou se atingir o fim do processo zero ou um dos processos atingir a coluna 0
	*/
	if (id == n_proc -1) {
		result = findLongestCommonSubsequence(lcs, lcs->cols);
		pr[PR_FOUND_INDEX] = result.found;
		pr[PR_NEXTID_INDEX] = result.nextId;
		pr[PR_NEXTCOLUMN_INDEX] = result.nextColumn;

		//printf("[%d] sending found:%d nextId:%d nextColumn:%d\n",id, pr[PR_FOUND_INDEX], pr[PR_NEXTID_INDEX], pr[PR_NEXTCOLUMN_INDEX]);
		//fflush(stdout);
		MPI_Send(&pr[0], PR_SIZE, MPI_INT, id-1, LCS_TAG, MPI_COMM_WORLD);
	
	} else {
		MPI_Recv(&pr[0], PR_SIZE, MPI_INT, id+1, LCS_TAG, MPI_COMM_WORLD, &status);
		//printf("[%d] received found:%d nextId:%d nextColumn:%d\n",id, pr[PR_FOUND_INDEX], pr[PR_NEXTID_INDEX], pr[PR_NEXTCOLUMN_INDEX]);
		//fflush(stdout);

		if (pr[PR_FOUND_INDEX] && id != 0) {
			// avisar procs que encontrou soluçao
			//MPI_Isend(&pr[0], PR_SIZE, MPI_INT, id-1, LCS_TAG, MPI_COMM_WORLD, &request);
			//printf("[%d] sending found:%d nextId:%d nextColumn:%d\n",id, pr[PR_FOUND_INDEX], pr[PR_NEXTID_INDEX], pr[PR_NEXTCOLUMN_INDEX]);
			//fflush(stdout);
			MPI_Send(&pr[0], PR_SIZE, MPI_INT, id-1, LCS_TAG, MPI_COMM_WORLD);

		} else {
			result = findLongestCommonSubsequence(lcs, pr[PR_NEXTCOLUMN_INDEX]);
			pr[PR_FOUND_INDEX] = result.found;
			pr[PR_NEXTID_INDEX] = result.nextId;
			pr[PR_NEXTCOLUMN_INDEX] = result.nextColumn;
			if (id != 0) {
				//printf("[%d] sending found:%d nextId:%d nextColumn:%d\n",id, pr[PR_FOUND_INDEX], pr[PR_NEXTID_INDEX], pr[PR_NEXTCOLUMN_INDEX]);
				//fflush(stdout);
				MPI_Send(&pr[0], PR_SIZE, MPI_INT, id-1, LCS_TAG, MPI_COMM_WORLD);
			}
		}
	} 
	
	//fflush(stdout);

	MPI_Barrier(MPI_COMM_WORLD);
	
	/*
	printf("[%d] result counter:%d sequence:%s\n", id, result.counter, result.sequence);
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);
	*/

	/*
	MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
               void *recvbuf, int recvcount, MPI_Datatype recvtype,
               int root, MPI_Comm comm)
	*/
	
	MPI_Gather(&result.counter, 1, MPI_INT, &sendcount[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if(id==0) {
		for (i=0;i<lcs->n_proc;i++) {
			//printf("[%d], sendcount[%d]=%d\n",id,i,sendcount[i]);
			counter += sendcount[i];
		}
		//fflush(stdout);

		displs[0] = 0;
		for (i=1;i<lcs->n_proc;i++) {
			displs[i] = displs[i-1] + sendcount[i-1];
			//printf("[%d], displs[%d]=%d\n",id,i,displs[i]);
		}
		//fflush(stdout);

		//printf("[%d] result string size:%d\n", id, counter);
		resultString = (char *)calloc(counter+1, sizeof(char));
	}

	/*
	MPI_Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                void *recvbuf, const int *recvcounts, const int *displs,
                MPI_Datatype recvtype, int root, MPI_Comm comm)
	*/

	MPI_Gatherv(&result.sequence[0], result.counter, MPI_CHAR, &resultString[0], &sendcount[0], &displs[0], MPI_CHAR, 0, MPI_COMM_WORLD);

	/*
	if (!id) {
		printf("[%d] RESULT:%s\n", id, resultString);
		fflush(stdout);
	}
	*/

	//depois de todos os processos terminarem, verificar o tempo que demorou a correr
	MPI_Barrier (MPI_COMM_WORLD);
    secs += MPI_Wtime();

	//Print dos dados e encerrar. [0] é o unico que leu o ficheiro, portanto é o unico que o encerra
	if (id == 0) {
		printf("Execution time:%f\n", secs);
		printf("%d\n",counter);
		printf("%s\n",resultString);
		fclose(fp);
	}
	
	//free_memory(lcs, result);
	
	MPI_Finalize();
	return 0;
}

