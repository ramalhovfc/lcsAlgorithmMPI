
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "structs.h"

short cost(int x)
{
    int i, n_iter = 20;
    double dcost = 0;
    for(i = 0; i < n_iter; i++)
		dcost += pow(sin((double) x),2) + pow(cos((double) x),2);
    return (short) (dcost / n_iter + 0.1);
}

/*
 
 c[i,j]
 
 0							if i=0 or j=0
 c[i-1,j-1]+1				if i,j>0 and xi=yj
 max(c[i,j-1], c[i-1,j])	if i,j>0 and xi!=yj
 
 */

void fillMatrix(LcsMatrix *lcsMtx) {
	
	int control;
	int n_proc, id;
	int cols, lines;
	int i, j;
	int **c;
	
	int *sender;
	int *receiver;
	int msgSize;

	
	int lines_start;
	int lines_finish;
	int cols_start;
	int cols_finish;
	
	char *seqLine;
	char *seqCol;

	MPI_Status status;
	
	//Para facilitar, guardar todas as informações de cada processo em variáveis próprias da função	
	n_proc=lcsMtx->n_proc;
	id=lcsMtx->id;
	cols=lcsMtx->cols;
	lines=lcsMtx->lines;
	c=lcsMtx->mtx;
	seqLine=lcsMtx->seqLine;
	seqCol=lcsMtx->seqColumn;
	
	/*
	printf("[%d] fillmatrix: n_proc:%d, cols:%d, lines:%d, seqLine:%s, seqCol:%s\n", id, n_proc, cols, lines, seqLine, seqCol);
	fflush(stdout);
	*/

	// tamanho maximo de mensagem enviada/recebida é quando tamos a tratar as ultimas colunas
	receiver = (int*) calloc((cols/n_proc) + (cols%n_proc), sizeof(int));
	sender = (int*) calloc((cols/n_proc) + (cols%n_proc), sizeof(int));
	
	/*
	printf("[%d] done allocing buffers of size:%d\n", id, (cols/n_proc) + (cols%n_proc));
	fflush(stdout);
	*/
	
	/*
	Depois de criar cada um dos blocos horizontais, vamos agora calcular as sub-matrizes em blocos verticais também de dimensão ncols/n_procc(explicado mais em detalhe no relatório).
	Assim, cada processo irá calcular um bloco de dados n vezes, e usa-se assim uma variavel control para controlar isso
	*/
	for(control=0; control<n_proc; control++) {

		/*
		printf("[%d] just started control %d\n", id, control);
		fflush(stdout);
		*/

		/*
		Para cada bloco a calcular, definir as posições inicial e final, tanto de linhas como de colunas, tendo em conta o valor actual de control.
		Agora, linhas serão constantes(pois cada bloco, excepto o final, tem o mesmo numero de linhas), mas colunas, vamos calcular em multiplos de cols/n_proc
		*/
		if(control == n_proc-1) {
			cols_finish = ((control+1) * ((cols-1)/n_proc)) + ((cols-1)%n_proc);
		} else {
			cols_finish = ((control+1) * ((cols-1)/n_proc));
		}

		if(control == 0) 
			cols_start = 1;
		else 
			cols_start = (control * ((cols-1)/n_proc))+1;
		
		lines_start=1;
		
		if(id == 0) 
			lines_finish = lines/n_proc;
		else 
			lines_finish=lines-1;
		
		/*
		printf("[%d] lines_start=%d lines_finish=%d cols_start=%d cols_finish=%d\n", id, lines_start, lines_finish, cols_start, cols_finish);
		fflush(stdout);
		*/

		//MPI_Barrier(MPI_COMM_WORLD);

		/*
		Depois de definir o tamanho da mensagem a receber e as dimensões, o processo espera que o anterior lhe envie as informações que precisa.
		Estas serão guardadas na linha auxiliar inicial de cada sub-matriz.
		Nota: id!=0 é devido ao facto de que [0], sendo o primeiro, não recebe de ninguém
		*/
		msgSize = cols_finish-cols_start+1;
		
		if(id != 0) {
			
			MPI_Recv(&receiver[0], msgSize, MPI_INT, id-1, LCS_TAG, MPI_COMM_WORLD, &status);
			
			memmove(&c[lines_start-1][cols_start], &receiver[0], msgSize*sizeof(int));

			/*
			printf("[%d] did copy %d received ints to [%d][%d]\n", id, msgSize,lines_start-1,cols_start);
			for(i=0; i<msgSize ;i++) {			
				printf("[%d] received[%d]=%d \n", id, i, receiver[i]);
			}
			fflush(stdout);
			*/
		}
		
		// Com os dados recebidos e posição determinada, calcular um bloco da sub-matriz
		for(i=lines_start; i<=lines_finish; i++) {
			for (j=cols_start; j<=cols_finish; j++) {
				
				if (seqLine[i] == seqCol[j]) {
					c[i][j] = c[i-1][j-1] + cost(i);
					
				} else {
					c[i][j] = c[i][j-1] >= c[i-1][j] ? c[i][j-1] : c[i-1][j];
				}
			}
		}

		/*
		if (!id) {
			printf("[%d] printing matrix ***************\n");
			printLcsMatrix(lcsMtx);
			fflush(stdout);
		}
		*/

		/*
		Depois de calculado, ele envia a ultima linha para o processo seguinte, que a guardará na primeira linha auxiliar e iniciará o que este processo acabou de terminar
		Nota: tal como o primeiro nao recebe, o ultimo processo também não envia informações a ninguem
		*/
		if (id < n_proc-1) {
			// printf("[%d] a copiar para o buffer a começar em linha:%d coluna:%d, %d elementos\n", id, lines_finish, cols_start, cols_finish-cols_start+1);
			memmove(&sender[0], &c[lines_finish][cols_start], (cols_finish-cols_start+1)*sizeof(int));

			/*
			for (j=0; j<=cols_finish-cols_start; j++) {
				printf("[%d] sender[%d]=%d\n", id, j, sender[j]);
			}
			*/
			
			MPI_Send(sender, msgSize, MPI_INT, id+1, LCS_TAG, MPI_COMM_WORLD);
		}

		/*
		if (id == n_proc-1) {
			printf("[%d] printing matrix\n", id);
			printLcsMatrix(lcsMtx);
			fflush(stdout);
		}
		*/
		
	}
}

LcsResult findLongestCommonSubsequence(LcsMatrix *lcsMtx, int column) {
	
	LcsResult result;
	
	int i,j;

	char *lcsStringInverted;
	int counter = 0;

	// começar no canto inferior direito
	if (lcsMtx->id == lcsMtx->n_proc-1) 
		j = lcsMtx->cols -1;
	else
		j = column;
	
	if (lcsMtx->id == 0)
		i = lcsMtx->lines/lcsMtx->n_proc;		
	else 
		i = lcsMtx->lines;

	//printf("[%d] finding LCS, starting line:%d, starting column:%d\n", lcsMtx->id, i, column);
	//fflush(stdout);

	// alocar espaço máximo da subsequencia
	lcsStringInverted = (char *)calloc(i<j ? i+1 : j+1, sizeof(char));

	while (i!=0 && j!=0) {
		// check match
		if (lcsMtx->seqLine[i] == lcsMtx->seqColumn[j]) {
			lcsStringInverted[counter] = lcsMtx->seqLine[i];
			counter = counter + 1;
			
			// move diagonally
			i--;
			j--;
			
		} else {
			// chech which is larger
			if (lcsMtx->mtx[i][j-1] >= lcsMtx->mtx[i-1][j]) {
				// move to largest
				j--;
			} else {
				i--;
			}
		}
	}

	if (i == 0) {
		
		// chegou ao topo do zero, acabou
		if (lcsMtx->id == 0) {
			result.found = 1;
			result.nextId = 0;
		} else {
			// chegou ao topo, mandar trabalho para outro
			result.found = 0;
			result.nextId = lcsMtx->id -1;
		}
		
		result.nextColumn = j;

		/*
		if (!result.found)
			printf("[%d] will send to:%d nextColumn:%d found:%d\n", lcsMtx->id, result.nextId, result.nextColumn, result.found);
		else
			printf("\t\t[%d] FOUND SOLUTION !!!\n", lcsMtx->id);
		fflush(stdout);
		*/

	} else if (j == 0) {
		// chegou á soluçao
		result.found = 1;
		//printf("[%d] found solution\n", lcsMtx->id);
	}	
	
	lcsStringInverted[counter] = '\0';

	// inverter subsequencia parcial
	result.sequence = (char *)calloc(counter+1, sizeof(char));
	for (i=counter-1; i>=0; i--) {
		result.sequence[counter-1-i] = lcsStringInverted[i];
	}
	
	result.counter = counter;

	/*	
	printf("[%d] partial sequence:%s\n", lcsMtx->id, result.sequence);
	fflush(stdout);
	*/
	
	return result;
}






