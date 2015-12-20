
#ifndef cpdProject2014_structs_h
#define cpdProject2014_structs_h

// group number
#define LCS_TAG 15

struct lcs_Matrix {
	int id;
	int n_proc;
	
	int lines;
	int cols;

	char *seqLine;
	char *seqColumn;
	int **mtx;
};
typedef struct lcs_Matrix LcsMatrix;

struct lcsResult {
	int counter;
	char *sequence;

	int found;
	int nextId;
	int nextColumn;
};
typedef struct lcsResult LcsResult;

#endif
