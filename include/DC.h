#ifndef DC_H
#define DC_H

#include <stdint.h>
#include <string>
#include <vector>
#include <iostream>

#define Forward 0
#define Backward 1
using namespace std;
// TDT tree structure
struct DCNode 
{
	int nbParts;
	int firstRow, lastRow;
    struct DCNode **son, *iso;

	~DCNode()
	{
		if (son != NULL)
			delete[] son;
		if (iso != NULL)
			delete iso;
	}
};

struct DCArgs 
{
    int firstRow, lastRow;
};

typedef struct 
{
    int *index, *value;
} index_t;

using namespace std;

class DC
{
	private:
		int *RowPerm, *RowRev;
		int nbParts, Level;
	
		int DC_partitioning (int **Row2Row, int *nRowPerRow, int *RowValue, int **local_Row2Row, int *local_nRowPerRow, int *local_RowValue, int firstRow, int lastRow, int *rowPart);

		int DC_create_normal (DCNode *tree, int **Row2Row, int *nRowPerRow, int *RowValue, int globalNbRow, int firstRow, int lastRow, int level);
		
		// Permute "tab" 2D array of int using "perm"
		void DC_permute (int **tab, int *ntab, int *val, int *rev, int *perm, int nbRow, int offset);

		void DC_permute_1D(int *ntab, int *perm, int nbRow, int offset);

		void DC_create_permutation (int *perm, int *part, int *nRowPerRow, int **Row2Row, int size);

	public:
	
		DCNode *treeRoot;

		int DC_creation (int **Row2Row, int *nRowPerRow, int *RowValue, int globalNbRow);

		void DC_traversal (void (*userSeqFctPtr)  (char **, DCArgs *), char **userArgs, int level = 1e9);		

		int* DC_get_RowPerm();

		int* DC_get_RowRev();
							
		DC(int globalNbRow, int nbparts, int lv) 
		{
			treeRoot = new DCNode();
			
			RowPerm = new int [globalNbRow];
			RowRev = new int [globalNbRow];

			nbParts = nbparts;
			Level = lv;
		}

		~DC()
		{
			delete[] RowPerm;
			delete[] RowRev;
			
			delete treeRoot;
		}
	
};

#endif
