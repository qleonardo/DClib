#ifdef CILK
    #include <cilk/cilk.h>
#endif
#include <limits.h>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <assert.h>
#include <unordered_map>
#include <vector>

#include "DC.h"
#include "omp.h"

// Initialize the content of DC tree nodes
void DC_node_init (DCNode *treePtr, int firstRow, int lastRow, int nbIsoRow, int *nbPartRow, int nbParts)
{
    treePtr->firstRow    = firstRow;
    treePtr->lastRow     = lastRow + 1;
    treePtr->iso          = nullptr;
	treePtr->son		  = nullptr;
    treePtr->nbParts      = nbParts;

    if (nbParts > 0) {
		
		treePtr->son = new DCNode*[nbParts];
		for (int i = 0; i < nbParts; i++)
			if (nbPartRow[i] > 0)
                treePtr->son[i] = new DCNode();
			else 
				treePtr->son[i] = nullptr;

        if (nbIsoRow > 0) 
            treePtr->iso = new DCNode();
    }
}

void DC_create_RowPart (int *RowPart, int **Row2Row, int *nRowPerRow, int *RowValue, int localNbRow, int *nbPartRow, int *nbIsoRow, int nbParts)
{

	for (int row = 0; row < localNbRow; row++) {
        int colorA = RowPart[row], colorB = -1;
		for (int j = 0; j < nRowPerRow[row]; j++) {
			int Row = Row2Row[row][j];
            if (RowPart[Row] == 1e9) continue;
			if (colorA != RowPart[Row]) 
                colorB = RowPart[Row];
		}
        if (colorB == -1) nbPartRow[colorA]++;
        else
        {
            (*nbIsoRow) ++;
            RowPart[row] = 1e9;
        }
	}
   
}

double I = 0;

// Create the DC tree and the element permutation, and compute the intervals of nodes
// and elements at each node of the tree-
int DC::DC_create_normal  (DCNode *tree, int **Row2Row, int *nRowPerRow, int *RowValue, int globalNbRow, int firstRow, int lastRow, int level)
{
    int localNbRow = lastRow - firstRow + 1;
    
    //cout << firstRow << ' ' << lastRow << endl << flush;
	if (level >= Level || localNbRow <= 1) {
        DC_node_init (tree, firstRow, lastRow, 0, 0, 0);
        if (localNbRow <= 1) return lastRow+1;
        return firstRow;
    }
     //cout << "tree_creation 0 " << endl << flush;

    int *local_nRowPerRow = new int[localNbRow]();
    int *local_RowValue = new int[localNbRow]();
    int **local_Row2Row = new int*[localNbRow];
    int *RowPart = new int [localNbRow];
   // cout << "tree_creation 0.1 " << endl << flush;
    int parts = DC_partitioning (Row2Row, nRowPerRow, RowValue, local_Row2Row, local_nRowPerRow, local_RowValue, firstRow, lastRow, RowPart);

    //cout << "tree_creation 0.5 " << endl << flush;
    int nbIsoRow = 0;
    int *nbPartRow = new int [parts]();

    DC_create_RowPart (RowPart, local_Row2Row, local_nRowPerRow, local_RowValue, localNbRow, nbPartRow, &nbIsoRow, nbParts);


    //cout << "tree_creation 2 " << endl << flush;

    // Create local element permutation
    int *localRowPerm = new int [localNbRow];
    DC_create_permutation (localRowPerm, RowPart, local_nRowPerRow, local_Row2Row, localNbRow);

    int *nnz = new int[parts]();
    for(int i = 0; i < localNbRow; i++)
    {
        if (RowPart[i] == 1e9) continue;
        nnz[RowPart[i]] += local_RowValue[i];
    }
    //cout << "tree_creation 3 " << endl << flush;
    // int mx = 0, mn = 1e9;
    // for(int i = 0; i < parts; i++)
    // {
    //     cout << nnz[i] << " ";
    //     mx = max(mx, nnz[i]+1);
    //     mn = min(mn, nnz[i]+1);
    // }
    // cout << endl << mx << " " << mn << endl;
    // I += 1.0 * mx / mn;

    delete[] RowPart;
#pragma omp parallel for
    for(int i = 0; i < localNbRow; i++)
        if (local_nRowPerRow[i] > 0)
            delete[] local_Row2Row[i];
    delete[] local_nRowPerRow;
    delete[] local_Row2Row;
    delete[] local_RowValue;

    //cout << "tree_creation 4 " << endl << flush;

    DC_permute (Row2Row, nRowPerRow, RowValue, RowRev, localRowPerm, localNbRow, firstRow);
#pragma omp parallel for
    for(int i = firstRow; i <= lastRow; i++)
        RowPerm[RowRev[i]] = i;

    //cout << "tree_creation 5 " << endl << flush;
    delete[] localRowPerm;

    // Initialize current node    
    DC_node_init (tree, firstRow, lastRow, nbIsoRow, nbPartRow, parts);
    // cout << tree->firstRow << ' ' << tree->lastRow << " " << tree->nbParts << endl << flush;

    // cout << "tree_creation 6 " << endl << flush;

    if (nbIsoRow == localNbRow) return firstRow;
    
    int *stRow = new int[parts];
    int *edRow = new int[parts];
    stRow[0] = firstRow;
    edRow[0] = stRow[0] + nbPartRow[0] - 1;
    for (int i = 1; i < parts; i++) {
        stRow[i] = edRow[i-1] + 1;
        edRow[i] = stRow[i] + nbPartRow[i] - 1;
    }
#ifndef FORKJOIN
#ifdef OMP
#pragma omp taskloop default(shared)
    for (int i = 0; i < parts; i++) {    
#elif CILK
    cilk_for (int i = 0; i < parts; i++) {
#endif
#else
// #pragma omp parallel for
    for (int i = 0; i < parts; i++) {    
#endif
        if (nbPartRow[i] > 0)
            DC_node_init (tree->son[i], stRow[i], edRow[i], 0, 0, 0);
    }

    // cout << "tree_creation 6.5 " << endl << flush;
    delete[] stRow, delete[] edRow;
    delete[] nbPartRow;
    
    // cout << "tree_creation 7 " << endl << flush;

    // DC partitioning of Isoarator elements
    if (nbIsoRow > 0)
		return DC_create_normal (tree->iso, Row2Row, nRowPerRow, RowValue, globalNbRow, lastRow-nbIsoRow+1, lastRow, level+1);
    else
        if (nbIsoRow > 0)
            return firstRow;
        else
            return lastRow+1;

    // cout << "tree_creation 8 " << endl << flush;

}

// Create the DC tree and the permutations
// qleonrdo: dimElem can be different
//void DC_create_tree (int *elemToNode, int nbElem, int dimElem, int nbNodes)
int DC::DC_creation (int **Row2Row, int *nRowPerRow, int *RowValue, int globalNbRow)
{
    double t1 = omp_get_wtime();
    //cout << "DC_Creation_1.0  "<< endl << flush;
#ifdef OMP
    #pragma omp parallel for
    for (int i = 0; i < globalNbRow; i++) 
        RowPerm[i] = RowRev[i] = i;
#elif CILK
    cilk_for (int i = 0; i < globalNbRow; i++) 
        RowPerm[i] = RowRev[i] = i; 
#endif

    //cout << "build Start!!!!!!!!!!  "<< endl << flush;

    int res = DC_create_normal (treeRoot, Row2Row, nRowPerRow, RowValue, globalNbRow, 0, globalNbRow-1, 0);

    //cout << "build finish!!!!!!!!!!  " << treeRoot->nbParts << endl << flush;

    double t2 = omp_get_wtime();

    // cout << "Build DC Times: " << t2 - t1 << endl << flush;

    return res;

}

int* DC::DC_get_RowPerm()
{
    return RowPerm;
}

int* DC::DC_get_RowRev()
{
    return RowRev;
}
