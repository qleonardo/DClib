#include <cstring>
#include <iostream>
#include <cmath>
#include <assert.h>
#include <omp.h>
#include <unordered_map>
#ifdef CILK
    #include <cilk/cilk.h>
#endif

#include "DC.h"

// Create a nodal graph from a item mesh for
void create_nodal_graph (int *graphIndex, int **graphValue, int **Row2Row, int *nRowPerRow, int localNbRow)
{
    graphIndex[0] = 0;
    for (int i = 0; i < localNbRow; i++) 
        graphIndex[i+1] = graphIndex[i] + nRowPerRow[i];
    if (graphIndex[localNbRow] == 0) return;
    (*graphValue) = new int[graphIndex[localNbRow]];
    int k = 0;
    for (int i = 0; i < localNbRow; i++) 
        for (int j = 0; j < nRowPerRow[i]; j++)
            (*graphValue)[k++] = Row2Row[i][j];
}
 
// Create local RowToNode array containing elements indexed contiguously from 0 to
// localNbItem and return the number of nodes
void create_local_Row2Row (int **local_Row2Row, int *local_nRowPerRow, int *local_RowValue, int **Row2Row, int *nRowPerRow, int *RowValue, int *RowPerm, int firstRow, int lastRow)
{
#pragma omp parallel for
    for (int i = firstRow; i <= lastRow; i++) {
        int newRow = i - firstRow, t = 0;
        local_nRowPerRow[newRow] = 0;
        local_Row2Row[newRow] = new int[nRowPerRow[i]];
        local_RowValue[newRow] = RowValue[i];
        for(int j = 0; j < nRowPerRow[i]; j++)
        {
            int row = RowPerm[Row2Row[i][j]];
            if (row < firstRow || row > lastRow) continue;
            Row2Row[i][t] = Row2Row[i][j];
            local_Row2Row[newRow][t++] = row - firstRow;
        }
        nRowPerRow[i] = local_nRowPerRow[newRow] = t;
    }
}

void PartGraphKway(int **Row2Row, int *nRowPerRow, int *RowValue, int *RowPart, int n, int k)
{
    int sum = 0;
    int *q = new int[n];
    for(int i = 0; i < n; i++) 
    {
        RowPart[i] = -1;
        sum += RowValue[i];
    }
    // cout << sum << endl;
    int maxRow = n / k;
    // int maxVal = sum / k;
    // int totRow = 0, totVal = 0;
    for(int i = 0, first = 0; i < k; i++, first++)
    {
        int nbRow = 0, now = 0, val = 0, flag;
        do
        {
            flag = 0;
            while(first < n)
            {
                if (RowPart[first] == -1)
                {
                    q[nbRow++] = first;
                    RowPart[first] = i;
                    // val += RowValue[first];
                    flag = 1;
                    break;
                }
                first++;
            }
            while(now < nbRow)
            {
                int Row = q[now++];
                // if ( val >= maxVal) break;
                if ( nbRow >= maxRow) break;
                for(int j = 0; j < nRowPerRow[Row]; j++)
                {
                    int nextRow = Row2Row[Row][j];
                    if (RowPart[nextRow] == -1)
                    {
                        q[nbRow++] = nextRow;
                        // val  += RowValue[nextRow];
                        RowPart[nextRow] = i;
                        // if ( val >= maxVal) break;
                        if ( nbRow >= maxRow) break;
                    }
                }
            }
        } while (nbRow <= maxRow && flag == 1);
        // } while (val <= maxVal && flag == 1);
        // totRow += nbRow;
        // totVal += val;
        // if (i < k - 1)
        // {
        //     maxRow = (n - totRow) / (k - 1 - i);
        //     maxVal = (sum - totVal) / (k - 1 - i);
        // }
    }
    delete[] q;
    for(int i = 0; i < n; i++)
        if (RowPart[i] == -1)
            RowPart[i] = k - 1;
}

// D&C partitioning of separators with more than MAX_ELEM_PER_PART elements
int DC::DC_partitioning (int **Row2Row, int *nRowPerRow, int *RowValue, int **local_Row2Row, int *local_nRowPerRow, int *local_RowValue, int firstRow, int lastRow, int *RowPart)
{
    //cout << "partitioning 0 " << endl << flush;

    int localNbRow = lastRow - firstRow + 1;

    create_local_Row2Row (local_Row2Row, local_nRowPerRow, local_RowValue, Row2Row, nRowPerRow, RowValue, RowPerm, firstRow, lastRow);
    
    //cout << "partitioning 1 " << firstRow << " " << lastRow << endl << flush;

    // int constraint = 1, objVal;
    // int *graphIndex = new int [localNbRow + 1]();
    // int *graphValue;
    
    // create_nodal_graph (graphIndex, &graphValue, local_Row2Row, local_nRowPerRow, localNbRow);
      
    int parts = min(nbParts, localNbRow);

    //cout << "partitioning 2 " << endl << flush;

    // Execution is correct without mutex although cilkscreen detects many race
    // conditions. Check if the problem is solved with future version of METIS (5.0)

    // if (localNbRow > 10000)
    //     METIS_PartGraphKway (&localNbRow, &constraint, graphIndex, graphValue,
    //                          local_RowValue, NULL, NULL, &parts, NULL, NULL,
    //                         NULL, &objVal, RowPart);
    // else
    //     METIS_PartGraphRecursive (&localNbRow, &constraint, graphIndex, graphValue,
    //                          local_RowValue, NULL, NULL, &parts, NULL, NULL,
    //                         NULL, &objVal, RowPart);

    PartGraphKway(local_Row2Row, local_nRowPerRow, local_RowValue, RowPart, localNbRow, parts);

    //cout << "partitioning 3 " << endl << flush;

    // if (graphIndex[localNbRow] > 0) 
    //     delete[] graphValue;
    // delete[] graphIndex; 
    
    // cout << "partitioning 4 " << endl << flush;

    return parts;
}
