#include <iostream>
#include <assert.h>
#include <algorithm>
#include "DC.h"

using namespace std;

void DC::DC_permute(int **tab, int *ntab, int *val, int *rev, int *perm, int nbRow, int offset)
{
    int **new_tab = new int*[nbRow];
    int *new_ntab = new int[nbRow];
    int *new_val = new int[nbRow];
    int *new_rev = new int[nbRow];
    
#pragma omp parallel for
    for(int i = 0; i < nbRow; i++)
    {
        int dst = perm[i] + offset;
        new_tab[i] = tab[dst];
        new_ntab[i] = ntab[dst];
        new_val[i] = val[dst];
        new_rev[i] = rev[dst];
    }
    
#pragma omp parallel for
    for(int i = 0; i < nbRow; i++)
    {
        int dst = i + offset;
        tab[dst] = new_tab[i];
        ntab[dst] = new_ntab[i];
        val[dst] = new_val[i];
        rev[dst] = new_rev[i];
    }

    delete[] new_tab, delete[] new_ntab, delete[] new_rev, delete[] new_val;
}

void DC::DC_permute_1D(int *ntab, int *perm, int nbRow, int offset)
{
    int *new_ntab = new int[nbRow];
    
#pragma omp parallel for
    for(int i = 0; i < nbRow; i++)
        new_ntab[i] = ntab[perm[i] + offset];
    
#pragma omp parallel for
    for(int i = 0; i < nbRow; i++)
        ntab[i + offset] = new_ntab[i];

    delete[] new_ntab;
}

#include <queue>
// Create permutation array from partition array
void DC::DC_create_permutation (int *perm, int *part, int *nRowPerRow, int **Row2Row, int size)
{
    vector<int> *rows = new vector<int>[nbParts+1];
#pragma omp parallel for
    for (int i = 0; i <= nbParts; i++)
    {
        int parts = i;
        if (parts == nbParts) parts = 1e9;
        for(int j = 0; j < size; j++)
            if (part[j] == parts) rows[i].push_back(j);
    }
    
    bool *vis = new bool[size]();
    int *offset = new int[nbParts+1]();
    for (int i = 0, j = 0; i <= nbParts; i++)
    {
        offset[i] = j;
        j += rows[i].size();
    }

#pragma omp parallel for
    for (int i = 0; i <= nbParts; i++)
    {
        int j = 0;
        queue<int> q;
        int sz = rows[i].size();
        for(int k = 0; k < sz; k++)
        {
            if (!vis[rows[i][k]])
            {
                q.push(rows[i][k]);
                vis[rows[i][k]] = 1;
                while(!q.empty())
                {
                    int now = q.front();
                    q.pop();
                    perm[offset[i] + j++] = now;
                    for(int jj = 0; jj < nRowPerRow[now]; jj++)
                    {
                        int nextRow = Row2Row[now][jj];
                        if (vis[nextRow] || part[nextRow] != part[rows[i][k]])
                            continue;
                        q.push(nextRow);
                        vis[nextRow] = 1;
                    }
                }
            }
        }
    }

    delete[] vis;
    delete[] offset;
}