#ifdef CILK
    #include <cilk/cilk.h>
#endif
#include <omp.h>
#include <bits/stdc++.h>
#include "DC.h"

using namespace std;

// Wrapper used to parallel user function
void DC::DC_traversal (void (*userSeqFctPtr)  (char **, DCArgs *), char **userArgs, int level)
{
        DCNode *treePtr = treeRoot;

        do
        {

#pragma omp parallel for
                for (int i = 0; i < treePtr->nbParts; i++)
                {
                        DCNode *ptr = treePtr->son[i];
                        if (ptr == nullptr) continue;
                        DCArgs treeArgs;
                        treeArgs.firstRow  = ptr->firstRow;
                        treeArgs.lastRow   = ptr->lastRow;
                
                        // Call user sequential function
                        userSeqFctPtr (userArgs, &treeArgs);   
                }

                if (treePtr->nbParts == 0)
                {
                        DCArgs treeArgs;
                        treeArgs.firstRow  = treePtr->firstRow;
                        treeArgs.lastRow   = treePtr->lastRow;
                
                        // Call user sequential function
                        userSeqFctPtr (userArgs, &treeArgs);   
                }

                treePtr = treePtr->iso;
                if (--level == 0) break;
        } while(treePtr != nullptr);
                
        
}
