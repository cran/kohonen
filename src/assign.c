/* assign.c: map objects to a trained SOM
   This version: Oct 13, 2005
   Author: Ron Wehrens
*/

#include <R.h>


void kohonenAssign(double *data, double *codes,
		   Sint *pncodes, Sint *pnd, Sint *pp,
		   Sint *classif, double *dists)
{
  int ncodes = *pncodes, nd = *pnd, p = *pp;
  int i, j, hit;
  double tmp;
  
  /* i is a counter over objects in data, j is a counter over variables. */
  for (i = 0; i < nd; i++) {
    hit = classif[i];
    
    dists[i] = 0.0;
    for (j = 0; j < p; j++) {
      tmp = data[i + j*nd] - codes[hit + j*ncodes];
      dists[i] += tmp * tmp;
    }
    
    dists[i] = sqrt(dists[i]);
  }    
}

