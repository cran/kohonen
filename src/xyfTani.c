/* xyfTani.c: supervised SOMs, started Sept 26, 2005
   This version: Nov 30, 2005
   Author: Ron Wehrens
*/

#include <R.h>

#define RANDIN  GetRNGstate()
#define RANDOUT PutRNGstate()
#define UNIF unif_rand()


void XYF_Tani(double *data, double *Ys, 
	      double *codes, double *codeYs,
	      double *nhbrdist,
	      double *alphas, double *radius,	
	      double *xweight,
	      double *changes,
	      double *xdists, double *ydists, /* working arrays */
	      Sint *pn, Sint *ppx, Sint *ppy, 
	      Sint *pncodes, Sint *prlen)
{
  int n = *pn, py = *ppy, px = *ppx, ncodes = *pncodes, rlen = *prlen;
  int cd, i, j, k, l, nearest, niter;
  double dm, xdist, ydist, dist, tmp, maxx, maxy, threshold, alpha, decay;

  RANDIN;

  niter = rlen * n;
  for (k = 0; k < niter; k++) {
    /* i is a counter over objects in data, cd is a counter over units
       in the map, and j is a counter over variables */
    i = (int)(n * UNIF);
    dm = DOUBLE_XMAX;

    maxx = maxy = 0.0;
    /* calculate distances in x and y spaces, and keep track of the
       largest values */
    for (cd = 0; cd < ncodes; cd++) {
      xdist = ydist = 0.0;
      for (j = 0; j < px; j++) {
 	tmp = data[i + j*n] - codes[cd + j*ncodes];
 	xdist += tmp * tmp;
      }
      xdists[cd] = sqrt(xdist);
      if (xdists[cd] > maxx) maxx = xdists[cd];
      
      /* Tanimoto distance */
      tmp = 0;
      for (j = 0; j < py; j++) {
	if ((Ys[i + j*n] > .5 && codeYs[cd + j*ncodes] < .5) ||
	    (Ys[i + j*n] <= .5 && codeYs[cd + j*ncodes] >= .5))
	  tmp += 1.0;
      }
      ydists[cd] = tmp/(double)py;
    }
    
    /* scale x and y distances so that largest value in both cases is
       1, and sum to take total distance. Find smallest distance. */ 
    dist = DOUBLE_XMAX;
    for (cd = 0; cd < ncodes; cd++) {
      xdists[cd] /= maxx;
      tmp = *xweight * xdists[cd] + (1.0 - *xweight) * ydists[cd];
      if (tmp < dist) {
 	dist = tmp;
 	nearest = cd;
      }
    }
    
    threshold = *radius - 
      (*radius - 1.0) * (3.0 * (double) k / (double) niter);
    if (threshold < 1.0) threshold = 0.5;
    alpha = alphas[0] - (alphas[0] - alphas[1]) * (double)k/(double)niter;

    l = (int)(k/n);

    for (cd = 0; cd < ncodes; cd++) {
      if(nhbrdist[cd + ncodes*nearest] > threshold) continue;

      for(j = 0; j < px; j++) {
	tmp = data[i + j*n] - codes[cd + j*ncodes];
	codes[cd + j*ncodes] += tmp * alpha;

	if (cd == nearest) changes[l] += tmp * tmp;
      }

      for(j = 0; j < py; j++) {
	tmp = Ys[i + j*n] - codeYs[cd + j*ncodes];
	codeYs[cd + j*ncodes] += tmp * alpha;

	if (cd == nearest) changes[l+rlen] += tmp * tmp;
      }
    }
  }
  
  for (l = 0; l < rlen; l++) {
    /* mean difference per variable per object */
    changes[l] = sqrt(changes[l]/px)/n;
    changes[l+rlen] = sqrt(changes[l+rlen]/py)/n;
  }

  RANDOUT;
}
