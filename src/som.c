#include <R.h>

#define RANDIN  GetRNGstate()
#define RANDOUT PutRNGstate()
#define UNIF unif_rand()

/* SOM function for the kohonen package, modelled after BDR's function
   in the class package. */

SOM_online(double *data, double *codes, double *nhbrdist,
	   double *alphas, double *pradius, double *changes,
	   Sint *pn, Sint *pp, Sint *pncodes, Sint *prlen)
{
  int n = *pn, p = *pp, ncodes = *pncodes, rlen = *prlen;
  double radius = *pradius;
  int cd, i, j, k, l, nearest, niter;
  double dm, dist, tmp, alpha, threshold;
  
  RANDIN;

  niter = rlen * n;
  for (k = 0; k < niter; k++) {
    /* pick a random data point */
    i = (int)(n * UNIF);

    /* find the nearest code 'near' */
    dm = DOUBLE_XMAX;
    for (cd = 0; cd < ncodes; cd++) {
      dist = 0.0;
      for (j = 0; j < p; j++) {
	tmp = data[i + j*n] - codes[cd + j*ncodes];
	dist += tmp * tmp;
      }

      if (dist < dm) {
	nearest = cd;
	dm = dist;
      }
    }

    /* update all codes within threshold of 'nearest'. Maximal radius
       at the start, minimal radius equals 1.0; maximal alpha at the
       start, minimal alpha at the end. Linear decrease for both. */
    threshold = radius - (radius-1.0) * (double)k/(double)niter;
    alpha = alphas[0] - (alphas[0] - alphas[1]) * (double)k/(double)niter;

    for (cd = 0; cd < ncodes; cd++) {
      if(nhbrdist[cd + ncodes*nearest] > threshold) continue;

      if (cd == nearest) l = (int)(k/n);

      for(j = 0; j < p; j++) {
	tmp = data[i + j*n] - codes[cd + j*ncodes];
	codes[cd + j*ncodes] += tmp * alpha;

	if (cd == nearest) changes[l] += tmp * tmp;
      }
    }
  }

  for (l = 0; l < rlen; l++)
    changes[l] /= n*p*p;  /* mean difference per variable per object */

  RANDOUT;
}
