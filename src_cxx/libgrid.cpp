#include <cmath>
#include <cstdio>
#include <algorithm>
#include <iostream>
#include <thread>
#include <vector>

#include "mathlib.hpp"
using namespace std;



// ******************************************************************************************* //

template<typename T>
inline void gwind0(T* __restrict__ gwx, T* __restrict__ gwy, T const gwid, int const nxa,
		   int const nxb, int const nya, int const nyb)
{
  double const wid = gwid*0.6005612;
  
  if (wid > 0) {
    double const xcen = (nxa + nxb)/2;

    double const ycen = (nya + nyb)/2; 
    for (int i = nxa; i <= nxb; ++i) {
      double const xq = (i - xcen)/wid;
      gwx[i] = exp(-(xq*xq));
    } 
    for (int i = nya; i <= nyb; ++i) {
      double const xq = (i - ycen)/wid;
      gwy[i] = exp(-(xq * xq));
    } 
  } else {
    for (int i = nxa; i <= nxb; ++i) 
      gwx[i] = 1.0;
    for (int i = nya; i <= nyb; ++i)
      gwy[i] = 1.0; 
  }
}

// ******************************************************************************************* //
template<typename T, typename U>
inline double averag(const T* const __restrict__ p, int const nxa, int const nxb, int const nya,
	     int const nyb, int const nxs, int const nys, int const idx,
	     int idy, const U* const __restrict__ gx, const U* const __restrict__ gy)
/* finds weighted average of array m over the block defined */
{
  double sum=0.0, sumg=0.0, sumgx=0.0;

  /* fix limits so sum doesn't run off edge of image */
  int const nxc = (nxa + idx < 0)? -idx: nxa;
  int const nyc = (nya + idy < 0)? -idy: nya;
  int const nxd = (nxb + idx > nxs)? nxs - idx: nxb;
  int const nyd = (nyb + idy > nys)? nys - idy: nyb;

  for (int i = nxc; i < nxd; i++)
    sumgx += gx[i];
  /* weighted sum in window */
  for (int j = nyc; j < nyd; j++) {
    double sumx = 0.0;
    int const jj = idx + nxs*(j + idy);
    for (int i = nxc; i < nxd; i++)
      sumx += gx[i]*p[i + jj];
    sum += gy[j]*sumx;
    sumg += gy[j]*sumgx;
  } /* end of j loop */
  sum /= sumg;
  return sum;
}

// ******************************************************************************************* //
template<typename T, typename U>
inline void unbias(const T* const m1, const U* const m2, int nxa, int nxb,
		   int nya, int nyb, int nxs, int nys,
		   U *gx, U *gy, U &av1, U &av2, U &cx,
		   U &cy, U &cxx, U &cxy, U &cyy,
		   int idelx, int idely)
{ 
  /*  find weighted means of m1 & m2 over the window 
      sets up quadratic fit to average of m2 as a fcn. of offsets */
  av1 = averag(m1, nxa, nxb, nya, nyb, nxs, nys, 0, 0, gx, gy); 
  double const t0 = averag(m2, nxa, nxb, nya, nyb, nxs, nys, idelx, idely, gx, gy); 
  double const t1 = averag(m2, nxa, nxb, nya, nyb, nxs, nys, idelx + 1, idely, gx, gy); 
  double const t2 = averag(m2, nxa, nxb, nya, nyb, nxs, nys, idelx - 1, idely, gx, gy); 
  double const t3 = averag(m2, nxa, nxb, nya, nyb, nxs, nys, idelx, idely + 1, gx, gy); 
  double const t4 = averag(m2, nxa, nxb, nya, nyb, nxs, nys, idelx, idely - 1, gx, gy); 
  double const t5 = averag(m2, nxa, nxb, nya, nyb, nxs, nys, idelx + 1, idely + 1, gx, gy); 
  av2 = t0; 
  cx = 0.5*(t1 - t2); 
  cy = 0.5*(t3 - t4); 
  cxx = 0.5*(t1 - 2*t0 + t2); 
  cyy = 0.5*(t3 - 2*t0 + t4); 
  cxy = t5 + t0 - t1 - t3;
}

// ******************************************************************************************* //

template<typename T>
inline int getmin(T *p, T &x0, T &y0)
{
  double	f11, f12, f13, f21, f22, f23, f31, f32, f33;
  double	fx, fy, t, fxx, fyy, fxy;

  /* find the min, p points to a 3x3 array */
  f11 = *p++;	f21 = *p++;	f31 = *p++;
  f12 = *p++;	f22 = *p++;	f32 = *p++;
  f13 = *p++;	f23 = *p++;	f33 = *p++;

  fx = 0.5 * ( f32 - f12 );	fy = 0.5 * ( f23 - f21 );
  t = 2.* ( f22 );
  fxx =  f32 + f12 - t;	fyy = f23 + f21 - t;
  /* find in which quadrant the minimum lies */
  if (f33 < f11) {
    if (f33 < f31) {
      if (f33 < f13)
	fxy = f33+f22-f32-f23;
      else
	fxy = f23+f12-f22-f13;
    } else {
      if (f31 < f13)
	fxy = f32+f21-f31-f22;
      else
	fxy = f23+f12-f22-f13;
    }
  } else { 
    if (f11 < f31) {
      if (f11 < f13) 
	fxy = f22+f11-f21-f12;
      else
	fxy = f23+f12-f22-f13;
    } else {
      if (f31 < f13)
	fxy = f32+f21-f31-f22;
      else
	fxy = f23+f12-f22-f13;
    }
  }
  t = -1./(fxx *fyy - fxy *fxy);
  x0 = t * (fx * fyy - fy * fxy);
  y0 = t * (fy * fxx - fx * fxy);
  if (std::abs(x0) >= 0.75 || std::abs(y0) >= 0.75) {
    x0 = -fx/fxx;
    y0 = -fy/fyy;
  }
  return 1;
}

// ******************************************************************************************* //

template <typename T, typename U>
inline double resid(T *m1, T *m2, int idx, int idy, int nxa,
	    int nxb, int nya, int nyb, int nxs,
	    int nys, int ndmx, U *gx, U *gy, U bs)
{
  int	nxc, nxd, nyc, nyd, nx, ny;
  U 	*p1, *p2, *ps;
  double	sum, sumx, t, ndmx2;
  int	i, j;
  double   sumg;
  //static int	mxc, mxd, myc, myd;
  //static double	gsum;

  /*set up limits */
  nxc = nxa;
  if (nxc + idx < 0)
    nxc = -idx;
  nyc = nya;
  if (nyc + idy < 0)
    nyc = - idy;
  nxd = nxb;
  if (nxd + idx >= nxs)
    nxd = nxs - idx - 1;
  nyd = nyb;
  if (nyd + idy >= nys)
    nyd = nys - idy - 1;
  sum = sumg = 0.0;

  nx = nxd - nxc +1;
  p2 = gy + nyc;
  ps = gx + nxc;

  // if (nxc != mxc || nxd != mxd || nyc != myc || nyd != myd) {
    /* sum gaussians over rectangle to get normalization */
    /* (only if limits change)*/
    j = nyd -nyc + 1;
    if (j <= 0 || nxd - nxc + 1 <= 0)
      return -1;		/* added 19feb95 LS */
    while (j) {
      i = nx;
      p1 = ps;
      while (i) {
	sumg += (*p1++) * (*p2);
	i--;
      }
      p2++;
      j--;
    }
    //gsum = sumg;
    //mxc = nxc;
    //mxd = nxd;
    //myc = nyc;
    //myd = nyd;
    // } else
    //  sumg = gsum;

  m1 += nyc*nxs + nxc;
  m2 += (nyc + idy)*nxs + nxc + idx;
  ny = nxs - nx;		/*residual increment after inner loop */
  p2 = gy + nyc;
  
  /* now the loop to compute the residual */
  j = nyd - nyc +1;
  ndmx2 = ndmx*ndmx;
  while (j) {
    i = nx;
    p1 = ps;
    sumx = 0.0;
    while (i) {
      t = *m1++ - *m2++;
      t = t + bs;
      t = t*t;
      t = std::min<double>(t, ndmx2);
      sumx += (*p1++) * t;
      i--;
    }
    sum += (*p2++) * sumx;
    m1 += ny;
    m2 += ny;
    j--;
  }
  /*return normalized residual */
  sum /= sumg;
  return sum;
}

// ******************************************************************************************* //

template<typename T, typename U>
inline void match_1(T *p1, T *p2, int nxa, int nxb,
		    int nya, int nyb, int nx, int ny,
		    U *gwx, U *gwy, U &xoffset, U &yoffset, int const& stretch_clip, int &badmatch)
{
  int idelx, idely, i, j, k, ndmx = 1000, done[9]={};
  int	di, dj, in, jn, iter, dd, badflag = 0;
  double	av1, av2, cx, cy, cxx, cxy, cyy, avdif, t, res[9], buf[9]={}, t1, t2;
  int const	itmax = 20;
  
  idelx = rint(xoffset);
  idely = rint(yoffset); 
  unbias(p1, p2, nxa, nxb, nya, nyb, nx, ny, gwx, gwy, av1, av2, cx, cy,
	 cxx, cxy, cyy, idelx, idely);
  /* look at a 3x3 matrix of residuals centered at 0 offset, find the location
     of the minimum, if not at center, then look at new 3x3 centered
     on the edge minimum; repeat until min at center */
  iter = itmax;
  badflag = 0;
  while (iter--) {
    for (k = 0; k < 9; k++) {
      if (done[k] == 0) {
	i = idelx + (k % 3) - 1;
	j = idely + (k / 3) - 1;
	avdif = av2 +  i*cx + j*cy + i*i*cxx + i*j*cxy + j*j*cyy - av1;
	res[k] = resid(p1, p2, i, j, nxa, nxb, nya, nyb, nx, ny, ndmx,
		       gwx, gwy, avdif);
      }
    }
    t = res[0];
    i = 0;
    for (k = 1; k < 9; k++) 
      if (res[k] < t) {
	t = res[k];
	i = k;
      }
    if (t < 0) {		/* added LS 19feb95 */
      printf("match - ran out of data at edge\n");
      badflag = 1;
      break;
    }
    idelx += (i % 3) - 1;
    idely += (i / 3) - 1;
    /* check if we have gone too far */
    if (std::abs(idelx) > stretch_clip || std::abs(idely) > stretch_clip) {
      badflag++;
      break;
    }
    if (i == 4)
      break;			/* done if it is the center one */
    /* not in center, shuffle what we have to put the edge min in center */
    di = (i % 3) - 1;
    dj = (i / 3) - 1;
    dd = dj * 3 + di;
    for (k = 0; k < 9; k++) {
      in = k%3 + di;
      jn = k/3 + dj;
      if (in >= 0 && jn >= 0 && in < 3 && jn < 3) { /* in range */
	done[k] = 1;
	buf[k] = res[k + dd];
      } else 
	done[k] = 0;		/* not in range, mark undone */
    }
    for (k = 0; k < 9; k++)
      res[k] = buf[k];		/* put back in res array */
  } /* end of iter while */
  /* done or reached itmax, which ? */
  if (iter <= 0) {
    badflag++;
  }
  if (badflag) {
    badmatch++;
    xoffset = yoffset = 0;
    return;
  }
				/* must have been OK so far */
  getmin(res, t1, t2);
  xoffset = idelx + t1;
  yoffset = idely + t2;
}

// ******************************************************************************************* //

template<typename T, typename U>
inline void do_one(int const ii, T* __restrict__ gwx, T* __restrict__ gwy, T gwid, 
		   const U* const __restrict__ p1, const U* const __restrict__ p2, int const stretch_clip, int const dx2, int const dy2,
		   int const nyg, int const nxg, int const nx, int const ny, U* __restrict__ out,
		   const int* const __restrict__ gx, const int* const __restrict__ gy)
{
  int badmatch = 0;

  int  const i1 = std::max<int>(0,gx[ii] - dx2) ;
  int  const i2 = std::min<int>(nx,gx[ii] + dx2) - 1;

  int  const j1 = std::max<int>(gy[ii] - dy2,0);
  int  const j2 = std::min<int>(gy[ii] + dy2, ny) - 1;

  U xoffset = 0, yoffset = 0;
  
  gwind0(gwx, gwy, gwid, i1, i2, j1, j2); /* get gaussian kernels */
  
  match_1(p1, p2, i1, i2, j1, j2, nx, ny, gwx, gwy,
	  xoffset, yoffset, stretch_clip, badmatch); /* get offsets */

  out[2*ii]   = xoffset;
  out[2*ii+1] = yoffset; 
}

// ******************************************************************************************* //

template<typename T, typename U>
void gridmatch(int const ny, int const nx, int const nyg, int const nxg, const T* const __restrict__ p1,
	       const T* const __restrict__ p2, const int* const gy, const int* const gx, int const dy,
	       int const dx, U const gwid, int stretch_clip, U* __restrict__ out, int const nthreads)
{
  
  /* the call is offsets = gridmatch(m1,m2,gx,gy,dx,dy,gwid,stretch_clip)
	where	m1 = reference input image
	m2 = image to compare with m1, m1 and m2 must be same size
	gx = array of x gridpoints
	gy = array of y gridpoints, gx and gy must have same size
	dx and dy are the window size, and gwid is the gaussian mask width
        stretch_clip = maximum allowed displacement before a bad match
          is declared

  Authors:  Richard A. Shine (original)
            Louis H. Strous (port from ANA to IDL)
            Lockheed-Martin Advanced Technology Center, Palo Alto, CA, USA

	    Parallel version adapted by J. de la Cruz Rodriguez (ISP-SU, 2021)

*/

  //double t0 = getTime();
  // 
  double	        /*xoffset, yoffset,*/ *gwx=NULL, *gwy=NULL;
  //int      i1, i2, j1, j2;
  //int  const dims[3] = {2,nxg,nxg};

  if (stretch_clip < 2)
    stretch_clip = 2;
  stretch_clip--;

  /* prepare the gaussian kernel */

  int const nc = nxg*nyg;			/* the number of subimages */
  int const dx2 = dx/2;
  int const dy2 = dy/2;
  int ii;
#pragma omp parallel default(shared) firstprivate(ii,gwx,gwy) num_threads(nthreads)
  {
    gwx = (double *) malloc((nx + ny)*sizeof(double));
    gwy = gwx + nx;
  
#pragma omp for schedule(static) 
    for( ii=0; ii<nc;++ii) {	
      do_one<double,T>(ii, gwx, gwy, gwid,  p1, p2, stretch_clip, dx2, dy2, nyg, nxg, nx, ny, out, gx, gy);
    }
    free(gwx);
  }
 
}

// ******************************************************************************************* //
template<typename T>
void stretch_matrix(int const ny, int const nx, int const npy, int const npx,
		    const T* const __restrict__ gr, T* const  __restrict__ out_y, T* const  __restrict__ out_x, int const nthreads)
{  
  using fp = double;
  

  int   jy(0), j1(0), j2(0);
  int   ix(0), iy(0), i1(0), i2(0), jx(0); //, intx, inty;
  fp    y(0), dy(0), dy1(0), dyp(0), dyp1(0);
  fp	x(0), dx(0), dx1(0), dxp(0), dxp1(0);
  fp	w1(0), w2(0), w3(0), w4(0), xl(0), yl(0);
  const T 	*jpbase(nullptr), *jbase(nullptr);

  int const n = nx;
  int const m = ny;

  const T* const __restrict__ xgbase = gr;
  int const nxg = npx;
  int const nyg = npy;
  int const nxgm = nxg - 1;
  int const nygm = nyg - 1;
  int const nx1 = nx-1;
  int const ny1 = ny-1;
  
  /* linearly interpolate the displacement grid values over array */
  /* similar to regrid3 in inner part */
  fp const xd = (double) n/nxg; // number of pixels in the image per grid cell
  fp const xinc = 1.0/xd; // maps the size of a pixel of the image in the matrix
  fp const xs   = xinc + (xd - 1.0)/(2.0*xd);
  fp const yd = (double) m/nyg;
  fp const yinc = 1.0/yd;
  fp const y0 = yinc + (yd - 1.0)/(2.0*yd);
  y = y0;
  
#pragma omp parallel default(shared) firstprivate(ix, iy, i1, i2, jx, jy, j1, j2, w1, w2, w3, w4, jpbase, jbase, x, y, dy, dy1, dx, dx1, xgbase, xl, yl, dyp, dyp1, dxp, dxp1) num_threads(nthreads)
    {
#pragma omp for schedule(static) 
      for (iy = 0; iy < m; ++iy) {
	y = y0 + iy * yinc;
	x = xs;
	jy = y;
	dy = y - jy;
	dy1 = 1.0 - dy;
	
	if (jy < 1)
	  j1 = j2 = 0;
	else if (jy >= nyg)
	  j1 = j2 = nygm;
	else {
	  j1 = jy - 1;
	  j2 = j1 + 1;
	}
	
	jbase  = xgbase + j1 * 2 * nxg;
	jpbase = xgbase + j2 * 2 * nxg;
	
	for (ix = 0; ix < n; ++ix) {
	  jx = x;
	  dx = x - jx;
	  dx1 = 1.0 - dx;
	  
	  if (jx < 1)
            i1 = i2 = 0;
	  else if (jx >= nxg)
            i1 = i2 = nxgm;
	  else {
            i1 = jx - 1;
            i2 = i1 + 1;
	  }
	  
	  // --- bilinear interpolation of the coordinates for a pixel in the image --- //
	  
	  w1 = dy1*dx1;
	  w2 = dy1*dx;
	  w3 = dy*dx1;
	  w4 = dy*dx;
	  
	  i1 = 2*i1;
	  i2 = 2*i2;
	  
	  xl = std::max<double>(std::min<double>(nx1, w1 * *(jbase+i1) + w2 * *(jbase+i2) + w3 * *(jpbase+i1) + w4 * *(jpbase+i2) + (double)ix), 0);
	  
	  i1 += 1;
	  i2 += 1;
	  
	  yl =  std::max<double>(std::min<double>(ny1, w1 * *(jbase+i1) + w2 * *(jbase+i2) + w3 * *(jpbase+i1) + w4 * *(jpbase+i2) + (double)iy), 0);
	  
	  out_x[iy*n+ix] = xl;
	  out_y[iy*n+ix] = yl;
	  
	  x += xinc;
	}
	// y += yinc;
      } // iy    
    } // pragma
}

void stretch_matrix_double(int const ny, int const nx,  int const npy, int const npx,
		    const double* const __restrict__ gr, double* const __restrict__ out_y, double* const __restrict__ out_x, int const nthreads)
{
  stretch_matrix<double>(ny,nx,npy,npx,gr,out_y,out_x,nthreads);
}

// ******************************************************************************************* //

template<typename T>
void stretch(int const ny, int const nx, const T* const __restrict__ im, int const npy, int const npx,
	   const T* const __restrict__ gr, T* const __restrict__ res, int const nthreads)
{  
  using fp = double;
  

  int   jy(0), j1(0), j2(0);
  int   ix(0), iy(0), i1(0), i2(0), jx(0); //, intx, inty;
  fp    y(0), dy(0), dy1(0), dyp(0), dyp1(0);
  fp	x(0), dx(0), dx1(0), dxp(0), dxp1(0);
  fp	w1(0), w2(0), w3(0), w4(0), xl(0), yl(0);
  const T 	*jpbase(nullptr), *jbase(nullptr);
  
  int const n = nx;
  int const m = ny;
  
  const T* const __restrict__ xgbase = gr;
  int const nxg = npx;
  int const nyg = npy;
  int const nxgm = nxg - 1;
  int const nygm = nyg - 1;
  int const nx1 = nx-1;
  int const ny1 = ny-1;
  
  /* linearly interpolate the displacement grid values over array */
  /* similar to regrid3 in inner part */
  fp const xd = (double) n/nxg; // number of pixels in the image per grid cell
  fp const xinc = 1.0/xd; // maps the size of a pixel of the image in the matrix
  fp const xs   = xinc + (xd - 1.0)/(2.0*xd);
  fp const yd = (double) m/nyg;
  fp const yinc = 1.0/yd;
  fp const y0 = yinc + (yd - 1.0)/(2.0*yd);
  y = y0;
  
#pragma omp parallel default(shared) firstprivate(ix, iy, i1, i2, jx, jy, j1, j2, w1, w2, w3, w4, jpbase, jbase, x, y, dy, dy1, dx, dx1, xgbase, xl, yl, dyp, dyp1, dxp, dxp1) num_threads(nthreads)
  {
#pragma omp for schedule(static) 
    for (iy = 0; iy < m; ++iy) {
      y = y0 + iy * yinc;
      x = xs;
      jy = y;
      dy = y - jy;
      dy1 = 1.0 - dy;
      
      if (jy < 1)
        j1 = j2 = 0;
      else if (jy >= nyg)
        j1 = j2 = nygm;
      else {
        j1 = jy - 1;
        j2 = j1 + 1;
      }
      
      jbase  = xgbase + j1 * 2 * nxg;
      jpbase = xgbase + j2 * 2 * nxg;
      
      for (ix = 0; ix < n; ++ix) {
        jx = x;
        dx = x - jx;
        dx1 = 1.0 - dx;
        
        if (jx < 1)
	  i1 = i2 = 0;
        else if (jx >= nxg)
	  i1 = i2 = nxgm;
        else {
	  i1 = jx - 1;
	  i2 = i1 + 1;
        }
	
        // --- bilinear interpolation of the coordinates for a pixel in the image --- //
        
        w1 = dy1*dx1;
        w2 = dy1*dx;
        w3 = dy*dx1;
        w4 = dy*dx;
        
        i1 = 2*i1;
        i2 = 2*i2;
        
        xl = std::max<double>(std::min<double>(nx1, w1 * *(jbase+i1) + w2 * *(jbase+i2) + w3 * *(jpbase+i1) + w4 * *(jpbase+i2) + (double)ix), 0);
        
        i1 += 1;
        i2 += 1;
        
        yl =  std::max<double>(std::min<double>(ny1, w1 * *(jbase+i1) + w2 * *(jbase+i2) + w3 * *(jpbase+i1) + w4 * *(jpbase+i2) + (double)iy), 0);
        
	
        // --- Now bilinear interpolation of the actual image --- //
        
        // -- bracket which pixels of the image must be taken --- //
        
        i2 = std::min<int>(std::max<int>(int(xl), 0), nx1-1);
        j2 = std::min<int>(std::max<int>(int(yl), 0), ny1-1);
        
        dxp  = xl - i2;
        dxp1 = 1.0 - dxp; 
        dyp  = yl - j2;
        dyp1 = 1.0  - dyp;

        w1 = *(im+j2*nx+i2);
        w2 = *(im+j2*nx+i2+1);
        w3 = *(im+(j2+1)*nx+i2);
        w4 = *(im+(j2+1)*nx+i2+1);
        
        res[iy*nx+ix] = w1*(dxp1*dyp1) + w2*(dxp*dyp1) + w3*(dxp1*dyp) + w4*(dxp*dyp);
        
        x += xinc;
      }
      // y += yinc;
    } // iy    
  } // pragma
  
}

// ******************************************************************************************* //

void stretch_double(int const ny, int const nx, const double* const __restrict__ im, int const npy, int const npx,
		       const double* const __restrict__ gr, double* const __restrict__ res, int const nthreads)
{  
  stretch<double>(ny, nx, im, npy, npx, gr, res, nthreads);
}

// ******************************************************************************************* //

void dsgridnest(int const ny, int const nx, const double* const __restrict__ im1, const double* const __restrict__ im2, int const nTile, const int* const __restrict__ tile, const int* const __restrict__ clip, double* const __restrict__ res, int const nthreads)
{
  using fp_t = double;
  
  int nprev = 0;
  std::vector<fp_t> displprev, prev, displ;
  std::vector<int> igx(nTile,0), igy(nTile,0);

  
    for(int k=0; k<nTile; ++k){
        
        int const n            = tile[k];
        int const stretch_clip = clip[k];
        int const ngw = int((2*nx)/float(n));
        int const nw  = int(1.25f * ngw);
        
        int const ngx = ((nx > ny) ? n : int(float(n*nx)/ny +0.5));
        int const ngy = ((nx > ny) ? int(float(n*ny)/nx +0.5) : n);
        
        igx[k] = ngx;
        igy[k] = ngy;
        
        fp_t const wx = fp_t(nx) / ngx;
        fp_t const wy = fp_t(ny) / ngy;
        
        int* const __restrict__ pgx = new int [ngx*ngy]();
        int* const __restrict__ pgy = new int [ngx*ngy]();
        
        // --- make subfields grid --- //
        
        for(int jj = 0; jj<ngy; ++jj){
            for(int ii = 0; ii<ngx; ++ii){
                pgy[jj*ngx+ii] = int(jj*wy+wy/2-1);
                pgx[jj*ngx+ii] = int(ii*wx+wx/2-1);
            }
        }

        int const dx = nw;
        int const dy = nw;
        fp_t const gwid = ngw;

        displ.resize(2*ngx*ngy,fp_t(0));
        
        if( k == 0 ){
            gridmatch<fp_t>(ny, nx, ngy, ngx, im1, im2, pgy, pgx, dy, dx, gwid, stretch_clip, &displ[0], nthreads);
        } else {
            if( n != nprev ){
                int const ngtotprev = igx[k-1]*igy[k-1];
                fp_t* const __restrict__ disx = new fp_t [ngtotprev];
                fp_t* const __restrict__ disy = new fp_t [ngtotprev];

                for(int ii=0; ii<ngtotprev; ++ii){
                    disx[ii] = displprev[2*ii];
                    disy[ii] = displprev[2*ii+1];
                }
                
                const fp_t* const __restrict__ fx = mth::congrid<fp_t,fp_t>(igy[k-1],igx[k-1], disx, ngy, ngx);
                const fp_t* const __restrict__ fy = mth::congrid<fp_t,fp_t>(igy[k-1],igx[k-1], disy, ngy, ngx);
		
                delete [] disx;
                delete [] disy;
                
                int const ngtot = ngx*ngy;
                prev.resize(2*ngtot,0);
                
                for(int ii =0; ii<ngtot; ++ii){
                    prev[2*ii] = fx[ii];
                    prev[2*ii+1] = fy[ii];
                }

                delete [] fx;
                delete [] fy;
                
            } else prev = displprev;

            {
                int const ngtot = 2*ngx*ngy;
                fp_t*  __restrict__ displnew = new fp_t [ngtot];
                
                fp_t* im3 = new fp_t [nx*ny]();

		stretch<fp_t>(ny, nx, im2, ngy, ngx, &prev[0],im3, nthreads);
                gridmatch<fp_t,fp_t>(ny, nx, ngy, ngx, im1, im3, pgy, pgx, dy, dx, gwid, stretch_clip, displnew, nthreads);
                delete [] im3;
                
                for(int ii=0; ii<ngtot; ++ii) displ[ii] = prev[ii] + displnew[ii];
                delete [] displnew;
            }
                    
            delete [] pgx;
            delete [] pgy;
        }
        if( k < (nTile-1) ){
            nprev = n;
            displprev = displ;
        }


    
    }
  
    // --- Now store result in a real IDL array --- //
    
    int const ngtot = igy[nTile-1]*igx[nTile-1]*2;
    for(int ii=0; ii<ngtot; ++ii){
      fp_t const tmp = (std::isnan(displ[ii])? 0.0 : displ[ii]);
      res[ii] = (std::isinf(displ[ii])? 0.0 : tmp);
    }
    
    
}

  
  
