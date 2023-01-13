#ifndef LIBGRIDHPP
#define LIBGRIDHPP


void dsgridnest(int const ny, int const nx, const double* const __restrict__ im1, const double* const __restrict__ im2,
		int const nTile, const int* const __restrict__ tile, const int* const __restrict__ clip,
		double* const __restrict__ res, int const nthreads);
  
void stretch_double(int const ny, int const nx, const double* const __restrict__ im, int const npy, int const npx,
		    const double* const __restrict__ gr, double* const __restrict__ res, int const nthreads);

void stretch_matrix_double(int const ny, int const nx, int const npy, int const npx,
			   const double* const __restrict__ gr, double*  __restrict__ out_y,
			   double*  __restrict__ out_x, int const nthreads);
#endif


