#ifndef IMTOOLS
#define IMTOOLS

void nearest2D(int const ny, int const nx, double* __restrict__ y, double* __restrict__ x, double* __restrict__ d_in, int const ny1, int const nx1, double* __restrict__ yy_in, double* __restrict__ xx_in, double* __restrict__ res_in, int const nthreads, double const missing);

void bilint2D(int const ny, int const nx, double* __restrict__ y, double* __restrict__ x, double* __restrict__ d_in, int const ny1, int const nx1, double* __restrict__ yy_in, double* __restrict__ xx_in, double* __restrict__ res_in, int const nthreads, double const missing);


void bilint_fast2D(int const ny, int const nx,
		   const double* const __restrict__ d_in, int const ny1, int const nx1, const double* const __restrict__ yy_in,
		   const double* const __restrict__ xx_in,  double* const __restrict__ res_in, int const nthreads);
#endif
