#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

int main (void)
{
  int i, n = 100;
  double data[n];
  double A[n/2+1]; 	// number of FFT modes

  gsl_fft_real_wavetable * real;
  gsl_fft_halfcomplex_wavetable * hc;
  gsl_fft_real_workspace * work;

  for (i = 0; i < n; i++)
    {
      data[i] = 0.0;
    }

  for (i = n / 3; i < 2 * n / 3; i++)
    {
      data[i] = 1.0;
    }
  for (i=0; i<n; i++) {
     //data[i] = sin(3.1415926/n*i);
  }
  for (i = 0; i < n; i++)
    {
      printf ("%d %e\n", i, data[i]);
    }
  printf ("\n");

  work = gsl_fft_real_workspace_alloc (n);
  real = gsl_fft_real_wavetable_alloc (n);

  gsl_fft_real_transform (data, 1, n, 
                          real, work);

  gsl_fft_real_wavetable_free (real);

  for (i = 0; i < n; i++)
    {
      printf ("%d %e\n", i, data[i]);
    }
  printf ("\n");


  A[0] = data[0];
  for (i=1; i<n-1; i+=2) {
     A[(i+1)/2] = data[i]*data[i]+ data[i+1]*data[i+1];
  }
  if (n%2==0) {    // n even
     A[n/2+1] = data[n-1]*data[n-1];
  }


  for (i = 0; i < n/2+1; i++)
    {
      printf ("%d %e\n", i, A[i]);
    }
  printf ("\n");

  for (i = 10; i < n; i++)
    {
      data[i] = 0;
    }

  hc = gsl_fft_halfcomplex_wavetable_alloc (n);

  gsl_fft_halfcomplex_inverse (data, 1, n, 
                               hc, work);
  gsl_fft_halfcomplex_wavetable_free (hc);

  for (i = 0; i < n; i++)
    {
      printf ("%d %e\n", i, data[i]);
    }

  gsl_fft_real_workspace_free (work);
  return 0;
}
