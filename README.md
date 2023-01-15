# bessel01

bessel function program with mpfr

static void bessel01(int fl, mpfr_t x, int l, mpfr_t eps, mpfr_t s, mpfr_t t, mpfr_t u, mpfr_t v, 
mpfr_t euler, mpfr_t pi, mp_rnd_t rmode)

//  comp zeroth and 1st order bessel function

//  original program is witten by Ooura.

//  x is argument of bessel function: given by user

//  s=j0, t=j1, u=y0, v=y1, when l=1

//  s=i0, t=i1, u=k0, v=k1, when l=-1

//  eps is convergence limit: given by user

//  euler is Euler's constant: given by user

//  pi is ratio of diameter vs. circumference: given by user

//  rmode is rounding mode: given by user



#ref 

Ooura's Mathematical Software Packages https://www.kurims.kyoto-u.ac.jp/~ooura/index.html

MPFR https://www.kurims.kyoto-u.ac.jp/~ooura/index.html
