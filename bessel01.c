#include "gmp.h"
#include "mpfr.h"
#include <stdlib.h>
#include <stdio.h>

static void bessel01(int fl, mpfr_t x, int l, mpfr_t eps, mpfr_t s, mpfr_t t, mpfr_t u, mpfr_t v, 
mpfr_t euler, mpfr_t pi, mp_rnd_t rmode)
{
//  comp zeroth and 1st order bessel function
//  original program is witten by Ooura.
//  x is argument of bessel function: given by user
//  s=j0, t=j1, u=y0, v=y1, when l=1
//  s=i0, t=i1, u=k0, v=k1, when l=-1
//  eps is convergence limit: given by user
//  euler is Euler's constant: given by user
//  pi is ratio of diameter vs. circumference: given by user
//  rmode is rounding mode: given by user
static mpfr_t w, y, z, h, prc;
static mpfr_t w1, w2, w3, w4, w5;
static mpfr_t d1, d2, d3, d4, d5, d6;
static mp_prec_t size;
static unsigned int k, n; int m;
//int flag=0;

if(fl==0){
//    get the precision of variables
size=mpfr_get_prec(euler);
mpfr_init2(w,size); mpfr_init2(y,size); mpfr_init2(z,size); mpfr_init2(prc,size); mpfr_init2(h,size);
mpfr_init2(w1,size); mpfr_init2(w2,size); mpfr_init2(w3,size); mpfr_init2(w4,size); mpfr_init2(w5,size);
mpfr_init2(d1,size); mpfr_init2(d2,size); mpfr_init2(d3,size); mpfr_init2(d4,size); mpfr_init2(d5,size);
mpfr_init2(d6,size);
}

mpfr_log(prc,eps,rmode);
mpfr_neg(prc,prc,rmode);
mpfr_set_d(w1,0.6,rmode);
mpfr_mul(w1,w1,prc,rmode);
//mpfr_out_str(stdout, 10, 0, eps, rmode); putchar('\n');
//mpfr_out_str(stdout, 10, 0, prc, rmode); putchar('\n');
//mpfr_out_str(stdout, 10, 0, w1, rmode); putchar('\n');
//mpfr_out_str(stdout, 10, 0, x, rmode); putchar('\n');
if(mpfr_cmp(w1,x)>0){
//printf(" #1 in bs.c\n");
// set values
//    s=eps; t=0.0; u=0.0; v=0.0
 mpfr_set(s,eps,rmode); mpfr_set_d(t,0.0,rmode);
 mpfr_set_d(u,0.0,rmode); mpfr_set_d(v,0.0,rmode);
//printf(" #1a in bs.c\n");
// w=2/x
 mpfr_set_d(w,2.0,rmode);
 mpfr_div(w,w,x,rmode);
// y=0.0
 mpfr_set_d(y,0.0,rmode);
//  z=euler-log(w)
 mpfr_log(w1,w,rmode);
 mpfr_sub(z,euler,w1,rmode);
 if(mpfr_cmp(x,eps)>0){
// k=int(prc/15+((2+l)/6*x+((prc**2/7*x)**(1/3))
//printf(" #3a\n");
  mpfr_set_d(d1,15.0,rmode);
  mpfr_set_d(d2,2.0,rmode);
  mpfr_set_d(d3,6.0,rmode);
  mpfr_set_d(d4,7.0,rmode);
  mpfr_set_d(d5,1.0,rmode);
  mpfr_set_d(d6,3.0,rmode);
//  w1=prc/15
  mpfr_div(w1,prc,d1,rmode);
//  w2=x*((2+l)/6)
  mpfr_set_si(w2,l,rmode);
  mpfr_add(w2,w2,d2,rmode); 
  mpfr_div(w2,w2,d3,rmode); 
  mpfr_mul(w2,w2,x,rmode); 
//  w3=x*prc*prc/7
  mpfr_mul(w3,prc,prc,rmode);
  mpfr_mul(w3,w3,x,rmode);
  mpfr_div(w3,w3,d4,rmode);
//   w3=w3**(1/3)
  mpfr_div(w4,d5,d6,rmode);
  mpfr_pow(w3,w3,w4,rmode);
//  w1=w1+w2
  mpfr_add(w1,w1,w2,rmode);
//  w1=w1+w3
  mpfr_add(w1,w1,w3,rmode);
//  k=int(w1)
  mpfr_floor(w1,w1);
  k=mpfr_get_d(w1,GMP_RNDD);
//printf(" #2\n");
 }
 else{
  k=1;
 }
 if(l>0){
  for(n=2*k;n>=2;n-=2){
//  y=y+s
   mpfr_add(y,y,s,rmode);
// u=s/n-u
   mpfr_div_ui(w1,s,n,rmode);
   mpfr_sub(u,w1,u,rmode);
//  v=t*(n+1)/(n*(n+2))-v
//   w1=t*(n+1)
   mpfr_set_d(w1,1.0,rmode);
   mpfr_add_ui(w1,w1,n,rmode);
   mpfr_mul(w1,w1,t,rmode);
//   w2=n*(n+2)
   mpfr_set_d(w2,2.0,rmode);
   mpfr_add_ui(w2,w2,n,rmode);
   mpfr_mul_ui(w2,w2,n,rmode);
//   w1=w1/w2
   mpfr_div(w1,w1,w2,rmode);
//  v=w1-v
   mpfr_sub(v,w1,v,rmode);
// t=s*n*w-t
   mpfr_mul_ui(w1,s,n,rmode);
   mpfr_mul(w1,w1,w,rmode);
   mpfr_sub(t,w1,t,rmode);
// s=t*(n-1)*w-s
   mpfr_set_d(w1,-1.0,rmode);
   mpfr_add_ui(w1,w1,n,rmode);
   mpfr_mul(w1,w1,t,rmode);
   mpfr_mul(w1,w1,w,rmode);
   mpfr_sub(s,w1,s,rmode);
  }
// y=1/(2*y+s)
  mpfr_set_d(w1,1.0,rmode);
  mpfr_set_d(w2,2.0,rmode);
  mpfr_mul(w2,w2,y,rmode);
  mpfr_add(w2,w2,s,rmode);
  mpfr_div(y,w1,w2,rmode);
// s=s*y
  mpfr_mul(s,s,y,rmode);
// t=t*y
  mpfr_mul(t,t,y,rmode);
//   w4=(2/pi)
  mpfr_set_d(w4,2.0,rmode);
  mpfr_div(w4,w4,pi,rmode);
// u=(4*y*u+s*z)*(2/pi)
  mpfr_set_d(w1,4.0,rmode);
  mpfr_mul(w1,w1,y,rmode);
  mpfr_mul(w1,w1,u,rmode);
  mpfr_mul(w2,s,z,rmode);
  mpfr_add(w1,w1,w2,rmode);
  mpfr_mul(u,w1,w4,rmode);
// v=(4*y*v+t*z-s*w/2-t)*(2/pi)
//   w1=4*y*v
  mpfr_set_d(w1,4.0,rmode);
  mpfr_mul(w1,w1,y,rmode);
  mpfr_mul(w1,w1,v,rmode);
//   w2=t*z
  mpfr_mul(w2,t,z,rmode);
//   w3=s*w/2
  mpfr_set_d(w3,0.5,rmode);
  mpfr_mul(w3,w3,w,rmode);
  mpfr_mul(w3,w3,s,rmode);
//   w1=w1+w2
  mpfr_add(w1,w1,w2,rmode);
//   w1=w1-w3
  mpfr_sub(w1,w1,w3,rmode);
//   w1=w1-t
  mpfr_sub(w1,w1,t,rmode);
//   v=w1*w4
  mpfr_mul(v,w1,w4,rmode);
 }
 else{
  for(n=2*k;n>=2;n-=2){
// y=s+y
   mpfr_add(y,y,s,rmode);
// u=s/n+u
   mpfr_div_ui(w1,s,n,rmode);
   mpfr_add(u,w1,u,rmode);
// t=s*n*w+t
   mpfr_mul_ui(w1,s,n,rmode);
   mpfr_mul(w1,w1,w,rmode);
   mpfr_add(t,w1,t,rmode);
// s=t*(n-1)*w+s
   mpfr_set_d(w1,-1.0,rmode);
   mpfr_add_ui(w1,w1,n,rmode);
   mpfr_mul(w1,w1,t,rmode);
   mpfr_mul(w1,w1,w,rmode);
   mpfr_add(s,s,w1,rmode);
  }
// y=cosh(x)/(2*y+s)
  mpfr_cosh(w1,x,rmode);
  mpfr_set_d(w2,2.0,rmode);
  mpfr_mul(w2,w2,y,rmode);
  mpfr_add(w2,w2,s,rmode);
  mpfr_div(y,w1,w2,rmode);
// s=s*y
  mpfr_mul(s,s,y,rmode);
// t=t*y
  mpfr_mul(t,t,y,rmode);
// d1=1.5
  mpfr_set_d(w1,1.5,rmode);
  if(mpfr_cmp(x,w1)<0){
// u=4*y*u-s*z
   mpfr_set_d(w1,4.0,rmode);
   mpfr_mul(w1,w1,y,rmode);
   mpfr_mul(w1,w1,u,rmode);
   mpfr_mul(w2,s,z,rmode);
   mpfr_sub(u,w1,w2,rmode);
  }
  else{
// v=-x/2
   mpfr_set_d(v,0.5,rmode);
   mpfr_mul(v,x,v,rmode);
   mpfr_neg(v,v,rmode);
// w=v
   mpfr_set(w,v,rmode);
// u=exp(-x)/2
   mpfr_neg(w1,x,rmode);
   mpfr_exp(w1,w1,rmode);
   mpfr_set_d(w2,0.5,rmode);
   mpfr_mul(u,w1,w2,rmode);
// h=9/(prc+x)
   mpfr_set_d(w1,9.0,rmode);
   mpfr_add(w2,prc,x,rmode);
   mpfr_div(h,w1,w2,rmode);
// y=exp(-h)
   mpfr_neg(w1,h,rmode);
   mpfr_exp(y,w1,rmode);
// k=int(log(2+2*prc/x)/h)
   mpfr_set_d(d2,2.0,rmode);
   mpfr_set_d(d1,1.0,rmode);
   mpfr_div(w1,prc,x,rmode);
   mpfr_add(w1,w1,d1,rmode);
   mpfr_mul(w1,w1,d2,rmode);
   mpfr_log(w1,w1,rmode);
   mpfr_div(w1,w1,h,rmode);
   mpfr_floor(w1,w1);
   k=mpfr_get_d(w1,GMP_RNDD);
   for(n=1;n<=k;n++){
//  v=v/y
    mpfr_div(v,v,y,rmode);
// w=w*y
    mpfr_mul(w,w,y,rmode);
// u=u+exp(v+w)
    mpfr_add(w1,v,w,rmode);
    mpfr_exp(w1,w1,rmode);
    mpfr_add(u,u,w1,rmode);
   }
// u=u*h
   mpfr_mul(u,u,h,rmode);
  }
// v=(1/x-t*u)/s
   mpfr_set_d(w1,1.0,rmode);
   mpfr_div(w1,w1,x,rmode);
   mpfr_mul(w2,t,u,rmode);
   mpfr_sub(v,w1,w2,rmode);
   mpfr_div(v,v,s,rmode);
  }
 }
 else{
// w=1/(4*x)
  mpfr_set_d(w1,0.25,rmode);
  mpfr_div(w,w1,x,rmode);
// k=int(prc*0.12)
  mpfr_set_d(w1,0.12,rmode);
  mpfr_mul(w1,w1,prc,rmode);
  mpfr_floor(w1,w1);
  k=mpfr_get_d(w1,GMP_RNDD);
  for(m=-1;m<=3;m+=4){
// v=sqrt(8/pi*w)
   mpfr_set_d(w1,8.0,rmode);
   mpfr_mul(w1,w1,w,rmode);
   mpfr_div(w1,w1,pi,rmode);
   mpfr_sqrt(v,w1,rmode);
// y=v
   mpfr_set(y,v,rmode);
// z=0
   mpfr_set_d(z,0.0,rmode);

   for(n=2;n<=8*k+6;n+=4){
// v=v*(w*(1.0*m/n-n+2))
    mpfr_set_si(w1,m,rmode);
    mpfr_div_ui(w2,w1,n,rmode);
    mpfr_set_d(d2,2.0,rmode);
    mpfr_sub_ui(w2,w2,n,rmode);
    mpfr_add(w2,w2,d2,rmode);
    mpfr_mul(w2,w2,w,rmode);
    mpfr_mul(v,w2,v,rmode);
// z=v-z*l
    mpfr_set_si(w2,l,rmode);
    mpfr_mul(w2,w2,z,rmode);
    mpfr_sub(z,v,w2,rmode);
// v=v*(w*(m*1.0*/(n+2.0)-n))
    mpfr_add_ui(w2,d2,n,rmode);
    mpfr_div(w2,w1,w2,rmode);
    mpfr_sub_ui(w2,w2,n,rmode);
    mpfr_mul(w2,w2,w,rmode);
    mpfr_mul(v,w2,v,rmode);
// y=v-y*l
    mpfr_set_si(w2,l,rmode);
    mpfr_mul(w2,w2,y,rmode);
    mpfr_sub(y,v,w2,rmode);
   }
   if(m==-1){
//  t=y
    mpfr_set(t,y,rmode);
//  u=z
    mpfr_set(u,z,rmode);
   }
  }
  if(l>0){
// v=sin(x-pi/4)
   mpfr_set_d(w1,0.25,rmode);
   mpfr_mul(w1,w1,pi,rmode);
   mpfr_sub(w1,x,w1,rmode);
   mpfr_sin(v,w1,rmode);
// w=cos(x-pi/4)
   mpfr_cos(w,w1,rmode);
// s=v*u+w*t
   mpfr_mul(w1,v,u,rmode);
   mpfr_mul(w2,w,t,rmode);
   mpfr_add(s,w1,w2,rmode);
// u=v*t-w*u
   mpfr_mul(w1,v,t,rmode);
   mpfr_mul(w2,w,u,rmode);
   mpfr_sub(u,w1,w2,rmode);
// t=v*y-w*z
   mpfr_mul(w1,v,y,rmode);
   mpfr_mul(w2,w,z,rmode);
   mpfr_sub(t,w1,w2,rmode);
// v=-v*z-w*y
   mpfr_mul(w1,v,z,rmode);
   mpfr_mul(w2,w,y,rmode);
   mpfr_add(v,w1,w2,rmode);
   mpfr_neg(v,v,rmode);
  }
  else{
// v=exp(x)/2
   mpfr_set_d(w1,0.5,rmode);
   mpfr_exp(v,x,rmode);
   mpfr_mul(v,v,w1,rmode);
// w=pi/(4*v)
   mpfr_set_d(w1,0.25,rmode);
   mpfr_mul(w1,w1,pi,rmode);
   mpfr_div(w,w1,v,rmode);
// s=v*(t-u)
   mpfr_sub(w1,t,u,rmode);
   mpfr_mul(s,v,w1,rmode);
// u=w*(t+u)
   mpfr_add(w1,t,u,rmode);
   mpfr_mul(u,w,w1,rmode);
// t=v*(y-z)
   mpfr_sub(w1,y,z,rmode);
   mpfr_mul(t,v,w1,rmode);
// v=w*(y+z)
   mpfr_add(w1,y,z,rmode);
   mpfr_mul(v,w,w1,rmode);
 }
}

//printf(" # end\n");
if(fl<0){
mpfr_clear(w);mpfr_clear(y);mpfr_clear(z);mpfr_clear(prc);mpfr_clear(h);
mpfr_clear(w1);mpfr_clear(w2);mpfr_clear(w3);mpfr_clear(w4);mpfr_clear(w5);
mpfr_clear(d1);mpfr_clear(d2);mpfr_clear(d3);mpfr_clear(d4);mpfr_clear(d5);
mpfr_clear(d6);
}

}
