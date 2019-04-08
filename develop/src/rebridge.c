/*
    program:	rebridge.c
    author:	Pieter J. in 't Veld for MIT Boston
    date:	May 4, 2001, January 28, 2002.
    purpose:	Trimer rebridging module, using Mavrantzas et al.,
    		Macromolecules 32, 5072 (1999).

    Notes:
      20010514	Check if spherical molecular notation corresponds with
		Mavrantas notation.
      20020128	Fixed angular incongruities in spherical coordinates.
      20090421  Mavrantas sphercal coordinates notation is different
		from the one used in Allen & Tiledesley, and the latter
		is adopted in our code.  Therefore, we need to adjust
		the implementation of Mavrantas's formulas in our coding.
		The difference is this, Mavrantas specified the bond angle
		theta_i as the angle at bead i, while A&T specified the bond
		angle theta_i as the angle at bead i-1. (PY)
		Therefore in function RebridgeSetup, A&T notation is used.
		Thus the Mavrantas's formula needs to be adjusted, e.g., the
		calculation of r_P in function Feasible().
      20091205  Added Jacobian calculation J(IV->I).

    Reference:  1. Macromolecules 1999, v32, 5072 
	        2. J.Chem.Phys. 2002, v117, 5465 
*/
#define __REBRIDGE_MODULE
#include "rebridge.h"


long Frame(
	vector *p0, vector *p1, vector *r_ref, 
	vector *x1, vector *x2, vector *x3)
{
  vector		t1, t2;
  double		f;
    
  t1			= V_Subtr(r_ref, p0);
  t2			= V_Subtr(p1, p0);
  if (!(f = sqrt(V_Dot(&t2, &t2)))) return 1;
  *x1			= V_Mult(1.0/f, &t2);
  *x2			= V_Mult(V_Dot(&t1, &t2), x1);
  *x3			= V_Cross(&t1, &t2);
  t2			= V_Mult(-f, &t1);
  *x2			= V_Add(x2, &t2);
  if (!(f = sqrt(V_Dot(x3, x3)))) return 1;
  *x2			= V_Mult(1.0/f, x2); 
  *x3			= V_Mult(1.0/f, x3);
  return 0;
}


double			R; // radius of C3
vector			u1, u2, u3,
			v1, v2, v3,
			w1, w2, w3;  // coordinate frame
vector			rN;

long Feasibility(
	vector *r, double l, double a, vector *u1, 
	double *phi_low, double *phi_high) // returns n_entries in phi interval
{
  static long		n_roots, n_real, i, j;
  static double		f, f2,
    			A, B, C, D, E, F,
  			c0, c1, c2, c3, c4, c5,
			c[5], phi[4];
  static double complex	z[4];
  static vector		t1;
    
  t1			= V_Subtr(r, &rN);
  c0			= V_Dot(&t1, &w2);
  c1			= a*V_Dot(&w2, u1);
  f			= 2/R;
  f			*= f;
  A			= f*(c0*c0+c1*c1);
  c2			= V_Dot(&t1, &w3);
  c3			= a*V_Dot(&w3, u1);
  B			= f*(c2*c2+c3*c3);
  C			= 2.0*f*(c1*c3+c0*c2);
  c4			= V_Dot(&t1, &t1);
  c5			= 2.0*a*V_Dot(&t1, u1);
  f			/= R;
  f2			= c4+R*R;
  F			= 0.25*f/R*(f2*(f2-2.0*(a*a+l*l))+
    				(a*a-l*l)*(a*a-l*l)+c5*c5);
  f2			-= a*a+l*l;
  D			= f*(f2*c0+c1*c5);
  E			= f*(f2*c2+c3*c5);
  
  c[0]			= -F-A+D;
  c[1]			= 2.0*(E-C);
  c[2]			= 2.0*(A-2.0*B-F);
  c[3]			= 2.0*(C+E);
  c[4]			= -(F+A+D);
  n_roots		= R_SolveFourthOrder(c, z);

  // Transcribe real roots
  
  n_real		= 0;
  for (i=0; i<n_roots; ++i)
  {
    f			= creal(z[i]);
    f2			= (((c[4]*f+c[3])*f+c[2])*f+c[1])*f+c[0];
    if (fabs(cimag(z[i]))<1e-9)
      phi[n_real++]	= 2.0*atan(creal(z[i]));
  }
  if (n_real&1)					// only even number of roots
    return -1;

  // Sort roots

  for (i=0; i<n_real-1; ++i)
    for (j=i+1; j<n_real; ++j)
      if (phi[j]<phi[i])
      {
	f		= phi[i];
	phi[i]		= phi[j];
	phi[j]		= f;
      }
	
  // Create feasibility intervals
    
  if (c[4]<0.0)
  {
    phi_low[0]		= phi[0];
    phi_high[0]		= phi[1];
    phi_low[1]		= phi[2];
    phi_high[1]		= phi[3];
    return (long) (n_real/2);
  }
  else
  {
    phi_low[0]		= -M_PI;
    phi_high[0]		= phi[0];
    phi_low[1]		= phi[1];
    phi_high[1]		= phi[2];
    phi_low[2]		= phi[3];
    phi_high[(long) (n_real/2)]		
    			= M_PI;
    return (long) (n_real/2)+1;
  }
}    


vector			rM,
 			rP, rQ,
  			r2, r3, r4;
double			a2, a4, phiL, phiR,
  			l3, cosa3, l4;
long			branch, R_FFAIL;
 
double Fx(double phi)
{ // revamped (vector manipulations written out)
  static double		f1, f2, fx, fy, fz, d;

  // r3
  
  r3.x			= (fx=R*cos(phi))*w2.x+(fy=R*sin(phi))*w3.x+rN.x;
  r3.y			= fx*w2.y+fy*w3.y+rN.y;
  r3.z			= fx*w2.z+fy*w3.z+rN.z;
  
  // r2
  
  f1			= (fx=rP.x-r3.x)*fx+(fy=rP.y-r3.y)*fy+(fz=rP.z-r3.z)*fz
  			    -l3*l3+a2*a2;
  f2			= 2.0*a2*(fx*u2.x+fy*u2.y+fz*u2.z);
  fy			= 4.0*a2*(fx*u3.x+fy*u3.y+fz*u3.z);
  fx			= f1-f2;
  fz			= f1+f2;
  d			= fy*fy-4.0*fx*fz;
  if ((d<0.0)||((fx==0.0)&&(fy==0.0)))
  {
    R_FFAIL		= TRUE;
    return 1e5;
  }
  fx			= a2*cos(phiL=2.0*atan(fx ? 0.5*(-fy+(branch&2 ?
  			  sqrt(d) : -sqrt(d)))/fx : -fz/fy));
  fy			= a2*sin(phiL);
  r2.x			= fx*u2.x+fy*u3.x+rP.x;
  r2.y			= fx*u2.y+fy*u3.y+rP.y;
  r2.z			= fx*u2.z+fy*u3.z+rP.z;

  // r4
  
  f1			= (fx=rQ.x-r3.x)*fx+(fy=rQ.y-r3.y)*fy+(fz=rQ.z-r3.z)*fz
  			    -l4*l4+a4*a4;
  f2			= 2.0*a4*(fx*v2.x+fy*v2.y+fz*v2.z);
  fy			= 4.0*a4*(fx*v3.x+fy*v3.y+fz*v3.z);
  fx			= f1-f2;
  fz			= f1+f2;
  d			= fy*fy-4.0*fx*fz;
  if ((d<0.0)||((fx==0.0)&&(fy==0.0)))
  {
    R_FFAIL		= TRUE;
    return 1e5;
  }
  fx			= a4*cos(phiR=2.0*atan(fx ? 0.5*(-fy+(branch&1 ?
  			  sqrt(d) : -sqrt(d)))/fx : -fz/fy));
  fy			= a4*sin(phiR);
  r4.x			= fx*v2.x+fy*v3.x+rQ.x;
  r4.y			= fx*v2.y+fy*v3.y+rQ.y;
  r4.z			= fx*v2.z+fy*v3.z+rQ.z;
  
  // Epilogue

  return ((fx=r4.x-r2.x)*fx+(fy=r4.y-r2.y)*fy+(fz=r4.z-r2.z)*fz
    -l3*l3-l4*(l4+2.0*l3*cosa3))/l3/l4;
}


double			phi_low[6], phi_high[6];

long Feasible(vector *p, sphere *s)
{
  long			i, j, k, nL, nR, ni;
  double		phiL_low[3], phiL_high[3],
  			phiR_low[3], phiR_high[3],
			f51, f31, f53;
  vector		rL, rR;
  register double	f;
  
  // Initialize frames and variables
  
  rM.x			= (p[0].x+p[6].x)/2.0;
  rM.y			= (p[0].y+p[6].y)/2.0;
  rM.z			= (p[0].z+p[6].z)/2.0;
  rL.x			= (p[0].x+p[1].x+p[2].x)/3.0;
  rL.y			= (p[0].y+p[1].y+p[2].y)/3.0;
  rL.z			= (p[0].z+p[1].z+p[2].z)/3.0;
  rR.x			= (p[4].x+p[5].x+p[6].x)/3.0;
  rR.y			= (p[4].y+p[5].y+p[6].y)/3.0;
  rR.z			= (p[4].z+p[5].z+p[6].z)/3.0;
 
  if (Frame(p, p+1, &rL, &u1, &u2, &u3)) return 0;
  if (Frame(p+6, p+5, &rR, &v1, &v2, &v3)) return 0;
  if (Frame(p+1, p+5, &rM, &w1, &w2, &w3)) return 0;
  w2.x			= -w2.x;
  w2.y			= -w2.y;
  w2.z			= -w2.z;
  w3.x			= -w3.x;
  w3.y			= -w3.y;
  w3.z			= -w3.z;
  
  rP.x			= (f=s[2].d*cos(s[2].alpha))*u1.x+p[1].x;
  rP.y			= f*u1.y+p[1].y;
  rP.z			= f*u1.z+p[1].z;
  rQ.x			= (f=s[5].d*cos(s[6].alpha))*v1.x+p[5].x;
  rQ.y			= f*v1.y+p[5].y;
  rQ.z			= f*v1.z+p[5].z;
  f31			= s[2].d*s[2].d+s[3].d*
  			    (s[3].d+2.0*s[2].d*cos(s[3].alpha));
  f53			= s[4].d*s[4].d+s[5].d*
  			    (s[5].d+2.0*s[4].d*cos(s[5].alpha));
  if (!(f51=(f=p[1].x-p[5].x)*f+(f=p[1].y-p[5].y)*f+(f=p[1].z-p[5].z)*f))
    return 0;
  f			= (f53-f31)/f51;
  rN.x			= 0.5*((1+f)*p[1].x+(1-f)*p[5].x);
  rN.y			= 0.5*((1+f)*p[1].y+(1-f)*p[5].y);
  rN.z			= 0.5*((1+f)*p[1].z+(1-f)*p[5].z);
  f			= 0.5*(f51+f31-f53);
  if ((f = floor((f31-f*f/f51)*1e15+0.5)/1e15)<=0.0)
    return 0;
  R			= sqrt(f);
 
  a2			= s[2].d*sin(s[2].alpha);
  a4			= s[5].d*sin(s[6].alpha);
  l3			= s[3].d;
  l4			= s[4].d;
  cosa3			= cos(s[4].alpha);
  
  // Feasibility intervals

  if ((nL = Feasibility(&rP, l3, a2, &u1, phiL_low, phiL_high))<0)
    return 0;
  if ((nR = Feasibility(&rQ, l4, a4, &v1, phiR_low, phiR_high))<0)
    return 0;
  
  // Construct phi interval from separate phiL and phiR intervals

  ni			= 0;
  for (i=0; i<nL; ++i)
    for (j=0; j<nR; ++j)
      if ((phiL_low[i]>=phiR_low[j])&&(phiL_low[i]<=phiR_high[j]))
      {
        phi_low[ni]	= phiL_low[i];
	if (phiL_high[i]<phiR_high[j])
	  phi_high[ni++]= phiL_high[i];
	else
	  phi_high[ni++]= phiR_high[j];
      }
  for (i=0; i<nR; ++i)
    for (j=0; j<nL; ++j)
      if ((phiR_low[i]>=phiL_low[j])&&(phiR_low[i]<=phiL_high[j]))
      {
        for (k=0; (k<ni)&&(phiR_low[i]!=phi_low[k]); ++k);
	if (k>=ni)
	{
	  phi_low[ni]	= phiR_low[i];
	  if (phiR_high[i]<phiL_high[j])
	    phi_high[ni++]= phiR_high[i];
	  else
	    phi_high[ni++]= phiL_high[j];
        }
      }

  return ni;
}


long Rebridge(vector *p, sphere *s)
{
  long			i, j, n, ni, nr, nr_old;
  double		phi[B_MAXNROOTS];
  vector		t1, r[B_MAXNROOTS][3];
  
  if (!(ni = Feasible(p, s)))
    return 0;

  // Regression:
  //   Search all feasible intervals over all branches
  
  nr			= 0;
  n			= 0;
  R_FFAIL		= FALSE;
  for (branch=0; branch<4; ++branch)
    for (i=0; (i<ni)&&(nr<B_MAXNROOTS); ++i)
    {
      nr_old		= nr;
      nr		+= R_SolveNumerical(Fx, phi_low[i], phi_high[i], 
      				20, 4, phi+nr, B_MAXNROOTS-nr);
      if (R_FFAIL)
        return 0;

      if (nr>=B_MAXNROOTS)
        return 0;
     
      // Chain segment reconstruction
      
      for (j=nr_old; (j<nr); ++j)
      {
	//fprintf(stderr, "F%d[%g] = %g\n", j, phi[j], Fx(phi[j]));
	Fx(phi[j]);
	t1		= V_Subtr(&r3, p+3);
        if (V_Dot(&t1, &t1)>1e-15)
	{
	  r[n][0]	= r2;
	  r[n][1]	= r3;
	  r[n][2]	= r4;
	  ++n;
	}
      }
    }
  
  // Wrap-up
 
  if (nr>1)
    if (nr%2) 					// Geometric fail
      return 0;
 
  i			= (long) ((double) n*ran1(seed));
  p[2]			= r[i][0];
  p[3]			= r[i][1];
  p[4]			= r[i][2];
  return n;
}


double Jacobian(vector *p, sphere *s)	// J(III->I), see ref.1 Appendix
{
  long			i;
  static double		l[7],
  			b1xb2b3, b1xb2b4, b1xb2b5,
			b2xb3b4, b2xb5b6,
			b3xb4b5, b3xb5b6,
			b4xb5b6,
  			J;
  static vector		b[7], v;

  for (i=1; i<7; ++i)
  {
    l[i]		= s[i].d;
    b[i]		= V_Subtr(p+i, p+i-1); 
    b[i]		= V_Mult(1.0/l[i], b+i);
  }
  v			= V_Cross(b+1, b+2);
  b1xb2b3		= V_Dot(&v, b+3);
  b1xb2b4		= V_Dot(&v, b+4);
  b1xb2b5		= V_Dot(&v, b+5);
  v			= V_Cross(b+2, b+3);
  b2xb3b4		= V_Dot(&v, b+4);
  v			= V_Cross(b+2, b+5);
  b2xb5b6		= V_Dot(&v, b+6);
  v			= V_Cross(b+3, b+4);
  b3xb4b5		= V_Dot(&v, b+5);
  v			= V_Cross(b+3, b+5);
  b3xb5b6		= V_Dot(&v, b+6);
  v			= V_Cross(b+4, b+5);
  b4xb5b6		= V_Dot(&v, b+6);

  J			= 1.0/(l[1]*l[6]*fabs(
  			    -l[3]*l[5]*b1xb2b3*b3xb4b5*(
			      l[2]*b2xb5b6+l[3]*b3xb5b6+l[4]*b4xb5b6)+
			    l[2]*l[4]*b2xb3b4*b4xb5b6*(
			      l[3]*b1xb2b3+l[4]*b1xb2b4+l[5]*b1xb2b5)));
  for (i=2; i<7; ++i)
    J			*= sqrt(l[i-1]*l[i-1]+l[i]*
			     (l[i]+2.0*l[i-1]*cos(s[i].alpha)));
  
  return J;
}

double new_Jacobian(vector *p, sphere *s)	// J(IV->I), see ref.2 Appendix
{
  long			i;
  static double		l[7],
  			b1xb2b3, b1xb2b4, b1xb2b5,
			b2xb3b4, b2xb5b6,
			b3xb4b5, b3xb5b6,
			b4xb5b6,
  			J;
  static vector		b[7], v;

  for (i=1; i<7; ++i)
  {
    l[i]		= s[i].d;
    b[i]		= V_Subtr(p+i, p+i-1); 
    b[i]		= V_Mult(1.0/l[i], b+i);
  }
  v			= V_Cross(b+1, b+2);
  b1xb2b3		= V_Dot(&v, b+3);
  b1xb2b4		= V_Dot(&v, b+4);
  b1xb2b5		= V_Dot(&v, b+5);
  v			= V_Cross(b+2, b+3);
  b2xb3b4		= V_Dot(&v, b+4);
  v			= V_Cross(b+2, b+5);
  b2xb5b6		= V_Dot(&v, b+6);
  v			= V_Cross(b+3, b+4);
  b3xb4b5		= V_Dot(&v, b+5);
  v			= V_Cross(b+3, b+5);
  b3xb5b6		= V_Dot(&v, b+6);
  v			= V_Cross(b+4, b+5);
  b4xb5b6		= V_Dot(&v, b+6);

  J			= 1.0/(l[1]*l[6]*fabs(
  			    -l[3]*l[5]*b1xb2b3*b3xb4b5*(
			      l[2]*b2xb5b6+l[3]*b3xb5b6+l[4]*b4xb5b6)+
			    l[2]*l[4]*b2xb3b4*b4xb5b6*(
			      l[3]*b1xb2b3+l[4]*b1xb2b4+l[5]*b1xb2b5)));
  for (i=1; i<=5; i++)
    J	*=	l[i];
  for (i=2; i<=6; i++)
    J	*=	l[i];
  for (i=2; i<=6; i++)
    J	*=	sin(s[i].alpha);
  
  return J;
}
// Setup input parameters for 7 vectors p[] and spherical coordinates s[]
// used by Feasible(), Rebridge(), and Jacobian() through the molecule
// parent list.  Used in conjunction with NeighborList(), NextEndBridge(),
// and NextRebridge().  
//
// Input:
//   molm:	Target molecule.
//   end:	Where to start in the parent list.
//   reverse:	Flag for transcription direction; 0: forwards, 1: reverse.
//
// Output:
//   p:		Transcribed positions.
//   s:		Calculated bond lengths and angles (no torsions!).
/*
long RebridgeSetup(
	molstruct *molm, long end, long reverse, vector *p, sphere *s)
{
   long			i, ip, flag = 0;
   sphere		*s1;
   vector		dr, drp;

   if (end >= molm->nsites) return 1;		// Exit on failure

   if (reverse)	{	p+=0;	s+=1;}		// end bead as p[0]
   else		{ 	p+=6;	s+=6;}		// end bead as p[6]

   ip	=	end;
   *p	=	molm->p[ip];

   for (i=0; (i<6)&&((ip = molm->parent[ip])>=0); ++i) {
      if (reverse)	p++;
      else		p--;

      dr	=	V_Subtr(p, molm->p+ip);		// bond vector
      *p	=	molm->p[ip];			// update p

      s->d	=	sqrt(V_Dot(&dr, &dr));		// bond length
      dr	=	V_Mult(1.0/s->d, &dr);		// unit bond vector

      if (i) {
         s1->alpha	=	acos(V_Dot(&dr, &drp));	// bond angle
      }
      drp	=	dr;			// previous dr

      if (reverse) {	s++;	s1=s;}		// or s1 = reverse ? ++s : s--;
      else	   {	s1=s;	s--;}
   }
   return (ip<0);	// return false if NO error occur
}
*/

long RebridgeSetup(
	molstruct *molm, long end, long reverse, vector *p, sphere *s)
{
  long			i, ip = end, flag = 0;
  sphere		*sl;
  vector		dr, drp, *pis, *pip;

  if (end>=molm->nsites) return 1;		// Exit on failure
  pis			= molm->p+ip;
  p			+= reverse ? 0 : 6;
  s			+= reverse ? 1 : 6;
  *p			= *pis;			// Transcribe positions
  for (i=0; (i<6)&&((ip = molm->parent[ip])>=0); ++i)
  {
    if (reverse) 
      ++p;
    else 
      --p;
    pip			= molm->p+ip;		// Parent position pointer
    *p			= *pip;			
    drp			= V_Subtr(pis, pip);	// Determine direction
    drp			= V_Mult(1.0/(s->d = sqrt(V_Dot(&drp, &drp))), &drp);
    if (i)
    {
      sl->alpha		= acos(V_Dot(&drp, &dr));
    }
    pis			= pip;
    dr			= drp;
    sl	 		= reverse ? ++s : s--;
  }
  return (ip<0);
}
