/*
    program:    sample.c
    author:     Peng Yi at MIT
    date:       October 22, 2006
		October 13, 2007, adding Pieter's functions
    purpose:    Sample system properties
    notes:	April 18,2011	Add structure factor calculation based on definition
*/

#define __SAMPLE_MODULE
#include "sample.h"


long mod (long numb, long divisor)
{
  while (numb<0) {
    numb+=divisor;
  }
  return numb%divisor;
}


long factorial(long number)			//calculate factorial n! 
{
   long	result;

   switch(number){
	case 12:		result	=	479001600;		break;
	case 11:		result	=	39916800;		break;
	case 10:		result	=	3628800;		break;
	case 9:		result	=	362880;		break;
	case 8:		result	=	40320;		break;
	case 7:		result	=	5040;		break;
	case 6:		result	=	720;		break;
	case 5:		result	=	120;		break;
	case 4:		result	=	24;		break;
	case 3:		result	=	6;		break;
	case 2:		result	=	2;		break;
	case 1: case 0:		result	=	1;		break;
   }
   return	result;
/*
  if (number<=1)
    return 1;
  else
    return (number*factorial(number-1));
*/
}

long intpow(long base, long power)		//calculate integer power, base^power
{
   long		i, result;
   result	=	1;
   for (i=0; i<power; i++) {
      result	*= base;
   }
   return 	result;
}

long intlog(long base, long value)		//calculate log_base^value
{
   long		power, result;
   power	=	0;
   result	=	1;
   while (result<value) {
      result	*=	base;
      power	++;
   }
   return	power;
}

double plgndr(int l, int m, double x)  		// associated Legendre polynomial
{                                      		// m=0, 1, ..., l, x = [-1,1]
  double fact, pll, pmm, pmmp1, somx2;
  int i, ll;

  if (m<-l || m>l || fabs(x)>1.0) {
    printf("Bad arguments in routine plgndr.\n");
    return 0;
  }

//  if (m<0)
//    return ((mod(-m,2)==0)?1.0:-1.0) * factorial(l+m) / factorial(l-m) * plgndr(l, -m, x);

  pmm = 1.0;                           		//compute P_m^m
  if (m>0) {
    somx2 = sqrt((1.0-x)*(1.0+x));
    fact = 1.0;
    for (i=1; i<=m; i++) {
      pmm *= -fact*somx2;
      fact += 2.0;
    }
  }
  if (l==m)
    return pmm;
  else {                              		//compute P_{m+1}^m
    pmmp1 = x*(2*m+1)*pmm;
    if (l==(m+1))
      return pmmp1;
    else {                            		//compute P_l^m, l>m+1
      for (ll=m+2; ll<=l; ll++) {
	pll = (x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
	pmm = pmmp1;
	pmmp1 = pll;
      }
      return pll;
    }
  }
}
/*
complex sfharmonics2(int l, int m, double costheta, double cosphi, double sinphi)	//spherical harmonics functions
{								//m = -l, -l+1, ..., l-1, l
  double 	Nlm, Plm;
  complex	Ylm;
  complex	eiphi, eimphi;
  int		i;

  if (abs(m)>l || fabs(costheta)>1.0) {
    printf("Bad arguments in routine sfharmonics.\n");
    return 0;
  }

  eiphi		=	cosphi + I*sinphi;
  eimphi	=	1.0;
  for (i=0; i<abs(m); i++) {
     eimphi	*=	(cosphi + I*sinphi);
  }

  Plm 	= 	plgndr(l, abs(m), costheta);
  Nlm	=	sqrt( (double) (2*l+1)/4/pi * factorial(l-abs(m)) / factorial(l+abs(m)));
//  Ylm	=	Nlm * Plm * cexp(I*abs(m)*phi);
  Ylm	=	Nlm * Plm * eimphi;

  if (m >= 0)
    return Ylm;
  else
//    return pow(-1,abs(m)) * conj(Ylm);
    return (mod(abs(m),2) ? -1 : 1) * conj(Ylm);
}
*/

#ifdef SAMPLE_CODE

double complex sfharmonics2(int l, int m, double costheta, double cosphi, double sinphi)
{
   double		SIN, SIN2, COS, COS2;
   double complex	eiphi, eimphi;		//exp(I*phi), exp(I*m*phi)
   int			i;
   double complex	Ylm;

   COS	=	costheta;
   COS2	=	COS * COS;
   SIN2	=	1-COS2;
   SIN	=	sqrt(SIN2);
   
   eiphi	=	cosphi + I*sinphi;
   eimphi	=	1.0;
   for (i=0; i<abs(m); i++) {
      eimphi	*=	eiphi;
   }

   switch (l) {
      case 6: {
         switch (abs(m)) {
            case 6:	Ylm	=	1.0/64 * sqrt(3003.0/pi) * eimphi * SIN2 * SIN2 * SIN2;
			break;
            case 5:	Ylm	=	-3.0/32 * sqrt(1001.0/pi) * eimphi * SIN2 * SIN2 * SIN * COS;
			break;
	    case 4:	Ylm	=	3.0/32 * sqrt(91.0/2/pi) * eimphi * SIN2 * SIN2 * (11*COS2-1);
			break;
            case 3:	Ylm	=	-1.0/32 * sqrt(1365.0/pi) * eimphi * SIN2 * SIN * (11*COS2-3) * COS;
			break;
	    case 2:	Ylm	=	1.0/64 * sqrt(1365.0/pi) * eimphi * SIN2 * (33*COS2*COS2 -18*COS2 +1);
			break;
	    case 1:	Ylm	=	-1.0/16 * sqrt(273.0/2/pi) * eimphi * SIN * (33*COS2*COS2-30*COS2+5) *COS;
			break;
	    case 0:	Ylm	=	1.0/32 * sqrt(13.0/pi) * (231*COS2*COS2*COS2-315*COS2*COS2+105*COS2-5);
			break;
	    default:	printf("l, m mismatch!\n");
			break;
	 }
         break;
      }
      case 4: {
         switch (abs(m)) {
            case 4:	Ylm	=	3.0/16 * sqrt(35.0/2/pi) * eimphi * SIN2 * SIN2;
			break;
	    case 3:	Ylm	=	-3.0/8 * sqrt(35.0/pi) * eimphi * SIN2 * SIN * COS;
			break;
	    case 2:	Ylm	=	3.0/8 * sqrt(5.0/2/pi) * eimphi * SIN2 * (7*COS2-1);
			break;
	    case 1:	Ylm	=	-3.0/8 * sqrt(5.0/pi) * eimphi * SIN * (7*COS2-3) * COS;
			break;
	    case 0:	Ylm	=	3.0/16 * sqrt(1.0/pi) * (35*COS2*COS2 - 30*COS2 +3);
			break;	
	    default:	printf("l, m mismatch!\n");
			break;
	 }
	 break;
      }
      default:	printf("could not find the l!\n");
		break;
   }

   if (m>=0)
      return	Ylm;
   else
      return	pow(-1, abs(m)) * conj(Ylm);
}

/*
complex sfharmonics(int l, int m, double costheta, double phi)	//spherical harmonics functions
{								//m = -l, -l+1, ..., l-1, l
  double 	Nlm, Plm;
  complex	Ylm;

  if (abs(m)>l || fabs(costheta)>1.0) {
    printf("Bad arguments in routine sfharmonics.\n");
    return 0;
  }

  Plm 	= 	plgndr(l, abs(m), costheta);
  Nlm	=	sqrt( (double) (2*l+1)/4/pi * factorial(l-abs(m)) / factorial(l+abs(m)));
  Ylm	=	Nlm * Plm * cexp(I*abs(m)*phi);

  if (m >= 0)
    return Ylm;
  else
//    return pow(-1,abs(m)) * conj(Ylm);
    return (mod(abs(m),2) ? -1 : 1) * conj(Ylm);
}
*/

double complex sfharmonics(int l, int m, double costheta, double phi)	//table from Wikipedia spherical harmonics entry
{
   double		SIN, SIN2, COS, COS2;
   double complex	Ylm;
   COS	=	costheta;
   COS2	=	COS * COS;
   SIN2	=	1-COS2;
   SIN	=	sqrt(SIN2);
   
   switch (l) {
      case 6: {
         switch (abs(m)) {
            case 6:	Ylm	=	1.0/64 * sqrt(3003.0/pi) * cexp(6*I*phi) * SIN2 * SIN2 * SIN2;
			break;
            case 5:	Ylm	=	-3.0/32 * sqrt(1001.0/pi) * cexp(5*I*phi) * SIN2 * SIN2 * SIN * COS;
			break;
	    case 4:	Ylm	=	3.0/32 * sqrt(91.0/2/pi) * cexp(4*I*phi) * SIN2 * SIN2 * (11*COS2-1);
			break;
            case 3:	Ylm	=	-1.0/32 * sqrt(1365.0/pi) * cexp(3*I*phi) * SIN2 * SIN * (11*COS2-3) * COS;
			break;
	    case 2:	Ylm	=	1.0/64 * sqrt(1365.0/pi) * cexp(2*I*phi) * SIN2 * (33*COS2*COS2 -18*COS2 +1);
			break;
	    case 1:	Ylm	=	-1.0/16 * sqrt(273.0/2/pi) * cexp(I*phi) * SIN * (33*COS2*COS2-30*COS2+5) *COS;
			break;
	    case 0:	Ylm	=	1.0/32 * sqrt(13.0/pi) * (231*COS2*COS2*COS2-315*COS2*COS2+105*COS2-5);
			break;
	    default:	printf("l, m mismatch!\n");
			break;
	 }
         break;
      }
      case 4: {
         switch (abs(m)) {
            case 4:	Ylm	=	3.0/16 * sqrt(35.0/2/pi) * cexp(4*I*phi) * SIN2 * SIN2;
			break;
	    case 3:	Ylm	=	-3.0/8 * sqrt(35.0/pi) * cexp(3*I*phi) * SIN2 * SIN * COS;
			break;
	    case 2:	Ylm	=	3.0/8 * sqrt(5.0/2/pi) * cexp(2*I*phi) * SIN2 * (7*COS2-1);
			break;
	    case 1:	Ylm	=	-3.0/8 * sqrt(5.0/pi) * cexp(I*phi) * SIN * (7*COS2-3) * COS;
			break;
	    case 0:	Ylm	=	3.0/16 * sqrt(1.0/pi) * (35*COS2*COS2 - 30*COS2 +3);
			break;	
	    default:	printf("l, m mismatch!\n");
			break;
	 }
	 break;
      }
      default:	printf("could not find the l!\n");
		break;
   }

   if (m>=0)
      return	Ylm;
   else
      return	pow(-1, abs(m)) * conj(Ylm);
}
#endif //SAMPLE_CODE
/*************************************************************************************/

#ifdef TEST
complex aveqlm (int l, int m, long i)	//calculate average q_lm value for particle i
{					//also update part[i].nbond & part[i].qlm !!
  long		jj, k, nbond;
  vector	dp;
  double	r2, alpha, alphasum;
  double	costheta, phi;    	//theta=[0,PI], phi=[0,2*PI]
  complex	aveqlm;

  nbond		=	0;  			//neighbor bond number
  aveqlm	=	0;
  alphasum	=	0;
 
  for (jj=0; jj<part[i].nverlet; jj++) {

    k		=	part[i].vlist[jj];
    dp.x	=	part[k].p.x	-	part[i].p.x;
    dp.y	=	part[k].p.y	-	part[i].p.y;
    dp.z	=	part[k].p.z	-	part[i].p.z;
    MapInBox(&dp);
    r2		=	dp.x*dp.x + dp.y*dp.y + dp.z*dp.z;
    if (r2 <= Rb*Rb) {
	nbond 	++;

	costheta =  dp.z/sqrt(r2);     	//theta=[0,PI]
//	if (dp.x >= 0)                 	//atan=[-PI/2, PI/2]
	  phi	=	atan2( dp.y, dp.x ); 
//	else
//	  phi	=	atan2( dp.y, dp.x ) + pi;

	alpha	=	(r2 + Rb*Rb - 2*sqrt(r2)*Rb) / 
			(sigma*sigma + Rb*Rb - 2*sigma*Rb);
	aveqlm 		+=	sfharmonics(l, m, costheta, phi) * alpha;
	alphasum	+=	alpha;
	//alpha is a damping factor, as explained in JCP v104 (1996) 9932
    }
  }
  if ( fabs(alphasum) < ZERO) {		//always be awared of divided-by-zero
    aveqlm		=	0;
  }
  else {
    aveqlm		=	aveqlm/alphasum;
  }
  part[i].nbond		=	nbond;		//update part[i].nbond
  part[i].qlm[m+l]	=	aveqlm;		//update part[i].qlm
  return aveqlm;
}

/*
complex tildeqlm(int l, int m, long i)	//calculate normalized q_lm for particle i
{
  int n;
  double sum;
  
  for (n=-l; n<=l; n++) {
    sum +=	pow(cabs(aveqlm(l,n,i)), 2); 
  }
  sum = sqrt(sum);
  
  if ( fabs(sum) < ZERO)
    return 0;
  else
    return aveqlm(l, m, i)/sum;
}
*/

/*
double ql(int l, long i)		//calculate q_l of particle i
{
  int 		m;
  double 	ql=0;
  complex 	term; 

  for (m=-l; m<=l; m++) {
    term	=	aveqlm(l, m, i);
    ql		+=	term * conj(term);
  }
  ql		*=	(double)4*pi/(2*l+1);
  return sqrt(ql);
}
*/

/******************************************************************************/

double ql(int l, long i)		//calculate q_l of particle i
{
  int		m;
  double	ql=0;

  for (m=-l; m<=l; m++) {
    ql	+=	part[i].qlm[m+l] * conj(part[i].qlm[m+l]);
  }
  ql	*=	(double)4*pi/(2*l+1);
  return	sqrt(ql);
}

/*
complex qlproduct(int l, long i, long j)	//calculate q_l(i) dot q_l(j)
{
  int		m;
  complex	termi, termj, term=0;
  double	sumi=0, sumj=0;

  for (m=-l; m<=l; m++) {
    termi	=	aveqlm(l, m, i);
    termj	=	aveqlm(l, m, j);
    term	+=	termi * conj(termj);
    sumi	+=	termi * conj(termi);
    sumj	+=	termj * conj(termj);
  }
  return	term/sqrt(sumi*sumj);		//need to check divided-by-zero
}
*/

/*********************************************************************************/

/**********************************************************************************/
/*
complex aveQlm1(int l, int m) 		//calculate average Q_lm
{
  long		i;
  long		sumbond	=	0;
  complex 	aveQlm	=	0;
  complex	term;

  for (i=0; i<NPARTS; i++)  {
    term	=	aveqlm(l, m, i);	//notice because aveqlm(l,m,i) updates part[i].nbond,
    aveQlm 	+=	part[i].nbond * term;	//then cannot appear in one expression at same time,
    sumbond	+=	part[i].nbond;		//e.g., aveQlm = part[i].bond * aveqlm(l,m,i).
  }
  
  if (sumbond == 0)
    return	0;
  else
    return aveQlm/sumbond;
}
*/

/****************************************************************************************/
/*
double CalcQl1(int l)    			//calculate Q_l
{
  int		m;
  double	Ql=0;
  complex	term;
 
  for (m=-l; m<=l; m++) {
    term	=	aveQlm1(l, m);
    Ql		+=	term * conj(term);    
  }
  Ql *= 4*pi/(2*l+1);
  return sqrt(Ql);
}
*/

/*************************************************************************************/

complex aveQlm(int l, int m) 		//calculate average Q_lm
{
  long		i;
  long		sumbond	=	0;
  complex 	aveQlm	=	0;

  for (i=0; i<NPARTS; i++)  {
    aveQlm 	+=	part[i].nbond * part[i].qlm[m+l];
    sumbond	+=	part[i].nbond;		
  }
  
  if (sumbond == 0)
    return	0;
  else
    return aveQlm/sumbond;
}

/******************************************************************************/

double CalcQl(int l)
{
  long		i;
  int		m;
  long		sumnbond;
  complex	aveQlm;
  double	Ql	=	0;
  
  for (m=-l; m<=l; m++) {
 
    sumnbond	=	0;	//notice where to set Ql=0 and where 
    aveQlm	=	0;	//to set aveQlm and sumnbond =0 !!

    for (i=0; i<NPARTS; i++) {		//calculate aveQlm
      aveQlm	+=	part[i].nbond * part[i].qlm[m+l];
      sumnbond	+=	part[i].nbond;
    }
    if (sumnbond	==	0)
      aveQlm	=	0;
    else
      aveQlm	=	aveQlm/sumnbond;

    Ql	+=	aveQlm * conj(aveQlm);
  }
  Ql	*=	(double) 4*pi/(2*l+1);
  return sqrt(Ql);
}

/****************************************************************************************/
/* void New_Qlm(int l)									*/
/*  											*/
/* Calculate aveqlm and ql of every particle, and also calculate			*/
/* aveQlm and Ql for the system.							*/
/* And since aveqlm = ylmalphasum/alphasum, we save two variables 			*/
/* ylmalphasum and alphasum, as well as aveqlm for each particle. 			*/
/* Same for Q.										*/
/*											*/
/* Input: configuration, l = l_of_Ylm, Verlet list.					*/
/****************************************************************************************/

void New_Qlmtemp(int l)
{
  int		m;
  long		i, jj, j;
  vector	dp;
  double	r2, alpha;	//alpha is a weighted # of bonds, see JCP V104 (1996) 9932
  double	costheta, phi;

   for (i=0; i<NPARTS; i++) {
      parttemp[i]	=	part[i];
   }

  AlphaSumtemp	=	0;			//initialize AlphaSum and YlmAlphaSum for
  for (m=-l; m<=l; m++) {			//the calculation of Qlm and Ql
    YlmAlphaSumtemp[m+l]	=	0;	
  }
  for (i=0; i<NPARTS; i++) {			//initialize alphasum and ylmalphasum for
    parttemp[i].alphasum	=	0;		//the calculation of qlm and ql
    for (m=-l; m<=l; m++) {			//for each particle
      parttemp[i].ylmalphasum[m+l]	=	0;
    }
  }

  for (i=0; i<NPARTS; i++) {			//Qlm is averaged among All particles
    for (jj=0; jj<part[i].nverlet; jj++) {	//search the Verlet list
			
      j		=	part[i].vlist[jj];
      if (j > i) {				//reduce theta and phi evaluation

        dp.x	=	part[j].p.x	-	part[i].p.x;
        dp.y	=	part[j].p.y	-	part[i].p.y;
        dp.z	=	part[j].p.z	-	part[i].p.z;
        MapInBox(&dp);
        r2		=	dp.x*dp.x + dp.y*dp.y + dp.z*dp.z;
    
        if (r2 < Rb*Rb) {
	  alpha	=	(r2 + Rb*Rb - 2*sqrt(r2)*Rb) / (sigma*sigma + Rb*Rb - 2*sigma*Rb);
	  parttemp[i].alphasum	+=	alpha;
	  parttemp[j].alphasum	+=	alpha;

          costheta	= dp.z/sqrt(r2);     	//theta=[0,PI]
	  phi		= atan2( dp.y, dp.x ); 
	  for (m=-l; m<=l; m++) {
	    parttemp[i].ylmalphasum[m+l]	+=	sfharmonics(l, m, costheta, phi) * alpha;
	    parttemp[j].ylmalphasum[m+l]	+=	sfharmonics(l, m, -costheta, pi+phi) * alpha;
	  }
        }
      }
    } 
  }

  for (i=0; i<NPARTS; i++) {			//calculate qlm and ql for each particle
    parttemp[i].ql		=	0; 
    if ( parttemp[i].alphasum < 1e-12) {	//always be aware of divided-by-zero
      if (parttemp[i].alphasum < 0)		//alphasum should always be positive
         printf("New_Qlm error, alphasum_i = %f\n", parttemp[i].alphasum);
      for (m=-l; m<=l; m++) {
        if (cabs(parttemp[i].ylmalphasum[m+l]) > ZERO) 
	  printf("Error, ylmalphasum[%ld] = %f+%f*i\n", m, creal(parttemp[i].ylmalphasum[m+l]), 
			cimag(parttemp[i].ylmalphasum[m+l]));
        parttemp[i].qlm[m+l] 	=	0;	//update qlm for each particle and each m
      }
      parttemp[i].ql		=	0;
    }
    else {
      for (m=-l; m<=l; m++) {
        parttemp[i].qlm[m+l] 	=	parttemp[i].ylmalphasum[m+l] / parttemp[i].alphasum;
	parttemp[i].ql		+=	parttemp[i].qlm[m+l] * conj(parttemp[i].qlm[m+l]);
      }
    }
    parttemp[i].ql	=	sqrt(parttemp[i].ql*4*pi/(2*l+1));
    
    AlphaSumtemp	+=	parttemp[i].alphasum;	//sum over all particles to calculate Qlm
    for (m=-l; m<=l; m++) {
      YlmAlphaSumtemp[m+l]	+=	parttemp[i].ylmalphasum[m+l];
    } 
  }
  
  Qltemp		=	0;			//global variable
  if (AlphaSumtemp < 1e-12) {			//calculate aveQlm and Ql
      if (AlphaSumtemp < 0)
         printf("New_Qlm error, AlphaSum = %f\n", AlphaSumtemp);
    for (m=-l; m<=l; m++) {
	Qlmtemp[m+l]	=	0;
    }
    Qltemp	=	0;
  }
  else {
    for (m=-l; m<=l; m++) {
	Qlmtemp[m+l]	=	YlmAlphaSumtemp[m+l]/AlphaSumtemp;
	Qltemp		+=	Qlmtemp[m+l] * conj(Qlmtemp[m+l]);
    }
    Qltemp	=	sqrt(Qltemp*4*pi/(2*l+1));
  }

  return;	//variables updated are global, so no need to have a return value
}

/****************************************************************************************/

void New_Qlm(int l)		//calculate q and Q using verlet list	
{
  int		m;
  long		i, jj, j, k, icell, jcell;
  vector	dp;
  double	r2, rxy, alpha;	//alpha is a weighted # of bonds, see JCP V104 (1996) 9932
  double	costheta, phi, cosphi, sinphi;

  AlphaSum	=	0;			//initialize AlphaSum and YlmAlphaSum for
  for (m=-l; m<=l; m++) {			//the calculation of Qlm and Ql
    YlmAlphaSum[m+l]	=	0;	
  }
  for (i=0; i<NPARTS; i++) {			//initialize alphasum and ylmalphasum for
    part[i].alphasum	=	0;		//the calculation of qlm and ql
    for (m=-l; m<=l; m++) {			//for each particle
      part[i].ylmalphasum[m+l]	=	0;
    }
  }

#ifdef VERLET_LIST
  for (i=0; i<NPARTS-1; i++) {			//Qlm is averaged among All particles
    for (jj=0; jj<part[i].nverlet; jj++) {	//search the Verlet list
			
      j		=	part[i].vlist[jj];

      if (j > i) {				//reduce theta and phi evaluation
      {
#elif CELL_LIST
  for (i=0; i<NPARTS-1; i++) {
     icell	=	part[i].icell;

     for (jj=0; jj<Cell[icell].nneigh; jj++) {
        jcell	=	Cell[icell].neigh[jj];

        for (k=0; k<Cell[jcell].sites; k++) {
           j	=	Cell[jcell].list[k]; 

           if (j > i) {
#else
  for (i=0; i<NPARTS-1; i++) {
     for (j=i+1; j<NPARTS; j++) {
        {{					//match the # of loops

#endif
        dp.x	=	part[j].p.x	-	part[i].p.x;
        dp.y	=	part[j].p.y	-	part[i].p.y;
        dp.z	=	part[j].p.z	-	part[i].p.z;
        MapInBox(&dp);
        r2		=	dp.x*dp.x + dp.y*dp.y + dp.z*dp.z;
    
        if (r2 < Rb*Rb) {
	  alpha	= 1.0;
	  //alpha	=	(r2 + Rb*Rb - 2*sqrt(r2)*Rb) / (sigma*sigma + Rb*Rb - 2*sigma*Rb);
	  part[i].alphasum	+=	alpha;
	  part[j].alphasum	+=	alpha;

          costheta	= dp.z/sqrt(r2);     	//theta=[0,PI]
//	  phi		= atan2( dp.y, dp.x ); 
	  rxy		=	sqrt(dp.x * dp.x + dp.y * dp.y);
	  cosphi	=	dp.x / rxy;
	  sinphi	=	dp.y / rxy;
	  if (rxy<ZERO) {
	     cosphi	=	0;
	     sinphi	=	0;
          }

	  for (m=-l; m<=l; m++) {
	    part[i].ylmalphasum[m+l]	+=	sfharmonics2(l, m, costheta, cosphi, sinphi) * alpha;
	    part[j].ylmalphasum[m+l]	+=	sfharmonics2(l, m, -costheta, -cosphi, -sinphi) * alpha;
	    
//	    part[i].ylmalphasum[m+l]	+=	sfharmonics(l, m, costheta, phi) * alpha;
//	    part[j].ylmalphasum[m+l]	+=	sfharmonics(l, m, -costheta, pi+phi) * alpha;
	  }
        }
      }
    } 
  }}

  for (i=0; i<NPARTS; i++) {			//calculate qlm and ql for each particle
    part[i].ql		=	0; 
    if ( part[i].alphasum < 1e-12) {	//always be aware of divided-by-zero
      if (part[i].alphasum < 0)		//alphasum should always be positive
         printf("New_Qlm error, alphasum_i = %f\n", part[i].alphasum);
      for (m=-l; m<=l; m++) {
        if (cabs(part[i].ylmalphasum[m+l]) > ZERO) 
	  printf("Error, ylmalphasum[%ld] = %f+%f*i\n", m, creal(part[i].ylmalphasum[m+l]), 
			cimag(part[i].ylmalphasum[m+l]));
        part[i].qlm[m+l] 	=	0;	//update qlm for each particle and each m
      }
      part[i].ql		=	0;
    }
    else {
      for (m=-l; m<=l; m++) {
        part[i].qlm[m+l] 	=	part[i].ylmalphasum[m+l] / part[i].alphasum;
	part[i].ql		+=	part[i].qlm[m+l] * conj(part[i].qlm[m+l]);
      }
    }
    part[i].ql	=	sqrt(part[i].ql*4*pi/(2*l+1));
    
    AlphaSum	+=	part[i].alphasum;	//sum over all particles to calculate Qlm
    for (m=-l; m<=l; m++) {
      YlmAlphaSum[m+l]	+=	part[i].ylmalphasum[m+l];
    } 
  }
  
  Ql		=	0;			//global variable
  if (AlphaSum < 1e-12) {			//calculate aveQlm and Ql
      if (AlphaSum < 0)
         printf("New_Qlm error, AlphaSum = %f\n", AlphaSum);
    for (m=-l; m<=l; m++) {
	Qlm[m+l]	=	0;
    }
    Ql	=	0;
  }
  else {
    for (m=-l; m<=l; m++) {
	Qlm[m+l]	=	YlmAlphaSum[m+l]/AlphaSum;
	Ql		+=	Qlm[m+l] * conj(Qlm[m+l]);
    }
    Ql	=	sqrt(Ql*4*pi/(2*l+1));
  }

  return;	//variables updated are global, so no need to have a return value
}

/****************************************************************************************/

void New_Qlm_NoVlist(int l)	//calculate New Qlm without using Verlet list	
{
  int		m;
  long		i, jj, j;
  vector	dp;
  double	r2, alpha;	//alpha is a weighted # of bonds, see JCP V104 (1996) 9932
  double	costheta, phi;

  AlphaSum	=	0;			//initialize AlphaSum and YlmAlphaSum for
  for (m=-l; m<=l; m++) {			//the calculation of Qlm and Ql
    YlmAlphaSum[m+l]	=	0;	
  }
  Ql		=	0;
  for (i=0; i<NPARTS; i++) {			//initialize alphasum and ylmalphasum for
    part[i].alphasum	=	0;		//the calculation of qlm and ql
    for (m=-l; m<=l; m++) {			//for each particle
      part[i].ylmalphasum[m+l]	=	0;
    }
    part[i].ql		=	0; 
  }

  for (i=0; i<NPARTS-1; i++) {			//Qlm is averaged among All particles
     for (j=i+1; j<NPARTS; j++) {		

        dp.x	=	part[j].p.x	-	part[i].p.x;
        dp.y	=	part[j].p.y	-	part[i].p.y;
        dp.z	=	part[j].p.z	-	part[i].p.z;
        MapInBox(&dp);
        r2		=	dp.x*dp.x + dp.y*dp.y + dp.z*dp.z;
    
        if (r2 < Rb*Rb) {
	  alpha	=	(r2 + Rb*Rb - 2*sqrt(r2)*Rb) / (sigma*sigma + Rb*Rb - 2*sigma*Rb);
	  part[i].alphasum	+=	alpha;
	  part[j].alphasum	+=	alpha;

          costheta	= dp.z/sqrt(r2);     	//theta=[0,PI]
	  phi		= atan2( dp.y, dp.x ); 
	  for (m=-l; m<=l; m++) {
	    part[i].ylmalphasum[m+l]	+=	sfharmonics(l, m, costheta, phi) * alpha;
	    part[j].ylmalphasum[m+l]	+=	sfharmonics(l, m, -costheta, pi+phi) * alpha;
	  }
        }
     } 
  }

  for (i=0; i<NPARTS; i++) {			//calculate qlm and ql for each particle
    if ( part[i].alphasum < 1e-12) {	//always be aware of divided-by-zero
      if (part[i].alphasum < 0)		//alphasum should always be positive
         printf("New_Qlm error, alphasum_i = %f\n", part[i].alphasum);
      for (m=-l; m<=l; m++) {
        if (cabs(part[i].ylmalphasum[m+l]) > ZERO) 
	  printf("Error, ylmalphasum[%ld] = %f+%f*i\n", m, creal(part[i].ylmalphasum[m+l]), 
			cimag(part[i].ylmalphasum[m+l]));
        part[i].qlm[m+l] 	=	0;	//update qlm for each particle and each m
      }
      part[i].ql		=	0;
    }
    else {
      for (m=-l; m<=l; m++) {
        part[i].qlm[m+l] 	=	part[i].ylmalphasum[m+l] / part[i].alphasum;
	part[i].ql		+=	part[i].qlm[m+l] * conj(part[i].qlm[m+l]);
      }
    }
    part[i].ql	=	sqrt(part[i].ql*4*pi/(2*l+1));
    
    AlphaSum	+=	part[i].alphasum;	//sum over all particles to calculate Qlm
    for (m=-l; m<=l; m++) {
      YlmAlphaSum[m+l]	+=	part[i].ylmalphasum[m+l];
    } 
  }
  
  if (AlphaSum < 1e-12) {			//calculate aveQlm and Ql
      if (AlphaSum < 0)
         printf("New_Qlm error, AlphaSum = %f\n", AlphaSum);
    for (m=-l; m<=l; m++) {
	Qlm[m+l]	=	0;
    }
    Ql	=	0;
  }
  else {
    for (m=-l; m<=l; m++) {
	Qlm[m+l]	=	YlmAlphaSum[m+l]/AlphaSum;
	Ql		+=	Qlm[m+l] * conj(Qlm[m+l]);
    }
    Ql	=	sqrt(Ql*4*pi/(2*l+1));
  }

  return;	//variables updated are global, so no need to have a return value
}
/****************************************************************************************/

/*
void local_q_update()	//q update of all particles affected by a move
{
  long	jj;
  int	m;
  
  for (jj=0; jj<nverletplus; jj++) {	
    for (m=-l_of_Ylm; m<=l_of_Ylm; m++) {
      part[vlistplus[jj]].qlm[m+l_of_Ylm]	=	aveqlm(l_of_Ylm, m, vlistplus[jj]);
    }
  }
  return;
}
*/

/****************************************************************************************/
/* void Update_Qlm(int l, long i, vector p_old, vector p_new)				*/
/*											*/
/* Update Qlm and Ql due to displacement of particle i.  Update qlm also.		*/
/* The update of qlm is through update of ylmalphasum and alphasum!			*/
/* And the same is for Qlm.								*/
/*											*/
/* Input: configuration, verlet listplus						*/
/*        correct Q and q from last step.						*/
/* Caution: when use this function to recover Q and q, we need to be careful of		*/
/*	    the role of oldmol and part[n].						*/
/****************************************************************************************/

void Update_Qlm(int l, long i, vector p_old, vector p_new)
{			
  int		m;
  long		jj, j, k, jcell;
  vector	dp;
  double	r2, rxy, alpha, costheta, phi, cosphi, sinphi;
  complex	termi, termj;

#ifdef VERLET_LIST
  for (jj=0; jj<nverletplus-1; jj++) {	//the last element in vlistplus is itself

    j	=	vlistplus[jj];
    {{
#elif CELL_LIST
  for (jj=0; jj<nneighcellplus; jj++){
     jcell	=	neighcellplus[jj];

     for (k=0; k<Cell[jcell].sites; k++) {
        j	=	Cell[jcell].list[k];

        if (j!=i) {
#else
    {{{
#endif

    dp.x	=	part[j].p.x	-	p_old.x;	//old i, j distance
    dp.y	=	part[j].p.y	-	p_old.y;
    dp.z	=	part[j].p.z	-	p_old.z;
    MapInBox(&dp);
    r2		=	dp.x*dp.x + dp.y*dp.y + dp.z*dp.z;

    if (r2 < Rb*Rb) {
      //alpha	=	1.0;
      alpha	=	(r2 + Rb*Rb - 2*sqrt(r2)*Rb) / (sigma*sigma + Rb*Rb - 2*sigma*Rb);
      part[i].alphasum	-=	alpha;		//substract contribution to qlm(i)
      part[j].alphasum	-=	alpha;		//substract contribution to qlm(j)
      AlphaSum		-=	2*alpha;	//substract contribution to Qlm

      costheta	=	dp.z/sqrt(r2);
//      phi	=	atan2( dp.y, dp.x ); 
      rxy		=	sqrt(dp.x*dp.x + dp.y*dp.y);
      cosphi	=	dp.x / rxy;
      sinphi	=	dp.y / rxy;
      if (rxy<ZERO) {
         cosphi	=	0;
         sinphi	=	0;
      }

      for (m=-l; m<=l; m++) {
        termi	=	sfharmonics2(l, m, costheta, cosphi, sinphi) * alpha;
	termj	=	sfharmonics2(l, m, -costheta, -cosphi, -sinphi) * alpha;
 
//	termi	=	sfharmonics(l, m, costheta, phi) * alpha;
//	termj	=	sfharmonics(l, m, -costheta, pi+phi) * alpha;	//reverse the vector direction
        part[i].ylmalphasum[m+l]	-=	termi;		//contribution to qlm(i)
	part[j].ylmalphasum[m+l]	-=	termj;		//contribution to qlm(j)
	YlmAlphaSum[m+l]	-=	(termi + termj);	//contribution to Qlm
      }
    }
   
    dp.x	=	part[j].p.x	-	p_new.x;	//new i, j distance
    dp.y	=	part[j].p.y	-	p_new.y;
    dp.z	=	part[j].p.z	-	p_new.z;
    MapInBox(&dp);
    r2		=	dp.x*dp.x + dp.y*dp.y + dp.z*dp.z;

    if (r2 < Rb*Rb) {
      //alpha	=	1.0;
      alpha	=	(r2 + Rb*Rb - 2*sqrt(r2)*Rb) / (sigma*sigma + Rb*Rb - 2*sigma*Rb);

      part[i].alphasum	+=	alpha;		//substract contribution to qlm(i)
      part[j].alphasum	+=	alpha;		//substract contribution to qlm(j)
      AlphaSum		+=	2*alpha;	//substract contribution to Qlm

      costheta	=	dp.z/sqrt(r2);
//      phi	=	atan2( dp.y, dp.x ); 
      rxy	=	sqrt(dp.x*dp.x + dp.y*dp.y);
      cosphi	=	dp.x / rxy;
      sinphi	=	dp.y / rxy;
      if (rxy<ZERO) {
         cosphi	=	0;
         sinphi	=	0;
      }
      for (m=-l; m<=l; m++) {
	termi	=	sfharmonics2(l, m, costheta, cosphi, sinphi) * alpha;
	termj	=	sfharmonics2(l, m, -costheta, -cosphi, -sinphi) * alpha;
//	termi	=	sfharmonics(l, m, costheta, phi) * alpha;
//	termj	=	sfharmonics(l, m, -costheta, pi+phi) * alpha;
        part[i].ylmalphasum[m+l]	+=	termi;		
	part[j].ylmalphasum[m+l]	+=	termj;	
	YlmAlphaSum[m+l]	+=	(termi + termj);
      }
    }

    if (fabs(part[j].alphasum) < 1e-14) {		//update qlm(j) 
      if (part[j].alphasum < 0)
         printf("alpha_j = %10.8f\n", part[j].alphasum);
      for (m=-l; m<=l; m++) {			
        part[j].qlm[m+l]	=	0;
      }
      part[j].ql	=	0;		//update ql(j)
    }
    else {
      part[j].ql	=	0;		//don't forget initialization
      for (m=-l; m<=l; m++) {
        part[j].qlm[m+l]	=	part[j].ylmalphasum[m+l] / part[j].alphasum;
	part[j].ql		+=	part[j].qlm[m+l] * conj(part[j].qlm[m+l]);
      }
      part[j].ql	=	sqrt(part[j].ql * 4 * pi/(2*l+1));
    }
  }    
  }}

  if (fabs(part[i].alphasum) < 1e-14) {		//now it is time to update qlm(i)
      if (part[i].alphasum < 0)
         printf("alpha_i = %10.8f\n", part[i].alphasum);
    for (m=-l; m<=l; m++) {
      part[i].qlm[m+l]	=	0;
    }
    part[i].ql	=	0;
  }
  else {
    part[i].ql	=	0;
    for (m=-l; m<=l; m++) {
      part[i].qlm[m+l]	=	part[i].ylmalphasum[m+l] / part[i].alphasum;
      part[i].ql	+=	part[i].qlm[m+l] * conj(part[i].qlm[m+l]);
    }
    part[i].ql	=	sqrt(part[i].ql * 4 * pi/(2*l+1));
  }	

  Ql	=	0;			//update Qlm and Ql
  if (fabs(AlphaSum) < 1e-14) {
      if (AlphaSum < 0)
         printf("Alpha = %10.8f\n", AlphaSum);
    for (m=-l; m<=l; m++) {
      Qlm[m+l]	=	0;
    }
    Ql	=	0;
  }
  else {
    for (m=-l; m<=l; m++) {
      Qlm[m+l]	=	YlmAlphaSum[m+l] / AlphaSum;
      Ql	+=	Qlm[m+l] * conj(Qlm[m+l]);
    }
  }
  Ql	=	sqrt(Ql * 4 * pi/ (2*l+1));

  return;

/*
  //this alternative doesn't rely on the vlistplus, but actually due to the fact that
  //oldmol.vlist is changed when vlist is updated (because we transfer only pointer
  //rather than the whole set of array from part[n] to oldmol.

  for (jj=0; jj<oldmol.nverlet; jj++) {		//search the old vlist of particle i, 
						//substract their contribution 
    j	=	oldmol.vlist[jj];

    dp.x	=	part[j].p.x	-	oldmol.p.x;	//old i, j distance
    dp.y	=	part[j].p.y	-	oldmol.p.y;
    dp.z	=	part[j].p.z	-	oldmol.p.z;
    MapInBox(&dp);
    r2		=	dp.x*dp.x + dp.y*dp.y + dp.z*dp.z;
    
    if (r2 < BONDCUTSQ) {
      alpha	=	(r2 + BONDCUTSQ - 2*sqrt(r2*BONDCUTSQ)) / 
			(sigma*sigma + BONDCUTSQ - 2*sigma*BONDCUT);
      part[i].alphasum	-=	alpha;		//substract contribution to qlm(i)
      part[j].alphasum	-=	alpha;		//substract contribution to qlm(j)
      AlphaSum		-=	2*alpha;	//substract contribution to Qlm

      costheta	=	dp.z/sqrt(r2);
      phi	=	atan2( dp.y, dp.x ); 
      for (m=-l; m<=l; m++) {
	termi	=	sfharmonics(l, m, costheta, phi) * alpha;
	termj	=	sfharmonics(l, m, -costheta, pi+phi) * alpha;	//note the vector direction change
        part[i].ylmalphasum[m+l]	-=	termi;			//contribution to qlm(i)
	part[j].ylmalphasum[m+l]	-=	termj;			//contribution to qlm(j)
	YlmAlphaSum[m+l]		-=	(termi + termj);	//contribution to Qlm
      }
    }
  }

  for (jj=0; jj<part[i].nverlet; jj++) {	//search the new vlist of particle i, 
						//plus their contribution
    j	=	part[i].vlist[jj];

    dp.x	=	part[j].p.x	-	part[i].p.x;
    dp.y	=	part[j].p.y	-	part[i].p.y;
    dp.z	=	part[j].p.z	-	part[i].p.z;
    MapInBox(&dp);
    r2		=	dp.x * dp.x + dp.y * dp.y + dp.z * dp.z;

    if (r2 < BONDCUTSQ) {
      alpha	=	(r2 + BONDCUTSQ - 2*sqrt(r2*BONDCUTSQ)) /
			(sigma * sigma + BONDCUTSQ - 2*sigma*BONDCUT);
      part[i].alphasum	+=	alpha;
      part[j].alphasum	+=	alpha;
      AlphaSum		+=	2*alpha;

      costheta	=	dp.z/sqrt(r2);
      phi	=	atan2( dp.y, dp.x);
      for (m=-l; m<=l; m++) {
	termi	=	sfharmonics(l, m, costheta, phi) * alpha;
	termj	=	sfharmonics(l, m, -costheta, pi+phi) * alpha;
	part[i].ylmalphasum[m+l]	+=	termi;
	part[j].ylmalphasum[m+l]	+=	termj;
	YlmAlphaSum[m+l]		+=	(termi + termj);
      }
    }
  }	

  for (jj=0; jj<oldmol.nverlet; jj++) {		//update qlm(j) in the old vlist of i
  						//do this after the update of new vlist 
    j	=	oldmol.vlist[jj];		//because some j's in the old vlist are 
    if (part[j].alphasum < ZERO) {		//also in the new vlist
      for (m=-l; m<=l; m++)
        part[j].qlm[m+l]	=	0;
    }
    else {
      for (m=-l; m<=l; m++) 
        part[j].qlm[m+l]	=	part[j].ylmalphasum[m+l] / part[j].alphasum;
    }
  }
  for (jj=0; jj<part[i].nverlet; jj++) {	//update qlm(j) in the new vlist of i
   
    j	=	part[i].vlist[jj];
    if (part[j].alphasum < ZERO) {
      for (m=-l; m<=l; m++)
        part[j].qlm[m+l]	=	0;
    }
    else {
      for (m=-l; m<=l; m++) 
        part[j].qlm[m+l]	=	part[j].ylmalphasum[m+l] / part[j].alphasum;
    }
  }

  if (part[i].alphasum < ZERO) {		//now it is time to update qlm(i)
    for (m=-l; m<=l; m++)
      part[i].qlm[m+l]	=	0;
  }
  else {
    for (m=-l; m<=l; m++)
      part[i].qlm[m+l]	=	part[i].ylmalphasum[m+l] / part[i].alphasum;
  }	

  Ql	=	0;				//update Qlm and Ql
  if (AlphaSum < ZERO) {
    for (m=-l; m<=l; m++) {
      Qlm[m+l]	=	0;
    }
    Ql	=	0;
  }
  else {
    for (m=-l; m<=l; m++) {
      Qlm[m+l]	=	YlmAlphaSum[m+l] / AlphaSum;
      Ql	+=	Qlm[m+l] * conj(Qlm[m+l]);
    }
  }
  Ql	=	sqrt(Ql * 4 * pi/ (2*l+1));

  return;
*/
}


int Qlbinfinder(double Ql)	//determine which Ql bin the sampled Ql value belongs to
{
   return	(int) ((Ql - Qlmiddle + Qlbins * Qlbinsize/2.0) / Qlbinsize); 
}

int NMAXbinfinder(int nmax)
{
   return	(int) ((nmax - NMAXmiddle + NMAXbins * NMAXbinsize/2.0) / NMAXbinsize);	
}

/*
double etaQl(double Ql)		//given etaQ[bins], find right eta using extrapolation and intrapolation
{
   int		i, bin;
   double 	eta;
   double	Qlmidpt[Qlbins];	//middle point of each Ql bin

   for (i=0; i<Qlbins; i++) {
      Qlmidpt[i]	=	(i+0.5-Qlbins/2.0) * Qlbinsize + Qlmiddle;
   }

   bin	=	Qlbinfinder(Ql);

   if ( fabs(Ql-Qlmidpt[bin]) < ZERO) {
      eta	=	etaQ[bin];
   }
   else if (Ql > Qlmidpt[bin]) {
      if (bin == Qlbins-1) {		//the rightmost bin, use extrapolation
	 eta	=	(Ql - Qlmidpt[bin-1])/Qlbinsize * (etaQ[bin] - etaQ[bin-1]) + etaQ[bin-1];
      }
      else {				//not the rightmost bin, use intrapolation
	 eta	=	(Ql - Qlmidpt[bin])/Qlbinsize * (etaQ[bin+1]-etaQ[bin]) + etaQ[bin];
      }
   }
   else {
      if (bin == 0) {			//the leftmost bin, use extrapolation
         eta	=	(Ql - Qlmidpt[bin+1])/Qlbinsize * (etaQ[bin+1]-etaQ[bin]) + etaQ[bin+1];
      }
      else {
	 eta	=	(Ql - Qlmidpt[bin-1])/Qlbinsize * (etaQ[bin]-etaQ[bin-1]) + etaQ[bin-1];
      }
   }
   return	eta;
}
*/

double etaQl(double Ql) 	//given Ql, use a fixed quadratic eta-Ql relation to calculate eta
{
   return	-0.5 * kQ * (Ql-Qlmiddle) * (Ql-Qlmiddle);	
}

double etaNMAX(long nmax)
{
   return	-0.5 * kN * ((double)nmax - NMAXmiddle) * ((double)nmax - NMAXmiddle);
}

double etaP2(double p2)
{
   return	-0.5 * kP2 * (p2 - P2middle) * (p2 - P2middle);
}

/*
double etaNMAX(int MAXSIZE)
{
   int		i, bin;
   double	etaN;
   double	NMAXmidpt[NMAXbins];

   for (i=0; i<NMAXbins; i++) {
      NMAXmidpt[i]	=	(i+0.5-NMAXbins/2.0) * NMAXbinsize + NMAXmiddle;
   }
   bin	=	NMAXbinfinder(MAXSIZE);   
   if ( fabs(MAXSIZE-NMAXmidpt[bin]) < ZERO) {
      etaN	=	eta[bin];
   }
   else if ( MAXSIZE > NMAXmidpt[bin]) {
      if (bin==NMAXbins-1) {
	 etaN	=	(MAXSIZE - NMAXmidpt[bin-1])/NMAXbinsize * (eta[bin] - eta[bin-1]) + eta[bin-1];	
      }
      else {
	 etaN	=	(MAXSIZE - NMAXmidpt[bin])/NMAXbinsize * (eta[bin+1] - eta[bin]) + eta[bin];	
      }
   }
   else {
      if (bin==0) {
	 etaN	=	(MAXSIZE - NMAXmidpt[bin+1])/NMAXbinsize * (eta[bin+1] - eta[bin]) + eta[bin+1];
      }
      else {
	 etaN	=	(MAXSIZE - NMAXmidpt[bin-1])/NMAXbinsize * (eta[bin] - eta[bin-1]) + eta[bin-1];
      }
   }
   return	etaN; 
}
*/
#endif	/* TEST */

#ifdef SAMPLE_CODE

void Calc_Qlm(long L)
{
   long		m;	// m=[-L, L];
   long		i, j, k, n, system;
   molstruct	*moli, *molj;
   vector	dp;
   double	r2, rxy, alpha, costheta, cosphi, sinphi,
  		AlphaSum[MAXNSYSTEMS];			// has to be MAXNSYSTEMS because NSYSTEMS unspecified
   double complex	YlmAlphaSum[MAXNSYSTEMS][2L+1], Qlm[MAXNSYSTEMS][2L+1];

#ifdef CELL_LIST
   cellstruct	*celli, *cellj, *celln;
#endif

   /* Initialization */

   for (system=0; system<NSYSTEMS; system++) {
      for (m=-L; m<=L; m++) {
         YlmAlphaSum[system][m+L]	=	0.0;
	 Qlm[system][m+L]		=	0.0;
      }
      AlphaSum[system]	=	0.0;
      Q6[system]	=	0.0;
   }
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
         for (m=-L; m<=L; m++) {
	    moli->qlm[i][m+L]		=	0.0;
	    moli->ylmalphasum[i][m+L]	=	0.0;
         }
         moli->alphasum[i]	=	0.0;
	 moli->q6[i]		=	0.0;
      }
   }
   /* Calculation */

   for (moli=mol; moli<mol+NMOLS; moli++) {
      system	=	moli->box;

      for (i=0; i<moli->nsites; i++) {

#ifdef CELL_LIST
	 celli	=	moli->cell[i];			// cell id of moli->i
	 for (n=0; n<celli->nneigh; n++) {
	    celln	=	celli->neigh[n];	// neighboring cell id

	    for (k=0; k<celln->nsites; k++) {		// sites in neighboring cell
	       if (molj=celln->mol[k]) {
		  j	=	celln->molsite[k];	// pick up one

                  if (moli!=molj || i!=j) {
#else
	 for (molj=mol; molj<mol+NMOLS; molj++) {	// note: molj starts from moli
	    if (molj->box == system) {
               for (j=0; j<molj->nsites; j++) {
                  if (molj>moli || j>i) {
#endif
                     dp	=	V_Subtr(molj->p+j, moli->p+i);
                     dp	=	MapInBox2(&dp, PBC, system);
		     r2	=	dp.x * dp.x + dp.y * dp.y + dp.z * dp.z;

		     if (r2 < Rb2) {
		        //alpha		=	1.0;
		        alpha	=	(r2 + Rb2 - 2*sqrt(r2)*Rb)/(1 + Rb2-2*Rb);
		        moli->alphasum[i]	+=	alpha;
		        molj->alphasum[j]	+=	alpha;

		        costheta	=	dp.z/sqrt(r2);
		        rxy	=	sqrt(dp.x*dp.x + dp.y*dp.y);
		        cosphi	=	(rxy<ZERO ? 0 : dp.x/rxy);
		        sinphi	=	(rxy<ZERO ? 0 : dp.y/rxy);
			
//			printf("costheta = %f cosphi = %f sinphi = %f\n", costheta, cosphi, sinphi);
                        for (m=-L; m<=L; m++) {
			   moli->ylmalphasum[i][m+L]	+=	
					sfharmonics2(L, m, costheta, cosphi, sinphi)*alpha;
			   molj->ylmalphasum[j][m+L]	+=	
					sfharmonics2(L, m, -costheta, -cosphi, -sinphi)*alpha;
		        }
                     }
		  }
         }  }  }
      }
   }

   /* NO need to correct double counting here because both numerator and denominator are double-counted */
   // only cell list case will result in double counting, no cell list will avoid double counting, see code

   /* Normalization */

   for (moli=mol; moli<mol+NMOLS; moli++) {
      system	=	moli->box;

      for (i=0; i<moli->nsites; i++) {
         if (moli->alphasum[i] < 0)			// alphasum should always > 0
	    Exit("sample.c", "Calc_Qlm", "alphasum < 0");
         else if (moli->alphasum[i] > ZERO) {		// always exclude divided-by-zero
	    for (m=-L; m<=L; m++) {
	       moli->qlm[i][m+L]	=	moli->ylmalphasum[i][m+L] / moli->alphasum[i];
       
               if (6==L)
                  moli->q6[i]		+=	moli->qlm[i][m+L] * conj(moli->qlm[i][m+L]);
               else if (4==L)
                  moli->q4[i]		+=	moli->qlm[i][m+L] * conj(moli->qlm[i][m+L]);
            }
            if (6==L)
               moli->q6[i]		=	sqrt(moli->q6[i] * 4 * pi/(2*L+1));
            else if (4==L)
               moli->q4[i]		=	sqrt(moli->q4[i] * 4 * pi/(2*L+1));

            AlphaSum[system]	+=	moli->alphasum[i];
            for (m=-L; m<=L; m++)
               YlmAlphaSum[system][m+L]	+=	moli->ylmalphasum[i][m+L];	// ylmalphasum=0 if alphasum=0
         }
      }
   }

   for (system=0; system<NSYSTEMS; system++) {
      if (AlphaSum[system] < 0)
	 Exit("sample.c", "Calc_Qlm", "AlphaSum < 0");
      else if (AlphaSum[system] > ZERO) {
         for (m=-L; m<=L; m++) {
            Qlm[system][m+L]	=	YlmAlphaSum[system][m+L] / AlphaSum[system];
	    
	    if (6==L)
	       Q6[system]	+=	Qlm[system][m+L]  * conj(Qlm[system][m+L]);
	    else if (4==L)
	       Q4[system]	+=	Qlm[system][m+L]  * conj(Qlm[system][m+L]);
         }
         if (6==L)
            Q6[system]		=	sqrt(Q6[system] * 4 * pi/(2*L+1));
         else if (4==L)
            Q4[system]		=	sqrt(Q4[system] * 4 * pi/(2*L+1));
      }
   }
}


double qlproductSQ(int l, molstruct *moli, long i, molstruct *molj, long j)	//calculate q_l(i) dot q_l(j) SQ
{						
  int 		m;
  float complex	qlproduct=0;
  double	sumi=0, sumj=0, sum;

  for (m=-l; m<=l; m++) {
    qlproduct	+=	moli->qlm[i][m+l] * conj(molj->qlm[j][m+l]);
    sumi	+=	moli->qlm[i][m+l] * conj(moli->qlm[i][m+l]);
    sumj	+=	molj->qlm[j][m+l] * conj(molj->qlm[j][m+l]);
  }
  sum		=	sumi * sumj;

  if (fabs(sum)<ZERO) 
    return	0;
  else 
    return 	qlproduct * qlproduct / sum;
}
#endif //SAMPLE_CODE

vector CenterofMass(molstruct *moli)		// center of mass of a chain
{
   long		i;
   double	mx, my, mz, m, M;
   vector	p;

   m	=	0.0;			// mass of one site
   mx	=	0.0;			// mass weighted x coordinate
   my	=	0.0;
   mz	=	0.0;
   M	=	0.0;			// total mass of one molecule

   for (i=0; i<moli->nsites; i++) {
      m		=	type[moli->type[i]].M;
      mx	+=	m * (moli->p[i].x);
      my	+=	m * (moli->p[i].y);
      mz	+=	m * (moli->p[i].z);
      M		+=	m;
	//printf("%f\t%f\t%f\t%f\t%f\n", m, mx, my, mz, M);
   }
   p.x	=	mx/M;
   p.y	=	my/M;
   p.z	=	mz/M;

   return	p;
}

vector groupCoM(beadstruct *group, long nsites)	// June 6,2009
{
   long		i, site; 
   double	mx, my, mz, m, M;
   vector	p;
   molstruct	*moli;

   m		=	0.0;
   mx = my = mz =	0.0;
   M		=	0.0;

   for (i=0; i<nsites; i++) {
      moli	=	group[i].moli;
      site	=	group[i].site;
      m		=	type[moli->type[site]].M;
      mx	+=	m * (moli->p[i].x);
      my	+=	m * (moli->p[i].y);
      mz	+=	m * (moli->p[i].z);
      M		+=	m;
   }
   p.x	=	mx/M;
   p.y	=	my/M;
   p.z	=	mz/M;
   return	p;
}


float R2_gyration(molstruct *molm)	// squared radius of gyration of one molecule
{					// chains need to be continuous in space
   int		i;
   vector	rmean;
   float	Rg2;

   rmean.x	=	0;
   rmean.y	=	0;
   rmean.z	=	0;
   
   for (i=0; i<molm->nsites; i++) {
      rmean.x	+=	molm->p[i].x;
      rmean.y	+=	molm->p[i].y;
      rmean.z	+=	molm->p[i].z;
   }
   rmean.x	/=	molm->nsites;
   rmean.y	/=	molm->nsites;
   rmean.z	/=	molm->nsites;

   Rg2	=	0;

   for (i=0; i<molm->nsites; i++) {
      Rg2	+=	(molm->p[i].x - rmean.x) * (molm->p[i].x - rmean.x);
      Rg2	+=	(molm->p[i].y - rmean.y) * (molm->p[i].y - rmean.y);
      Rg2	+=	(molm->p[i].z - rmean.z) * (molm->p[i].z - rmean.z);
   }
   Rg2	/=	molm->nsites;

   return	Rg2;
}


float rx2_gyration(molstruct *molm)
{
   int		i;
   float	xmean;
   float	rx2;

   xmean	=	0.0;
   for (i=0; i<molm->nsites; i++) {
      xmean	+=	molm->p[i].x;
   }
   xmean	/=	molm->nsites;

   rx2	=	0.0;
   for (i=0; i<molm->nsites; i++) {
      rx2	+=	(molm->p[i].x - xmean) * (molm->p[i].x - xmean);
   }
   rx2	/=	molm->nsites;
   return	rx2;
}

float ry2_gyration(molstruct *molm)
{
   int		i;
   float	ymean;
   float	ry2;

   ymean	=	0.0;
   for (i=0; i<molm->nsites; i++) {
      ymean	+=	molm->p[i].y;
   }
   ymean	/=	molm->nsites;

   ry2	=	0.0;
   for (i=0; i<molm->nsites; i++) {
      ry2	+=	(molm->p[i].y - ymean) * (molm->p[i].y - ymean);
   }
   ry2	/=	molm->nsites;
   return	ry2;
}

float rz2_gyration(molstruct *molm)
{
   int		i;
   float	zmean;
   float	rz2;

   zmean	=	0.0;
   for (i=0; i<molm->nsites; i++) {
      zmean	+=	molm->p[i].z;
   }
   zmean	/=	molm->nsites;

   rz2	=	0.0;
   for (i=0; i<molm->nsites; i++) {
      rz2	+=	(molm->p[i].z - zmean) * (molm->p[i].z - zmean);
   }
   rz2	/=	molm->nsites;
   return	rz2;
}

float group_Rg2(vector *rbead, int nsites)	// R2 gyration of one group
{
   int		i;
   vector	rmean;
   float	R2;

   rmean.x	=	0;
   rmean.y	=	0;
   rmean.z	=	0;
   
   for (i=0; i<nsites; i++) {
      rmean.x	+=	rbead[i].x;
      rmean.y	+=	rbead[i].y;
      rmean.z	+=	rbead[i].z;
   }
   rmean.x	/=	nsites;
   rmean.y	/=	nsites;
   rmean.z	/=	nsites;

   R2	=	0;

   for (i=0; i<nsites; i++) {
      R2	+=	(rbead[i].x - rmean.x) * (rbead[i].x - rmean.x);
      R2	+=	(rbead[i].y - rmean.y) * (rbead[i].y - rmean.y);
      R2	+=	(rbead[i].z - rmean.z) * (rbead[i].z - rmean.z);
   }
   R2	/=	nsites;

   return	R2;
}

//-------------------------------------------------------------------
// geometrical center of a group of atoms, i.e., not weighted by mass
// added 6/12/2016
//
vector group_com(vector *rbead, long nsites)	
{
   int		i;
   vector	com;

   com.x = 0;
   com.y = 0;
   com.z = 0;
   for (i=0; i<nsites; i++) {
      com.x += rbead[i].x;
      com.y += rbead[i].y;
      com.z += rbead[i].z;
   }
   com.x /= nsites;
   com.y /= nsites;
   com.z /= nsites;

   return com;
}
// end: group_com()
//--------------------------------------------------------------------


matrix GyrationTensor(molstruct *molm)		// Gyration tensor of one molecule
{
   long         i;
   matrix       tensor;
   double       x, y, z;
   vector       rmean;

   M_Null(&tensor);

   rmean.x	=	0;
   rmean.y	=	0;
   rmean.z	=	0;
   
   for (i=0; i<molm->nsites; i++) {
      rmean.x	+=	molm->p[i].x;
      rmean.y	+=	molm->p[i].y;
      rmean.z	+=	molm->p[i].z;
   }
   rmean.x	/=	molm->nsites;
   rmean.y	/=	molm->nsites;
   rmean.z	/=	molm->nsites;

   for (i=0; i<molm->nsites; i++) {
      x         =       molm->p[i].x;
      y         =       molm->p[i].y;
      z         =       molm->p[i].z;

      tensor.x.x        +=      x * x;
      tensor.x.y        +=      x * y;
      tensor.x.z        +=      x * z;
      tensor.y.y        +=      y * y;
      tensor.y.z        +=      y * z;
      tensor.z.z        +=      z * z;
   }
   tensor	=	M_Mult(1.0/molm->nsites, &tensor);

   tensor.x.x   -=      rmean.x * rmean.x; 	// substract the center of mass contribution
   tensor.x.y   -=      rmean.x * rmean.y;
   tensor.x.z   -=      rmean.x * rmean.z;
   tensor.y.y   -=      rmean.y * rmean.y;
   tensor.y.z   -=      rmean.y * rmean.z;
   tensor.z.z   -=      rmean.z * rmean.z;
   tensor.y.x   =       tensor.x.y;             // a symmetric matrix by definition
   tensor.z.x   =       tensor.x.z;
   tensor.z.y   =       tensor.y.z;

   return       tensor;
}


matrix groupGyraTensor(beadstruct *group, long nsites)		// Gyration tensor of 
								// a group of particles
{
   long         i, site;
   matrix       tensor;
   double       x, y, z;
   vector       rmean;
   molstruct	*molm;

   M_Null(&tensor);

   rmean.x	=	0;
   rmean.y	=	0;
   rmean.z	=	0;

   for (i=0; i<nsites; i++) {		// center of mass, NO mass weighted
      molm	=	group[i].moli;   
      site	=	group[i].site;

      rmean.x	+=	molm->p[site].x;
      rmean.y	+=	molm->p[site].y;
      rmean.z	+=	molm->p[site].z;
   }
   rmean.x	/=	nsites;
   rmean.y	/=	nsites;
   rmean.z	/=	nsites;

   for (i=0; i<nsites; i++) {
      molm	=	group[i].moli;
      site	=	group[i].site;

      x         =       molm->p[site].x;
      y         =       molm->p[site].y;
      z         =       molm->p[site].z;

      tensor.x.x        +=      x * x;
      tensor.x.y        +=      x * y;
      tensor.x.z        +=      x * z;
      tensor.y.y        +=      y * y;
      tensor.y.z        +=      y * z;
      tensor.z.z        +=      z * z;
   }
   tensor	=	M_Mult(1.0/nsites, &tensor);

   tensor.x.x   -=      rmean.x * rmean.x; 	// substract the center of mass contribution
   tensor.x.y   -=      rmean.x * rmean.y;
   tensor.x.z   -=      rmean.x * rmean.z;
   tensor.y.y   -=      rmean.y * rmean.y;
   tensor.y.z   -=      rmean.y * rmean.z;
   tensor.z.z   -=      rmean.z * rmean.z;
   tensor.y.x   =       tensor.x.y;             // a symmetric matrix by definition
   tensor.z.x   =       tensor.x.z;
   tensor.z.y   =       tensor.y.z;

   return       tensor;
}


vector fshape(vector *eigvalue)		// shape descriptors, given eigvalues of gyration tensor
{
   vector	fshape;
   double	a, b, c, temp;

   a	=	eigvalue->x;
   b	=	eigvalue->y;
   c	=	eigvalue->z;

   if (a>b) {				// arrange the ordering so that a < b < c
      temp=b;	b=a;	a=temp;
   }
   if (b>c) {
      temp=c; 	c=b; 	b=temp;
   }
   if (a>b) {
      temp=b; 	b=a; 	a=temp;
   }
   fshape.x	=	c - 0.5* (a+b);		// asphericity
   fshape.y	=	b - a;			// acylindricity
   fshape.z	=	((fshape.x * fshape.x) + 0.75*(fshape.y * fshape.y)) / ((a+b+c)*(a+b+c));
						// relative shape anisotropy between 0 and 1
   return	fshape;
}

matrix InertiaTensor(molstruct *molm)	// added 1/30/08, calculate moment of inertia tensor
{					// inertia tensor is a mass weighted Gyration tensor
   long         i;
   matrix       tensor;
   double       x, y, z, m, mass;
   vector       com;

   M_Null(&tensor);
   com  =       CenterofMass(molm);

   mass	=	0;
   for (i=0; i<molm->nsites; i++) {
      mass	+=	type[molm->type[i]].M;
   }

   for (i=0; i<molm->nsites; i++) {
      x         =       molm->p[i].x - com.x;
      y         =       molm->p[i].y - com.y;
      z         =       molm->p[i].z - com.z;
      m         =       type[molm->type[i]].M;

      tensor.x.x        +=       m * (y*y + z*z);
      tensor.x.y        +=      -m * (x*y);
      tensor.x.z        +=      -m * (x*z);
      tensor.y.y        +=       m * (x*x + z*z);
      tensor.y.z        +=      -m * (y*z);
      tensor.z.z        +=       m * (x*x + y*y);
   }
   tensor.y.x   =       tensor.x.y;             //it is a symmetric matrix by definition
   tensor.z.x   =       tensor.x.z;
   tensor.z.y   =       tensor.y.z;

   return       tensor;
}


matrix groupInerTensor(beadstruct *group, long nsites)	// added June 6,2009, 
							// calculate moment of inertia tensor
							// of a group of beads		
{
   long         i, site;
   matrix       tensor;
   double       x, y, z, m, mass;
   vector       com;
   molstruct	*molm;

   M_Null(&tensor);
   com  =       groupCoM(group, nsites);

   mass	=	0;
   for (i=0; i<nsites; i++) {
      molm	=	group[i].moli;
      site	=	group[i].site;
      mass	+=	type[molm->type[site]].M;
   }

   for (i=0; i<nsites; i++) {
      molm	=	group[i].moli;
      site	=	group[i].site;

      x         =       molm->p[site].x;
      y         =       molm->p[site].y;
      z         =       molm->p[site].z;
      m         =       type[molm->type[site]].M;

      tensor.x.x        +=      m * ((y-com.y)*(y-com.y) + (z-com.z)*(z-com.z));
      tensor.x.y        +=      -m * (x-com.x) * (y-com.y);
      tensor.x.z        +=      -m * (x-com.x) * (z-com.z);
      tensor.y.y        +=      m * ((x-com.x)*(x-com.x) + (z-com.z)*(z-com.z));
      tensor.y.z        +=      -m * (y-com.y) * (z-com.z);
      tensor.z.z        +=      m * ((x-com.x)*(x-com.x) + (y-com.y)*(y-com.y));
   }
   tensor.y.x   =       tensor.x.y;             //it is a symmetric matrix by definition
   tensor.z.x   =       tensor.x.z;
   tensor.z.y   =       tensor.y.z;

   return       tensor;
}

//==============================================================//
//	grpgytensor(): gyration tensor of a group of beads	//
//		       added on 5/9/2013			//
//==============================================================//
matrix	grpgytensor(int size, vector *rbead)
{
   int		i;
   float	m = 1.0;
   float	x, y, z;
   vector	com;
   matrix	tensor;

   V_Null(&com);
   for (i=0; i<size; i++) {
      com	=	V_Add(&com, rbead+i);
   }
   com		=	V_Mult(1.0/size, &com);		// center of mass
   
   M_Null(&tensor);
   for (i=0; i<size; i++) {
      x		=	rbead[i].x - com.x;
      y		=	rbead[i].y - com.y;
      z		=	rbead[i].z - com.z;

      tensor.x.x        +=       m * (y*y + z*z);
      tensor.x.y        +=      -m * (x*y);
      tensor.x.z        +=      -m * (x*z);
      tensor.y.y        +=       m * (x*x + z*z);
      tensor.y.z        +=      -m * (y*z);
      tensor.z.z        +=       m * (x*x + y*y);
   }
   tensor.y.x	=	tensor.x.y;             //inertia tensor is symmetric
   tensor.z.x	=	tensor.x.z;
   tensor.z.y	=	tensor.y.z;

   return	tensor;
}

//==============================================================//
//	R2_n2n(): squared end-to-end distance of one chain	//
//==============================================================//
float R2_n2n(molstruct *moli)
{
   int		i, ib=moli->box;
   vector	p, dp;

   V_Null(&p);
   dp	=	V_Subtr(moli->p, moli->p+(moli->nsites)-1);

   return	V_Dot(&dp, &dp);
}


vector	CenterofNucleus(long nuclid, molstruct *molm) 	// center of specified nucleus, which contains molm
{
   vector	com, 
		rA,			// center of one molecule in nucleus
		rO,			// center of nucleus 
		rBA, 			// rB-rA, B and A belong to same nucleus
		rOA;			// rO-rA
   molstruct	*moli;
   long		system=0, n;

   // Step 1: find one molecule A belonging to the nucleus as a reference

   rA	=	CenterofMass(molm);
   system	=	molm->box;

   // Step 2: calc. the shift of center of nucleus to this molecule

   n	=	0;					// # of molecules in this nucleus
   V_Null(&rBA);
   V_Null(&rOA);

   for (moli=mol; moli<mol+NMOLS; moli++) {
      if (moli->box == system && moli->nuclid[0] == nuclid) {
	 n	++;
	 com	=	CenterofMass(moli);
	 rBA	=	V_Subtr(&com, &rA);
	 rBA	=	MapInBox2(&rBA, PBC, system);
	 rOA	=	V_Add(&rOA, &rBA);		// assuming all molecules of same weight
      }
   }
   rOA	=	V_Mult(1.0/n, &rOA);
   rO	=	V_Add(&rA, &rOA);			// center of nucleus
   return	rO;					// rO might be out of central box, but it is fine
}

//==============================================//
//	Sample_Pressure(): calculate pressure	//
//==============================================//
void Sample_Pressure()
{
   long		i, j, k, n, nsites[2];
   double	vol, temp, virial, rc3, rc9;
   double	pres, dpres;
   double	Sigma, Epsilon;

   for (i=0; i<NBOX; i++) {
      n		=	NSites[i];			// # of atoms (not chains)
      vol	=	BOX[i].vol;
      temp	=	BOX[i].temp;
      virial	=	vir[i].tot;

      pres	=	( n*temp - virial/3) / vol;

      if (V_LJLRC) {
         if (2==NTYPES) {
            nsites[0]	=	2 * NMols[i];
            nsites[1]	=	NSites[i] - nsites[0];
            
            for (j=0; j<NTYPES; j++) {
               for (k=0; k<NTYPES; k++) {
                  Sigma		=	type[j].mix[k].SIGMA;
                  Epsilon	=	type[j].mix[k].EPSILON;

                  rc3	=	BOX[i].rc/Sigma;
                  rc3	=	rc3 * rc3 * rc3;
                  rc9	=	rc3 * rc3 * rc3;

                  dpres	=	2.0/(9*rc9) - 1.0/(3*rc3);
		  dpres	*=	nsites[j] * nsites[k] / (vol*vol);
                  dpres	*=	16.0 * pi * Epsilon * Sigma * Sigma * Sigma;

                  pres		+=	dpres;
               }
            }
         }
         if (1==NTYPES) {
            rc3	=	BOX[i].rc/type[0].SIGMA;
            rc3	=	rc3 * rc3 * rc3;
            rc9	=	rc3 * rc3 * rc3;

            dpres	=	16.0 * pi * n * n / (vol*vol) * (2.0/(9 * rc9) - 1.0/(3 * rc3));
            dpres	*=	type[0].EPSILON * type[0].SIGMA * type[0].SIGMA * type[0].SIGMA;
            pres	+=	dpres;
         }
      }
      BOX[i].pres	=	pres;
   }
}

//===========================================================================//

void SampleDrift()
{
   molstruct	*moli;
   long		system;
   vector	p1, dp;

   for (system=0; system<NSYSTEMS; system++)
      drift2[system]	=	0.0;

   for (moli=mol; moli<mol+NMOLS; moli++) 
      if ( (system=moli->box)>= 0) {
         p1		=	CenterofMass(moli);
         dp		=	V_Subtr(&p1, &(moli->origin));
         drift2[system]	+=	V_Dot( &dp, &dp );
      }

   for (system=0; system<NSYSTEMS; system++)
      drift2[system]	/=	NMols[system];

/*
   if (D_DRIFT)				// temporary
      for (system=0; system<NSYSTEMS; system++) 
         for (moli=mol; moli<mol+NMOLS; moli++)
            if (moli->box == system) {
               drift2	=	V_Dot( &(moli->drift), &(moli->drift) );
               PutInDistribution(D_Drift+system, drift2, 1.0, 1.0);
            }
*/
}


void SampleSpherical()		// sample spherical coordinates distribution, 11/12/07
{				// given spherical coordinates updated
   molstruct	*moli;
   long		system, i, j, k, n[MAXNSYSTEMS];
   long		t1, t2, t3;		// value list 0: trans, 1: g+, 2: g-
   double	phi, phi2, phi3;
   double	HB=2.0*M_PI/3, LB=-HB;

   for (system=0; system<NSYSTEMS; system++) {
      n[system]		=	0;
      transfrac[system]	=	0.0;
      for (i=0; i<3; i++) {
         for (j=0; j<3; j++) {
            Ptors2[system][i][j]	=	0.0;
            for (k=0; k<3; k++) {
               Ptors3[system][i][j][k]	=	0.0;
      }  }  }
   }

   AllSpherical();		// calculate spherical coordinates

   for (system=0; system<NSYSTEMS; system++)
      for (moli=mol; moli<mol+NMOLS; moli++)
         if (moli->box == system) {
            if (D_TORSION) {
               for (i=3; i<moli->nsites; i++) {
                  phi	=	AdjustAngle(moli->s[i].beta);		// torsion angle
                  //PutInDistribution(D_Torsion+system, phi, 1.0, 1.0);
		  n[system]	++;

                  if (phi >=HB || phi <LB)		t1	=	0;
                  else if (phi >=0 && phi <HB)		t1	=	1;
                  else					t1	=	2;
 
                  if (t1==0)	transfrac[system]	+=	1.0;

                  if (i+1 <moli->nsites) {
                     phi2	=	AdjustAngle(moli->s[i+1].beta);
                     if (phi2 >=HB || phi2 <LB)		t2	=	0;
                     else if (phi2 >=0 && phi2 <HB)	t2	=	1;
                     else				t2	=	2;
                     Ptors2[system][t1][t2]		+=	1.0;

                     if (i+2 <moli->nsites) {
                        phi3	=	AdjustAngle(moli->s[i+2].beta);
                        if (phi3 >=HB || phi3 <LB)		t3	=	0;
                        else if (phi3 >=0 && phi3 <HB)		t3	=	1;
                        else					t3	=	2;
                        Ptors3[system][t1][t2][t3]	+=	1.0;
                     }
                  }
	       }
/*
            if (D_BONDA)
               for (i=2; i<moli->nsites; i++)
                  PutInDistribution(D_Bonda +system, AdjustAngle(moli->s[i].alpha), 1.0, 1.0);
            if (D_BONDL)
               for (i=1; i<moli->nsites; i++)
                  PutInDistribution(D_Bondl +system, moli->s[i].d, 1.0, 1.0);
*/
         }  }

   for (system=0; system<NSYSTEMS; system++) {
      transfrac[system]	/=	n[system];
   }
   return;
}


void Dist_Spherical()
{
   long		system, i;
   molstruct	*moli;

   for (system=0; system <NSYSTEMS; system++) {
      for (moli=mol; moli<mol+NMOLS; moli++) {
         if (moli->box == system) {
            if (D_TORSION) {
               for (i=3; i<moli->nsites; i++) 
                  PutInDistribution(D_Torsion + system, AdjustAngle(moli->s[i].beta), 1.0, 1.0);
            }
            if (D_BONDA) {
               for (i=2; i<moli->nsites; i++)
                  PutInDistribution(D_Bonda + system, AdjustAngle(moli->s[i].alpha), 1.0, 1.0);
            }
            if (D_BONDL) {
               for (i=1; i<moli->nsites; i++)
                  PutInDistribution(D_Bondl + system, moli->s[i].d, 1.0, 1.0);
            }
         }
   }  }
   return;
}


void SampleRadial()
{
   long		i, j, system;
   molstruct	*moli, *molj;

   if (D_RADIAL)
      for (moli=mol; moli<mol+NMOLS; moli++)
         if ( (system=moli->box) >=0)
            for (molj=moli; molj<mol+NMOLS; molj++)
               if ( molj->box == system )
                  for (i=0; i<moli->nsites; i++)
                     for (j=0; j<molj->nsites; j++)
                        if ( moli==molj ? (j>i+2) : 1) {
                           PutInDistribution(D_Radial +system, sqrt(DistSQ(moli->p[i], molj->p[j], system)), 1.0, 1.0);
			}
}


void SampleP2()				// sample global orientational order
{
   molstruct	*moli, *molj;
   long		i, j, system, n[MAXNSYSTEMS], nz[MAXNSYSTEMS];
   vector	vi, vj;
   double	idotj;

   for (i=0; i<NSYSTEMS; i++) {		// global orientational order
      P2[i]	=	0.0;
      P2z[i]	=	0.0;
      n[i]	=	0;
      nz[i]	=	0;
   }
   for (moli=mol; moli<mol+NMOLS; moli++)		// calculate P2
      if ( (system=moli->box)>=0 )
         for (i=1; i<moli->nsites-1; i++) {
            vi	=	V_Subtr(moli->p+i+1, moli->p+i-1);

            for (molj=moli; molj<mol+NMOLS; molj++)
               if ( molj->box == system )
                  for (j=1; j<molj->nsites-1; j++)
                     if ( (moli!=molj) ? 1 : j>i ) {
			vj		=	V_Subtr(molj->p+j+1, molj->p+j-1);			
                        idotj		=	V_Dot(&vi, &vj);
			idotj		*=	idotj;
			idotj		/=	V_Dot(&vi, &vi) * V_Dot(&vj, &vj);
                        P2[system]	+=	idotj;
			n[system]	++;
                     }
         }

   for (moli=mol; moli<mol+NMOLS; moli++)		// calculate P2z
      if ( (system=moli->box)>= 0)
         for (i=1; i<moli->nsites-1; i++) {
            vi		=	V_Subtr(moli->p+i+1, moli->p+i-1);
            P2z[system]	+=	(vi.z * vi.z) / V_Dot(&vi, &vi); 
            nz[system]	++;
         }
            
   for (i=0; i<NSYSTEMS; i++) {				// normalize
      P2[i]	/=	n[i];
      P2[i]	*=	1.5;
      P2[i]	-=	0.5;
      P2z[i]	/=	nz[i];
      P2z[i]	*=	1.5;
      P2z[i]	-=	0.5;
   }
   return;
}

//**********************************************************************//
//	Sample orientation order parameter, global P2 and local p2	//
//**********************************************************************//

void SampleP2All()
{
   molstruct	*moli, *molj;
   long		i, j, system, n[MAXNSYSTEMS], nz[MAXNSYSTEMS];
   vector	vi, vj;
   double	idotj, r2;
#ifdef CELL_LIST
   long		k, l;
   static cellstruct	*celli, *cellk;
#endif

time_t	start, end;

   // Initialization

   for (i=0; i<NSYSTEMS; i++) {		// global orientational order
      P2[i]	=	0.0;
      P2z[i]	=	0.0;
      P2M[i]	=	0.0;
      n[i]	=	0;
      nz[i]	=	0;
   }

   for (moli=mol; moli<mol+NMOLS; moli++) {	// local orientational order
      for (i=0; i<moli->nsites; i++) {
         moli->p2[i]	=	0.0;
         moli->np2[i]	=	0;
      }
   }

   // Calculation

   for (moli=mol; moli<mol+NMOLS; moli++)
      if ( (system=moli->box)>=0 )
         for (i=1; i<moli->nsites-1; i++) {			// exclude 2 end sites
            vi	=	V_Subtr(moli->p+i+1, moli->p+i-1);
            vi	=	MapInBox2(&vi, PBC, system);

	    // NOTE (8/27/2012): P2 calculation will not use cell list because it is global
/*
#ifdef CELL_LIST
            celli	=	moli->cell[i];
            for (l=0; l<celli->nneigh; l++) {
               cellk	=	celli->neigh[l];
               for (k=0; k<cellk->nsites; k++) {
                  if (molj=cellk->mol[k]) {
		     j	=	cellk->molsite[k];
                     if (j!=0 && j!=molj->nsites-1 && (moli!=molj || i!=j)) {
#else
*/
            for (molj=moli; molj<mol+NMOLS; molj++) {
               if ( molj->box == system ) {
                  for (j=1; j<molj->nsites-1; j++) {
                     if ( (moli!=molj) ? 1 : j>i ) {
/*
#endif	// CELL_LIST
*/
			vj		=	V_Subtr(molj->p+j+1, molj->p+j-1);			
            		vj		=	MapInBox2(&vj, PBC, system);
                        idotj		=	V_Dot(&vi, &vj);
			idotj		*=	idotj;
			idotj		/=	V_Dot(&vi, &vi) * V_Dot(&vj, &vj);
                        P2[system]	+=	idotj;
			n[system]	++;
                        
                        r2	=	DistSQ(moli->p[i], molj->p[j], system);
                        if (r2 < Rp2*type[0].SIGMA*type[0].SIGMA) {	// calc. local p2
			   moli->p2[i]	+=	idotj;
			   molj->p2[j]	+=	idotj;
			   moli->np2[i]	++;
			   molj->np2[j]	++;
			}
            }  }  }  }
         }

   // Calc. the orientation along z-axis P2z
   for (moli=mol; moli<mol+NMOLS; moli++)
      if ( (system=moli->box)>= 0)
         for (i=1; i<moli->nsites-1; i++) {
            vi		=	V_Subtr(moli->p+i+1, moli->p+i-1);
            vi		=	MapInBox2(&vi, PBC, system);
            P2z[system]	+=	(vi.z * vi.z) / V_Dot(&vi, &vi); 
            nz[system]	++;
         }
            
   // Normalization

   for (i=0; i<NSYSTEMS; i++) {
      if (n[i]) {
         P2[i]	/=	n[i];
         P2[i]	*=	1.5;
         P2[i]	-=	0.5;
      }
      else	P2[i]	=	0;
      if (nz[i]) {
         P2z[i]	/=	nz[i];
         P2z[i]	*=	1.5;
         P2z[i]	-=	0.5;
      }
      else	P2z[i]	=	0;
   }

   // Calculate modified P2M using local p2(i)
   for (i=0; i<NSYSTEMS; i++) {
      n[i]	=	0;
   }
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=1; i<moli->nsites-1; i++) {
         if (moli->np2[i]) {
            moli->p2[i]	/=	moli->np2[i];
            moli->p2[i]	*=	1.5;
	    moli->p2[i]	-=	0.5;
         }
	 else	moli->p2[i]	=	0;

         P2M[moli->box]	+=	moli->p2[i];
   	 n[moli->box]	++;
      }
      moli->p2[0]		=	moli->p2[1];
      moli->p2[moli->nsites-1]	=	moli->p2[moli->nsites-2];
   }
   for (i=0; i<NSYSTEMS; i++) {
      P2M[i]	/=	n[i];
   }
   return;
}


//==============================================================//
//	SampleP2All_MPI(): MPI version of SampleP2All()		//
//			   5/6/2013				//
//==============================================================//
void SampleP2All_MPI()
{
   molstruct	*moli, *molj;
   long		i, j, k, l, system, n[MAXNSYSTEMS], nz[MAXNSYSTEMS];
   vector	vi, vj;
   double	idotj, r2;
#ifdef CELL_LIST
   static cellstruct	*celli, *cellk;
#endif

#ifdef	myMPI
   int		proc;
   MPI_Status	status;
   int		data_nP2[MAXNSYSTEMS];				// for global P2
   double	data_P2[MAXNSYSTEMS];				// for global P2
   int		data_np2[MAXNMOLS*MAXNMOLSITES];		// for local p2's
   double	data_p2[MAXNMOLS*MAXNMOLSITES];			// for local p2's
   
   static int	molid[MAXNMOLS];
   static int	init = 1;
   int		ave, extra;

   if (init) {
      init	=	0;

      //------Assign molecules to processes, try to leave proc #0 has small load------//
      ave	=	NMOLS/MPI_numprocs;
      extra	=	NMOLS%MPI_numprocs;
      molid[0]	=	0;
      for (i=1; i<MPI_numprocs-extra+1; i++) {
	 molid[i]	=	molid[i-1] + ave;
      }
      for (i=MPI_numprocs-extra+1; i<MPI_numprocs+1; i++) {
	 molid[i]	=	molid[i-1] + ave + 1;
      }
      
      if (MPI_myid==0) {
	 printf("# Molecules assignment\n");
	 for (i=0; i<MPI_numprocs; i++) {
	    printf("# Process # %d: molecules # %d - %d\n", i, molid[i], molid[i+1]-1);
	 }
      }
      
   }

   for (i=0; i<MAXNSYSTEMS; i++) {
      data_nP2[i]	=	0;
      data_P2[i]	=	0.0;
   }
   for (i=0; i<MAXNMOLS*MAXNMOLSITES; i++) {
      data_np2[i]	=	0;
      data_p2[i]	=	0.0;
   }
   
#endif	//myMPI

   // Initialization

   for (i=0; i<NSYSTEMS; i++) {		// global orientational order
      P2[i]	=	0.0;
      P2z[i]	=	0.0;
      P2M[i]	=	0.0;
      n[i]	=	0;
      nz[i]	=	0;
   }

   for (moli=mol; moli<mol+NMOLS; moli++) {	// local orientational order
      for (i=0; i<moli->nsites; i++) {
         moli->p2[i]	=	0.0;
         moli->np2[i]	=	0;
      }
   }

   // Calculation
   
#ifdef	myMPI
   for (moli=mol+molid[MPI_myid]; moli<mol+molid[MPI_myid+1]; moli++) {
#else
   for (moli=mol; moli<mol+NMOLS; moli++) {
#endif
      if ( (system=moli->box)>=0 ) {
         for (i=1; i<moli->nsites-1; i++) {			// exclude 2 end sites
            vi	=	V_Subtr(moli->p+i+1, moli->p+i-1);
            vi	=	MapInBox2(&vi, PBC, system);

            for (molj=moli; molj<mol+NMOLS; molj++) {
               if ( molj->box == system ) {
                  for (j=1; j<molj->nsites-1; j++) {
                     if ( (moli!=molj) ? 1 : j>i ) {

			vj		=	V_Subtr(molj->p+j+1, molj->p+j-1);			
            		vj		=	MapInBox2(&vj, PBC, system);
                        idotj		=	V_Dot(&vi, &vj);
			idotj		*=	idotj;
			idotj		/=	V_Dot(&vi, &vi) * V_Dot(&vj, &vj);
                        P2[system]	+=	idotj;
			n[system]	++;
                        
                        r2	=	DistSQ(moli->p[i], molj->p[j], system);
                        if (r2 < Rp2*type[0].SIGMA*type[0].SIGMA) {	// calc. local p2
			   moli->p2[i]	+=	idotj;
			   molj->p2[j]	+=	idotj;
			   moli->np2[i]	++;
			   molj->np2[j]	++;
			}
            }  }  }  }
         }
      }
   }

   //------Transfer data to process #0 and process data------//
#ifdef	myMPI
   if (MPI_myid >0) {
      // Send data to process #0
      MPI_Send(n, NSYSTEMS, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(P2, NSYSTEMS, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

      k	=	0;
      for (moli=mol; moli<mol+NMOLS; moli++) {
	 for (i=0; i<moli->nsites; i++) {
	    data_np2[k]	=	moli->np2[i];
	    data_p2[k]	=	moli->p2[i];
	    k		++;
	 }
      }
      if (k!=NSITES)	printf("SampleP2All_MPI(): Error k!=NSITES\n");

      MPI_Send(data_np2, NSITES, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(data_p2, NSITES, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
   }
   else {
      // Receive data from other processes
      for (proc=1; proc<MPI_numprocs; proc++) {
	 MPI_Recv(data_nP2, NSYSTEMS, MPI_INT, proc, 0, MPI_COMM_WORLD, &status);
	 MPI_Recv(data_P2, NSYSTEMS, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &status);

	 for (i=0; i<NSYSTEMS; i++) {
	    P2[i]	+=	data_P2[i];
	    n[i]	+=	data_nP2[i];
	 }

	 MPI_Recv(data_np2, NSITES, MPI_INT, proc, 0, MPI_COMM_WORLD, &status);
	 MPI_Recv(data_p2, NSITES, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &status);

	 k	=	0;
	 for (moli=mol; moli<mol+NMOLS; moli++) {
	    for (i=0; i<moli->nsites; i++) {
	       moli->np2[i]	+=	data_np2[k];
	       moli->p2[i]	+=	data_p2[k];
	       //printf("%f %d\n", data_p2[i], data_np2[i]);
	       
	       k		++;
	    }
	 }
      }
   }
#endif

   // Calc. the orientation along z-axis P2z
   for (moli=mol; moli<mol+NMOLS; moli++)
      if ( (system=moli->box)>= 0)
         for (i=1; i<moli->nsites-1; i++) {
            vi		=	V_Subtr(moli->p+i+1, moli->p+i-1);
            vi		=	MapInBox2(&vi, PBC, system);
            P2z[system]	+=	(vi.z * vi.z) / V_Dot(&vi, &vi); 
            nz[system]	++;
         }
            
   // Normalization

   for (i=0; i<NSYSTEMS; i++) {
      if (n[i]) {
         P2[i]	/=	n[i];
         P2[i]	*=	1.5;
         P2[i]	-=	0.5;
      }
      else	P2[i]	=	0;
      if (nz[i]) {
         P2z[i]	/=	nz[i];
         P2z[i]	*=	1.5;
         P2z[i]	-=	0.5;
      }
      else	P2z[i]	=	0;
   }

   // Calculate modified P2M using local p2(i)

   for (i=0; i<NSYSTEMS; i++) {
      n[i]	=	0;
   }
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=1; i<moli->nsites-1; i++) {
         if (moli->np2[i]) {
            moli->p2[i]	/=	moli->np2[i];
            moli->p2[i]	*=	1.5;
	    moli->p2[i]	-=	0.5;
         }
	 else	moli->p2[i]	=	0;

//printf("%d->p2[%d]=%f\n", moli-mol, i, moli->p2[i]);

         P2M[moli->box]	+=	moli->p2[i];
   	 n[moli->box]	++;
      }
      moli->p2[0]		=	moli->p2[1];
      moli->p2[moli->nsites-1]	=	moli->p2[moli->nsites-2];
   }
   for (i=0; i<NSYSTEMS; i++) {
      P2M[i]	/=	n[i];
   }
   return;
}

//////////////////////////////////////////////////////////
/* Sample number of connections based on chord vectors	*/
//////////////////////////////////////////////////////////
void SampleConnection()	
{
   molstruct	*moli, *molj;
   long		i, j, system, n[MAXNSYSTEMS], nz[MAXNSYSTEMS];
   vector	vi, vj;
   double	idotj, r2;
   static long		nconnect[1000];		// maximum 999 connections for one site
   static long		init=1;
#ifdef CELL_LIST
   static cellstruct	*cellm, *celli;
#endif
   if (init) {
      for (i=0; i<1000; i++)
          nconnect[i]	=	0;
      init	=	0;
   }

   for (i=0; i<NSYSTEMS; i++) {		// global orientational order
      P2[i]	=	0.0;
      P2z[i]	=	0.0;
      n[i]	=	0;
      nz[i]	=	0;
   }

   for (moli=mol; moli<mol+NMOLS; moli++)	// local orientational order
      for (i=0; i<moli->nsites; i++) {		// include end beads
         moli->p2[i]	=	0.0;
         moli->np2[i]	=	0;
      }

   for (moli=mol; moli<mol+NMOLS; moli++)
      if ( (system=moli->box)>=0 )
         for (i=0; i<moli->nsites; i++) {
            if (i==0)
               vi	=	V_Subtr(moli->p+2, moli->p);
            else if (i==moli->nsites-1)
               vi	=	V_Subtr(moli->p+moli->nsites-1, moli->p+moli->nsites-3);
            else
	       vi	=	V_Subtr(moli->p+i+1, moli->p+i-1);

            for (molj=moli; molj<mol+NMOLS; molj++)
               if ( molj->box == system )
                  for (j=0; j<molj->nsites; j++)
                     if ( (moli!=molj) ? 0 : j>i ) {
                        if (j==0)
                           vj	=	V_Subtr(molj->p+2, molj->p);
                        else if (j==molj->nsites-1)
                           vj	=	V_Subtr(molj->p+molj->nsites-1, molj->p+molj->nsites-3);
                        else
			   vj	=	V_Subtr(molj->p+j+1, molj->p+j-1);

                        idotj		=	V_Dot(&vi, &vj);
			idotj		*=	idotj;
			idotj		/=	V_Dot(&vi, &vi) * V_Dot(&vj, &vj);	//cos^2
			
                        if (idotj >= 0.933) {	// 0.933=cos^2 15degrees
//                        if (idotj >= 0.75) {	// 0.933=cos^2 30degrees
                        
                        if ((r2=DistSQ(moli->p[i], molj->p[j], system)) < Rconn2*type[0].SIGMA*type[0].SIGMA) {	// calc. distance between i and j 
			   moli->np2[i]	++;
			   molj->np2[j]	++;
			}

			}
                     }
         }
   for (moli=mol; moli<mol+NMOLS; moli++)
      for (i=0; i<moli->nsites; i++)
         nconnect[moli->np2[i]]	++;

   for (i=0; i<200; i++)
      printf("%ld  %ld\n", i, nconnect[i]);
/*
   // calc. the orientation along z-axis
   for (moli=mol; moli<mol+NMOLS; moli++)
      if ( (system=moli->box)>= 0)
         for (i=1; i<moli->nsites-1; i++) {
            vi		=	V_Subtr(moli->p+i+1, moli->p+i-1);
            P2z[system]	+=	(vi.z * vi.z) / V_Dot(&vi, &vi); 
            nz[system]	++;
         }
            
   // normalize
   for (i=0; i<NSYSTEMS; i++) {
      if (n[i]) {
         P2[i]	/=	n[i];
         P2[i]	*=	1.5;
         P2[i]	-=	0.5;
      }
      else	P2[i]	=	0;
      if (nz[i]) {
         P2z[i]	/=	nz[i];
         P2z[i]	*=	1.5;
         P2z[i]	-=	0.5;
      }
      else	P2z[i]	=	0;
   }

   // calc. modified P2M using local p2(i)
   for (i=0; i<NSYSTEMS; i++) {
      P2M[i]	=	0.0;
      n[i]	=	0;
   }
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=1; i<moli->nsites-1; i++) {
         if (moli->np2[i]) {
            moli->p2[i]	/=	moli->np2[i];
            moli->p2[i]	*=	1.5;
	    moli->p2[i]	-=	0.5;
         }
	 else	moli->p2[i]	=	0;

         P2M[moli->box]	+=	moli->p2[i];
   	 n[moli->box]	++;
      }
      moli->p2[0]		=	moli->p2[1];
      moli->p2[moli->nsites-1]	=	moli->p2[moli->nsites-2];
   }
   for (i=0; i<NSYSTEMS; i++)
      P2M[i]	/=	n[i];
*/
}

//////////////////////////////////////////
/* Collect local p2 distribution	*/
//////////////////////////////////////////
void Dist_p2()		
{
   molstruct	*moli;
   long		i;

   if (D_LOCALP2)
      for (moli=mol; moli<mol+NMOLS; moli++)
         //for (i=1; i<moli->nsites-1; i++)
         for (i=0; i<moli->nsites; i++)
            PutInDistribution(D_Localp2+moli->box, moli->p2[i], 1.0, 1.0); 
}


//////////////////////////////////////////////////////////
/* Sample LC orientation matrix for all sites		*/
/* Principles of condensed matter physics, page 168	*/
//////////////////////////////////////////////////////////
void SampleM_Q()
{	
   molstruct		*moli, *molj, *molk;
   long			i, j, k, n, system, sitek;
   double		r2, ri;
   vector		evalue, evector, vi, vj;
   matrix		Q;
   static long		nconnect[1000];		// maximum 999 connections for one site
   static long		init=1;
   static double	cosangle;		// cosine of critical angle
#ifdef CELL_LIST
   cellstruct		*cellj, *cellk;
#endif

   if (init) {
      critangle	=	15.0/180*M_PI;
      cosangle	=	cos(critangle);
      for (i=0; i<200; i++)
         nconnect[i]	=	0;
      init	=	0;
   }

   for (moli=mol; moli<mol+NMOLS; moli++) {
      //system=moli->box;
      system	=	0;				// as of Aug11,09, only one box
      for (i=1; i<moli->nsites-1; i++) {		// not for end sites

         moli->nconn[i]	=	0;

         for (molj=mol; molj<mol+NMOLS; molj++)		// list all chord vectors near i
            if ( molj->box == system )
               for (j=1; j<molj->nsites-1; j++)
                  if ((moli!=molj || i!=j) && (r2=DistSQ(moli->p[i], molj->p[j], system)) < Rp2*type[0].SIGMA*type[0].SIGMA) {

                     n		=	moli->nconn[i];
		     moli->connsite[i][n]	=	j;
		     moli->connmol[i][n]	=	molj;
		     moli->nconn[i]		++;
                  }

	 M_Null(&Q);

         vi	=	V_Subtr(moli->p+i+1, moli->p+i-1);	// self
         r2	=	V_Dot(&vi, &vi);
         vi	=	V_Mult(1.0/sqrt(r2), &vi);	// normalize
         Q.x.x	+=	vi.x * vi.x;
         Q.y.y	+=	vi.y * vi.y;
         Q.z.z	+=	vi.z * vi.z;
         Q.x.y	+=	vi.x * vi.y;
         Q.x.z	+=	vi.x * vi.z;
         Q.y.z	+=	vi.y * vi.z;

         for (k=0; k<moli->nconn[i]; k++) {		// neighbors
            molk	=	moli->connmol[i][k];
            sitek	=	moli->connsite[i][k];

            vi	=	V_Subtr(molk->p+sitek+1, molk->p+sitek-1);
            r2	=	V_Dot(&vi, &vi);
            vi	=	V_Mult(1.0/sqrt(r2), &vi);	// normalize
            Q.x.x	+=	vi.x * vi.x;
            Q.y.y	+=	vi.y * vi.y;
            Q.z.z	+=	vi.z * vi.z;
            Q.x.y	+=	vi.x * vi.y;
            Q.x.z	+=	vi.x * vi.z;
            Q.y.z	+=	vi.y * vi.z;
         }
         Q	=	M_Mult(1.0/(moli->nconn[i]+1), &Q);	// normalize
         Q.x.x	-=	0.333333;
         Q.y.y	-=	0.333333;
         Q.z.z	-=	0.333333;
         Q.y.x	=	Q.x.y;
         Q.z.x	=	Q.x.z;
	 Q.z.y	=	Q.y.z;

         evalue	=	M_eig(Q);				// find eigenvalue
         evector=	V_eig(Q, MAX(MAX(evalue.x, evalue.y),evalue.z));// find eigenvector
         r2	=	V_Dot(&evector, &evector);
         evector=	V_Mult(1.0/sqrt(r2), &evector);		// normalize eigenvector

         moli->vp2[i]	=	evector;		// assign value
//r2	=	V_Dot(&evector, &evector); 
//printf("%f %f %f %f ", evalue.x, evalue.y, evalue.z, evalue.x+evalue.y+evalue.z);
//printf("%f %f %f %f\n", evector.x, evector.y, evector.z, r2);
      }
      moli->vp2[0]	=	V_Add(moli->vp2+1, moli->vp2+2);	// end beads
      moli->vp2[0]	=	V_Mult(0.5, moli->vp2+0);

      moli->vp2[moli->nsites-1]	=	V_Add(moli->vp2+moli->nsites-2, moli->vp2+moli->nsites-3);
      moli->vp2[moli->nsites-1]	=	V_Mult(0.5, moli->vp2+moli->nsites-1);
   }

   // Now determine the number of connection for each site
   // Note, connectivity cutoff should be smaller than p2 cutoff !!!
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
         moli->np2[i]	=	0;			// NOTE: used to store connection #

         for (molj=mol; molj<mol+NMOLS; molj++) {
            for (j=0; j<molj->nsites; j++) {
               if (molj==moli && j==i)	break;
               r2	=	 DistSQ(moli->p[i], molj->p[j], moli->box);
               if (r2<Rconn2*type[0].SIGMA*type[0].SIGMA && fabs(V_Dot(moli->vp2+i, molj->vp2+j)) > cosangle)
                  moli->np2[i]	++;
/*
         for (k=0; k<moli->nconn[i]; k++) {		// search the neighbors
            molk	=	moli->connmol[i][k];
            sitek	=	moli->connsite[i][k];
	    r2		=	DistSQ(molk->p[sitek], moli->p[i], moli->box);
            if (r2 < Rconn2 && fabs(V_Dot(moli->vp2+i, molk->vp2+sitek)) > cosangle)
               moli->np2[i]	++;
         }            
*/	}}
         nconnect[moli->np2[i]]	++;
      }	    
   }
   printf("population distribution of number of connection\n");
   for (i=0; i<200; i++)	printf("%ld\n", nconnect[i]);
}


//////////////////////////////////////////////////////////
/* sample the mean squared of end-to-end distance	*/
//////////////////////////////////////////////////////////
void SampleN2NSQ()	
{
   molstruct	*moli;
   double	R2 = 0.0;
   long		system;

   if (D_DENSITY)			// temporary
      for (system=0; system<NSYSTEMS; system++) 
         for (moli=mol; moli<mol+NMOLS; moli++)
            if (moli->box == system)
               PutInDistribution(D_Density+system, R2_n2n(moli), 1.0, 1.0);
}


//////////////////////////////////////////////////////////
/* sample energy					*/
//////////////////////////////////////////////////////////
void SampleEnergy()
{
   long		system;
   if (D_ENERGY)
      for (system=0; system<NSYSTEMS; system++)
         PutInDistribution(D_Energy+system, v[system].tot, 1.0, 1.0);
}

//////////////////////////////////////////////////////////
/* test function for 2D distribution			*/
//////////////////////////////////////////////////////////
void Sample2D()			
{
   long		system;
   double	dyn[2], y=1.0, weight=1.0;

   if (D_ENERGY)
      for (system=0; system<NSYSTEMS; system++) {
         dyn[0]	=	v[system].tot;
         dyn[1]	=	BOX[system].pres;
         D_Submit(D_Energy+system, dyn, &y, &weight);
      }
}

//////////////////////////////////////////////////////////
/* Print out distributions				*/
//////////////////////////////////////////////////////////
void S_PrintDistribution(diststruct **d)
{
   double		L = 1.0;				// change data units
   long			system, n = d-D_Distributions;

   if (*d) {
      for (system=0; system<NSYSTEMS; system++) {
         printf("%s\n", (*d+system)->header);
         D_Print(*d+system, 0, &L);
      }
   }
}

void S_PrintAll()
{
   long		i;
   for (i=0; i<D_NDIST; i++)
      S_PrintDistribution(D_Distributions+i);
}

void InitAllDistributions(diststruct **d, double binsize)	// add new distribution only to the end
{								// increase D_NDIST when adding
   long		i;
   double	dn[1], datom[1], denergy[1], ddensity[1], dpressure[1], dx[1];
   double	dp2[1];
   double	dtorsa[1], dbonda[1];

   *dn		=	1.0;
   *datom	=	binsize;		// D_BINSIZE
/*   denergy[0]	=	2.000;
   denergy[1]	=	2.000;
*/
   *denergy	=	0.1;
   *ddensity	=	4.00;
   *dpressure	=	0.002;
   *dx		=	0.1;
   *dtorsa	=	M_PI / 45;
   *dbonda	=	M_PI / 180;
   *dp2		=	0.01;

   for (i=0; i<D_NDIST; i++) {
      switch(i) {
         case 0: InitDist(D_DENSITY, d+i, "Density", 1, ddensity);	break;	// 1D distribution
	 case 1: InitDist(D_ENERGY, d+i, "Energy", 1, denergy);		break;	// 2D distribution
         case 2: InitDist(D_PRESSURE, d+i, "Pressure", 1, dpressure);	break;
         case 3: InitDist(D_DRIFT, d+i, "Drift", 1, dx);		break;
         case 4: InitDist(D_TORSION, d+i, "torsion angle", 1, dtorsa);	break;
         case 5: InitDist(D_BONDA, d+i, "bond angle", 1, dbonda);	break;
         case 6: InitDist(D_BONDL, d+i, "bond length", 1, dpressure);	break;
         case 7: InitDist(D_RADIAL, d+i, "radial dist", 1, dx);		break;
         case 8: InitDist(D_LOCALP2, d+i, "local p2", 1, dp2);		break;
         case 9: InitDist(D_XTALSIZE, d+i, "Xtal size", 1, dn);		break;
	 default:	break;
      }
   }
}

void LinkDistributions(diststruct **d)
{
   D_Density	=	*d++;			// i.e., D_Density = *d; d++;
   D_Energy	=	*d++;
   D_Pressure	=	*d++;
   D_Drift	=	*d++;
   D_Torsion	=	*d++;
   D_Bonda	=	*d++;
   D_Bondl	=	*d++;
   D_Radial	=	*d++;
   D_Localp2	=	*d++;
   D_Xtalsize	=	*d++;
}

void InitSample(char *argv[])
{
   static long		init = 1;
   char			name[256];
   long			n;
   double		binsize, d[1];

   if (init) {
      for (n=0; n<D_NDIST; n++) 
          D_Distributions[n]	=	NULL;
      init	=	0;
   }
   n		=	0;
   binsize	=	-1.0;

   if (binsize<=0.0) 
      binsize	=	D_BINSIZE;
   InitAllDistributions(D_Distributions, binsize);
   LinkDistributions(D_Distributions);
}

void ReinitSample(diststruct **d)	// only use when setting binsize afterwards
{
   InitAllDistributions(d, D_BINSIZE);
   if (d==D_Distributions)
      LinkDistributions(d);
}

void Sample_Histogram()
{
/*
   if (dynvar==1) {		//update the histogram
      p[MAXSIZE]	++;
      if (MAXSIZE == sizeright && mod(samplecycle, 2)==0)		samplecycle ++;
      else if (MAXSIZE == sizeleft && mod(samplecycle, 2)==1)	samplecycle ++;
   }
   else if (dynvar==2) {
      pQ[Qstatefinder(Ql)]	++;
      if (Qstatefinder(Ql) == Qrightid && mod(samplecycle, 2)==0)		samplecycle ++;
      else if (Qstatefinder(Ql) == Qleftid && mod(samplecycle, 2)==1)	samplecycle ++;
   }

   if (mod(counter*5, NCYCLE) == 0) {
      Print_Histogram();
   }
*/
}

vector Center_of_Mass(long *list, long number)		// center of mass of the molecules in a list
{						// need to be fixed to include different mass
   long         i, j;
   vector       p;

   V_Null(&p);

   for (i=0; i<number; i++) {
      for (j=0; j<Nsites; j++) {
         p.x       +=      mol[list[i]].p[j].x;
         p.y       +=      mol[list[i]].p[j].y;
         p.z       +=      mol[list[i]].p[j].z;
      }
   }
   V_Mult(1.0/(number*Nsites), &p);
   return       p;
}


matrix Gyration_Tensor(long *list, long number)
{
   long         i, j;
   matrix       tensor;
   double	x, y, z;
   vector	CoM;

   M_Null(&tensor);

   for (i=0; i<number; i++) {
      for (j=0; j<Nsites; j++) {

	 x	=	mol[list[i]].p[j].x;
	 y	=	mol[list[i]].p[j].y;
	 z	=	mol[list[i]].p[j].z;

         tensor.x.x        +=      x * x;
         tensor.x.y        +=      x * y;
         tensor.x.z        +=      x * z;
         tensor.y.y        +=      y * y;
         tensor.y.z        +=      y * z;
         tensor.z.z        +=      z * z;
      }
   }

   M_Mult(1.0/(number*Nsites), &tensor);

   CoM	=	Center_of_Mass(list, number);

   tensor.x.x	-=	CoM.x * CoM.x;		// substract the center of mass contribution
   tensor.x.y	-=	CoM.x * CoM.y;
   tensor.x.z	-=	CoM.x * CoM.z;
   tensor.y.y	-=	CoM.y * CoM.y;
   tensor.y.z	-=	CoM.y * CoM.z;
   tensor.z.z	-=	CoM.z * CoM.z;

   tensor.y.x   =       tensor.x.y;             //it is a symmetric matrix by definition
   tensor.z.x   =       tensor.x.z;
   tensor.z.y   =       tensor.y.z;

   return       tensor;
}


void Init_Sample()				// initialize sampling
{
   long		i;

   for (i=0; i<NBOX; i++) {
      chp[i]	=	0.0;
      cchp[i]	=	0.0;
      ichp[i]	=	0;
   }
   return;
}


void Sample_All()				// instant sampling during simulation
{
   if (V_VIRIAL)
      Sample_Pressure();
   SampleDrift();
   //SampleP2();
}


void Sample_Done()
{
   long		i;
/*
   for (i=0; i<NBOX; i++) {
      fprintf(foutput, "\n\tBOX[%ld]:\n\n", i);
      fprintf(foutput, "\tChemical potential:\t%f\n", -log(cchp[i]/ichp[i]) * kT);
      fprintf(foutput, "\tTotal sampling:\t%ld\n", ichp[i]);
   }
*/
   fflush(foutput);
   return;
}

/*==============================================================*/
/*	Realtime composition gradient in the system, 1/26/17	*/
/*==============================================================*/
void composition_profile(long timestep, char direction, vector vnucl, double* local_comp, double* local_grad, double presski)
{
#define shifting 0			// 1: do profile shifting; 0: do not
#define roll   4			// perform average every roll snapshots processed, except for t=0
   static int init = 1;
   static FILE *fPtr, *fPtr2, *fPtr3;
   static float type1[2048]; 	// number of type1 atoms in each bin
   static float total[2048];	// number of total atoms in each bin
   static float ch_type1[2048];	// for chosen atoms only
   static float ch_total[2048];	// for chosen atoms only
   static float stress[2048][6];
   // cna value range 1-5, 1-fcc 2-hcp 3-bcc, 4-icosohedral 5-unknown
   static float cna[2048][6];  
   static float cna2[2048][6];
   static float cna3[2048][6];

   static int count =0;
   static int index;		// how many averaging performed

   int i, j, k, n, bin;
   molstruct *moli;
   int system = 0;
   int itmp;
   int itmp3 = 0;
   int itmp4 = 0;
   int itmp5 = 0;
   int itmp6 = 0;
   
   vector vtmp0, vtmp;
   double coord, Ldim, dr;
   double vol, binvol;
   double presske; 		// kinetic contribution to pressure
   vol = BOX[system].lx * BOX[system].ly * BOX[system].lz;

   int 		nbins;    				// number of bins
   double 	binwidth;  				// binwidth in Angstrom
   binwidth = 2.5; 			    		// approximate binwidth
   //binwidth = 10.0;

   switch(direction) {
      case 'x':
		Ldim = BOX[system].lx;
		break;
      case 'y': 
		Ldim = BOX[system].ly;
      		break;
      case 'z': 
		Ldim = BOX[system].lz;
      		break;
      default:  
                printf("Error: No direction is defined in composition_profile()\n");
		exit(1);
   }
   nbins = (int) (Ldim / binwidth);
   binvol = vol/nbins;
   binwidth = Ldim / nbins;  				// slightly adjust binwidth

   if (nbins > 2048) {
      printf("composition_profile error: nbins > 2048\n");
      exit(1);
   }

   double 	samplewidth;				// width when collect particle type info
   samplewidth = binwidth; 				// samplewidth might be different than binwidth
   //samplewidth = 3*binwidth;
   //samplewidth = 4.0; 	

   if (init) {
      init = 0;
      fPtr = fopen("composition_profile", "w");
      fPtr2 = fopen("composition_gradient", "w");
      fPtr3 = fopen("2dprofile", "w");
      index = 0;
   }

   if (count==0 || (count-1) % roll == 0) {  		// for average over snapshots
   //if ( count % roll == 0) {  			// for average over snapshots
      for (i=0; i<2048; i++) {
	 type1[i] = 0;
	 total[i] = 0;
	 ch_type1[i] = 0;
	 ch_total[i] = 0;
	 for (j=0; j<6; j++) {
	    stress[i][j] = 0.0;
	 }
	 for (j=1; j<6; j++) {
	    cna[i][j] 	= 0.0;
	    cna2[i][j]	= 0.0;
	    cna3[i][j]	= 0.0;
	 }
      }
   }

   for (moli=mol; moli < mol + NMOLS; moli++) {
      for (i=0; i < moli->nsites; i++) {

	 // Move all the atoms into the original box

	 // -- step 1: move the origin (0,0,0) to the center of the box
	 moli->p[i].x -= (BOX[system].xlo + 0.5 * BOX[system].lx);
	 moli->p[i].y -= (BOX[system].ylo + 0.5 * BOX[system].ly);
	 moli->p[i].z -= (BOX[system].zlo + 0.5 * BOX[system].lz);
	 
	 // -- step 2: map in atoms
	 vtmp = MapInBox2(moli->p+i, PBC, system);

	 // Adjust coordinate as if the origin (0,0,0) is at the  bottom-left corner of the box
	 //
	 vtmp.x += BOX[system].lx * 0.5;		// x coordinates start at 0
	 vtmp.y += BOX[system].ly * 0.5;		// y coordinates start at 0
	 vtmp.z += BOX[system].lz * 0.5;		// z coordinates start at 0

	 // keep origin(0,0,0) at the bottom-left corner of the box
	 moli->p[i] = vtmp;				

	 for (n=0; n<nbins; n++) {  			// loop over n to find the right bin
	    switch (direction) {
	       case 'x': 
	                 coord = vtmp.x;
	                 break;
	       case 'y': 
	                 coord = vtmp.y;
	                 break;
	       case 'z': 
	                 coord = vtmp.z;
	                 break;
	       default:  printf("Error: No direction is defined in composition_profile()\n");
	                 exit(1);
	    }

	    if ( fabs(coord  - (n+0.5)*binwidth) <= 0.5*samplewidth ) {
	       bin = n;
	    }
	    else if ( (0+0.5)*binwidth - (coord - Ldim) <= 0.5*samplewidth ) {
	       bin = 0;
	    }
	    else if ( coord + Ldim - (nbins-0.5)*binwidth <= 0.5*samplewidth ) {
	       bin = nbins-1;
	    }
	    else {   		// found the right bin already
	       continue;
	    }

	    //bin = (int) (coord/binwidth);
	    if (bin<0 || bin>=nbins) {
	       printf("bin not defined, bin = %d %f %f\n", bin, vtmp.x, BOX[system].lx);
	       exit(1);
	    }

	    total[bin] ++;
	    if (moli->cna2[i]==3)
	       ch_total[bin] ++;

	    //if (moli->type[i] == 0) 			// LAMMPS type index starts with 1, rather than 0
	    if (moli->type[i] %2 == 0) 			// type 0, 2, 4, ...
	       type1[bin] ++;

	    if (moli->cna2[i]==3)
	       if (moli->type[i] %2 == 0) 			// type 0, 2, 4, ...
		  ch_type1[bin] ++;

	    cna[bin][moli->cna[i]]	+= 1;
	    cna2[bin][moli->cna2[i]]	+= 1;
	    cna3[bin][moli->cna3[i]]	+= 1;

	    for (j=0; j<6; j++) {
	       stress[bin][j] += moli->stress[i][j];
	    }
	 }
      }
   }

   float grad[21];			// variable to record gradient as a function of composition
   int   gradcount[21];

   if (count % roll == 0) {  		// do averaging and output
   //if ((count+1) % roll == 0) {  		// for average

      //---- loop over bins to do average
      for (i=0; i<nbins; i++) {
	 if (count!=0) {  		// not the first configuration
	    type1[i] /= roll;
	    total[i] /= roll;
	    ch_type1[i] /= roll;
	    ch_total[i] /= roll;
	    for (j=0; j<6; j++) {
	       stress[i][j] /= roll;
	       stress[i][j] /= binvol*1e4; 	// virial convert to pressure, then unit from bar to GPa
	    }
	    for (j=1; j<6; j++) {
	       cna[i][j]   /= roll;
	       cna2[i][j]  /= roll;
	       cna3[i][j]  /= roll;
	    }
	 }
      }
      
      //---- find the middle of Ni block for shifting output, such that Ni blocks on the sides
      if (shifting) {
	 itmp = 0;
	 n = 0;
	 for (i=0; i<nbins; i++) {
	    if ( 1-type1[i]/total[i] < 1e-6) {
	       itmp += i;
	       n ++;
	    }
	 }
	 itmp = (int) ((float)itmp/n + 0.5);
      }
      else {
         itmp = 0; 	// no shifting
      }

      //---- Output profiles
      fprintf(fPtr, "\'timestep=%d\' \'index=%d\' \'time=%-6.3f(ns)\'\n", timestep, index, timestep/2.5e5);

      for (j=0; j<nbins; j++) {
	 i = (j+itmp)%nbins;

	 // plus presski because post-analysis kinetic energy not calculated, i.e., 0K.
	 // so need to add it back. (pressure = virial stress + kinetic contribution)

	 presske =  presski * total[i] / binvol / 1e4; 	// each bin has different total
	 //presske = 0;

	 fprintf(fPtr, "%d %f %5.2f %5.2f %5.4f %5.4f %5.4f %5.4f %8.3f %8.3f %8.3f %5.4f %8.3f %8.3f %5.4f\n", 
	         j, (j+0.5)*binwidth, type1[i], total[i], type1[i]/total[i], 
		 cna[i][2]/total[i], cna2[i][3]/total[i], cna3[i][3]/total[i],
		 -stress[i][0]+presske, -stress[i][1]+presske, -stress[i][2]+presske,
		 ch_type1[i]/ch_total[i], ch_type1[i], ch_total[i], (type1[i]-ch_type1[i])/(total[i]-ch_total[i]));
      }
      fprintf(fPtr, "\n\n");
      fflush(fPtr);

      // Output 2D profile of selected slice (bin)
      //
      //float _2dtype1[16][16];
      //float _2dtotal[16][16];
      //float _2dB2[16][16];
      float _2dtype1[8][8];
      float _2dtotal[8][8];
      float _2dB2[8][8];
      int bin1, bin2;
      int Nbins12 = 8;

      for (i=0; i<Nbins12; i++) {
	 for (j=0; j<Nbins12; j++) {
	    _2dtype1[i][j] = 0;
	    _2dtotal[i][j] = 0;
	    _2dB2[i][j] = 0;
	 }
      }

      //---- find the slice in the 1st dimension
      j = 0;
      for (i=1; i<nbins; i++) {
         if (cna2[i][3]/total[i] > cna2[j][3]/total[j]) {
	    j=i;
	 }
      }
      fprintf(fPtr3, "timestep = %d, itmp=%d, j=%d, %5.4f, %5.4f, %5.4f, %5.4f\n", 
                                 timestep, itmp, j, cna2[j][3]/total[j], type1[j]/total[j],
				 ch_type1[j]/ch_total[j], (type1[j]-ch_type1[j])/(total[j]-ch_total[j]));

      for (moli=mol; moli < mol + NMOLS; moli++) {
	 for (i=0; i < moli->nsites; i++) {
	    // remember by default the origin (0,0,0) is at the bottom-left corner of the box
	    vtmp.x = moli->p[i].x;		// x coordinates start at 0
	    vtmp.y = moli->p[i].y;		// y coordinates start at 0
	    vtmp.z = moli->p[i].z;		// z coordinates start at 0

            if ( j == (int) (nbins * vtmp.x/BOX[system].lx) ) {   // in the chosen bin
	       bin1 = (int) (Nbins12 * vtmp.y/BOX[system].ly);
	       bin2 = (int) (Nbins12 * vtmp.z/BOX[system].lz);

	       _2dtotal[bin1][bin2] ++;

	       if (moli->type[i] %2 == 0) 			// type 0, 2, 4, ...
		  _2dtype1[bin1][bin2] ++;

	       if (moli->cna2[i] == 3) 
		  _2dB2[bin1][bin2] ++;
	    }
	 }
      }

      for (j=0; j<Nbins12; j++) {
	 for (k=0; k<Nbins12; k++) {
	    fprintf(fPtr3, "%5.2f", _2dtype1[j][k]/_2dtotal[j][k]);
	 }
	 fprintf(fPtr3, "\n");
      }
      fprintf(fPtr3, "\n");

      for (j=0; j<Nbins12; j++) {
	 for (k=0; k<Nbins12; k++) {
	    fprintf(fPtr3, "%5.2f", _2dB2[j][k]/_2dtotal[j][k]);
	 }
	 fprintf(fPtr3, "\n");
      }
      fprintf(fPtr3, "\n");
      fflush(fPtr3);

      //---- Calculate composition gradient as a function of composition
      
      for (i=0; i<21; i++) {
	 grad[i] = 0.0;
	 gradcount[i] = 0;
      }

      for (i=0; i<nbins; i++) {
	 itmp3 = ((i-1)>=0 ? (i-1) : (i-1+nbins));
	 itmp4 = ((i+1)<nbins ? (i+1) : (i+1-nbins));

	 if (itmp3<0 || itmp3>=nbins) {
	    printf("itmp3 error\n");
	    exit(1);
	 }
	 if (itmp4<0 || itmp4>=nbins) {
	    printf("itmp4 error\n");
	    exit(1);
	 }

	 if (total[i] >0) {
	    j = (int) ((1.0*type1[i]/total[i])/(1.0/20));
	    if (j<0 || j>=21) {
	       printf("j error j = %d i = %d type1[i] = %d total[i] = %d\n", j, i, type1[i], total[i]);
	       exit(1);
	    }
	 
	    if (total[itmp3] >0 && total[itmp4]>0) {
	       grad[j] += fabs( 1.0*type1[itmp3]/total[itmp3]-1.0*type1[itmp4]/total[itmp4] );
	       gradcount[j] ++;
	    }
	 }
      }

      fprintf(fPtr2, "\'timestep=%d\' \'index=%d\' \'time=%04g(ns)\'\n", timestep, index, timestep/2.5e5);

      for (i=0; i<21; i++) {
	 grad[i] /= gradcount[i];
	 fprintf(fPtr2, "%f\t%f\t%d\n", i*1.0/20, grad[i], gradcount[i]);
      }

      fprintf(fPtr2, "\n\n");
      fflush(fPtr2);

      index ++; 
   }
   count ++;	// increment count
   
   //======== Calculate local composition and gradient at vnucl
   //
   int total_left = 0;
   int total_right = 0;
   int total_middle = 0;
   int type1_left = 0;
   int type1_right = 0;
   int type1_middle = 0;

   // vnucl is the coordinate of the center of the nucleus
   vnucl.x -= (BOX[system].xlo + 0.5 * BOX[system].lx);   	// move (0,0,0) to the center of the box
   vnucl.y -= (BOX[system].ylo + 0.5 * BOX[system].ly);
   vnucl.z -= (BOX[system].zlo + 0.5 * BOX[system].lz);

   vtmp0 = MapInBox2(&vnucl, PBC, system);
   vtmp0.x += BOX[system].lx * 0.5;		// x coordinates start at 0
   vtmp0.y += BOX[system].ly * 0.5;		// y coordinates start at 0
   vtmp0.z += BOX[system].lz * 0.5;		// z coordinates start at 0

   for (moli=mol; moli < mol + NMOLS; moli++) {
      for (i=0; i < moli->nsites; i++) {

	 //vtmp = MapInBox2(moli->p+i, PBC, system);
	 vtmp = moli->p[i];

	 switch(direction) {
	    case 'x': dr = vtmp0.x - vtmp.x; 	break;
	    case 'y': dr = vtmp0.y - vtmp.y; 	break;
	    case 'z': dr = vtmp0.z - vtmp.z; 	break;
	    default : break;
	 }

         if ( dr > 0.5*samplewidth && dr < 1.5*samplewidth ) {
	    total_left ++;
	    if (moli->type[i] == 0)
	       type1_left ++;
	 }
	 else if ( dr < -0.5*samplewidth && dr > -1.5*samplewidth ) {
	    total_right ++;
	    if (moli->type[i] == 0)
	       type1_right ++;
	 }
	 else if ( fabs(dr) < 0.5*samplewidth) {
	    total_middle ++;
	    if (moli->type[i] == 0)
	       type1_middle ++;
	 }
      }
   }
   *local_comp = 1.0 * type1_middle/total_middle;				// local composition
   *local_grad = (1.0*type1_left/total_left - 1.0*type1_right/total_right);	// local gradient

   return;
}

/*======================================*/
/*	Realtime 2D profile, 5/14/18	*/
/*======================================*/
void profile_2D(long timestep, char* direction1, char* direction2, double presski)
{
#define nbins1 40
#define nbins2 20
#define roll   1 			// perform average every roll snapshots
   static int init = 1;
   static FILE *fPtr;
   static double type1[nbins1][nbins2]; // number of type1 atoms in each bin
   static double total[nbins1][nbins2];	// number of total atoms in each bin
   static float pxx[nbins1][nbins2];
   static float pyy[nbins1][nbins2];
   static float pzz[nbins1][nbins2];
   static float pxy[nbins1][nbins2];
   static float pxz[nbins1][nbins2];
   static float pyz[nbins1][nbins2];
   static float cna[nbins1][nbins2];
   static float pe[nbins1][nbins2];

   static int count =0;
   static int index;			// how many averaging performed

   int i, j, n;
   molstruct *moli;
   int system = 0;
   double squaredroot2 = sqrt(2);   

   vector vtmp0, vtmp;
   double coord, dr;
   double vol;
   vol = BOX[system].lx * BOX[system].ly * BOX[system].lz;
   double binvol;
   binvol = vol/nbins1/nbins2;
   double presske;

   int dir1, dir2;
   int bin1, bin2;
   double Ldim1, Ldim2, binwidth1, binwidth2;

   if (strcmp(direction1, "x")==0)
      dir1	=	1;
   else if (strcmp(direction1, "y")==0)
      dir1	=	2;
   else if (strcmp(direction1, "z")==0)
      dir1	=	3;
   else if (strcmp(direction1, "xy")==0)
      dir1	=	4;
   else if (strcmp(direction1, "xz")==0)
      dir1	=	5;
   else if (strcmp(direction1, "yz")==0)
      dir1	=	6;
   else {
      printf("direction1 error!\n");
      exit(1);
   }

   if (strcmp(direction2, "x")==0)
      dir2	=	1;
   else if (strcmp(direction2, "y")==0)
      dir2	=	2;
   else if (strcmp(direction2, "z")==0)
      dir2	=	3;
   else if (strcmp(direction2, "xy")==0)
      dir2	=	4;
   else if (strcmp(direction2, "xz")==0)
      dir2	=	5;
   else if (strcmp(direction2, "yz")==0)
      dir2	=	6;
   else {
      printf("direction2 error!\n");
      exit(1);
   }

   switch(dir1) {
      case 1:   Ldim1 = BOX[system].lx;
		break;
      case 2:   Ldim1 = BOX[system].ly;
      		break;
      case 3:   Ldim1 = BOX[system].lz;
      		break;
      case 4:   Ldim1 = MAX(BOX[system].lx, BOX[system].ly) * squaredroot2;
      		break;
      case 5:   
      		break;
      case 6:   
      		break;
      default:  
                printf("Error: direction1 in profile_2D()\n");
		exit(1);
   }

   switch(dir2) {
      case 1:   Ldim2 = BOX[system].lx;
		break;
      case 2:   Ldim2 = BOX[system].ly;
      		break;
      case 3:   Ldim2 = BOX[system].lz;
      		break;
      case 4:   Ldim2 = MAX(BOX[system].lx, BOX[system].ly) * squaredroot2;
      		break;
      case 5:   
      		break;
      case 6:   
      		break;
      default:  
                printf("Error: direction2 in profile_2D()\n");
		exit(1);
   }

   binwidth1 = Ldim1/nbins1;
   binwidth2 = Ldim2/nbins2;

   if (init) {
      init = 0;
      fPtr = fopen("profile_2D", "w");
      index = 0;
   }

   if (count==0 || (count-1) % roll == 0) {  		// for average over snapshots
   //if ( count % roll == 0) {  			// for average over snapshots
      printf("initialize type1 and total\n");
      for (i=0; i<nbins1; i++) {
	 for (j=0; j<nbins2; j++) {
	    type1[i][j] = 0;
	    total[i][j] = 0;
	    pxx[i][j] = 0;
	    pyy[i][j] = 0;
	    pzz[i][j] = 0;
	    pxy[i][j] = 0;
	    pxz[i][j] = 0;
	    pyz[i][j] = 0;
	    cna[i][j] = 0;
	    pe[i][j] = 0;
	 }
      }
   }

   for (moli=mol; moli < mol + NMOLS; moli++) {
      for (i=0; i < moli->nsites; i++) {

	 // Move all the atoms into the original box

	 // -- step 1: move the origin (0,0,0) to the center of the box
	 
	 //if (moli==mol && i==0)
	 //   printf("atom #0: %f %f %f\n", moli->p[i].x, moli->p[i].y, moli->p[i].z); 
	 moli->p[i].x -= (BOX[system].xlo + 0.5 * BOX[system].lx);
	 moli->p[i].y -= (BOX[system].ylo + 0.5 * BOX[system].ly);
	 moli->p[i].z -= (BOX[system].zlo + 0.5 * BOX[system].lz);
	 //if (moli==mol && i==0)
	 //   printf("atom #0: %f %f %f\n", moli->p[i].x, moli->p[i].y, moli->p[i].z); 
	 
	 // -- step 2: map in atoms
	 vtmp = MapInBox2(moli->p+i, PBC, system);
	 //if (moli==mol && i==0)
	 //   printf("atom #0: %f %f %f\n", moli->p[i].x, moli->p[i].y, moli->p[i].z); 
	 //vtmp = moli->p[i];
	 //if (moli==mol && i==0)
	 //   printf("atom #0: %f %f %f\n", vtmp.x, vtmp.y, vtmp.z); 

	 // Adjust coordinate as if the original (0,0,0 is at the  bottom-left corner of the box
	 //
	 
	 vtmp.x += BOX[system].lx * 0.5;		// x coordinates start at 0
	 vtmp.y += BOX[system].ly * 0.5;		// y coordinates start at 0
	 vtmp.z += BOX[system].lz * 0.5;		// z coordinates start at 0

	 moli->p[i] = vtmp;

	 //if (moli==mol && i==0)
	 //   printf("atom #0: %f %f %f\n", vtmp.x, vtmp.y, vtmp.z); 

	 switch (dir1) {
	    case 1: coord = vtmp.x;
		      break;
	    case 2: coord = vtmp.y;
		      break;
	    case 3: coord = vtmp.z;
		      break;
	    case 4: coord = (vtmp.x + vtmp.y) * squaredroot2 * 0.5;
		      break;
	    case 5: 
		      break;
	    case 6: 
		      break;
	    default:  printf("Error: No direction is defined in composition_profile()\n");
		      exit(1);
	 }
	 bin1 = (int) coord/binwidth1;

	 switch (dir2) {
	    case 1: coord = vtmp.x;
		      break;
	    case 2: coord = vtmp.y;
		      break;
	    case 3: coord = vtmp.z;
		      break;
	    case 4: coord = (vtmp.x + vtmp.y) * squaredroot2*0.5;
		      break;
	    case 5: 
		      break;
	    case 6: 
		      break;
	    default:  printf("Error: No direction is defined in composition_profile()\n");
		      exit(1);
	 }
	 bin2 = (int) coord/binwidth2;

	 if (bin1<0 || bin1>=nbins1 || bin2<0 || bin2>=nbins2) {
	    printf("bin not defined, bin = %d %d %f %f\n", bin1, bin2, vtmp.x, BOX[system].lx);
	    fflush(stdout);
	    exit(1);
	 }

	 if (fabs(vtmp.y-0.5*BOX[system].ly)<10) { 	// analyzing only one slice
	 //if (fabs(vtmp.y-vtmp.x) < 10) { 	 	// analyzing only one slice
	    total[bin1][bin2] += 1;

	    //if (moli->type[i] == 0) 			// LAMMPS type index starts with 1, rather than 0
	    //if (moli->type[i] % 2 == 0) 			// type 0, 2, 4, ...
	    if (moli->type[i] == 1) 			// type 0, 2, 4, ...
	       type1[bin1][bin2] += 1;

	    pxx[bin1][bin2] += moli->stress[i][0];
	    pyy[bin1][bin2] += moli->stress[i][1];
	    pzz[bin1][bin2] += moli->stress[i][2];
	    pxy[bin1][bin2] += moli->stress[i][3];
	    pxz[bin1][bin2] += moli->stress[i][4];
	    pyz[bin1][bin2] += moli->stress[i][5];

	    if (moli->cna[0]==1) { 			// fcc crystal
	       cna[bin1][bin2] += 1;
	    }
	    pe[bin1][bin2] += moli->pe[0];
	 }
      }
   }

   if (count % roll == 0) {  			// for average over snapshots
   //if ((count+1) % roll == 0) {  		// for average over snapshots

      // loop over bins to do average
      for (i=0; i<nbins1; i++) {
	 for (j=0; j<nbins2; j++) {
	    if (count!=0) {  		// not the first configuration
	       type1[i][j] /= roll;
	       total[i][j] /= roll;
	       pxx[i][j] /= roll;
	       pyy[i][j] /= roll;
	       pzz[i][j] /= roll;
	       pxy[i][j] /= roll;
	       pxz[i][j] /= roll;
	       pyz[i][j] /= roll;
	       cna[i][j] /= roll;
	       pe[i][j] /= roll;
	    }
	    pxx[i][j] /= binvol*1e4; 	// convert unit of pressure from bar to GPa
	    pyy[i][j] /= binvol*1e4;
	    pzz[i][j] /= binvol*1e4;
	    pxy[i][j] /= binvol*1e4;
	    pxz[i][j] /= binvol*1e4;
	    pyz[i][j] /= binvol*1e4;
	    cna[i][j] /= total[i][j];
	    pe[i][j] /= total[i][j];
	 }
      }
      presske = presski * total[i][j] / binvol / 1e4;  		// unit GPa

      //---- Output
      fprintf(fPtr, "\'timestep=%d\' \'index=%d\' \'time=%-6.3f(ns)\'\n", timestep, index, timestep/2.5e5);

      /*
      for (i=0; i<nbins1; i++) {
	 for (j=0; j<nbins2; j++) {
      */
      for (j=0; j<nbins2; j++) {
	 for (i=0; i<nbins1; i++) {
	    // plus presski because post-analysis kinetic energy not calculated, 
	    // so need to add it back.
	    // (pressure = - virial stress + kinetic contribution)
	    //
	    //fprintf(fPtr, "%d %d %f %f %f\n", i, j, type1[i][j], total[i][j], 1-type1[i][j]/total[i][j]);
	    //fprintf(fPtr, "%f %f %f\n", (i+0.5)*binwidth1, (j+0.5)*binwidth2, type1[i][j]/total[i][j]);
	    fprintf(fPtr, "%f %f %f %f %f %f %f %f %f %f\n", (i+0.5)*binwidth1, (j+0.5)*binwidth2, 
	             -pxx[i][j]+presske, -pyy[i][j]+presske, -pzz[i][j]+presske, 
		     -pxy[i][j], -pxz[i][j], -pyz[i][j], pe[i][j], cna[i][j]);
	             //(pxx[i][j]+pyy[i][j]+pzz[i][j])/total[i][j]);
	 }
      }
      /*
      for (j=0; j<nbins2; j++) {
	 for (i=0; i<nbins1; i++) {
      */
      fprintf(fPtr, "\n\n");
      fflush(fPtr);
      index ++;
   }
   count ++;	// increment count

   /*******************/
   
   double min1, min2;
   int intface1[nbins2], intface2[nbins2];

   for (j=0; j<nbins2; j++) {
      min1 = 1e6;
      min2 = 1e6;
      for (i=0; i<nbins1/2; i++) {
	 if (fabs(pe[i][j]+4.1) < min1) {
	    min1 = fabs(pe[i][j]+4.1);
	    intface1[j] = i;
	 }
      }
      for (i=nbins1/2; i<nbins1; i++) {
	 if (fabs(pe[i][j]+4.1) < min2) {
	    min2 = fabs(pe[i][j]+4.1);
	    intface2[j] = i;
	 }
      }
   }
   printf("interface:\n");
   for (j=0; j<nbins2; j++) {
      printf("%f %f %f\n", (j+0.5)*binwidth2, (intface1[j]+0.5)*binwidth1, (intface2[j]+0.5)*binwidth1);
   }
   printf("\n");

   /*******************/

   float itmp1=0, itmp2=0;
   for (i=0; i<nbins1; i++) {
      for (j=0; j<nbins2; j++) {
	 printf("%f %f %f %f\n", (i+0.5)*binwidth1, (j+0.5)*binwidth2, total[i][j], type1[i][j]);
	 itmp1 += total[i][j];
	 itmp2 += type1[i][j];
      }
   }
   printf("PBC = %d\n", PBC);
   printf("total: %f %f\n", itmp1, itmp2);
   printf("binwidths: %f %f\n", binwidth1, binwidth2);
   return;
}

/////////////////////////////////////////////////////////////////
/* Calculate radial distribution function and structure factor */
/////////////////////////////////////////////////////////////////

void radial(char *grswitch)		// calculate and print radial distribution function
{
   long		i, j, system, igr, i_q;
   float	r, vb, nid, q, temp;
   vector	rij, comi, comj, rij0;
   molstruct	*moli, *molj;
   static long		init = 1, ngr=0;	// ngr: # of sampling
   static float		binsize[MAXNSYSTEMS], 
			Lupper[MAXNSYSTEMS], 
			*gr[MAXNSYSTEMS], 	// g(r) of all beads
			*grtot[MAXNSYSTEMS],	// for averaging
                        *grcom[MAXNSYSTEMS],	// g(r) of com of chains
                        *grcomtot[MAXNSYSTEMS],
			*grcomxy[MAXNSYSTEMS], 	// g(r) of com on xy plane
			*grcomxytot[MAXNSYSTEMS],
			*grcomz[MAXNSYSTEMS],	// g(r) of com in z direction
			*grcomztot[MAXNSYSTEMS];
   static FILE		*frdf;
   static float		ndensity[MAXNSYSTEMS];	// average number density

   // variables below for structure factor calculation, 06/25/09

   static long		nsq=0;			// nsq: # of sampling
   static float		*sq[MAXNSYSTEMS], 
			*sqtot[MAXNSYSTEMS], 
			dq, qmax;

   qmax	=	20.0;
   //dq	=	qmax/NGRBINS;
   dq	=	0.02;

   // Initialization

   if (init) {
      init	=	0;
      if (!(frdf=fopen("radial.dat", "w")))
         Exit("sample", "radial", "rdf output file failed to open");

      for (i=0; i<NSYSTEMS; i++) {
         //Lupper[i]	=	0.45*MIN(BOX[i].lx, MIN(BOX[i].ly, BOX[i].lz));
         Lupper[i]	=	4.0 * BOX[i].rc;
	 binsize[i]	=	0.05;			// system unit
         if (Lupper[i]/binsize[i] > NGRBINS) {	// make sure NGRBINS big enough
            printf("NGRBINS too small!\n");
            exit(0);
         }
         ndensity[i]		=	0.0;

         gr[i]		=	(float *) calloc(NGRBINS, sizeof(float));
         grtot[i]	=	(float *) calloc(NGRBINS, sizeof(float));
         grcom[i]	=	(float *) calloc(NGRBINS, sizeof(float));
         grcomtot[i]	=	(float *) calloc(NGRBINS, sizeof(float));
         grcomz[i]	=	(float *) calloc(NGRBINS, sizeof(float));
         grcomztot[i]	=	(float *) calloc(NGRBINS, sizeof(float));
         grcomxy[i]	=	(float *) calloc(NGRBINS, sizeof(float));
         grcomxytot[i]	=	(float *) calloc(NGRBINS, sizeof(float));
	 sq[i]		=	(float *) calloc(NGRBINS, sizeof(float));
	 sqtot[i]	=	(float *) calloc(NGRBINS, sizeof(float));
      }
   }
   // Note: halfL and binsize are calculated once for all, thus we need to 
   //       make sure that the system is in equilibrium from the beginning!

   if (!strcmp(grswitch,"sample")) {		// do sampling
      for (i=0; i<NSYSTEMS; i++) {
         for (j=0; j<NGRBINS; j++) {
            gr[i][j]		=	0.0;
	    grcom[i][j]		=	0.0;
            grcomz[i][j]	=	0.0;
            grcomxy[i][j]	=	0.0;
	 }
      }

      // Calculate g(r) for all beads

      for (moli=mol; moli<mol+NMOLS; moli++) {
         if ( (system=moli->box)>=0 ) {
            for (i=0; i<moli->nsites; i++) {

               for (molj=moli; molj<mol+NMOLS; molj++) {
                  if ( molj->box == system ) {
                     for (j=0; j<molj->nsites; j++) {
                        if ( (moli!=molj) ? 1 : j>=i+DLJ ) {
                        //if ( (moli!=molj) ? 1 : j>=i+1 ) {
	
                           rij0	=	V_Subtr(molj->p+j, moli->p+i);
			   rij	=	MapInBox2(&rij0, PBC, system);
			   r	=	sqrt(V_Dot(&rij, &rij));	// pair distance

			   if (r < Lupper[system]) {	// no need, because of MapInBox2
			      igr	=	(int) (r/binsize[system]);
			      gr[system][igr]	+=	2;	// contribution from i and j
			   }
                        }
               }  }  }
      }  }  } 

      // Calculate g(r) for center of mass of chains

      for (moli=mol; moli<mol+NMOLS-1; moli++) {
         if ( (system=moli->box)>=0 ) {
            comi	=	CenterofMass(moli); 

	    for (molj=moli+1; molj<mol+NMOLS; molj++) {
	       if (molj->box == system) {
		  comj	=	CenterofMass(molj);

		  rij0	=	V_Subtr(&comi, &comj);
		  rij	=	MapInBox2(&rij0, PBC, system);
		  r	=	sqrt(V_Dot(&rij, &rij));

		  if (r<Lupper[system]) {
	             igr	=	(int) (r/binsize[system]);
		     grcom[system][igr]	+=	2;

		     if (fabs(rij.z)*unit.LENGTH > 8.0)		// 8 angstrom apart in z direction
			grcomz[system][igr]	+=	2; 
		     else
			grcomxy[system][igr]	+=	2;
		  }
            }  }
      }  }

      // normalize g(r), do it every time because BOX.vol changes
      // In fact, it is not necessary as we assume the system in equilibrium

      for (i=0; i<NSYSTEMS; i++) {
         //fprintf(frdf, "RDF in System %ld (%ld mols %ld sites and vol = %f)\n", 
	//	i, NMols[i], NSites[i], BOX[i].vol);

         //fprintf(frdf, "r \t gr \t grcom \t grcomz \t grcomxy\n");
         for (j=0; j<NGRBINS; j++) {
            vb	=	4.0*M_PI/3 
			* ((j+1)*(j+1)*(j+1)-j*j*j) * binsize[i] * binsize[i] * binsize[i];
            nid	=	vb * NSites[i]/BOX[i].vol; 
            gr[i][j]	/=	(NSites[i] * nid);
            grtot[i][j]	+=	gr[i][j];

	    nid =	vb * NMols[i]/BOX[i].vol;
            grcom[i][j]		/=	(NMols[i] * nid);
	    grcomz[i][j]	/=	(NMols[i] * nid);
	    grcomxy[i][j]	/=	(NMols[i] * nid);

	    grcomtot[i][j]	+=	grcom[i][j];
	    grcomztot[i][j]	+=	grcomz[i][j];
	    grcomxytot[i][j]	+=	grcomxy[i][j];

            //fprintf(frdf, "%f %f %f %f %f\n", binsize[i]*(j+0.5)*unit.LENGTH, 
		//gr[i][j], grcom[i][j], grcomz[i][j], grcomxy[i][j]);
         }
         ndensity[i]		+=	NSites[i]/BOX[i].vol;
      }
      ngr	++;				// # of sampling increase by 1
   }
   else if (!strcmp(grswitch, "print")) {

      // Calculate ensemble average

      for (i=0; i<NSYSTEMS; i++) {
         for (j=0; j<NGRBINS; j++) {
            grtot[i][j]		/=	ngr;
            grcomtot[i][j]	/=	ngr;
	    grcomztot[i][j]	/=	ngr;
	    grcomxytot[i][j]	/=	ngr;
         }
         ndensity[i]		/=	ngr;
      }

      // Calculate structure factor S(q) by doing FFT to rdf

      for (system=0; system<NSYSTEMS; system++) {
         for (i_q=0; i_q<NGRBINS; i_q++) {
            sq[system][i_q]	=	0.0;
            q	=	(i_q+0.5) * dq;

            for (i=0; i<NGRBINS; i++) {
               r	=	(i+0.5) * binsize[system];

               if (r>0.45*Lupper[system])	break;
		
               temp	=	(grtot[system][i]-1) * sin(q*r) / (q*r);
               temp	*=	ndensity[system];
               temp	*=	4*M_PI/3 * ((i+1)*(i+1)*(i+1)-i*i*i) *
               			binsize[system]*binsize[system]*binsize[system];
               sq[system][i_q]	+=	temp;
            }
            sq[system][i_q]	+=	1.0;
         }
      }
/*
      for (system=0; system<NSYSTEMS; system++) {
         for (i_q=0; i_q<NGRBINS; i_q++) {
            sq[system][i_q]	=	0.0;
  
            for (i=0; i<NGRBINS; i++) {
               if ((i+0.5)*binsize[system] >= 0.625*Lupper[system]) {
                  break;
               }			// upper limit of integral set to half box size
               sq[system][i_q]	+=	(gr[system][i]-1) * (i+0.5)*binsize[system]
					* sin((i_q+0.5)*dq*(i+0.5)*binsize[system])
					* binsize[system];
            }
            sq[system][i_q]	*=	4*M_PI*NSites[system]/BOX[system].vol/((i_q+0.5)*dq);
            sq[system][i_q]	+=	1.0;
            sqtot[system][i_q]	+=	sq[system][i_q];
         }
      }
*/
      // Output radial distribution function g(r) and structure factor S(q)

      for (i=0; i<NSYSTEMS; i++) {
         fprintf(frdf, "RDF and structure factor in system %ld. # of sample = %ld\n", system, ngr);
         fprintf(frdf, "r \t grtot \t grcomtot \t grcomztot \t grcomxytot \t q \t S(q)\n");
         for (j=0; j<NGRBINS; j++) {
            fprintf(frdf, "%f %f %f %f %f %f %f\n", binsize[i]*(j+0.5)*unit.LENGTH, 
		grtot[i][j], grcomtot[i][j], grcomztot[i][j], grcomxytot[i][j],
		(j+0.5)*dq/unit.LENGTH, sq[i][j]);
         }
      }
      fclose(frdf);
      for (i=0; i<NSYSTEMS; i++) {
         free(gr[i]);
         free(grtot[i]);
         free(grcom[i]);
         free(grcomtot[i]);
         free(grcomz[i]);
         free(grcomztot[i]);
         free(grcomxy[i]);
         free(grcomxytot[i]);
	 free(sq[i]);
	 free(sqtot[i]);
      }
   }
   else {
      Exit("sample.c", "radial", "parameter not found");
   }
   return;
}

//////////////////////////////////////////
/* 	S(k) calculated by definition	*/
/*	Add: 4/17/2012	- 4/20/2012	*/
//////////////////////////////////////////

// sq2 and sq3 are the correct ones that calculate S(q) from definition 
// Both both used advanced memory allocation method
// sq3 need smaller memory due to separation of x, y, z component
// for sq2, when Mx, My, Mz >20 memory is not enough

// both use discrete q's, rather than ranor() with continuous q's
// due to the discrete nature of q's, both applied smart method
// to calculate next q from previous q.

// results of sq2 and sq3 are same, using Mx=My=Mz=10 as test
// sq3 is better because it starts from q=0, where real=1 and imag=0 and 
// it is very precise as a starting point.  And sq3 uses the contugate
// relations to reduce computation burden.

void sq4(char *sqswitch)
{
   molstruct		*moli;
   vector		vq, r;
   long			i, j, k, n, i_q, id;
   long			nx, ny, nz, ix, iy, iz;
   static long		NSQBINS, Mx, My, Mz;
   static float		qmax, dq;
   float		x, y, z, inv_Lx, inv_Ly, inv_Lz;
   float		sumre, sumim, q;
   float		coskx, cosky, coskz, sinkx, sinky, sinkz;
   static float		*cosx, *cosy, *cosz, *sinx, *siny, *sinz;
   static float		**realx, **realy, **realz, **imagx, **imagy, **imagz;
   static float		*sq, *sqtot, *sq2tot;
   static long		*nsq, ntot;
   static long		init = 1;

   static FILE 		*fPtr;

   if (init) {
      init	=	0;
      fPtr = fopen("sq4.dat", "w");
      fprintf(fPtr, "sq4 file opened.\n");
      fflush(fPtr);

      dq	=	0.1;		// in system unit, sigma^-1
      //dq	=	0.02;		// in system unit, sigma^-1
      qmax	=	4*M_PI;		// in system unit, sigma^-1
      NSQBINS	=	(long) (qmax/dq);
      Mx	=	(long) (qmax*BOX[0].lx/(2*M_PI));
      My	=	(long) (qmax*BOX[0].ly/(2*M_PI));
      Mz	=	(long) (qmax*BOX[0].lz/(2*M_PI));

      cosx	=	(float *) calloc((size_t)NSITES, sizeof(float));
      cosy	=	(float *) calloc((size_t)NSITES, sizeof(float));
      cosz	=	(float *) calloc((size_t)NSITES, sizeof(float));
      sinx	=	(float *) calloc((size_t)NSITES, sizeof(float));
      siny	=	(float *) calloc((size_t)NSITES, sizeof(float));
      sinz	=	(float *) calloc((size_t)NSITES, sizeof(float));

      realx	=	(float **) calloc((size_t) (2*Mx+1), sizeof(float *));
      realy	=	(float **) calloc((size_t) (2*My+1), sizeof(float *));
      realz	=	(float **) calloc((size_t) (2*Mz+1), sizeof(float *));
      imagx	=	(float **) calloc((size_t) (2*Mx+1), sizeof(float *));
      imagy	=	(float **) calloc((size_t) (2*My+1), sizeof(float *));
      imagz	=	(float **) calloc((size_t) (2*Mz+1), sizeof(float *));
      realx[0]	=	(float *) calloc((size_t) (2*Mx+1)*NSITES, sizeof(float));
      realy[0]	=	(float *) calloc((size_t) (2*My+1)*NSITES, sizeof(float));
      realz[0]	=	(float *) calloc((size_t) (2*Mz+1)*NSITES, sizeof(float));
      imagx[0]	=	(float *) calloc((size_t) (2*Mx+1)*NSITES, sizeof(float));
      imagy[0]	=	(float *) calloc((size_t) (2*My+1)*NSITES, sizeof(float));
      imagz[0]	=	(float *) calloc((size_t) (2*Mz+1)*NSITES, sizeof(float));

      for (n=1; n<2*Mx+1; n++) {
         realx[n]	=	realx[n-1] + NSITES;
         imagx[n]	=	imagx[n-1] + NSITES;
      }
      for (n=1; n<2*My+1; n++) {
         realy[n]	=	realy[n-1] + NSITES;
         imagy[n]	=	imagy[n-1] + NSITES;
      }
      for (n=1; n<2*Mz+1; n++) {
         realz[n]	=	realz[n-1] + NSITES;
         imagz[n]	=	imagz[n-1] + NSITES;
      }

      sq	=	(float *) calloc((size_t)NSQBINS, sizeof(float));
      sqtot	=	(float *) calloc((size_t)NSQBINS, sizeof(float));
      sq2tot	=	(float *) calloc((size_t)NSQBINS, sizeof(float));
      nsq	=	(long *) calloc((size_t)NSQBINS, sizeof(long));
      ntot	=	0;

      fprintf(fPtr,"sq4 initialization done\n");
      fprintf(fPtr,"%d %d %d %d\n", Mx, My, Mz, NSQBINS);
      fflush(fPtr);
   }
   
   if (!strcmp(sqswitch,"sample")) {		// calculate S(q) for one conf.

      inv_Lx	=	2*M_PI/BOX[0].lx;
      inv_Ly	=	2*M_PI/BOX[0].ly;
      inv_Lz	=	2*M_PI/BOX[0].lz;

      id	=	0;
      for (moli=mol; moli<mol+NMOLS; moli++) {
         for (i=0; i<moli->nsites; i++) { 
            r	=	moli->p[i];
            x	=	r.x * inv_Lx;
            y	=	r.y * inv_Ly;
            z	=	r.z * inv_Lz;

            cosx[id]	=	cos(x);           sinx[id]	=	sin(x);
            cosy[id]	=	cos(y);           siny[id]	=	sin(y);
            cosz[id]	=	cos(z);           sinz[id]	=	sin(z);

            realx[Mx][id]	=	1;	imagx[Mx][id]	=	0;	// qx=0
            realy[My][id]	=	1;	imagy[My][id]	=	0;	// qy=0
            realz[Mz][id]	=	1;	imagz[Mz][id]	=	0;	// qz=0

            for (nx=1; nx<=Mx; nx++) {
               realx[nx+Mx][id]	=	realx[nx-1+Mx][id] * cosx[id]
				-	imagx[nx-1+Mx][id] * sinx[id];
               imagx[nx+Mx][id]	=	imagx[nx-1+Mx][id] * cosx[id]
				+	realx[nx-1+Mx][id] * sinx[id];
               realx[-nx+Mx][id]=	realx[nx+Mx][id];
               imagx[-nx+Mx][id]=	-imagx[nx+Mx][id];
            }
            for (ny=1; ny<=My; ny++) {
               realy[ny+My][id]	=	realy[ny-1+My][id] * cosy[id]
				-	imagy[ny-1+My][id] * siny[id];
               imagy[ny+My][id]	=	imagy[ny-1+My][id] * cosy[id]
				+	realy[ny-1+My][id] * siny[id];
               realy[-ny+My][id]=	realy[ny+My][id];
               imagy[-ny+My][id]=	-imagy[ny+My][id];
            }
            for (nz=1; nz<=Mz; nz++) {
               realz[nz+Mz][id]	=	realz[nz-1+Mz][id] * cosz[id]
				-	imagz[nz-1+Mz][id] * sinz[id];
               imagz[nz+Mz][id]	=	imagz[nz-1+Mz][id] * cosz[id]
				+	realz[nz-1+Mz][id] * sinz[id];
               realz[-nz+Mz][id]=	realz[nz+Mz][id];
               imagz[-nz+Mz][id]=	-imagz[nz+Mz][id];
            }
            id	++;
      }  }	// loop over all particles

      for (i_q=0; i_q<NSQBINS; i_q++) {		// initialize
        nsq[i_q]	=	0;
        sq[i_q]		=	0.0;
      }

      for (nx=-Mx; nx<=Mx; nx++) {		// for all q's
         ix	=	nx+Mx;
         for (ny=-My; ny<=My; ny++) {
            iy	=	ny+My;
            for (nz=-Mz; nz<=Mz; nz++) {
               iz	=	nz+Mz;

               sumre	=	0.0;		// sum over all particles
               sumim	=	0.0;
               for (id=0; id<NSITES; id++) {	// sum over all particles
                  coskx	=	realx[ix][id];	sinkx	=	imagx[ix][id];
                  cosky	=	realy[iy][id];	sinky	=	imagy[iy][id];
                  coskz	=	realz[iz][id];	sinkz	=	imagz[iz][id];

                  sumre	+=	(coskx * cosky - sinkx * sinky) * coskz
			 -	(sinkx * cosky + coskx * sinky) * sinkz; 
                  sumim	+=	(sinkx * cosky + coskx * sinky) * coskz
			 +	(coskx * cosky - sinkx * sinky) * sinkz;
               }

               vq.x	=	nx * inv_Lx;	// find out the q bin that fit
               vq.y	=	ny * inv_Ly;	// this {nx, ny, nz}
               vq.z	=	nz * inv_Lz;
               q	=	V_Dot(&vq, &vq);
	       q	=	sqrt(q);
               i_q	=	(int) (q/dq);
               if (i_q<NSQBINS) { 
                  sq[i_q]	+=	(sumre * sumre + sumim * sumim)/NSITES;
                  nsq[i_q]	++;
               }
      }  }  } 	//	for all q's
      for (i_q=0; i_q<NSQBINS; i_q++) {		// for ensemble average
         if (nsq[i_q]>0) {
            sqtot[i_q]	+=	sq[i_q]/nsq[i_q];
            sq2tot[i_q]	+=	sq[i_q]*sq[i_q]/(nsq[i_q]*nsq[i_q]);
         }
      }
      ntot	++;
      fprintf(fPtr,"sq4 %d\n", ntot);
   }
   else if (!strcmp(sqswitch, "print")) {	// note "print" after unit conversion to SI
      printf("\n### SQ4: STRUCTURE FACTOR S(q) CALCULATED BASED ON DEFINITION\n\n");

      printf("System Units:\n");
      printf("qmax = %f\tdq = %f\tNSQBINS = %ld\n", qmax, dq, NSQBINS);
      printf("(Lx, Ly, Lz) = (%f, %f, %f)\n", BOX[0].lx/unit.LENGTH, BOX[0].ly/unit.LENGTH, BOX[0].lz/unit.LENGTH);
      printf("(Mx, My, Mz) = (%ld, %ld, %ld)\n", Mx, My, Mz);

      printf("\nAngstrom:\n");
      printf("qmax = %f\tdq = %f\tNSQBINS = %ld\n", qmax/unit.LENGTH, dq/unit.LENGTH, NSQBINS);
      printf("(Lx, Ly, Lz) = (%f, %f, %f)\n", BOX[0].lx, BOX[0].ly, BOX[0].lz);
      printf("2PI/(Lx, Ly, Lz) = (%f, %f, %f)\n", 2*M_PI/BOX[0].lx, 2*M_PI/BOX[0].ly, 2*M_PI/BOX[0].lz);
      printf("q_cut = 2PI*sqrt(1/Lx^2+1/Ly^2+1/Lz^2) = %f\n",
              2*M_PI*sqrt( 1/(BOX[0].lx*BOX[0].lx) + 1/(BOX[0].ly*BOX[0].ly) + 1/(BOX[0].lz*BOX[0].lz) ) );
      
      printf("\nq(Ang)\t\tS(q)\t\tErr\t\tnsq_lastframe\n");
      for (i_q=0; i_q<NSQBINS; i_q++) {
         sqtot[i_q]	/=	ntot;
         sq2tot[i_q]	/=	ntot;
         sq2tot[i_q]	=	sqrt(sq2tot[i_q] - sqtot[i_q]*sqtot[i_q]);
         printf("%f\t%f\t%f\t%ld\n", (i_q+0.5)*dq/unit.LENGTH, sqtot[i_q], sq2tot[i_q], nsq[i_q]);
      }

      free(nsq); 	free(sq);
      free(sqtot);	free(sq2tot);
      free(cosx);	free(sinx);
      free(cosy);	free(siny);
      free(cosz);	free(sinz);
      free(realx[0]);	free(imagx[0]);
      free(realy[0]);	free(imagy[0]);	
      free(realz[0]);	free(imagz[0]);	
      free(realx);	free(imagx);
      free(realy);	free(imagy);	
      free(realz);	free(imagz);	

      fprintf(fPtr,"sq4 done\n");
      fclose(fPtr);
   }
   return;   
}

void sq3(char *sqswitch)
{
   molstruct		*moli;
   vector		vq, r;
   long			i, j, k, n, i_q, id;
   long			M, Mx, My, Mz, nx, ny, nz, ix, iy, iz;
   float		x, y, z, inv_Lx, inv_Ly, inv_Lz;
   float		sumre, sumim, q, dq=0.2;
   float		coskx, cosky, coskz, sinkx, sinky, sinkz;
   static float		*cosx, *cosy, *cosz, *sinx, *siny, *sinz;
   static float		**realx, **realy, **realz, **imagx, **imagy, **imagz;
   static float		*sq, *sqtot;
   static long		*nsq, ntot;
   static long		init = 1;

   Mx	=	(long) (BOX[0].lx*2);
   My	=	(long) (BOX[0].ly*2);
   Mz	=	(long) (BOX[0].lz*2);

   if (init) {
      init	=	0;
      cosx	=	(float *) calloc((size_t)NSITES, sizeof(float));
      cosy	=	(float *) calloc((size_t)NSITES, sizeof(float));
      cosz	=	(float *) calloc((size_t)NSITES, sizeof(float));
      sinx	=	(float *) calloc((size_t)NSITES, sizeof(float));
      siny	=	(float *) calloc((size_t)NSITES, sizeof(float));
      sinz	=	(float *) calloc((size_t)NSITES, sizeof(float));

      realx	=	(float **) calloc((size_t) (2*Mx+1), sizeof(float *));
      realy	=	(float **) calloc((size_t) (2*My+1), sizeof(float *));
      realz	=	(float **) calloc((size_t) (2*Mz+1), sizeof(float *));
      imagx	=	(float **) calloc((size_t) (2*Mx+1), sizeof(float *));
      imagy	=	(float **) calloc((size_t) (2*My+1), sizeof(float *));
      imagz	=	(float **) calloc((size_t) (2*Mz+1), sizeof(float *));
      realx[0]	=	(float *) calloc((size_t) (2*Mx+1)*NSITES, sizeof(float));
      realy[0]	=	(float *) calloc((size_t) (2*My+1)*NSITES, sizeof(float));
      realz[0]	=	(float *) calloc((size_t) (2*Mz+1)*NSITES, sizeof(float));
      imagx[0]	=	(float *) calloc((size_t) (2*Mx+1)*NSITES, sizeof(float));
      imagy[0]	=	(float *) calloc((size_t) (2*My+1)*NSITES, sizeof(float));
      imagz[0]	=	(float *) calloc((size_t) (2*Mz+1)*NSITES, sizeof(float));

      for (n=1; n<2*Mx+1; n++) {
         realx[n]	=	realx[n-1] + NSITES;
         imagx[n]	=	imagx[n-1] + NSITES;
      }
      for (n=1; n<2*My+1; n++) {
         realy[n]	=	realy[n-1] + NSITES;
         imagy[n]	=	imagy[n-1] + NSITES;
      }
      for (n=1; n<2*Mz+1; n++) {
         realz[n]	=	realz[n-1] + NSITES;
         imagz[n]	=	imagz[n-1] + NSITES;
      }

      sq	=	(float *) calloc((size_t)NGRBINS, sizeof(float));
      sqtot	=	(float *) calloc((size_t)NGRBINS, sizeof(float));
      nsq	=	(long *) calloc((size_t)NGRBINS, sizeof(long));
      ntot	=	0;
   }
   
   if (!strcmp(sqswitch,"sample")) {		// calculate S(q) for one conf.

      inv_Lx	=	2*M_PI/BOX[0].lx;
      inv_Ly	=	2*M_PI/BOX[0].ly;
      inv_Lz	=	2*M_PI/BOX[0].lz;

      id	=	0;
      for (moli=mol; moli<mol+NMOLS; moli++) {
         for (i=0; i<moli->nsites; i++) { 
            r	=	moli->p[i];
            x	=	r.x * inv_Lx;
            y	=	r.y * inv_Ly;
            z	=	r.z * inv_Lz;

            cosx[id]	=	cos(x);           sinx[id]	=	sin(x);
            cosy[id]	=	cos(y);           siny[id]	=	sin(y);
            cosz[id]	=	cos(z);           sinz[id]	=	sin(z);

            realx[Mx][id]	=	1;	imagx[Mx][id]	=	0;	// qx=0
            realy[My][id]	=	1;	imagy[My][id]	=	0;	// qy=0
            realz[Mz][id]	=	1;	imagz[Mz][id]	=	0;	// qz=0

            for (nx=1; nx<=Mx; nx++) {
               realx[nx+Mx][id]	=	realx[nx-1+Mx][id] * cosx[id]
				-	imagx[nx-1+Mx][id] * sinx[id];
               imagx[nx+Mx][id]	=	imagx[nx-1+Mx][id] * cosx[id]
				+	realx[nx-1+Mx][id] * sinx[id];
               realx[-nx+Mx][id]=	realx[nx+Mx][id];
               imagx[-nx+Mx][id]=	-imagx[nx+Mx][id];
            }
            for (ny=1; ny<=My; ny++) {
               realy[ny+My][id]	=	realy[ny-1+My][id] * cosy[id]
				-	imagy[ny-1+My][id] * siny[id];
               imagy[ny+My][id]	=	imagy[ny-1+My][id] * cosy[id]
				+	realy[ny-1+My][id] * siny[id];
               realy[-ny+My][id]=	realy[ny+My][id];
               imagy[-ny+My][id]=	-imagy[ny+My][id];
            }
            for (nz=1; nz<=Mz; nz++) {
               realz[nz+Mz][id]	=	realz[nz-1+Mz][id] * cosz[id]
				-	imagz[nz-1+Mz][id] * sinz[id];
               imagz[nz+Mz][id]	=	imagz[nz-1+Mz][id] * cosz[id]
				+	realz[nz-1+Mz][id] * sinz[id];
               realz[-nz+Mz][id]=	realz[nz+Mz][id];
               imagz[-nz+Mz][id]=	-imagz[nz+Mz][id];
            }
            id	++;
      }  }	// loop over all particles

      for (i_q=0; i_q<NGRBINS; i_q++) {		// initialize
        nsq[i_q]	=	0;
        sq[i_q]		=	0.0;
      }

      for (nx=-Mx; nx<=Mx; nx++) {		// for all q's
         ix	=	nx+Mx;
         for (ny=-My; ny<=My; ny++) {
            iy	=	ny+My;
            for (nz=-Mz; nz<=Mz; nz++) {
               iz	=	nz+Mz;

               sumre	=	0.0;		// sum over all particles
               sumim	=	0.0;
               for (id=0; id<NSITES; id++) {	// sum over all particles
                  coskx	=	realx[ix][id];	sinkx	=	imagx[ix][id];
                  cosky	=	realy[iy][id];	sinky	=	imagy[iy][id];
                  coskz	=	realz[iz][id];	sinkz	=	imagz[iz][id];

                  sumre	+=	(coskx * cosky - sinkx * sinky) * coskz
			 -	(sinkx * cosky + coskx * sinky) * sinkz; 
                  sumim	+=	(sinkx * cosky + coskx * sinky) * coskz
			 +	(coskx * cosky - sinkx * sinky) * sinkz;
               }

               vq.x	=	nx * inv_Lx;	// find out the q bin that fit
               vq.y	=	ny * inv_Ly;	// this {nx, ny, nz}
               vq.z	=	nz * inv_Lz;
               q	=	V_Dot(&vq, &vq);
	       q	=	sqrt(q);
               i_q	=	(int) (q/dq);
               if (i_q<NGRBINS) { 
                  sq[i_q]	+=	(sumre * sumre + sumim * sumim)/NSITES;
                  nsq[i_q]	++;
               }
      }  }  } 	//	for all q's
      for (i_q=0; i_q<NGRBINS; i_q++) {		// for ensemble average
         if (nsq[i_q]>0)
            sqtot[i_q]	+=	sq[i_q]/nsq[i_q];
      }
      ntot	++;
   }
   else if (!strcmp(sqswitch, "print")) {
      printf("\n### SQ3: STRUCTURE FACTOR S(q) CALCULATED BASED ON DEFINITION\n\n");
      printf("(Mx, My, Mz) = (%ld, %ld, %ld)\n", Mx, My, Mz);
      printf("ntot = %ld\n", ntot);
      printf("maxL = %f\n", MAX(MAX(BOX[0].lx, BOX[0].ly), BOX[0].lz)*unit.LENGTH);
      printf("minL = %f\n", MIN(MIN(BOX[0].lx, BOX[0].ly), BOX[0].lz)*unit.LENGTH);
      printf("min(minq) = %f\n", 2*M_PI/(MAX(MAX(BOX[0].lx, BOX[0].ly), BOX[0].lz)*unit.LENGTH));
      printf("max(minq) = %f\n", 2*M_PI/(MIN(MIN(BOX[0].lx, BOX[0].ly), BOX[0].lz)*unit.LENGTH));

      for (i_q=0; i_q<NGRBINS; i_q++) {
         sqtot[i_q]	/=	ntot;
         printf("%f\t%f\t%ld\n", (i_q+0.5)*dq/unit.LENGTH, sqtot[i_q], nsq[i_q]);
      }
   }
   return;   
}


void sq2(char *sqswitch)
{
   molstruct		*moli;
   vector		vq, r;
   long			i, j, k, n, i_q, id;
   long			M, Mx, My, Mz, nx, ny, nz, ix, iy, iz;
   long			pre_nx, pre_ny, pre_nz;
   float		x, y, z, inv_Lx, inv_Ly, inv_Lz;
   float		sumre, sumim, q, dq=0.2;
   //static int		***grid, *ncount;
   static float		*cosx, *cosy, *cosz, *sinx, *siny, *sinz;
   static float		****re, ****im;
   //static float		re[21][21][21][9000], im[21][21][21][9000];
   //static float		***real, ***imag, *realbase, *imagbase;
   static float		*sq, *sqtot;		// S(q) and ensemble average
   static long		*nsq, ntot;		// # of samplings
   static long		init = 1;

   //M	=	10;
   Mx	=	(long) (BOX[0].lx*2);
   My	=	(long) (BOX[0].ly*2);
   Mz	=	(long) (BOX[0].lz*2);
   Mx=10;	
   My=10;
   Mz=10;

   if (init) {
      init	=	0;
      cosx	=	(float *) calloc((size_t)NSITES, sizeof(float));
      cosy	=	(float *) calloc((size_t)NSITES, sizeof(float));
      cosz	=	(float *) calloc((size_t)NSITES, sizeof(float));
      sinx	=	(float *) calloc((size_t)NSITES, sizeof(float));
      siny	=	(float *) calloc((size_t)NSITES, sizeof(float));
      sinz	=	(float *) calloc((size_t)NSITES, sizeof(float));
/*     
      re	=	(float ****) calloc(2*Mx+1, sizeof(float ***));
      im	=	(float ****) calloc(2*Mx+1, sizeof(float ***));

      for (i=0; i<=2*Mx; i++) {
         re[i]		=	(float ***) calloc(2*My+1, sizeof(float **));
         im[i]		=	(float ***) calloc(2*My+1, sizeof(float **));
         for (j=0; j<=2*My; j++) {
            re[i][j]	=	(float **) calloc(2*Mz+1, sizeof(float *));
            im[i][j]	=	(float **) calloc(2*Mz+1, sizeof(float *));
            for (k=0; k<NSITES; k++) {
               re[i][j][k]	=	(float *) calloc(NSITES, sizeof(float));
               im[i][j][k]	=	(float *) calloc(NSITES, sizeof(float));
            }
         }
      }
*/

      re	=	(float ****) calloc((size_t) (2*Mx+1), sizeof(float ***));
      re[0]	=	(float ***)  calloc((size_t) ((2*Mx+1)*(2*My+1)), sizeof(float **));
      re[0][0]	=	(float **)   calloc((size_t) ((2*Mx+1)*(2*My+1)*(2*Mz+1)), sizeof(float *));
      re[0][0][0]	=	(float *) calloc((size_t) ((2*Mx+1)*(2*My+1)*(2*Mz+1)*NSITES), sizeof(float));
      if (!re[0][0][0]) {
         printf("re allocation failed\n");
         exit(1);
      }
      im	=	(float ****) calloc((size_t) (2*Mx+1), sizeof(float ***));
      im[0]	=	(float ***)  calloc((size_t) ((2*Mx+1)*(2*My+1)), sizeof(float **));
      im[0][0]	=	(float **)   calloc((size_t) ((2*Mx+1)*(2*My+1)*(2*Mz+1)), sizeof(float *));
      im[0][0][0]	=	(float *) calloc((size_t) ((2*Mx+1)*(2*My+1)*(2*Mz+1)*NSITES), sizeof(float));
      if (!im[0][0][0]) {
         printf("im allocation failed\n");
         exit(1);
      }
      for (n=1; n<2*Mx+1; n++) {
         re[n]	=	re[n-1]  + (2*My+1);
         im[n]	=	im[n-1]  + (2*My+1);
      }
      for (n=1; n<(2*Mx+1)*(2*My+1); n++) {
         re[0][n]	=	re[0][n-1] + (2*Mz+1);
         im[0][n]	=	im[0][n-1] + (2*Mz+1);
      }
      for (n=1; n<(2*Mx+1)*(2*My+1)*(2*Mz+1); n++) {
         re[0][0][n]	=	re[0][0][n-1] + NSITES;
         im[0][0][n]	=	im[0][0][n-1] + NSITES;
      }

      sq	=	(float *) calloc((size_t)NGRBINS, sizeof(float));
      sqtot	=	(float *) calloc((size_t)NGRBINS, sizeof(float));
      nsq	=	(long *) calloc((size_t)NGRBINS, sizeof(long));
      ntot	=	0;
//printf("memory allocation done\n");
//exit(1);
   }

   if (!strcmp(sqswitch,"sample")) {		// calculate S(q) for one conf.

      for (i_q=0; i_q<NGRBINS; i_q++) {		// initialize
        nsq[i_q]	=	0;
        sq[i_q]		=	0.0;
      }

      inv_Lx	=	2*M_PI/BOX[0].lx;
      inv_Ly	=	2*M_PI/BOX[0].ly;
      inv_Lz	=	2*M_PI/BOX[0].lz;

      id	=	0;
      for (moli=mol; moli<mol+NMOLS; moli++) {
        for (i=0; i<moli->nsites; i++) { 
           r	=	moli->p[i];
           x	=	r.x * inv_Lx;
           y	=	r.y * inv_Ly;
           z	=	r.z * inv_Lz;

           cosx[id]	=	cos(x);           sinx[id]	=	sin(x);
           cosy[id]	=	cos(y);           siny[id]	=	sin(y);
           cosz[id]	=	cos(z);           sinz[id]	=	sin(z);

           re[0][0][0][id]	=	cos(Mx*x+My*y+Mz*z);	// real part
           im[0][0][0][id]	=	-sin(Mx*x+My*y+Mz*z);	// imaginary part
           id	++;
      } }

      pre_nx	=	-Mx;
      pre_ny	=	-My;	
      pre_nz	=	-Mz;
      for (nx=-Mx; nx<=Mx; nx++) {
        ix	=	nx+Mx;
        for (ny=-My; ny<=My; ny++) {
          iy	=	ny+My;
          for (nz=-Mz; nz<=Mz; nz++) {
            iz	=	nz+Mz;

            sumre	=	0.0;	// over all particles
            sumim	=	0.0;

            id		=	0;
            for (moli=mol; moli<mol+NMOLS; moli++) {
              for (i=0; i<moli->nsites; i++) {

                if (nz==(pre_nz+1)) {
                  re[ix][iy][iz][id]	=	re[ix][iy][iz-1][id] * cosz[id]
					-	im[ix][iy][iz-1][id] * sinz[id];	
                  im[ix][iy][iz][id]	=	im[ix][iy][iz-1][id] * cosz[id]
					+	re[ix][iy][iz-1][id] * sinz[id];
                }
                else if (ny==(pre_ny+1)) {
                  re[ix][iy][iz][id]	=	re[ix][iy-1][iz][id] * cosy[id]
					-	im[ix][iy-1][iz][id] * siny[id];	
                  im[ix][iy][iz][id]	=	im[ix][iy-1][iz][id] * cosy[id]
					+	re[ix][iy-1][iz][id] * siny[id];
                }
                else if (nx==(pre_nx+1)) {
                  re[ix][iy][iz][id]	=	re[ix-1][iy][iz][id] * cosx[id]
					-	im[ix-1][iy][iz][id] * sinx[id];	
                  im[ix][iy][iz][id]	=	im[ix-1][iy][iz][id] * cosx[id]
					+	re[ix-1][iy][iz][id] * sinx[id];
                }
                sumre	+=	re[ix][iy][iz][id];
                sumim	+=	im[ix][iy][iz][id];
		id	++;

            } }	// loop over particles
            vq.x	=	nx * inv_Lx;
            vq.y	=	ny * inv_Ly;
            vq.z	=	nz * inv_Lz;
            q		=	V_Dot(&vq, &vq);
            q		=	sqrt(q);
            i_q		=	(int) (q/dq);
            if (i_q<NGRBINS) {
               sq[i_q]	+=	(sumre*sumre + sumim*sumim)/NSITES;	// one q contribution
               nsq[i_q]	++; 
            }
            pre_nx	=	nx;
            pre_ny	=	ny;
            pre_nz	=	nz;
      } } }	// loop over nx, ny, nz

      for (i_q=0; i_q<NGRBINS; i_q++) {		// for ensemble average
         if (nsq[i_q]>0)
            sqtot[i_q]	+=	sq[i_q]/nsq[i_q]; 
      }
      ntot	++;
   }
   else if (!strcmp(sqswitch,"print")) {
      printf("\n### SQ2: STRUCTURE FACTOR S(q) CALCULATED BASED ON DEFINITION\n\n");
      printf("(Mx, My, Mz) = (%ld, %ld, %ld)\n", Mx, My, Mz);
      printf("ntot = %ld\n", ntot);
      printf("maxL = %f\n", MAX(MAX(BOX[0].lx, BOX[0].ly), BOX[0].lz)*unit.LENGTH);
      printf("minL = %f\n", MIN(MIN(BOX[0].lx, BOX[0].ly), BOX[0].lz)*unit.LENGTH);
      printf("min(minq) = %f\n", 2*M_PI/(MAX(MAX(BOX[0].lx, BOX[0].ly), BOX[0].lz)*unit.LENGTH));
      printf("max(minq) = %f\n", 2*M_PI/(MIN(MIN(BOX[0].lx, BOX[0].ly), BOX[0].lz)*unit.LENGTH));

      for (i_q=0; i_q<NGRBINS; i_q++) {
         sqtot[i_q]	/=	ntot;
         printf("%f\t%f\t%ld\n", (i_q+0.5)*dq/unit.LENGTH, sqtot[i_q], nsq[i_q]);
      }
   }
   return;
}

void sq1(char *sqswitch)
{
   molstruct		*moli;
   long			i, j, k, n, i_q, nx, ny, nz, M;
   double		q, dq=0.2, real, imag, arg;
   vector		vq, r;
   static long		*nsq, ntot;			// nsq: # of sampling
   static double	*sqtot, *sqre, *sqim; 
   static int		init = 1;

   if (init) {
      nsq	=	(long *) calloc(NGRBINS, sizeof(long));
      sqre	=	(double *) calloc(NGRBINS, sizeof(double));
      sqim	=	(double *) calloc(NGRBINS, sizeof(double));
      sqtot	=	(double *) calloc(NGRBINS, sizeof(double));
      ntot	=	0;
      init	=	0;
   }

   if (!strcmp(sqswitch,"sample")) {		// do sampling
      //M	=	NGRBINS*dq/(2*M_PI/MAX(MAX(BOX[0].lx, BOX[0].ly), BOX[0].lz));
      M	=	10;

      for (i_q=0; i_q<NGRBINS; i_q++) {	
         sqre[i_q]	=	0.0;
         sqim[i_q]	=	0.0;
         nsq[i_q]	=	0;
      }

      for (nx=-M; nx<=M; nx++) {		// sample q space
         for (ny=-M; ny<=M; ny++) {
            for (nz=-M; nz<=M; nz++) {
               vq.x	=	nx/BOX[0].lx*2*M_PI;
               vq.y	=	ny/BOX[0].ly*2*M_PI;
               vq.z	=	nz/BOX[0].lz*2*M_PI;
               q	=	V_Dot(&vq, &vq);
               q	=	sqrt(q);

               real	=	0.0;
	       imag	=	0.0;
               for (moli=mol; moli<mol+NMOLS; moli++) {
                  for (i=0; i<moli->nsites; i++) {
                     r		=	moli->p[i];
                     arg	=	V_Dot(&vq, &r);
                     real	+=	cos(arg);
                     imag	+=	sin(arg);
               }  }
               i_q		=	(int) (q/dq);
               if (i_q<NGRBINS) {
                  sqre[i_q]	+=	real;
                  sqim[i_q]	+=	imag;
                  nsq[i_q]		++;
               }
      }  }  }	// nx, ny, nz loops
      for (i_q=0; i_q<NGRBINS; i_q++) {	
         sqre[i_q]	/=	nsq[i_q];
         sqim[i_q]	/=	nsq[i_q];
      }

      for (i_q=0; i_q<NGRBINS; i_q++) {		// for ensemble average
         sqtot[i_q]	+=	(sqre[i_q]*sqre[i_q] + sqim[i_q]*sqim[i_q])/NSITES;
      }
      ntot	++;
   }
   else if(!strcmp(sqswitch, "print")) {
      printf("\n### STRUCTURE FACTOR S(q) CALCULATED BASED ON DEFINITION\n\n");
      printf("M = %ld\n", M);
      printf("ntot = %ld\n", ntot);
      for (i_q=0; i_q<NGRBINS; i_q++) {
         sqtot[i_q]	/=	ntot;
         printf("%f\t%f\t%ld\n", (i_q+0.5)*dq/unit.LENGTH, sqtot[i_q], nsq[i_q]);
      }
   }
   return;
}
 
/*
void sq1(char *sqswitch) 
{
   molstruct		*moli;
   long			i, n, i_q;
   double		arg, real, imag;
   double		q, dq=0.02;
   vector		vq, r;			// q vector
   static FILE		*fPtr;
   static long		nsq;			// nsq: # of sampling
   static double	*sqtot, *sq; 
   static int		init = 1;

   if (init) {
      fPtr	=	fopen("sq.dat", "w");
      nsq	=	0;
      sq	=	(double *) calloc(NGRBINS, sizeof(double));
      sqtot	=	(double *) calloc(NGRBINS, sizeof(double));
      init	=	0;
   }

   if (!strcmp(sqswitch,"sample")) {		// do sampling
      for (i_q=0; i_q<NGRBINS; i_q++) {
         q	=	(i_q+0.5) * dq;

         real	=	0.0;
         imag	=	0.0;

         for (n=0; n<10; n++) {		// 10 random q vectors
            vq	=	ranor();
            vq	=	V_Mult(q, &vq);
	
            for (moli=mol; moli<mol+NMOLS; moli++) {
               for (i=0; i<moli->nsites; i++) {
                  r	=	moli->p[i];
                  arg	=	V_Dot(&vq, &r);
                  real	+=	cos(arg);
                  imag	+=	sin(arg);
               }
            }
         }	
         real	/=	n;
         imag	/=	n;

         sq[i_q]	=	(real*real+imag*imag)/NSITES;
         sqtot[i_q]	+=	sq[i_q];
      }
      nsq		++;
   }
   else if(!strcmp(sqswitch, "print")) {
      for (i_q=0; i_q<NGRBINS; i_q++) {
         sqtot[i_q]	/=	nsq;
         fprintf(fPtr, "%f\t%f\n", (i_q+0.5)*dq, sqtot[i_q]);
      }      
   }
   fclose(fPtr);
   return;
}
*/

//------------------------------------------------------//
// 	atom_neighbor() build nearest neighbor list	//
//	Add: 12/31/2013					//
//------------------------------------------------------//
void	atom_neighbor()
{
   int		i, j;
   int		ni, nj;
   molstruct *	moli;
   molstruct *	molj;
   float	r2;
   //float	Rnn2 = 14.0652;		// 3.75 ^2 where 3.75 first minimum from g(r)
   //float	Rnn2 = 11.56;		// 3.4 ^2
   float	Rnn = 3.6;
   float	Rnn2 = Rnn*Rnn;		// 3.8 ^2
   float	min, max, ave;
#ifdef CELL_LIST
   static int		k, l;
   static cellstruct	*celli, *cellk;
#endif

   // initialization
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
	 moli->nneigh[i]	=	0;
      }
   }

   // find nearest neighbors
   for (moli=mol; moli<mol+NMOLS; moli++) {
   for (i=0; i<moli->nsites; i++) {

      min=1e6;
      max=0;
      ave=0;

#ifdef CELL_LIST
      celli	=	moli->cell[i];		// identify the local cell id
      for (l=0; l<celli->nneigh; l++) {		// loop over all neighboring cells (self included)
	 cellk	=	celli->neigh[l];

	 for (k=0; k<cellk->nsites; k++) {	// loop over neighboring sites
	    molj	=	cellk->mol[k];
	    j		=	cellk->molsite[k];

	    if ( (molj>moli) || (molj==moli && j>i) ) {
#else
      for (molj=moli; molj<mol+NMOLS; molj++) {
      for (j=0; j<molj->nsites; j++) {{
#endif
      if ( (moli!=molj) ? 1 : j>i ) {

	 r2	=	DistSQ(moli->p[i], molj->p[j], moli->box);


	 if (r2 < Rnn2) {
	    ni	=	moli->nneigh[i];
	    nj	=	molj->nneigh[j];
	    
	    moli->neigh[i][ni]	=	(molj-mol)*(MAXNMOLSITES+1) + j;
	    molj->neigh[j][nj]	=	(moli-mol)*(MAXNMOLSITES+1) + i;

	    moli->nneigh[i]	++;
	    molj->nneigh[j]	++;

	    min	=	MIN(r2, min);
	    max	=	MAX(r2, max);
	    ave	+=	r2;
	 }
      }
      }}
      }

      moli->rmin[i]	=	sqrt(min);
      moli->rmax[i]	=	sqrt(max);
      moli->rave[i]	=	sqrt(ave/moli->nneigh[i]);
   }
   }
   return;
}

//--------------------------------------------------------------//
// 	comm_neigh_para(): common neighborhood parameter	//
//	Add: 12/31/2013						//
//--------------------------------------------------------------//
void	comm_neigh_para()
{
   int		i, j, k, ll, m, n;
   molstruct *	moli;
   molstruct *	molj;
   molstruct *  molk;
   int		ni, nj;
   vector	rik, rjk;
   vector	r;
   float	r2;
   int		ncommon;

   for (moli=mol; moli<mol+NMOLS; moli++) {
   for (i=0; i<moli->nsites; i++) {

      r2	=	0.0;
      ni	=	moli->nneigh[i];		// # of n.n. for i

      // search thru i's  n.n.
      for (m=0; m<ni; m++) {
	 //printf("m = %d\n", m);

	 V_Null(&r);

	 molj	=	moli->neigh[i][m] / (MAXNMOLSITES+1) + mol;
	 j	=	moli->neigh[i][m] % (MAXNMOLSITES+1);

         nj	=	molj->nneigh[j];		// # of n.n. for j

	 //ncommon	=	0;

	 // search thru j's n.n.
	 for (n=0; n<nj; n++) {

            for (ll=0; ll<ni; ll++) {
	       if (molj->neigh[j][n] == moli->neigh[i][ll]) {		// if k is common neighbor for i and j

		  molk	=	molj->neigh[j][n] / (MAXNMOLSITES+1) + mol;
		  k	=	molj->neigh[j][n] % (MAXNMOLSITES+1);

		  rik	=	V_Subtr(molk->p+k, moli->p+i);
		  rjk	=	V_Subtr(molk->p+k, molj->p+j);
		    
		  rik	=	MapInBox2(&rik, PBC, moli->box);
		  rjk	=	MapInBox2(&rjk, PBC, moli->box);

		  r	=	V_Add(&r, &rik);
		  r	=	V_Add(&r, &rjk);

		  //ncommon++;
		  
		  break;
	       }
	    }
	    //printf("n = %d\n", n);
	    //printf("ncommon = %d\n", ncommon);
	 }
	 //printf("%f\n", V_Dot(&r, &r));

	 r2	+=	V_Dot(&r, &r);
      }
      r2	/=	ni;
      
      moli->cnp[i]	=	r2;
   }
   }
   return;
}

//------------------------------------------------------//
// 	hyper_distance(): calculate hyper-distance	//
//	Add: 3/24/2014					//
//------------------------------------------------------//
float	hyper_distance(int N, float *r_ref, float *r_trj)
{
   molstruct	*moli;
   int		i, j;
   int		n;
   int		system = 0;
   float	l = 0.0;
   vector	p, q, dr;
   float	par[MAXNMOLS*MAXNMOLSITES];
   float	max_par, min_par;

   /*
   n	=	0;
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
	 par[n]	=	moli->cnp[i];
	 n	++;
      }
   }
   if (3*n!=N) {
      printf("Error, hyper_distance()\n");
      exit(0);
   }

   max_par	=	6.9;
   min_par	=	6.7;

   // calc. average shift of normal atoms b/w r_trj and r_ref 
   n	=	0;					// # of normal atoms
   V_Null(&p);
   V_Null(&q);

   for (j=0; j<N; j++){
      if (par[j] < max_par && par[j] > min_par) {	// normal atoms
	 i	=	3*j;
	 p.x	+=	r_trj[i];
	 p.y	+=	r_trj[i+1];
         p.z	+=	r_trj[i+2];

	 q.x	+=	r_ref[i];
	 q.y	+=	r_ref[i+1];
         q.z	+=	r_ref[i+2];
      	 n	++;		
      }
   }
   dr	=	V_Subtr(&p, &q); 
   dr	=	MapInBox2(&dr, PBC, system);
   dr	=	V_Mult(1.0/n, &dr);

   // fix the shift
   for (i=0; i<3*N; i+=3) {
      r_trj[i]		-=	dr.x;
      r_trj[i+1]	-=	dr.y;
      r_trj[i+2]	-=	dr.z;
   }


   l	=	0.0;
   for (i=0; i<3*N; i+=3) {
      if (par[i/3] >= 16.0) {			// atoms used to calc. hyper-distance
	 p.x	=	r_trj[i];
	 p.y	=	r_trj[i+1];
	 p.z	=	r_trj[i+2];

	 q.x	=	r_ref[i];
	 q.y	=	r_ref[i+1];
	 q.z	=	r_ref[i+2];

	 l	+=	DistSQ(p, q, system);
      }
   }
   */
   l	=	0.0;
   for (i=0; i<3*N; i+=3) {
      p.x	=	r_trj[i];
      //p.y	=	r_trj[i+1];
      //p.z	=	r_trj[i+2];

      q.x	=	r_ref[i];
      //q.y	=	r_ref[i+1];
      //q.z	=	r_ref[i+2];

      l	+=	DistSQ(p, q, system);
   }
   l	=	sqrt(l);
   return	l;
}

//------------------------------------------------------//
// 	slip_distance(): calculate hyper-distance	//
//	Add: 4/6/2014					//
//------------------------------------------------------//
float	slip_distance(int N, float *r_ref, float *r_trj)
{
   int		i;
   int		system = 0;
   float	l = 0.0;
   vector	p, q;

   V_Null(&p);
   V_Null(&q);
   for (i=0; i<3*N; i+=3) {
      p.x	+=	r_trj[i];
      //p.y	+=	r_trj[i+1];
      //p.z	+=	r_trj[i+2];

      q.x	+=	r_ref[i];
      //q.y	+=	r_ref[i+1];
      //q.z	+=	r_ref[i+2];
   }
   l	=	DistSQ(p, q, system);
   l	=	sqrt(l);
   return	l;
}

//------------------------------------------------------//
//	Calculate non-cubic factor of simulation box	//
//	Add: 2/21/2015					//
//------------------------------------------------------//
float noncubic()
{
   float	tmp1, tmp2, tmp3;
   float	asphericity;
   int		system=0;

   tmp1	=	MAX( MAX(BOX[system].lx, BOX[system].ly), BOX[system].lz);
   tmp3	=	MIN( MIN(BOX[system].lx, BOX[system].ly), BOX[system].lz);
   if ( fabs(BOX[system].lx-tmp1) > ZERO && fabs(BOX[system].lx-tmp3) > ZERO) {
      tmp2	=	BOX[system].lx;
   }
   else if ( fabs(BOX[system].ly-tmp1) > ZERO && fabs(BOX[system].ly-tmp3) > ZERO) {
      tmp2	=	BOX[system].ly;
   }
   else if ( fabs(BOX[system].lz-tmp1) > ZERO && fabs(BOX[system].lz-tmp3) > ZERO) {
      tmp2	=	BOX[system].lz;
   }
   asphericity	=	tmp1 * tmp1 - 0.5*(tmp2*tmp2 + tmp3*tmp3);
   asphericity	/=	(tmp1 + tmp2 + tmp3) * (tmp1 + tmp2 + tmp3) / 9;

   return	asphericity;
}

//----------------------------------------------//
//	function: get_property()		//
//						//
//	extract value for atomic property  	//
//	Add: 1/17/2018				//
//----------------------------------------------//

void *get_atom_property(molstruct *moli, int i, char *property)
{
   if (!strcmp(property, "cna2")) {
      return (void *)&moli[i].cna2;
   }
   return NULL;
}

//----------------------------------------------//
//	function: ave_property()		//
//						//
//	calculate local average of property	//
//	Add: 1/17/2018				//
//----------------------------------------------//

void ave_property(char *property, float rcutoff)
{
   molstruct *moli, *molj;
   int 	 i, j, n, nsites, system;
   int   nave;
   float rcutoff2 = rcutoff * rcutoff, r2;
   short *var_short;
   float var;

   static int  init=1;
   static float *L;
#ifdef CELL_LIST
   long		k, l;
   static cellstruct	*celli, *cellk;
#endif

   if (init) {
      init	=	0;

      nsites = 0;
      for (moli=mol; moli<mol+NMOLS; moli++) {
	 for (i=0; i<moli->nsites; i++) {
	    nsites ++;
	 }
      }

      L = (float *) calloc(nsites, sizeof(float));	// only for one system now
      if (L==NULL)
         Exit("sample", "ave_property", "out of memory");
   }

   // Perform local averaging
   
   n = 0;

   for (moli=mol; moli<mol+NMOLS; moli++) {
      if ( (system=moli->box) >= 0 ) {
         for (i=0; i<moli->nsites; i++) {

	    if (!strcmp(property, "cna2")) {
	       //var_short = (short *)get_atom_property(molj, j, property);
	       var = *(short *)get_atom_property(moli, i, property);
	    }
	    else {
	       var = 0.0;
	    }
	    L[n] = var;
	    nave = 1;

//#ifdef CELL_LIST
	    /*
	    celli	=	moli->cell[i];
	    for (l=0; l<celli->nneigh; l++) {
	       cellk	=	celli->neigh[l];
	       for (k=0; k<cellk->nsites; k++) {
		  if (molj=cellk->mol[k]) {
		     j	=	cellk->molsite[k];
		     if (j!=0 && j!=molj->nsites-1 && (moli!=molj || i!=j)) {
	    */
//#else
            for (molj=mol; molj<mol+NMOLS; molj++) {
               if ( molj->box == system ) {
                  for (j=0; j<molj->nsites; j++) {
                     if ( (moli!=molj) || j!=i ) {
//#endif
                        r2 = DistSQ(moli->p[i], molj->p[j], system);
                        if (r2 < rcutoff2) {
			   if (!strcmp(property, "cna2")) {
			      var = *(short *)get_atom_property(molj, j, property);
			   }
			   else {
			      var = 0.0;
			   }
			   L[n] += var;
			   nave ++;
			}
		     }
		  }
	       }
	    }
	    L[n] /= nave;
	    n ++;

		     }
		  }
	       }

   n = 0;
   for (moli=mol; moli<mol+NMOLS; moli++) {
      if ( (system=moli->box) >= 0 ) {
	 for (i=0; i<moli->nsites; i++) {
	    moli->cna2[i] = (short)L[n];
	 }
      }
   }
}


