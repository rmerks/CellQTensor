#include <time.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

#ifndef TYPES_DECLARATION

typedef struct {
	
	long index;
        long boundary;
        long adhesion;
        double s;
        double q1;
        double q2;
	
	double x;
	double y;
        double lam;
        double kappa;
        double theta;

} Vertex;

typedef struct {
	
	long index;
        long type;
	
	double x;
	double y;
        double fx;
        double fy;

} Forces;

typedef struct {
	
	long index;
        long bulk;
        long sharedbulk;
	
	double x;
	double y;
        double phi;
        double q1;
        double q2;
        double s;
        long loworder;
       

} Grid;

typedef struct {
	
	double phi;
	double phi_exp;
	long phi_exp_degrees;
	long counts;
	double coherency;
        long bulk;
        long bulk_exp;
        double q1;
        double q2;
        double s;
        long loworder;
        
} Grid_exp;

typedef struct {
	
        long bulk;
        double q1;
        double q2;
        
} Grid_old;

typedef struct {
	
	long border;
	double x;
	double y;

} Border_exp;

typedef struct {
	
	long index;
        long edgecross;
	
	double x;
	double y;
        double psi;

} Lines;

typedef struct {

  long d;
  long h;
  long m;
  long s;

} CPU_Time;

#endif
#define TYPES_DECLARATION
