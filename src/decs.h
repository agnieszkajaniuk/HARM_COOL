/***********************************************************************************
    Copyright 2006 Charles F. Gammie, Jonathan C. McKinney, Scott C. Noble, 
                   Gabor Toth, and Luca Del Zanna

                        HARM  version 1.0   (released May 1, 2006)

    This file is part of HARM.  HARM is a program that solves hyperbolic 
    partial differential equations in conservative form using high-resolution
    shock-capturing techniques.  This version of HARM has been configured to 
    solve the relativistic magnetohydrodynamic equations of motion on a 
    stationary black hole spacetime in Kerr-Schild coordinates to evolve
    an accretion disk model. 

    You are morally obligated to cite the following two papers in his/her 
    scientific literature that results from use of any part of HARM:

    [1] Gammie, C. F., McKinney, J. C., \& Toth, G.\ 2003, 
        Astrophysical Journal, 589, 444.

    [2] Noble, S. C., Gammie, C. F., McKinney, J. C., \& Del Zanna, L. \ 2006, 
        Astrophysical Journal, 641, 626.

   
    Further, we strongly encourage you to obtain the latest version of 
    HARM directly from our distribution website:
    http://rainman.astro.uiuc.edu/codelib/


    HARM is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    HARM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HARM; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

***********************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifdef MPI_USED
#include "mpi.h"
#endif

/*************************************************************************
      COMPILE-TIME PARAMETERS : 
*************************************************************************/
/** here are the few things that we change frequently **/

#define N1       (256)        /* number of physical zones in X1-direction */
#define N2       (256)        /* number of physical zones in X2-direction */
#define N3       (1)        /* number of physical zones in X3-direction */

#define BL       (1)          /* whether or not to use BL coords */

#define NG       (2)


#define COORDS_GAMMIE            0
#define COORDS_CYLINDRIFY_GAMMIE 2

#define COORDS         COORDS_GAMMIE //Use original coords
//#define COORDS         COORDS_CYLINDRIFY_GAMMIE



//add mass in the drift frame (=1) instead of fluid frame (=0)
#define DRIFT_FLOOR (1)
#define JBOUNDAXISYMMETRIC (0)     //Enable for 3D
//#define JBOUNDAXISYMMETRIC (N3>1)     //Enable for 3D
#define THETAZEROPI                 (0) /* Theta [0, PI] instead of [delta, PI - delta] */
#define DOKTOT   0

/// whether to flip gdet sign over coordinate singularities
/// completely generally, this should be 1 so that \detg is smooth across the axis.  
/// So then standard boundary conditions on primitives give correct non-kinked behavior through polar axis (including for ram pressure flux term).
#define FLIPGDETAXIS 0

#define GHOSTZONESIGNCHANGEINMETRIC (0) /*Sign change for V and B components included in metric */ 
#define FIXUPAVG                    (1) /*whether or not average in fixup(), zero if THETAZEROPI=1 */
/* how many cells near the poles to stabilize, choose 0 for no stabilization */
#if N3>1
#define POLEFIX (2)                         
#else
#define POLEFIX (0)
#endif

#if(DOKTOT)
#  define NPR        (9)
#else
#  define NPR        (8)      /* number of primitive variables */
#endif

#define NDIM       (4)        /* number of total dimensions.  Never changes */
#define NPG        (5)        /* number of positions on grid for grid functions */
#define COMPDIM    (2)        /* number of non-trivial spatial dimensions used in computation */

#define NIMG       (4)        /* Number of types of images to make, kind of */

/* whether or not to use Cool EOS */
#define COOL_EOS  (1)
#define COOL_EOS_SIMPLECS  (0)


#define TRACERSDUMP (1)

/* HDF5 dumping selection */
#define HDFDUMP (1)
#if HDFDUMP
	#define HDFDUMP_Rho        (1)
	#define HDFDUMP_Energy     (1)
	#define HDFDUMP_DivB       (1)
	#define HDFDUMP_U1U2U3gdet (1)
	#define HDFDUMP_B1B2B3gdet (1)
	#define HDFDUMP_Ucon       (1)
	#define HDFDUMP_Ucov       (1)
	#define HDFDUMP_Bcon       (1)
	#define HDFDUMP_Bcov       (1)
	#define HDFDUMP_Metric     (1)
    #if HDFDUMP_Metric
        #define HDFDUMP_gcov      (1)
        #define HDFDUMP_gcon      (1)
    #endif
        #define HDFDUMP_mtrc_trns  (1)
#endif

#undef  PHSLICE 
#undef  THSLICE


/* whether or not to use Font's  adiabatic/isothermal prim. var. inversion method: */
#define DO_FONT_FIX (1)

/* whether or not to rescale primitive variables before interpolating them for flux/BC's */
#define RESCALE     (0)

/** FIXUP PARAMETERS, magnitudes of rho and u, respectively, in the floor : **/
#define RHOMIN	(1.e-6)
#define UUMIN	(1.e-8)
#define RHOMINLIMIT (1.e-20)
#define UUMINLIMIT  (1.e-20)
#define POWRHO (2.0)

#define FLOORFACTOR (1.)
#define BSQORHOMAX (50.*FLOORFACTOR)
#define BSQOUMAX (2500.*FLOORFACTOR)
#define UORHOMAX (50.*FLOORFACTOR)

/* A numerical convenience to represent a small non-zero quantity compared to unity:*/
#define SMALL	(1.e-20)

/* Max. value of gamma, the lorentz factor */
#define GAMMAMAX (50.)

/* maximum fractional increase in timestep per timestep */
#define SAFE	(1.3)


#define COORDSINGFIX 1
// whether to move polar axis to a bit larger theta
// theta value where singularity is displaced to
#define SINGSMALL (1.E-20)


/* I/O format strings used herein : */
#define FMT_DBL_OUT "%32.22e"
#define FMT_INT_OUT "%10d"



/*************************************************************************
    MNEMONICS SECTION 
*************************************************************************/
/* boundary condition mnemonics */
#define OUTFLOW	(0)
#define SYMM	(1)
#define ASYMM	(2)
#define FIXED	(3)

/* mnemonics for primitive vars; conserved vars */
#define RHO	(0)	
#define UU	(1)
#define U1	(2)
#define U2	(3)
#define U3	(4)
#define B1	(5)
#define B2	(6)
#define B3	(7)
#define KTOT    (8)

/* mnemonics for dimensional indices */
#define TT	(0)     
#define RR	(1)
#define TH	(2)
#define PH	(3)

/* mnemonics for centering of grid functions */
#define FACE1	(0)	
#define FACE2	(1)
#define CORN	(2)
#define CENT	(3)
#define FACE3	(4)

/* mnemonics for slope limiter */
#define MC	(0)
#define VANL	(1)
#define MINM	(2)

/* mnemonics for diagnostic calls */
#define INIT_OUT	(0)
#define DUMP_OUT	(1)
#define IMAGE_OUT	(2)
#define LOG_OUT		(3)
#define FINAL_OUT	(4)

/* Directional Mnemonics */
// -------------> r
// |         3    
// |        1-0   
// |         2    
// v            
// theta      
#define X1UP    (0)
#define X1DN    (1)
#define X2UP    (2)
#define X2DN    (3)


/* failure modes */
#define FAIL_UTOPRIM        (1)
#define FAIL_VCHAR_DISCR    (2)
#define FAIL_COEFF_NEG	    (3)
#define FAIL_COEFF_SUP	    (4)
#define FAIL_GAMMA          (5)
#define FAIL_METRIC         (6)


/* For rescale() operations: */
#define FORWARD 1
#define REVERSE 2

#define NCOOLVAL        (7)        /* Number of cooling valuses: Qnu, P, Hd*/
#define COOL_LAMBDA_SIM (0)
#define COOL_P          (1)
#define COOL_Hd         (2)
#define COOL_Qnu        (3)
#define COOL_Tau        (4)
#define COOL_T          (5)
#define COOL_Yee        (6)



/*************************************************************************
    GLOBAL ARRAY SECTION 
*************************************************************************/

extern int SN1;
extern int SN2;
extern int SN3;

#if 0
extern double CoolVal[N3+2][N1+4][N2+4][NCOOLVAL]; /*array with cooling*/
extern double  a_p[N3+2][N1+4][N2+4][NPR] ;	/* space for primitive vars */
extern double  a_dq[N3+2][N1+4][N2+4][NPR] ;  /* slopes */
extern double  a_F1[N3+2][N1+4][N2+4][NPR] ;	/* fluxes */
extern double  a_F2[N3+2][N1+4][N2+4][NPR] ;	/* fluxes */
extern double  a_F3[N3+2][N1+4][N2+4][NPR] ;	/* fluxes */
extern double  a_ph[N3+2][N1+4][N2+4][NPR] ;	/* half-step primitives */
extern char    a_pflag[N3+2][N1+4][N2+4];	/* identifies failure points */
#endif
/* for debug */
//extern double psave[N3][N1][N2][NPR] ;
extern double fimage[NIMG][N1*N2];
extern int    failimage[5][N1*N2];

/* grid functions */
#if 0
extern double a_conn[N3+2][N1+4][N2+4][NDIM][NDIM][NDIM] ;
extern double a_gcon[N3+2][N1+4][N2+4][NPG][NDIM][NDIM] ;
extern double a_gcov[N3+2][N1+4][N2+4][NPG][NDIM][NDIM] ;
extern double a_gdet[N3+2][N1+4][N2+4][NPG] ;
#endif
/*
extern double (*   p)[N1+4][N2+4][NPR] ;
extern double (*  dq)[N1+4][N2+4][NPR] ;
extern double (*  F1)[N1+4][N2+4][NPR] ;
extern double (*  F2)[N1+4][N2+4][NPR] ;
extern double (*  F3)[N1+4][N2+4][NPR] ;
extern double (*  ph)[N1+4][N2+4][NPR] ;
extern char   (*  pflag)[N1+4][N2+4];

extern double (* conn)[N1+4][N2+4][NDIM][NDIM][NDIM] ;
extern double (* gcon)[N1+4][N2+4][NPG][NDIM][NDIM] ;
extern double (* gcov)[N1+4][N2+4][NPG][NDIM][NDIM] ;
extern double (* gdet)[N1+4][N2+4][NPG] ;

extern double (* CoolVal)[N1+4][N2+4][NCOOLVAL] ;
extern double (*   A)[N1+4][N2+4] ;
*/

extern double *   _p ;
extern double *  _dq ;
extern double *  _F1 ;
extern double *  _F2 ;
extern double *  _F3 ;
extern double *  _ph ;
extern char   *  _pflag ;

extern double * _conn ;
extern double * _gcon ;
extern double * _gcov ;
extern double * _gdet ;
extern double * _coords ;
extern double * _sf ;  /*Sign Flip */

extern double * _CoolVal ;
extern double * _A ;
extern double * _A2 ;
extern double * _A3 ;

#if(DO_FONT_FIX)
extern double Katm[N1+4];
#endif


/*************************************************************************
    GLOBAL VARIABLES SECTION 
*************************************************************************/
/* physics parameters */
extern double a ;
extern double gam ;

/* numerical parameters */
extern double Rin,Rout,hslope,R0 ;
extern double Rbr, cyl_R0, cyl_Th0;
extern double cour ;
extern double dV,dx[NDIM],startx[NDIM] ;
extern double dt ;
extern double t,tf ;
extern double x1curr,x2curr ;
extern int nstep ;

/* output parameters */
extern double DTd ;
extern double DTl ;
extern double DTi ;
extern int    DTr ;
extern int    dump_cnt ;
extern int    image_cnt ;
extern int    rdump_cnt ;
extern int    nstroke ;

/* global flags */
extern int failed ;
extern int lim ;
extern double defcon ;

/* diagnostics */
extern double mdot ;
extern double edot ;
extern double ldot ;

extern double mdot_in ;
extern double mdot_out ;
extern double mass_start ;
extern double mass_curr ;
extern double mass_in ;
extern double mass_out ;

/* set global variables that indicate current local metric, etc. */
extern int hcurr,icurr,jcurr,pcurr ;
struct of_geom {
	double gcon[NDIM][NDIM] ;
	double gcov[NDIM][NDIM] ;
	double g ;
} ;

struct of_state {
	double ucon[NDIM] ;
	double ucov[NDIM] ;
	double bcon[NDIM] ;
	double bcov[NDIM] ;
} ;





/*************************************************************************
    MACROS
*************************************************************************/
/* Active zone definition */
extern int hstart,hstop, istart,istop;
const int jstart=2,jstop=N2+1;
extern int rank, nprocs, hrank,hsize, irank,isize;
const int jrank=0,jsize=1;
extern int rank_prevh, rank_nexth, rank_previ, rank_nexti;
extern int dh_rank_start,di_rank_start,dj_rank_start;

#ifdef MPI_USED
extern MPI_Status MPIStatus;

extern MPI_Datatype slice_iNPR1_start;
extern MPI_Datatype slice_iNPR1_startP;
extern MPI_Datatype slice_iNPR1_stop;
extern MPI_Datatype slice_iNPR1_stopN;

/*
extern MPI_Datatype slice_iNPR2_start;
extern MPI_Datatype slice_iNPR2_startP;
extern MPI_Datatype slice_iNPR2_stop;
extern MPI_Datatype slice_iNPR2_stopN;
*/

extern MPI_Datatype slice_i1_start;
extern MPI_Datatype slice_i1_startP;
extern MPI_Datatype slice_i1_stop;
extern MPI_Datatype slice_i1_stopN;

extern MPI_Datatype slice_i1c_start;
extern MPI_Datatype slice_i1c_startP;
extern MPI_Datatype slice_i1c_stop;
extern MPI_Datatype slice_i1c_stopN;

extern MPI_Datatype slice_hNPR1_start;
extern MPI_Datatype slice_hNPR1_startP;
extern MPI_Datatype slice_hNPR1_stop;
extern MPI_Datatype slice_hNPR1_stopN;

extern MPI_Datatype slice_h1_start;
extern MPI_Datatype slice_h1_startP;
extern MPI_Datatype slice_h1_stop;
extern MPI_Datatype slice_h1_stopN;

#if N3 > 1
extern MPI_Datatype slice_hopNPR2_m2[N3];
extern MPI_Datatype slice_hopNPR2_n1[N3];

extern int slice_h_rankop[N3];
extern int slice_h_hop[N3];
#endif

#endif

/* loop over all active zones */
#if N3==1
#define ZLOOP for(h=0,i=istart;i<=istop;i++)for(j=jstart;j<=jstop;j++)
#else
#define ZLOOP for(h=hstart;h<=hstop;h++)for(i=istart;i<=istop;i++)for(j=jstart;j<=jstop;j++)
#endif
/* loop over all active zones */
#define IMAGELOOP for(j=jstart;j<=jstop;j++)for(i=istart;i<=istop;i++)

/* specialty loop */
/*
#define ZSLOOP(istart,istop,jstart,jstop) \
        for(i=istart;i<=istop;i++)\
	for(j=jstart;j<=jstop;j++)
*/

#if N3==1
#define ZSLOOP(hstartM,hstopM, istartM,istopM,jstartM,jstopM) \
        for(h=0,i=(irank==0?istart+istartM:istart);i<=(irank==isize-1?istop+istopM:istop);i++)\
        for(j=(jrank==0?jstart+jstartM:jstart);j<=(jrank==jsize-1?jstop+jstopM:jstop);j++)
#else
#define ZSLOOP(hstartM,hstopM, istartM,istopM,jstartM,jstopM) \
        for(h=(hrank==0?hstart+hstartM:hstart);h<=(hrank==hsize-1?hstop+hstopM:hstop);h++)\
        for(i=(irank==0?istart+istartM:istart);i<=(irank==isize-1?istop+istopM:istop);i++)\
        for(j=(jrank==0?jstart+jstartM:jstart);j<=(jrank==jsize-1?jstop+jstopM:jstop);j++)
#endif


#if N3==1
#define ZSLOOPH(hstartM,hstopM) \
        for(h=0;h<=0;h++)
#else
#define ZSLOOPH(hstartM,hstopM) \
        for(h=(hrank==0?hstart+hstartM:hstart);h<=(hrank==hsize-1?hstop+hstopM:hstop);h++)
#endif

#define ZSLOOPI(istartM,istopM) \
        for(i=(irank==0?istart+istartM:istart);i<=(irank==isize-1?istop+istopM:istop);i++)

#define ZSLOOPJ(jstartM,jstopM) \
        for(j=(jrank==0?jstart+jstartM:jstart);j<=(jrank==jsize-1?jstop+jstopM:jstop);j++)

    	
/* loop over Primitive variables */
#define PLOOP  for(k=0;k<NPR;k++)
/* loop over all Dimensions; second rank loop */
#define DLOOP  for(j=0;j<NDIM;j++) for(k=0;k<NDIM;k++)
/* loop over all Dimensions; first rank loop */
#define DLOOPA for(j=0;j<NDIM;j++)
/* loop over all Space dimensions; second rank loop */
#define SLOOP  for(j=1;j<NDIM;j++) for(k=1;k<NDIM;k++)
/* loop over all Space dimensions; first rank loop */
#define SLOOPA for(j=1;j<NDIM;j++)
/* loop over Cooling variables */
#define CLOOP  for(k=0;k<NCOOLVAL;k++)


//extern double fval1,fval2;
#define MY_MIN(fval1,fval2) ( ((fval1) < (fval2)) ? (fval1) : (fval2))
#define MY_MAX(fval1,fval2) ( ((fval1) > (fval2)) ? (fval1) : (fval2))

#define delta(i,j) ( (i == j) ? 1. : 0.)
#define dot(a,b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3]) 


/*************************************************************************
    FUNCTION DECLARATIONS 
*************************************************************************/

//neutrino functions
//void read_neutrino_table(void);
//void neutrino_cooling_func(double rho, double u, double r, double u_eq, double ro_eq, double *val);
//void cool_down(double prim[][N2 + 4][NPR], double Dt);



double bl_gdet_func(double r, double th) ;
double bsq_calc(double *pr, struct of_geom *geom) ;
int    gamma_calc(double *pr, struct of_geom *geom, double *gamma) ;
double gdet_func(double lgcov[NDIM][NDIM]) ;
double mink(int j, int k) ; 
double ranc(int seed) ;
double slope_lim(double y1, double y2, double y3) ;
void ut_calc_3vel(double *vcon, struct of_geom *geom, double *ut);

int restart_init(void) ;

void area_map(int h, int i, int j, double *_prim) ;
void bcon_calc(double *pr, double *ucon, double *ucov, double *bcon) ;
void blgset(int h, int i, int j, struct of_geom *geom);
void bl_coord(double *X, double *r, double *th, double *ph) ;
void bl_coord_vec(double *X, double *V);
void bl_gcon_func(double r, double th, double gcov[NDIM][NDIM]) ;
void bl_gcov_func(double r, double th, double gcov[NDIM][NDIM]) ;
void bound_prim(double *_pr) ;
void dxdxp_func(double *X, double dxdxp[NDIM][NDIM]);
void conn_func(double *X, struct of_geom *geom, double lconn[NDIM][NDIM][NDIM]) ;
void coord(int h, int i, int j, int loc, double *X) ;
void diag(int call_code) ;
void diag_flux(double *_F1, double *_F2, double *_F3) ;
void diag_mass(double Dt);

void get_coord_vec(int h, int i, int j, double *V);
void get_coord_r(int h, int i, int j, double *r);
void get_coord(int h, int i, int j, double *r, double *theta, double *phi);

#if HDFDUMP
void dumph5(int dump_cnt) ;
#else
void dump(FILE *fp) ;
#endif

void fail(int fail_type) ;
void fixup(double *_pv) ;
void fixup1zone( int h, int i, int j, double prim[NPR] ) ;
void fixup_utoprim( double *_pv )  ;
void set_Katm(void);
int  get_G_ATM( double *g_tmp );
void fix_flux(double *_F1, double *_F2, double *_F3) ;
void flux_ct(double *_F1,double *_F2, double *_F3) ;
void gaussj(double **tmp, int n, double **b, int m) ;
void gcon_func(double lgcov[NDIM][NDIM], double lgcon[NDIM][NDIM]) ;
void gcov_func(double *X, double lgcov[NDIM][NDIM]) ;
void get_geometry(const int h, const int i, const int j, const int loc, struct of_geom *geom) ;
void get_state(double *pr, struct of_geom *geom, struct of_state *q) ;
void image_all( int image_count ) ;
void init(void) ;
void init_bondi(void) ;
void lower(double *a, struct of_geom *geom, double *b) ;
void ludcmp(double **a, int n, int *indx, double *d) ;
void mhd_calc(double *pr, int dir, struct of_state *q, double *mhd)  ;
void mhd_calc_CoolEOS(double r, double u, double P, int dir, struct of_state *q, double *mhd) ;

void misc_source(double *ph, const int hh, const int ii, const int jj, struct of_geom *geom,
		                struct of_state *q, double *dU, double Dt) ;
void primtoflux(double *pa, struct of_state *q, const int dir, struct of_geom *geom,
			double *fl) ;
void primtoU(double *p, struct of_state *q, struct of_geom *geom, double *U);
void raise(double *v1, struct of_geom *geom, double *v2) ;
void rescale(double *pr, int which, int dir, int hh, int ii, int jj, int face, 
			struct of_geom *geom) ;
void restart_write(void) ;
void restart_read(FILE *fp) ;
void set_arrays(void) ;
void set_grid(void) ;
void set_points(void) ;
void step_ch(void) ;
void source(double *pa, struct of_geom *geom, int hh, int ii, int jj, double *Ua,double Dt) ;
void timestep(void) ;
void u_to_v(double *pr, int i, int j) ;
void ucon_calc(double *pr, struct of_geom *geom, double *ucon) ;
void ucon_to_utcon(double *ucon,struct of_geom *geom, double *utcon);
void usrfun(double *pr,int n,double *beta,double **alpha) ;
void Utoprim(double *Ua, struct of_geom *geom, double *pa) ;
int Utoprim_2d(double U[NPR], double gcov[NDIM][NDIM], double gcon[NDIM][NDIM], 
               double gdet, double prim[NPR]);
int Utoprim_5d(double U[NPR], double gcov[NDIM][NDIM], double gcon[NDIM][NDIM], 
               double gdet, double prim[NPR]);
int Utoprim_1dvsq2fix1(double U[NPR], double gcov[NDIM][NDIM], double gcon[NDIM][NDIM], 
               double gdet, double prim[NPR], double K );
int Utoprim_1dfix1(double U[NPR], double gcov[NDIM][NDIM], double gcon[NDIM][NDIM], 
               double gdet, double prim[NPR], double K );

void vchar(double *pr, struct of_state *q, struct of_geom *geom,
		int dir, double *cmax, double *cmin) ;

int invert_matrix( double A[NDIM][NDIM], double Ainv[NDIM][NDIM] );
int LU_decompose( double A[NDIM][NDIM], int permute[] );
void LU_substitution( double A[NDIM][NDIM], double B[], int permute[] );

void init_entropy();
void compute_ktot(double *_pi,double *_prh, double *_pr, int h,int i, int j, double Dt);

