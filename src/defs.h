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

/*************************************************************************
    GLOBAL ARRAYS SECTION 
*************************************************************************/

int SN1 = N1+4;
int SN2 = N2+4;
int SN3 = N3+2;

#if 0
double CoolVal[N3+2][N1+4][N2+4][NCOOLVAL]; /*array with neutrino cooling*/

double   a_p[N3+2][N1+4][N2+4][NPR] ;  /* space for primitive vars */
double  a_dq[N3+2][N1+4][N2+4][NPR] ;  /* slopes */
double  a_F1[N3+2][N1+4][N2+4][NPR] ;  /* fluxes */
double  a_F2[N3+2][N1+4][N2+4][NPR] ;  /* fluxes */
double  a_F3[N3+2][N1+4][N2+4][NPR] ;  /* fluxes */
double  a_ph[N3+2][N1+4][N2+4][NPR] ;  /* half-step primitives */
char    a_pflag[N3+2][N1+4][N2+4];	/* identifies failure points */
#endif

/* for debug */
//double psave[N3][N1][N2][NPR] ;
double fimage[NIMG][N1*N2];
int    failimage[5][N1*N2];

/* grid functions */
#if 0
double a_conn[N3+2][N1+4][N2+4][NDIM][NDIM][NDIM] ;
double a_gcon[N3+2][N1+4][N2+4][NPG][NDIM][NDIM] ;
double a_gcov[N3+2][N1+4][N2+4][NPG][NDIM][NDIM] ;
double a_gdet[N3+2][N1+4][N2+4][NPG] ;
#endif
/*
double (*   p)[N1+4][N2+4][NPR] ;
double (*  dq)[N1+4][N2+4][NPR] ;
double (*  F1)[N1+4][N2+4][NPR] ;
double (*  F2)[N1+4][N2+4][NPR] ;
double (*  F3)[N1+4][N2+4][NPR] ;
double (*  ph)[N1+4][N2+4][NPR] ;
char   (*  pflag)[N1+4][N2+4];

double (* conn)[N1+4][N2+4][NDIM][NDIM][NDIM] ;
double (* gcon)[N1+4][N2+4][NPG][NDIM][NDIM] ;
double (* gcov)[N1+4][N2+4][NPG][NDIM][NDIM] ;
double (* gdet)[N1+4][N2+4][NPG] ;

double (* CoolVal)[N1+4][N2+4][NCOOLVAL] ;
double (*   A)[N1+4][N2+4] ;
*/

double *   _p ;
double *  _dq ;
double *  _F1 ;
double *  _F2 ;
double *  _F3 ;
double *  _ph ;
char   *  _pflag ;

double * _conn ;
double * _gcon ;
double * _gcov ;
double * _gdet ;
double * _coords ;
double * _sf ;

double * _CoolVal ;
double * _A ;
double * _A2 ;
double * _A3 ;


#if(DO_FONT_FIX)
double Katm[N1+4];
#endif

/*************************************************************************
    GLOBAL VARIABLES SECTION 
*************************************************************************/
/* physics parameters */
double a ;
double gam ;

/* numerical parameters */
double Rin,Rout,hslope,R0 ;
double Rbr=-1.;
double cyl_R0,cyl_Th0;
double cour ;
double dV,dx[NDIM],startx[NDIM] ;
double dt ;
double t,tf ;
int hstart=1,hstop=N3,istart=2,istop=N1+1;
//const int jstart=2,jstop=N2+1;
int rank=0,nprocs=1, hrank=0,hsize=1, irank=0,isize=1;
//const int jrank=0,jsize=1;
int rank_prevh=-1, rank_nexth=-1, rank_previ=-1, rank_nexti=-1;
int dh_rank_start=0,di_rank_start=0,dj_rank_start=0;

int hcurr,icurr,jcurr,pcurr,ihere,jhere,phere ;
double dminarg1,dminarg2 ;
int nstep ;
//double fval1,fval2;

/* output parameters */
double DTd ;
double DTl ;
double DTi ;
int    DTr ;
int    dump_cnt ;
int    image_cnt ;
int    rdump_cnt ;
int    nstroke ;

/* global flags */
int    failed ;
int    lim ;
double defcon ;

/* diagnostics */
double mdot = 0. ;
double edot = 0. ;
double ldot = 0. ;

double mdot_in = 0. ;
double mdot_out = 0. ;
double mass_start = 0. ;
double mass_curr = 0. ;
double mass_in = 0. ;
double mass_out = 0. ;

#ifdef MPI_USED
	MPI_Status MPIStatus;

	MPI_Datatype slice_iNPR1_start;
	MPI_Datatype slice_iNPR1_startP;
	MPI_Datatype slice_iNPR1_stop;
	MPI_Datatype slice_iNPR1_stopN;
/*
	MPI_Datatype slice_iNPR2_start;
	MPI_Datatype slice_iNPR2_startP;
	MPI_Datatype slice_iNPR2_stop;
	MPI_Datatype slice_iNPR2_stopN;
*/
	MPI_Datatype slice_i1_start;
	MPI_Datatype slice_i1_startP;
	MPI_Datatype slice_i1_stop;
	MPI_Datatype slice_i1_stopN;

	MPI_Datatype slice_i1c_start;
	MPI_Datatype slice_i1c_startP;
	MPI_Datatype slice_i1c_stop;
	MPI_Datatype slice_i1c_stopN;


	MPI_Datatype slice_hNPR1_start;
	MPI_Datatype slice_hNPR1_startP;
	MPI_Datatype slice_hNPR1_stop;
	MPI_Datatype slice_hNPR1_stopN;

	MPI_Datatype slice_h1_start;
	MPI_Datatype slice_h1_startP;
	MPI_Datatype slice_h1_stop;
	MPI_Datatype slice_h1_stopN;

#if N3 > 1
	MPI_Datatype slice_hopNPR2_m2[N3];
	MPI_Datatype slice_hopNPR2_n1[N3];

	int slice_h_rankop[N3];
	int slice_h_hop[N3];
#endif

#endif
