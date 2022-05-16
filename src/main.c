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
#include <sys/time.h>
#include <fcntl.h>
#include <errno.h>
#include <time.h>
                  
#include "decs.h"
#include "defs.h"

#ifdef MPI_USED
#include "mpi.h"
#endif

#if( COOL_EOS )
#include "cooleos.h"
#endif

#if TRACERSDUMP
#include "tracers.h"
#endif

inline unsigned long long getTicks()
{
	struct timeval tv;
	unsigned long long ticks;
		
	if (gettimeofday(&tv, 0) == 0)
	{
		ticks = tv.tv_sec;
		ticks *= 1000;
		ticks += tv.tv_usec / 1000;
		return ticks; 
	}
	else return 0;
}

#ifdef MPI_USED
void set_MPI(int argc,char *argv[]);
void set_MPI_Slices(void);
#endif

double tdump,timage,tlog=0. ;
int nfailed = 0 ;

/*****************************************************************/
/*****************************************************************
   main():
   ------

     -- Initializes, time-steps, and concludes the simulation. 
     -- Handles timing of output routines;
     -- Main is main, what more can you say.  

-*****************************************************************/
int main(int argc,char *argv[])
{
	time_t trestore ;
	unsigned long long lastClock = getTicks();
	char stdbuf[128];
    int  restart_save;

	nstep = 0 ;

#if( COOL_EOS )
	SetupThreads_CoolEOS ();
#endif

#ifdef MPI_USED
    MPI_Init(&argc, &argv);          //Start MPI
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );  //Pobranie numeru procesu, 0 - proces glowny
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs );  //Pobranie liczby procesow
#endif

	if(rank == 0)  
    {
		freopen ("stdout.txt","w",stdout);
        freopen ("stderr.txt","w",stderr);
	}
    else
	{
/*
		sprintf (stdbuf, "stdout.txt.%d",rank);
		freopen (stdbuf,"w",stdout);
		sprintf (stdbuf, "stderr.txt.%d",rank);
		freopen (stdbuf,"w",stderr);
*/
		freopen ("stdout.txt.all","w",stdout);
        freopen ("stderr.txt.all","w",stderr);
	}

    //Ustawienie zerowej dlugosci buforow std
//    if(nprocs == 1)
    	setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

	nstep = 0 ;
	
#ifdef MPI_USED
	set_MPI(argc, argv);
#endif



	if(rank == 0)	
		printf("Harm code A&I Janiuk mod v. 4.8.6 - 16.05.2022\n");

	printf("rank %d RHO=[%d, %d]\n", rank, istart - (irank==0 ? 2 : 1) + di_rank_start, istop - (irank==0 ? 2 : 1) + di_rank_start);
	printf("rank %d TH=[%d, %d]\n", rank, jstart - (jrank==0 ? 2 : 1) + dj_rank_start, jstop - (jrank==0 ? 2 : 1) + dj_rank_start);
	printf("rank %d PH=[%d, %d]\n", rank, N3 == 1 ? 0 : hstart - 1 + dh_rank_start, N3 == 1 ? 0 : hstop - 1 + dh_rank_start);
	printf("rank %d irank=%d, hrank=%d\n", rank, irank, hrank);
	printf("rank %d rank_prevh=%d, rank_nexth=%d\n", rank, rank_prevh, rank_nexth);
	printf("rank %d rank_previ=%d, rank_nexti=%d\n", rank, rank_previ, rank_nexti);


	/* Perform Initializations, either directly or via checkpoint */
	system("mkdir -p dumps images cache tracers/in tracers/out tracers/jet");  

    int restart_status = restart_init(); 
	if(!restart_status) { 
		int i,j,h;

	    init() ;
#if( DOKTOT )
        init_entropy();
#endif

		{
		    double (*   p)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_p);
		    double (*   gdet)[SN2][NPG] = (double (*) [SN2][NPG])(_gdet);

			mass_start = 0.;         
			ZLOOP {
				mass_start += p[h][i][j][RHO]*gdet[i][j][CENT]*dx[RR]*dx[TH]*dx[PH];
    		}
		}
#ifdef MPI_USED
		MPI_Allreduce(MPI_IN_PLACE,&mass_start,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif

		diag_mass(0.);

#if TRACERSDUMP
		init_tracers();
#endif

	} 

	/* do initial diagnostics */
	diag(INIT_OUT) ;

    if(!restart_status) {
#if TRACERSDUMP
   	    dump_tracers();
#endif
	    tdump = t+DTd ;
	    timage = t+DTi ;
	    tlog = t+DTl ;
    }
	trestore = time(NULL)+DTr ;

	defcon = 1. ;
	
    	
	while(t < tf) {

		/* step variables forward in time */
		nstroke = 0 ;
		step_ch() ;

		if(rank == 0) {
			fprintf(stdout,"t=%10.5g dt=%10.5g nstep=%8d %10.5g calcTime=%.3lf val=" FMT_DBL_OUT "\n",t,dt,nstep,
					nstroke/(2.*N1*N2), ((double)(getTicks()-lastClock))/1000., mass_curr) ;
			fflush(stdout);
			lastClock = getTicks();
		}

		/* Handle output frequencies: */
		
		if(t >= tdump) {
			diag(DUMP_OUT) ;
			tdump += DTd ;
		}
		if(t >= timage) {
			diag(IMAGE_OUT) ;
			timage += DTi ;
		}
		if(t >= tlog) {
			diag(LOG_OUT) ;
#if TRACERSDUMP
			dump_tracers();
#endif
			tlog += DTl ;
		}

	
		/* restart dump */
		nstep++ ;

        restart_save = 0;

		if(rank == 0) {
			restart_save = (time(NULL) >= trestore);
		}

#ifdef MPI_USED
		MPI_Allreduce(MPI_IN_PLACE,&restart_save,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif

		if(restart_save) {
			restart_write() ;
			trestore += DTr;
		}


		/* deal with failed timestep, though we usually exit upon failure */
		if(failed) {
			restart_init() ;
			failed = 0 ;
			nfailed = nstep ;
			defcon = 0.3 ;
		}
		if(nstep > nfailed + DTr*4.*(1 + 1./defcon)) defcon = 1. ;


	}
//	fprintf(stderr,"ns,ts: %d %d\n",nstep,nstep*N1*N2) ;

	/* do final diagnostics */
	diag(FINAL_OUT) ;

	return(0) ;
}


/*****************************************************************/
/*****************************************************************
  set_arrays():
  ----------

       -- sets to zero all arrays, plus performs pointer trick 
          so that grid arrays can legitimately refer to ghost 
          zone quantities at positions  i = -2, -1, N1, N1+1 and 
          j = -2, -1, N2, N2+1

 *****************************************************************/
void set_arrays()
{
#ifndef MPI_USED
	int i,j,k ;
#endif

#if 0
	p     =  (double (*) [N2+4][NPR])(&(  p[2][2][0])) ;
	dq    =  (double (*) [N2+4][NPR])(&( dq[2][2][0])) ;
	F1    =  (double (*) [N2+4][NPR])(&( F1[2][2][0])) ;
	F2    =  (double (*) [N2+4][NPR])(&( F2[2][2][0])) ;
	ph    =  (double (*) [N2+4][NPR])(&( ph[2][2][0])) ;

	pflag =  (int    (*) [N2+4])(     &( pflag[2][2] )) ;
#endif
/*
	p      =  new double[N3+2][N1+4][N2+4][NPR] ;
	dq     =  new double[N3+2][N1+4][N2+4][NPR] ;
	F1     =  new double[N3+2][N1+4][N2+4][NPR] ;
	F2     =  new double[N3+2][N1+4][N2+4][NPR] ;
	F3     =  new double[N3+2][N1+4][N2+4][NPR] ;
	ph     =  new double[N3+2][N1+4][N2+4][NPR] ;

	pflag  =  new char[N3+2][N1+4][N2+4] ;

	CoolVal = new double[N3+2][N1+4][N2+4][NCOOLVAL] ;

	A       = new double[N3+2][N1+4][N2+4] ;
*/

	_p      =  (double *)calloc(SN3*SN1*SN2*NPR,sizeof(double));
	_dq     =  (double *)calloc(SN3*SN1*SN2*NPR,sizeof(double));
	_F1     =  (double *)calloc(SN3*SN1*SN2*NPR,sizeof(double));
	_F2     =  (double *)calloc(SN3*SN1*SN2*NPR,sizeof(double));
	_F3     =  (double *)calloc(SN3*SN1*SN2*NPR,sizeof(double));
	_ph     =  (double *)calloc(SN3*SN1*SN2*NPR,sizeof(double));
	_pflag  =  (char *)calloc(SN3*SN1*SN2,sizeof(char));
	_A      =  (double *)calloc(SN3*SN1*SN2,sizeof(double));
	_A2      =  (double *)calloc(SN3*SN1*SN2,sizeof(double));
	_A3      =  (double *)calloc(SN3*SN1*SN2,sizeof(double));

	if(!_p || !_dq || !_F1 || !_F2 || !_F3 || !_ph || !_pflag || !_A || !_A2 || !_A3)
	{
		printf("Invalid memory alloc!\n");
		exit(-1);
	}

#if( COOL_EOS )
	_CoolVal=  (double *)calloc(SN3*SN1*SN2*NCOOLVAL,sizeof(double));
	if(!_CoolVal)
	{
		printf("Invalid memory alloc!\n");
		exit(-1);
	}
#endif

	/* everything must be initialized to zero */
/*
//	ZSLOOP(-2,N1+1,-2,N2+1) {
	ZSLOOP(-1,1,-2,2,-2,2) {
		PLOOP {
			p[h][i][j][k]   = 0. ;
			ph[h][i][j][k]  = 0. ;
			dq[h][i][j][k]  = 0. ;
			F1[h][i][j][k]  = 0. ;
			F2[h][i][j][k]  = 0. ;
			F3[h][i][j][k]  = 0. ;
		}
		pflag[h][i][j] = 0 ;
	}
*/

#ifndef MPI_USED
	k = 0;
	IMAGELOOP {  
	    failimage[0][k] = failimage[1][k] = failimage[2][k] = failimage[3][k] = failimage[4][k++] = 0 ; 
	}
#endif

#if 0
	/* grid functions */
	conn = (double (*) [N2+4][NDIM][NDIM][NDIM])
		(& ( conn[2][2][0][0][0])) ;
	gcon = (double (*) [N2+4][NPG][NDIM][NDIM])
		(& ( gcon[2][2][0][0][0])) ;
	gcov = (double (*) [N2+4][NPG][NDIM][NDIM])
		(& ( gcov[2][2][0][0][0])) ;
	gdet = (double (*) [N2+4][NPG])
		(& ( gdet[2][2][0])) ;
#endif
/*
	conn   =  new double[N3+2][N1+4][N2+4][NDIM][NDIM][NDIM] ;
	gcon   =  new double[N3+2][N1+4][N2+4][NPG][NDIM][NDIM] ;
	gcov   =  new double[N3+2][N1+4][N2+4][NPG][NDIM][NDIM] ;
	gdet   =  new double[N3+2][N1+4][N2+4][NPG] ;
*/
	_conn  =  (double *)calloc(SN1*SN2*NDIM*NDIM*NDIM,sizeof(double));
	_gcon  =  (double *)calloc(SN1*SN2*NPG*NDIM*NDIM,sizeof(double));
	_gcov  =  (double *)calloc(SN1*SN2*NPG*NDIM*NDIM,sizeof(double));
	_gdet  =  (double *)calloc(SN1*SN2*NPG,sizeof(double));
	_coords=  (double *)calloc(SN1*SN2*NDIM,sizeof(double));
#if GHOSTZONESIGNCHANGEINMETRIC
	_sf    =  (double *)calloc(SN1*SN2*NPG*NPR,sizeof(double));
#endif

	if(!_conn || !_gcon || !_gcov || !_gdet)
	{
		printf("Invalid memory alloc 2!\n");
		exit(-1);
	}

}


/*****************************************************************/
/*****************************************************************
  set_grid():
  ----------

       -- calculates all grid functions that remain constant 
          over time, such as the metric (gcov), inverse metric 
          (gcon), connection coefficients (conn), and sqrt of 
          the metric's determinant (gdet).

 *****************************************************************/
void set_grid()
{
  double (*   conn)[SN2][NDIM][NDIM][NDIM] = (double (*) [SN2][NDIM][NDIM][NDIM])(_conn);
  double (*   gcon)[SN2][NPG][NDIM][NDIM] = (double (*) [SN2][NPG][NDIM][NDIM])(_gcon);
  double (*   gcov)[SN2][NPG][NDIM][NDIM] = (double (*) [SN2][NPG][NDIM][NDIM])(_gcov);
  double (*   gdet)[SN2][NPG] = (double (*) [SN2][NPG])(_gdet);
  double (*   coords)[SN2][NDIM] = (double (*) [SN2][NDIM])(_coords);

	int h,i,j,k;
	double X[NDIM], V[NDIM] ;
	struct of_geom geom ;

	/* set up boundaries, steps in coordinate grid */
	set_points() ;
	dV = dx[RR]*dx[TH]*dx[PH] ;

	DLOOPA X[j] = 0. ;

//	ZSLOOP(-2,N1+1,-2,N2+1) {
	h = (N3>1 ? 1 : 0);

	for(i=0;i<SN1;++i)for(j=0;j<SN2;++j)
	{
//		printf("h=%3d i=%3d j=%3d\n",h,i,j);  
		/* zone-centered */
		coord(h,i,j,CENT,X) ;
		gcov_func(X, gcov[i][j][CENT]) ;
		gdet[i][j][CENT] = gdet_func(gcov[i][j][CENT]) ;
		gcon_func(gcov[i][j][CENT], gcon[i][j][CENT]) ;

		get_geometry(h,i,j,CENT,&geom) ;
		conn_func(X, &geom, conn[i][j]) ;

		bl_coord_vec(X, V);
		for(k=0;k<NDIM;k++) coords[i][j][k] = V[k];

		/* corner-centered */
		coord(h,i,j,CORN,X) ;
		gcov_func(X,gcov[i][j][CORN]) ;
		gdet[i][j][CORN] = gdet_func(gcov[i][j][CORN]) ;

		gcon_func(gcov[i][j][CORN],gcon[i][j][CORN]) ;

		/* r-face-centered */
		coord(h,i,j,FACE1,X) ;
		gcov_func(X,gcov[i][j][FACE1]) ;
		gdet[i][j][FACE1] = gdet_func(gcov[i][j][FACE1]) ;
		gcon_func(gcov[i][j][FACE1],gcon[i][j][FACE1]) ;

		/* theta-face-centered */
		coord(h,i,j,FACE2,X) ;
		gcov_func(X,gcov[i][j][FACE2]) ;
		gdet[i][j][FACE2] = gdet_func(gcov[i][j][FACE2]) ;
		gcon_func(gcov[i][j][FACE2],gcon[i][j][FACE2]) ;

		/* phi-face-centered */
		coord(h,i,j,FACE3,X) ;
		gcov_func(X,gcov[i][j][FACE3]) ;
		gdet[i][j][FACE3] = gdet_func(gcov[i][j][FACE3]) ;
		gcon_func(gcov[i][j][FACE3],gcon[i][j][FACE3]) ;

#if (FLIPGDETAXIS == 1)
    	{
			double 	r,th,ph;
            int ll;

            for(ll=0;ll<NPG;ll++){
				coord(h,i,j,ll,X); 
				bl_coord(X,&r,&th,&ph) ;
				
                if(th < -SINGSMALL || th > M_PI + SINGSMALL)
					gdet[i][j][ll] *= -1.;
			}
		}
#endif
#if GHOSTZONESIGNCHANGEINMETRIC
    	{
		    double (*   sf)[SN2][NPG][NPR] = (double (*) [SN2][NPG][NPR])(_sf);
			double 	r,th,ph;
            int     k,ll;

            for(ll=0;ll<NPG;ll++){
                PLOOP sf[i][j][ll][k] = 1.;

			    coord(h,i,j,ll,X); 
                bl_coord(X,&r,&th,&ph) ;
			    if(th < -SINGSMALL || th > M_PI + SINGSMALL) {
					sf[i][j][ll][U2] = -1.;
					sf[i][j][ll][U3] = -1.;
					sf[i][j][ll][B2] = -1.;
					sf[i][j][ll][B3] = -1.;
                }
         	}
		}
#endif
	}

	/* done! */

#if( COOL_EOS )
	SetupCache_CoolEOS();
#endif

#ifdef MPI_USED
	set_MPI_Slices();
#endif

}

