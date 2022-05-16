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

#include "decs.h"

#ifdef MPI_USED
#include "mpi.h"
#endif

#if( COOL_EOS )
#include "cooleos.h"
#endif

extern double tlog ;

/* all diagnostics subroutine */

void diag(int call_code)
{
	double (*   p)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_p);
    double (*   gdet)[SN2][NPG] = (double (*) [SN2][NPG])(_gdet);
#if( COOL_EOS )
    double (*   CoolVal)[SN1][SN2][NCOOLVAL] = (double (*) [SN1][SN2][NCOOLVAL])(_CoolVal);
#endif

	char dfnam[100];
	int h,i,j ;
	FILE *dump_file;
	double U[NPR],pp,e,rmed,divb,divbmax,e_fin,m_fin, divbmaxglob;
	struct of_geom geom ;
	struct of_state q ;
	int imax[NDIM];


	double lnu = 0.0;
	double lnu2 = 0.0;
#if( COOL_EOS )
	double lnu_rho_1=0.0;
#endif
	double lnu_rho_2=0.0;
	double lnu_rho=0.0;
    double mdot_=mdot;
    double edot_=edot;
    double ldot_=ldot;

	static double e_init,m_init ;
	static FILE *ener_file ;

	if(call_code==INIT_OUT) {
		/* set things up */
		ener_file = fopen("ener.out","a") ;
		if(ener_file==NULL) {
			fprintf(stderr,"error opening energy output file\n") ;
			exit(1) ;
		}
	}

#ifdef MPI_USED
		MPI_Allreduce(MPI_IN_PLACE,&mdot_,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE,&edot_,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE,&ldot_,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif

	/* calculate conserved quantities */
	if((call_code==INIT_OUT || 
	   call_code==LOG_OUT ||
	   call_code==FINAL_OUT) &&
	   !failed) {
		pp = 0. ;
		e = 0. ;
		rmed = 0. ;
		divbmax = 0. ;


#ifdef MPI_USED
	MPI_Sendrecv(_p, 1, slice_iNPR1_stop, rank_nexti, 124,
                _p, 1, slice_iNPR1_startP, rank_previ, 124,
                MPI_COMM_WORLD, &MPIStatus);
//	MPI_Sendrecv(_p, 1, slice_iNPR1_start, rank_previ, 125,
//                _p, 1, slice_iNPR1_stopN, rank_nexti, 125,
//                MPI_COMM_WORLD, &MPIStatus);
	if(N3>2)
	{
		MPI_Sendrecv(_p, 1, slice_hNPR1_stop, rank_nexth, 125,
                _p, 1, slice_hNPR1_startP, rank_prevh, 125,
                MPI_COMM_WORLD, &MPIStatus);
//		MPI_Sendrecv(_p, 1, slice_hNPR1_start, rank_prevh, 125,
//                _p, 1, slice_hNPR1_stopN, rank_nexth, 125,
//                MPI_COMM_WORLD, &MPIStatus);
	}

#endif

		
//		ZSLOOP(0,N1-1,0,N2-1) {
		ZLOOP {
			get_geometry(h,i,j,CENT,&geom) ;

#if( COOL_EOS )
			InitCache_CoolEOS(h, i, j);
#endif

//			fprintf(stdout, "rank=%d, i=%d, j=%d\n", rank, i, j);
			get_state(p[h][i][j],&geom,&q) ;
			primtoU(p[h][i][j],&q,&geom,U) ;

#if( COOL_EOS )
			CoolVal_rho0_u_CoolEOS(p[h][i][j][RHO], p[h][i][j][UU], h, i, j, CoolVal[h][i][j]);
#endif

			rmed += U[RHO]*dV ;
			pp += U[U3]*dV ;
			e += U[UU]*dV ;

#if N3>2
			/* flux-ct defn */
			divb = fabs( 0.25*(
			 	  p[h][i][j][B1]*gdet[i][j][CENT] 
			 	+ p[h-1][i][j][B1]*gdet[i][j][CENT] 
				+ p[h][i][j-1][B1]*gdet[i][j-1][CENT]
				+ p[h-1][i][j-1][B1]*gdet[i][j-1][CENT]
				- p[h][i-1][j][B1]*gdet[i-1][j][CENT] 
				- p[h-1][i-1][j][B1]*gdet[i-1][j][CENT]
				- p[h][i-1][j-1][B1]*gdet[i-1][j-1][CENT]
				- p[h-1][i-1][j-1][B1]*gdet[i-1][j-1][CENT]
				)/dx[RR] +
				0.25*(
				  p[h][i][j][B2]*gdet[i][j][CENT] 
				+ p[h-1][i][j][B2]*gdet[i][j][CENT] 
				+ p[h][i-1][j][B2]*gdet[i-1][j][CENT]
				+ p[h-1][i-1][j][B2]*gdet[i-1][j][CENT]
				- p[h][i][j-1][B2]*gdet[i][j-1][CENT] 
				- p[h-1][i][j-1][B2]*gdet[i][j-1][CENT] 
				- p[h][i-1][j-1][B2]*gdet[i-1][j-1][CENT]
				- p[h-1][i-1][j-1][B2]*gdet[i-1][j-1][CENT]
				)/dx[TH] +
				0.25*(
				  p[h][i][j][B3]*gdet[i][j][CENT] 
				+ p[h][i-1][j][B3]*gdet[i-1][j][CENT] 
				+ p[h][i][j-1][B3]*gdet[i][j-1][CENT]
				+ p[h][i-1][j-1][B3]*gdet[i-1][j-1][CENT]
				- p[h-1][i][j][B3]*gdet[i][j][CENT] 
				- p[h-1][i-1][j][B3]*gdet[i-1][j][CENT] 
				- p[h-1][i][j-1][B3]*gdet[i][j-1][CENT]
				- p[h-1][i-1][j-1][B3]*gdet[i-1][j-1][CENT]
				)/dx[PH]) ;
#else
			/* flux-ct defn */
			divb = fabs( 0.5*(
			 	  p[h][i][j][B1]*gdet[i][j][CENT] 
				+ p[h][i][j-1][B1]*gdet[i][j-1][CENT]
				- p[h][i-1][j][B1]*gdet[i-1][j][CENT] 
				- p[h][i-1][j-1][B1]*gdet[i-1][j-1][CENT]
				)/dx[RR] +
				0.5*(
				  p[h][i][j][B2]*gdet[i][j][CENT] 
				+ p[h][i-1][j][B2]*gdet[i-1][j][CENT]
				- p[h][i][j-1][B2]*gdet[i][j-1][CENT] 
				- p[h][i-1][j-1][B2]*gdet[i-1][j-1][CENT]
				)/dx[TH]) ;
#endif

			if(divb > divbmax) {
				imax[RR] = i - (irank==0 ? 2 : 1) + di_rank_start ;
				imax[TH] = j - (jrank==0 ? 2 : 1) + dj_rank_start ;
                imax[PH] = N3 == 1 ? 0 : h - 1 + dh_rank_start ;
				divbmax = divb ;
			}
			
#if( COOL_EOS )
			lnu += CoolVal[h][i][j][COOL_LAMBDA_SIM]*gdet[i][j][CENT];

			lnu_rho_1 += CoolVal[h][i][j][COOL_LAMBDA_SIM]*gdet[i][j][CENT]*p[h][i][j][RHO];
			lnu_rho_2 += gdet[i][j][CENT]*p[h][i][j][RHO];
#endif			
			
		}

	}


#if( COOL_EOS )
	ZSLOOPH(0,0) ZSLOOPI(0,0) {
		double Hd_avg = 0.0;
//        double X[NDIM];
        double r, th, ph;
		double hmax = 0.0;

		ZSLOOPJ(N2/2,0) {
			Hd_avg += CoolVal[h][i][j][COOL_Hd];
		}
		Hd_avg /= N2/2;
		
		//Find j equal to Hd_avg;
		ZSLOOPJ(N2/2,0) {
//	        coord (h, i, j, CENT, X);
//    	    bl_coord (X, &r, &th, &ph);
    	    get_coord (h, i, j, &r, &th, &ph);
    	    
			hmax = fabs(r*tan(th));
	    	if(Hd_avg <= fabs(r*tan(th))) {
				break;
			}
		}
		if(j > jstop) j = jstop;

		lnu2 += CoolVal[h][i][j][COOL_Qnu]*gdet[i][j][CENT];

//		printf("lnu ir=%d, j=%d, Hd=%E, hmax=%E\n", i, j, Hd_avg, hmax);
		
		
	}
#endif


//	fprintf(stdout, "rank=%d finished\n", rank);
	
#ifdef MPI_USED
	MPI_Allreduce(MPI_IN_PLACE,&pp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE,&e,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE,&rmed,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&divbmax,&divbmaxglob,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#else
    divbmaxglob = divbmax;
#endif


	if(call_code == INIT_OUT) {
		e_init = e ;
		m_init = rmed ;
	}

	if(call_code == FINAL_OUT) {
		e_fin = e ;
		m_fin = rmed ;
		if(rank == 0)
		{
			fprintf(stderr,"\n\nEnergy: ini,fin,del: %g %g %g\n",
				e_init,e_fin,(e_fin-e_init)/e_init) ;
			fprintf(stderr,"mass: ini,fin,del: %g %g %g\n",
				m_init,m_fin,(m_fin-m_init)/m_init) ;
		}
	}

	if((call_code == INIT_OUT && tlog == 0.)|| 
	   call_code == LOG_OUT ||
	   call_code == FINAL_OUT) {

#ifdef MPI_USED
		MPI_Allreduce(MPI_IN_PLACE,&lnu,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE,&lnu2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#if( COOL_EOS )
		MPI_Allreduce(MPI_IN_PLACE,&lnu_rho_1,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE,&lnu_rho_2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
#endif

#if( COOL_EOS )
		lnu_rho=lnu_rho_1/lnu_rho_2;
#endif			
        if(rank != 0 && divbmax == divbmaxglob)
        {
	        MPI_Send(imax, NDIM, MPI_INT, 0, 126, MPI_COMM_WORLD);
        }
        else if(rank == 0 && divbmax != divbmaxglob)
        {
            divbmax = divbmaxglob;
	        MPI_Recv(imax, NDIM, MPI_INT, MPI_ANY_SOURCE, 126, MPI_COMM_WORLD, &MPIStatus);
        }

		if(rank == 0)
		{
			fprintf(stderr,"LOG      t=%g \t divbmax: %d %d %d %g\n",
				t,imax[RR],imax[TH],imax[PH],divbmax) ;

			fprintf(ener_file,"%10.5g %10.5g %10.5g %10.5g %15.8g %15.8g ",
				t,rmed,pp,e,0.0, 0.0 /*p[N1/2][N2/2][UU]*pow(p[N1/2][N2/2][RHO],-gam),
				p[N1/2][N2/2][UU]*/) ;
			fprintf(ener_file,"%15.8g %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g",mdot,edot,ldot,lnu,lnu_rho,lnu_rho_2, lnu2) ;

			fprintf(ener_file," -0 %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g",mdot_in, mdot_out, mass_start, mass_curr, mass_in, 
								mass_out, (mass_curr-mass_in-mass_out)/mass_start*100.0) ;

			fprintf(ener_file,"\n") ;
			fflush(ener_file) ;
		}
	}


	/* dump at regular intervals */
	if(call_code == INIT_OUT || 
	   call_code == DUMP_OUT ||
	   call_code == FINAL_OUT) {
#if HDFDUMP	
		/* make regular dump file */
		dumph5(dump_cnt) ;
#else

		/* make regular dump file */
		sprintf(dfnam,"dumps/dump%03d",dump_cnt) ;
		if (rank == 0) {
			fprintf(stderr,"DUMP     file=%s\n",dfnam) ;
        }

		remove(dfnam);

#ifdef MPI_USED
		MPI_Barrier(MPI_COMM_WORLD);
#endif

		dump_file = fopen(dfnam,"a") ;

		if(dump_file==NULL) {
			fprintf(stderr,"error opening dump file\n") ;
			exit(2) ;
		}

		dump(dump_file) ;
		fclose(dump_file) ;

#endif
		dump_cnt++ ;
	}

	/* image dump at regular intervals */
	if(call_code == IMAGE_OUT || 
	   call_code == INIT_OUT ||
	   call_code == FINAL_OUT) {

		image_all( image_cnt );

		image_cnt++ ;
	}
}


/** some diagnostic routines **/


void fail(int fail_type)
{
	double (*   p)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_p);

	failed = 1 ;

	fprintf(stderr,"\n\nfail: %d %d %d %d\n",hcurr, icurr,jcurr,fail_type) ;

	area_map(hcurr,icurr,jcurr,(double*)p) ;
	
	fprintf(stderr,"fail\n") ;

	diag(FINAL_OUT) ;

	/* for diagnostic purposes */
	exit(0) ;
}



/* map out region around failure point */
void area_map(int h, int i, int j, double *_prim)
{
	double (*   prim)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_prim);
	int k ;

	fprintf(stderr,"area map\n") ;

	PLOOP {
		fprintf(stderr,"variable %d \n",k) ;
		fprintf(stderr,"h = \t %12d %12d %12d\n",h-1,h,h+1) ;
		fprintf(stderr,"i = \t %12d %12d %12d\n",i-1,i,i+1) ;
		fprintf(stderr,"j = %d \t %12.5g %12.5g %12.5g\n",
				j+1,
				prim[h][i-1][j+1][k],
				prim[h][i][j+1][k],
				prim[h][i+1][j+1][k]) ;
		fprintf(stderr,"j = %d \t %12.5g %12.5g %12.5g\n",
				j,
				prim[h][i-1][j][k],
				prim[h][i][j][k],
				prim[h][i+1][j][k]) ;
		fprintf(stderr,"j = %d \t %12.5g %12.5g %12.5g\n",
				j-1,
				prim[h][i-1][j-1][k],
				prim[h][i][j-1][k],
				prim[h][i+1][j-1][k]) ;
	}

	/* print out other diagnostics here */

}

/* evaluate fluxed based diagnostics; put results in
 * global variables */

void diag_flux(double *_F1, double *_F2, double *_F3)
{
	double (*   F1)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_F1);
//	double (*   F2)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_F2);
//	double (*   F3)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_F3);
	struct of_geom geom ;
	int h, j ;

    mdot = edot = ldot = 0. ;
    mdot_in = mdot_out = 0.;

//        for(j=0;j<N2;j++) {
	if(irank==0) ZSLOOPH(0,0) ZSLOOPJ(0,0) {
//		    get_geometry(h,istart,j,CENT,&geom) ;
//			double alpha = 1./sqrt(-geom.gcon[0][0]) ;

            mdot += F1[h][istart][j][RHO]*dx[TH]*dx[PH] ;
//            mdot_in += F1[h][istart][j][RHO]*dx[TH]*dx[PH]*alpha ;
            edot -= (F1[h][istart][j][UU] - F1[h][istart][j][RHO])*dx[TH]*dx[PH] ;
            ldot += F1[h][istart][j][U3]*dx[TH]*dx[PH] ;
    }
/*
	if(irank==isize-1) ZSLOOPH(0,0) ZSLOOPJ(0,0) {
//		    get_geometry(h,istop,j,CENT,&geom) ;
//			double alpha = 1./sqrt(-geom.gcon[0][0]) ;
            mdot_out -= F1[h][istop][j][RHO]*dx[TH]*dx[PH] ;
    }
*/
}

#define NEU_AGA           1
#define NEU_AGA_GHOST1    2
#define NEU_AGA_AVG010    3
#define NEU_AGA_AVG1020   4

//#define MDOTTYPE   NEU_AGA
#define MDOTTYPE   NEU_AGA_GHOST1
//#define MDOTTYPE   NEU_AGA_AVG010
//#define MDOTTYPE   NEU_AGA_AVG1020

void diag_mass(double Dt)
{
	double (*   p)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_p);
    double (*   gdet)[SN2][NPG] = (double (*) [SN2][NPG])(_gdet);
	struct of_geom geom ;
    double ucon[NDIM] ;
	int h, i, j ;

    mass_curr = 0.;
    mdot_in = mdot_out = 0.;

	ZLOOP {
    	geom.g = gdet[i][j][CENT] ;
		mass_curr += p[h][i][j][RHO]*geom.g*dx[RR]*dx[TH]*dx[PH];
   	}

#if MDOTTYPE == NEU_AGA
	i = istart;
	if(irank==0) ZSLOOPH(0,0) ZSLOOPJ(0,0) {
        get_geometry(h,i,j,CENT,&geom) ;
        ucon_calc(p[h][i][j], &geom, ucon) ;

        mdot_in += p[h][i][j][RHO]*ucon[1]*geom.g*dx[TH]*dx[PH] ;
    }

	i = istop;
	if(irank==isize-1) ZSLOOPH(0,0) ZSLOOPJ(0,0) {
        get_geometry(h,i,j,CENT,&geom) ;
        ucon_calc(p[h][i][j], &geom, ucon) ;

        mdot_out -= p[h][i][j][RHO]*ucon[1]*geom.g*dx[TH]*dx[PH] ;
    }
#elif MDOTTYPE == NEU_AGA_GHOST1
	i = istart-1;
	if(irank==0) ZSLOOPH(0,0) ZSLOOPJ(0,0) {
        get_geometry(h,i,j,CENT,&geom) ;
        ucon_calc(p[h][i][j], &geom, ucon) ;

        mdot_in += p[h][i][j][RHO]*ucon[1]*geom.g*dx[TH]*dx[PH] ;
    }

	i = istop+1;
	if(irank==isize-1) ZSLOOPH(0,0) ZSLOOPJ(0,0) {
        get_geometry(h,i,j,CENT,&geom) ;
        ucon_calc(p[h][i][j], &geom, ucon) ;

        mdot_out -= p[h][i][j][RHO]*ucon[1]*geom.g*dx[TH]*dx[PH] ;
    }
#elif MDOTTYPE == NEU_AGA_AVG010
	if(irank==0) for(i=istart;i<=istart+10;i++) ZSLOOPH(0,0) ZSLOOPJ(0,0) {
        get_geometry(h,i,j,CENT,&geom) ;
        ucon_calc(p[h][i][j], &geom, ucon) ;

        mdot_in += p[h][i][j][RHO]*ucon[1]*geom.g*dx[TH]*dx[PH]/11.0 ;
    }

	if(irank==isize-1) for(i=istop-10;i<=istop;i++) ZSLOOPH(0,0) ZSLOOPJ(0,0) {
        get_geometry(h,i,j,CENT,&geom) ;
        ucon_calc(p[h][i][j], &geom, ucon) ;

        mdot_out -= p[h][i][j][RHO]*ucon[1]*geom.g*dx[TH]*dx[PH]/11.0 ;
    }
#elif MDOTTYPE == NEU_AGA_AVG1020
	if(irank==0) for(i=istart+10;i<=istart+20;i++) ZSLOOPH(0,0) ZSLOOPJ(0,0) {
        get_geometry(h,i,j,CENT,&geom) ;
        ucon_calc(p[h][i][j], &geom, ucon) ;

        mdot_in += p[h][i][j][RHO]*ucon[1]*geom.g*dx[TH]*dx[PH]/11.0 ;
    }

	if(irank==isize-1) for(i=istop-20;i<=istop-10;i++) ZSLOOPH(0,0) ZSLOOPJ(0,0) {
        get_geometry(h,i,j,CENT,&geom) ;
        ucon_calc(p[h][i][j], &geom, ucon) ;

        mdot_out -= p[h][i][j][RHO]*ucon[1]*geom.g*dx[TH]*dx[PH]/11.0 ;
    }
#else
Ququ
#endif


#ifdef MPI_USED
/*
		MPI_Allreduce(MPI_IN_PLACE,&mdot_in,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE,&mdot_out,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE,&mass_curr,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
*/
	{	
		double reduceAll[3];
		reduceAll[0] = mdot_in; 
		reduceAll[1] = mdot_out; 
		reduceAll[2] = mass_curr; 

		MPI_Allreduce(MPI_IN_PLACE,reduceAll,sizeof(reduceAll)/sizeof(reduceAll[0]),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

		mdot_in = reduceAll[0]; 
		mdot_out = reduceAll[1]; 
		mass_curr = reduceAll[2]; 
	}
#endif

   	mass_in += mdot_in*Dt;   
   	mass_out += mdot_out*Dt;   
}
