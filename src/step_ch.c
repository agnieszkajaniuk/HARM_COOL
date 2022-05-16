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

/**
 *
 * this contains the generic piece of code for advancing
 * the primitive variables 
 *
**/


#include "decs.h"

#if( COOL_EOS )
#include "cooleos.h"
#endif

#if TRACERSDUMP
#include "tracers.h"
#endif

/** algorithmic choices **/

/* use local lax-friedrichs or HLL flux:  these are relative weights on each numerical flux */
#define HLLF  (0.0)
#define LAXF  (1.0)

/** end algorithmic choices **/


double advance( double *_pi, double *_pb, 
	double Dt, double *_pf, int fullstep) ;
double fluxcalc( double *_pr, double *_F, 
	const int dir, const int fullstep ) ;
void   flux_cd(double *_F1, double *_F2) ;

void debug2(double *pr, const char *where);

/***********************************************************************************************/
/***********************************************************************************************
  step_ch():
  ---------
     -- handles the sequence of making the time step, the fixup of unphysical values, 
        and the setting of boundary conditions;

     -- also sets the dynamically changing time step size;

***********************************************************************************************/
void step_ch()
{
//	double (*   p)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_p);
//	double (*  ph)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_ph);
	double ndt ;
//	int h,i,j,k ;

//	debug2(_p, "step_ch 1");

	if(rank == 0)	
		fprintf(stdout,"h") ;
	ndt = advance(_p, _p, 0.5*dt, _ph, 0) ;   /* time step primitive variables to the half step */

	fixup(_ph) ;         /* Set updated densities to floor, set limit for gamma */
	bound_prim(_ph) ;    /* Set boundary conditions for primitive variables, flag bad ghost zones */
	fixup_utoprim(_ph);  /* Fix the failure points using interpolation and updated ghost zone values */
	bound_prim(_ph) ;    /* Reset boundary conditions with fixed up points */

	/* Repeat and rinse for the full time (aka corrector) step:  */
	if(rank == 0)	
		fprintf(stdout,"f") ;
//	ZLOOP PLOOP psave[h][i][j][k] = p[h][i][j][k] ;
	ndt = advance(_p, _ph, dt,    _p, 1) ;

	fixup(_p) ;
    bound_prim(_p) ;
	fixup_utoprim(_p);
	bound_prim(_p) ;

	diag_mass(dt);

#if TRACERSDUMP
	update_tracers(dt);
#endif

	/* Determine next time increment based on current characteristic speeds: */
        if(dt < 1.e-9) {
                fprintf(stderr,"^timestep too small\n") ;
                exit(11) ;
        }

        /* increment time */
        t += dt ;

        /* set next timestep */
        if(ndt > SAFE*dt) ndt = SAFE*dt ;
        dt = ndt ;
        if(t + dt > tf) dt = tf - t ;  /* but don't step beyond end of run */

        /* done! */
}

/***********************************************************************************************/
/***********************************************************************************************
  advance():
  ---------
     -- responsible for what happens during a time step update, including the flux calculation, 
         the constrained transport calculation (aka flux_ct()), the finite difference 
         form of the time integral, and the calculation of the primitive variables from the 
         update conserved variables;
     -- also handles the "fix_flux()" call that sets the boundary condition on the fluxes;

***********************************************************************************************/
double advance(
	double *_pi, 
	double *_pb, 
	double Dt,
	double *_pf,
	int fullstep
	)
{
  double (*   F1)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_F1);
  double (*   F2)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_F2);
  double (*   F3)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_F3);
  double (*   pi)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_pi);
  double (*   pb)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_pb);
  double (*   pf)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_pf);
  char   (*pflag)[SN1][SN2] = (char (*) [SN1][SN2])(_pflag);

  int h,i,j,k;
  double ndt,ndt1,ndt2,ndt3,U[NPR],dU[NPR] ;
  struct of_geom geom ;
  struct of_state q ;
  
//  double Uo[NPR],po[NPR] ;
  
  if(pf != pi)	
	  ZLOOP PLOOP pf[h][i][j][k] = pi[h][i][j][k] ;        /* needed for Utoprim */
  
  if(rank == 0)	
	  fprintf(stdout,"0") ;

#ifdef MPI_USED
	MPI_Sendrecv(pb, 1, slice_iNPR1_stop, rank_nexti, 122,
                pb, 1, slice_iNPR1_startP, rank_previ, 122,
                MPI_COMM_WORLD, &MPIStatus);
	MPI_Sendrecv(pb, 1, slice_iNPR1_start, rank_previ, 123,
                pb, 1, slice_iNPR1_stopN, rank_nexti, 123,
                MPI_COMM_WORLD, &MPIStatus);
#if	N3>2
	{
		MPI_Sendrecv(pb, 1, slice_hNPR1_stop, rank_nexth, 124,
                pb, 1, slice_hNPR1_startP, rank_prevh, 124,
                MPI_COMM_WORLD, &MPIStatus);
		MPI_Sendrecv(pb, 1, slice_hNPR1_start, rank_prevh, 125,
                pb, 1, slice_hNPR1_stopN, rank_nexth, 125,
                MPI_COMM_WORLD, &MPIStatus);
	}
#endif

#endif

  ndt1 = fluxcalc(_pb, _F1, 1, fullstep) ;
  ndt2 = fluxcalc(_pb, _F2, 2, fullstep) ;

  if(N3>2) 
      ndt3 = fluxcalc(_pb, _F3, 3, fullstep) ;

  fix_flux(_F1,_F2,_F3) ;

  flux_ct(_F1,_F2,_F3) ;

  /* evaluate diagnostics based on fluxes */
  if(fullstep)
  {
       diag_flux(_F1,_F2,_F3) ;
  }

#ifdef MPI_USED
//	MPI_Sendrecv(F1, 1, slice_iNPR1_stop, rank_nexti, 126,
//                F1, 1, slice_iNPR1_startP, rank_previ, 126,
//                MPI_COMM_WORLD, &MPIStatus);
	MPI_Sendrecv(F1, 1, slice_iNPR1_start, rank_previ, 127,
                F1, 1, slice_iNPR1_stopN, rank_nexti, 127,
                MPI_COMM_WORLD, &MPIStatus);
	if(N3>1)
	{
		MPI_Sendrecv(F3, 1, slice_hNPR1_start, rank_prevh, 128,
                F3, 1, slice_hNPR1_stopN, rank_nexth, 128,
                MPI_COMM_WORLD, &MPIStatus);
	}
#endif


  if(rank == 0)	
  	fprintf(stdout,"1") ;
  /** now update pi to pf **/
  ZLOOP {

    get_geometry(h,i,j,CENT,&geom) ;

#if( COOL_EOS )
	InitCache_CoolEOS(h, i, j);
#endif

    source(pb[h][i][j],&geom,h,i,j,dU,Dt) ;

    get_state(pi[h][i][j],&geom,&q) ;
    primtoU(pi[h][i][j],&q,&geom,U) ;

   	PLOOP {
   	U[k] += Dt*(
#if( N1 != 1 )
		  - (F1[h][i+1][j][k] - F1[h][i][j][k])/dx[RR]
#endif
#if( N2 != 1 )
		  - (F2[h][i][j+1][k] - F2[h][i][j][k])/dx[TH]
#endif
#if N3>2
		  - (F3[h+1][i][j][k] - F3[h][i][j][k])/dx[PH]
#endif
		  + dU[k]
		  ) ;
    }

#if( COOL_EOS )
    pflag[h][i][j] = Utoprim_5d(U, geom.gcov, geom.gcon, geom.g, pf[h][i][j]);
#else
    pflag[h][i][j] = Utoprim_2d(U, geom.gcov, geom.gcon, geom.g, pf[h][i][j]);
#endif
//	if(pflag[h][i][j])
//		fprintf(stdout, "pflag[%d][%d][%d] = %d\n", h, i, j, pflag[h][i][j]);

#ifndef MPI_USED
    if( pflag[h][i][j] ) failimage[0][i-istart+(j-jstart)*N1]++ ;
#endif

#if( !COOL_EOS )
#if( DO_FONT_FIX  || DOKTOT) 
    if( pflag[h][i][j] ) { 
      double kappa;

      if (DOKTOT) {
        kappa = pf[h][i][j][KTOT];
      } else {
        kappa = Katm[i];
      }

      pflag[h][i][j] = Utoprim_1dvsq2fix1(U, geom.gcov, geom.gcon, geom.g, pf[h][i][j], kappa );
      if( pflag[h][i][j] ) { 
#ifndef MPI_USED
	failimage[1][i-istart+(j-jstart)*N1]++ ;
#endif
	pflag[h][i][j] = Utoprim_1dfix1(U, geom.gcov, geom.gcon, geom.g, pf[h][i][j], kappa );
#ifndef MPI_USED
	if( pflag[h][i][j] ) failimage[2][i-istart+(j-jstart)*N1]++ ;
#endif
      }
    }
#endif
#endif //#if( !COOL_EOS )

#if(DOKTOT)
    compute_ktot(_pi,_pb,_pf,h,i,j,Dt);
#endif
		
  }

#if N3>2
	  ndt = defcon * 1./(1./ndt1 + 1./ndt2 + 1./ndt3) ;
#else
	  ndt = defcon * 1./(1./ndt1 + 1./ndt2) ;
#endif

  if(rank == 0)	
    fprintf(stdout,"2 ") ;

  return(ndt) ;
}


/***********************************************************************************************/
/***********************************************************************************************
  fluxcalc():
  ---------
     -- sets the numerical fluxes, avaluated at the cell boundaries using the slope limiter
        slope_lim();

     -- only has HLL and Lax-Friedrichs  approximate Riemann solvers implemented;
        
***********************************************************************************************/
double fluxcalc(
	double *_pr, 
	double *_F, 
	const int dir,
	const int fullstep
	)
{
    double (*   pr)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_pr);
    double (*   F)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_F);
	double (*   dq)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_dq);
    double (*   sf)[SN2][NPG][NPR] = (double (*) [SN2][NPG][NPR])(_sf);

	int h,i,j,k,hdel,idel,jdel,face ;
	double p_l[NPR],p_r[NPR],F_l[NPR],F_r[NPR],U_l[NPR],U_r[NPR] ;
	double cmax_l,cmax_r,cmin_l,cmin_r,cmax,cmin,ndt,dtij ;
	double ctop ;
	struct of_geom geom ;
	struct of_state state_l,state_r ;
	void rescale(double *pr, int which, int dir, int hh, int ii, int jj, int face, struct of_geom *geom) ;
//	double bsq ;

    if     (dir == 1) {hdel = 0; idel = 1; jdel = 0; face = FACE1;}
	else if(dir == 2) {hdel = 0; idel = 0; jdel = 1; face = FACE2;}
	else if(dir == 3) {hdel = 1; idel = 0; jdel = 0; face = FACE3;}
	else { exit(10); }

#if(RESCALE)
	/** evaluate slopes of primitive variables **/
	/* first rescale */
//	ZSLOOP(-2,N1+1,-2,N2+1) {
	ZSLOOP(0,0,-2,2,-2,2) {
		get_geometry(h,i,j,CENT,&geom) ;
		rescale(pr[h][i][j],FORWARD, dir, h,i,j,CENT,&geom) ;
	}
#endif
	/* then evaluate slopes */
//	ZSLOOP(-1,N1,-1,N2) PLOOP {
	ZSLOOP(0,0,-1,1,-1,1) PLOOP {
   	        //-new get_geometry(h,i,j,CENT,&geom) ;
		//-new bsq = bsq_calc(pr[i][j],&geom) ;
		//-new if(bsq/pr[i][j][RHO] > 10. ||
		//-new    bsq/pr[i][j][UU]  > 1.e3) lim = MINM ;
		//-new else lim = MC ;

//#if GHOSTZONESIGNCHANGEINMETRIC
#if 0
		dq[h][i][j][k] = slope_lim(
				pr[h-hdel][i-idel][j-jdel][k]*sf[i-idel][j-jdel][face][k],
				pr[h][i][j][k]*sf[i][j][face][k],
				pr[h+hdel][i+idel][j+jdel][k]*sf[i+idel][j+jdel][face][k]
				) ;
#else
		dq[h][i][j][k] = slope_lim(
				pr[h-hdel][i-idel][j-jdel][k],
				pr[h][i][j][k],
				pr[h+hdel][i+idel][j+jdel][k]
				) ;
#endif
	}

#ifdef MPI_USED
	if(idel>0)
		MPI_Sendrecv(dq, 1, slice_iNPR1_stop, rank_nexti, 129,
                dq, 1, slice_iNPR1_startP, rank_previ, 129,
                MPI_COMM_WORLD, &MPIStatus);

	if(hdel>0)
		MPI_Sendrecv(dq, 1, slice_hNPR1_stop, rank_nexth, 130,
                dq, 1, slice_hNPR1_startP, rank_prevh, 130,
                MPI_COMM_WORLD, &MPIStatus);
#endif


	ndt = 1.e9 ;
//        ZSLOOP(-jdel,N1,-idel,N2) {
//        ZSLOOP(0,0,-jdel,1,-idel,1) {
        ZSLOOP(0,0,(-hdel-jdel),1,(-hdel-idel),1) {

                /* this avoids problems on the pole */
//                if(dir == 2 && (j == 0 || j == N2)) {

#if !JBOUNDAXISYMMETRIC
                if(dir == 2 && ((jrank == 0 && j == jstart) || (jrank == jsize-1 && j == jstop+1))) {
                        PLOOP F[h][i][j][k] = 0. ;
                }
		else 
#endif
		{

//#if GHOSTZONESIGNCHANGEINMETRIC
#if 0
                PLOOP {
                        p_l[k] = pr[h-hdel][i-idel][j-jdel][k]*sf[i-idel][j-jdel][face][k] 
					+ 0.5*dq[h-hdel][i-idel][j-jdel][k] ;
                        p_r[k] = pr[h][i][j][k]*sf[i][j][face][k]   
					- 0.5*dq[h][i][j][k] ;
                }
#else
                PLOOP {
                        p_l[k] = pr[h-hdel][i-idel][j-jdel][k] 
					+ 0.5*dq[h-hdel][i-idel][j-jdel][k] ;
                        p_r[k] = pr[h][i][j][k]   
					- 0.5*dq[h][i][j][k] ;
                }
#endif

		get_geometry(h,i,j,face,&geom) ;

#if( COOL_EOS )
		InitCache_CoolEOS(h, i, j);
#endif


#if(RESCALE)
		rescale(p_l,REVERSE,dir,h,i,j,face,&geom) ;
		rescale(p_r,REVERSE,dir,h,i,j,face,&geom) ;
#endif
		get_state(p_l,&geom,&state_l) ;
		get_state(p_r,&geom,&state_r) ;
		
		primtoflux(p_l,&state_l,dir,&geom,F_l) ;
		primtoflux(p_l,&state_l,TT, &geom,U_l) ;
		vchar(p_l,&state_l,&geom,dir,&cmax_l,&cmin_l) ;

		primtoflux(p_r,&state_r,dir,&geom,F_r) ;
		primtoflux(p_r,&state_r,TT, &geom,U_r) ;
		vchar(p_r,&state_r,&geom,dir,&cmax_r,&cmin_r) ;

		cmax = fabs(MY_MAX(MY_MAX(0., cmax_l),  cmax_r)) ;
		cmin = fabs(MY_MAX(MY_MAX(0.,-cmin_l), -cmin_r)) ;
		ctop = MY_MAX(cmax,cmin) ;


		PLOOP F[h][i][j][k] = 
			HLLF*(
			(cmax*F_l[k] + cmin*F_r[k] 
				- cmax*cmin*(U_r[k] - U_l[k]))/
					(cmax + cmin + SMALL) 
			) +
			LAXF*(
			0.5*(F_l[k] + F_r[k] 
				- ctop*(U_r[k] - U_l[k])) 
			) ;

            /* evaluate restriction on timestep */

            if(fullstep && i>=istart && i<=istop && j>=jstart && j<=jstop)
		    {
//            	cmax = MY_MAX(cmax,cmin) ;
//            	dtij = cour*dx[dir]/cmax ;
            	dtij = cour*dx[dir]/ctop ;

				if(dtij < ndt) ndt = dtij ;
			}

		}
	}

#if(RESCALE)
//	ZSLOOP(-2,N1+1,-2,N2+1) {
	ZSLOOP(0,0,-2,2,-2,2) {
		get_geometry(h,i,j,CENT,&geom) ;
		rescale(pr[h][i][j],REVERSE,dir,h,i,j,CENT,&geom) ;
	}
#endif

#ifdef MPI_USED
	if(fullstep)
		MPI_Allreduce(MPI_IN_PLACE,&ndt,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
#endif

	return(ndt) ;

}

#if N3 <= 2
/***********************************************************************************************/
/***********************************************************************************************
  flux_ct():
  ---------
     -- performs the flux-averaging used to preserve the del.B = 0 constraint (see Toth 2000);
        
***********************************************************************************************/
void flux_ct(double *_F1, double *_F2, double *_F3)
{
  double (*   F1)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_F1);
  double (*   F2)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_F2);
  double (*   F3)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_F3);
  
  double (*   emf)[SN1][SN2] = (double (*) [SN1][SN2])(_A);
	int h,i,j ;
//	static	double emf[N3+2][N1+4][N2+4] ;

//Flux-CT 3D is consistent with: "Dynamic Fisheye Grids for Binary Black Holes Simulations" Miguel Zilhao, Scott C. Noble


#ifdef MPI_USED
	MPI_Sendrecv(F2, 1, slice_iNPR1_stop, rank_nexti, 131,
                F2, 1, slice_iNPR1_startP, rank_previ, 131,
                MPI_COMM_WORLD, &MPIStatus);
#endif

	//EMF F1,F2,B1,B2
	/* calculate EMFs */
	/* Toth approach: just average */
	ZSLOOP(0,0,0,1,0,1) emf[h][i][j] = 0.25*(F1[h][i][j][B2] + F1[h][i][j-1][B2]
				  - F2[h][i][j][B1] - F2[h][i-1][j][B1]) ;
	

#ifdef MPI_USED
//	MPI_Sendrecv(emf, 1, slice_i1_stop, rank_nexti, 128,
//                emf, 1, slice_i1_startP, rank_previ, 128,
//                MPI_COMM_WORLD, &MPIStatus);
	MPI_Sendrecv(emf, 1, slice_i1_start, rank_previ, 132,
                emf, 1, slice_i1_stopN, rank_nexti, 132,
                MPI_COMM_WORLD, &MPIStatus);
#endif

	/* rewrite EMFs as fluxes, after Toth */
//        ZSLOOP(0,N1,0,N2-1) {
        ZSLOOP(0,0,0,1,0,0) {
                F1[h][i][j][B1] = 0. ;
                F1[h][i][j][B2] =  0.5*(emf[h][i][j] + emf[h][i][j+1]) ;
        }
//        ZSLOOP(0,N1-1,0,N2) {
        ZSLOOP(0,0,0,0,0,1) {
                F2[h][i][j][B1] = -0.5*(emf[h][i][j] + emf[h][i+1][j]) ;
                F2[h][i][j][B2] = 0. ;
	    }

		if(N3>2)
		{

#ifdef MPI_USED
			MPI_Sendrecv(F3, 1, slice_iNPR1_stop, rank_nexti, 133,
                F3, 1, slice_iNPR1_startP, rank_previ, 133,
                MPI_COMM_WORLD, &MPIStatus);

			MPI_Sendrecv(F1, 1, slice_hNPR1_stop, rank_nexth, 134,
                F1, 1, slice_hNPR1_startP, rank_prevh, 134,
                MPI_COMM_WORLD, &MPIStatus);

#endif
			//EMF F1,F3,B1,B3
			ZSLOOP(0,0,0,1,0,0) emf[h][i][j] = 0.25*(F3[h][i][j][B1] + F3[h][i-1][j][B1]
                                 - F1[h][i][j][B3] - F1[h-1][i][j][B3]) ;
	

#ifdef MPI_USED
			MPI_Sendrecv(emf, 1, slice_i1_start, rank_previ, 135,
                emf, 1, slice_i1_stopN, rank_nexti, 135,
                MPI_COMM_WORLD, &MPIStatus);

			MPI_Sendrecv(emf, 1, slice_h1_start, rank_prevh, 136,
                emf, 1, slice_h1_stopN, rank_nexth, 136,
                MPI_COMM_WORLD, &MPIStatus);
#endif

			/* rewrite EMFs as fluxes, after Toth */
    	    ZSLOOP(0,0,0,1,0,0) {
            	    F1[h][i][j][B1] = 0. ;
                	F1[h][i][j][B3] = -0.5*(emf[h][i][j] + emf[h+1][i][j]) ;
        	}
        	ZSLOOP(0,0,0,0,0,0) {
            	    F3[h][i][j][B1] = 0.5*(emf[h][i][j] + emf[h][i+1][j]) ;
                	F3[h][i][j][B3] = 0. ;
	    	}



#ifdef MPI_USED
			MPI_Sendrecv(F2, 1, slice_hNPR1_stop, rank_nexth, 137,
                F2, 1, slice_hNPR1_startP, rank_prevh, 137,
                MPI_COMM_WORLD, &MPIStatus);

#endif

			//EMF F2,F3,B2,B3
			ZSLOOP(0,0,0,0,0,1) emf[h][i][j] = 0.25*(F2[h][i][j][B3] + F2[h-1][i][j][B3]
				  - F3[h][i][j][B2] - F3[h][i][j-1][B2]) ;
	

#ifdef MPI_USED
			MPI_Sendrecv(emf, 1, slice_h1_start, rank_prevh, 138,
                emf, 1, slice_h1_stopN, rank_nexth, 138,
                MPI_COMM_WORLD, &MPIStatus);
#endif

			/* rewrite EMFs as fluxes, after Toth */
    	    ZSLOOP(0,0,0,0,0,1) {
            	    F2[h][i][j][B2] = 0. ;
                	F2[h][i][j][B3] = 0.5*(emf[h][i][j] + emf[h+1][i][j]) ;
        	}
        	ZSLOOP(0,0,0,0,0,0) {
            	    F3[h][i][j][B2] = -0.5*(emf[h][i][j] + emf[h][i][j+1]) ;
                	F3[h][i][j][B3] = 0. ;
	    	}
		}
}
#endif

#if 0
/***********************************************************************************************/
/***********************************************************************************************
  flux_ct():
  ---------
     -- performs the flux-averaging used to preserve the del.B = 0 constraint (see Toth 2000);
        
***********************************************************************************************/
void flux_ct(double *_F1, double *_F2, double *_F3)
{
  double (*   F1)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_F1);
  double (*   F2)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_F2);
  double (*   F3)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_F3);
  
  double (*   e1)[SN1][SN2] = (double (*) [SN1][SN2])(_A);
  double (*   e2)[SN1][SN2] = (double (*) [SN1][SN2])(_A2);
  double (*   e3)[SN1][SN2] = (double (*) [SN1][SN2])(_A3);
  int h,i,j ;

#ifdef MPI_USED
	MPI_Sendrecv(F1, 1, slice_iNPR1_stop, rank_nexti, 122,
                F1, 1, slice_iNPR1_startP, rank_previ, 122,
                MPI_COMM_WORLD, &MPIStatus);
	MPI_Sendrecv(F2, 1, slice_iNPR1_stop, rank_nexti, 123,
                F2, 1, slice_iNPR1_startP, rank_previ, 123,
                MPI_COMM_WORLD, &MPIStatus);
	MPI_Sendrecv(F3, 1, slice_iNPR1_stop, rank_nexti, 124,
                F3, 1, slice_iNPR1_startP, rank_previ, 124,
                MPI_COMM_WORLD, &MPIStatus);

	MPI_Sendrecv(F1, 1, slice_hNPR1_stop, rank_nexth, 125,
                F1, 1, slice_hNPR1_startP, rank_prevh, 125,
                MPI_COMM_WORLD, &MPIStatus);
	MPI_Sendrecv(F2, 1, slice_hNPR1_stop, rank_nexth, 126,
                F2, 1, slice_hNPR1_startP, rank_prevh, 126,
                MPI_COMM_WORLD, &MPIStatus);
	MPI_Sendrecv(F3, 1, slice_hNPR1_stop, rank_nexth, 127,
                F3, 1, slice_hNPR1_startP, rank_prevh, 127,
                MPI_COMM_WORLD, &MPIStatus);
#endif


#define E1(h,i,j)  (-F3[h][i][j][B2] + F2[h][i][j][B3])
#define E2(h,i,j)  (-F1[h][i][j][B3] + F3[h][i][j][B1])
#define E3(h,i,j)  (-F2[h][i][j][B1] + F1[h][i][j][B2])

#define E1j(h,i,j) (0.5*(E1(h,i,j) + E1(h,i,j-1)))
#define E1h(h,i,j) (0.5*(E1(h,i,j) + E1(h-1,i,j)))

#define E2i(h,i,j) (0.5*(E2(h,i,j) + E2(h,i-1,j)))
#define E2h(h,i,j) (0.5*(E2(h,i,j) + E2(h-1,i,j)))

#define E3i(h,i,j) (0.5*(E3(h,i,j) + E3(h,i-1,j)))
#define E3j(h,i,j) (0.5*(E3(h,i,j) + E3(h,i,j-1)))

   ZSLOOP(0,0,0,1,0,1)
   {
   	e1[h][i][j] = 0.25*( E1h(h,i,j) + E1h(h,i,j-1) 
                          + E1j(h,i,j) + E1j(h-1,i,j) );
   	e2[h][i][j] = 0.25*( E2i(h,i,j) + E2i(h-1,i,j) 
                          + E2h(h,i,j) + E2h(h,i-1,j) );
   	e3[h][i][j] = 0.25*( E3i(h,i,j) + E3i(h,i,j-1) 
                          + E3j(h,i,j) + E3j(h,i-1,j) );
   }


        ZSLOOP(0,0,0,1,0,1) {
                F1[h][i][j][B1] = 0. ;
                F2[h][i][j][B1] = -e3[h][i][j] ;
                F3[h][i][j][B1] = e2[h][i][j] ;
        }
        ZSLOOP(0,0,0,1,0,1) {
                F1[h][i][j][B2] = e3[h][i][j] ;
                F2[h][i][j][B2] = 0. ;
                F3[h][i][j][B2] = -e1[h][i][j] ;
        }
        ZSLOOP(0,0,0,1,0,1) {
                F1[h][i][j][B3] = -e2[h][i][j] ;
                F2[h][i][j][B3] = e1[h][i][j] ;
                F3[h][i][j][B3] = 0. ;
        }
  
}

#endif

#if N3 > 2
/***********************************************************************************************/
/***********************************************************************************************
  flux_ct():
  ---------
     -- performs the flux-averaging used to preserve the del.B = 0 constraint (see Toth 2000);
        
***********************************************************************************************/
void flux_ct(double *_F1, double *_F2, double *_F3)
{
  double (*   F1)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_F1);
  double (*   F2)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_F2);
  double (*   F3)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_F3);
  
  double (*   e1)[SN1][SN2] = (double (*) [SN1][SN2])(_A);
  double (*   e2)[SN1][SN2] = (double (*) [SN1][SN2])(_A2);
  double (*   e3)[SN1][SN2] = (double (*) [SN1][SN2])(_A3);
  int h,i,j ;

#ifdef MPI_USED
    MPI_Sendrecv(F2, 1, slice_iNPR1_stop, rank_nexti, 132,
                F2, 1, slice_iNPR1_startP, rank_previ, 132,
                MPI_COMM_WORLD, &MPIStatus);
    MPI_Sendrecv(F3, 1, slice_iNPR1_stop, rank_nexti, 133,
                F3, 1, slice_iNPR1_startP, rank_previ, 133,
                MPI_COMM_WORLD, &MPIStatus);

    MPI_Sendrecv(F1, 1, slice_hNPR1_stop, rank_nexth, 134,
                F1, 1, slice_hNPR1_startP, rank_prevh, 134,
                MPI_COMM_WORLD, &MPIStatus);
    MPI_Sendrecv(F2, 1, slice_hNPR1_stop, rank_nexth, 135,
                F2, 1, slice_hNPR1_startP, rank_prevh, 135,
                MPI_COMM_WORLD, &MPIStatus);
#endif


   ZSLOOP(0,0,0,1,0,1) {
     e1[h][i][j] = 0.25*((F2[h][i][j][B3] + F2[h-1][i][j][B3])
                          - (F3[h][i][j][B2] + F3[h][i][j-1][B2])
                            ) ;
     e2[h][i][j] = 0.25*((F3[h][i][j][B1] + F3[h][i-1][j][B1])
                          - (F1[h][i][j][B3] + F1[h-1][i][j][B3])
                            ) ;
     e3[h][i][j] = 0.25*((F1[h][i][j][B2] + F1[h][i][j-1][B2])
                          - (F2[h][i][j][B1] + F2[h][i-1][j][B1])
                           ) ;
   }

#ifdef MPI_USED
    MPI_Sendrecv(e2, 1, slice_i1_start, rank_previ, 136,
                e2, 1, slice_i1_stopN, rank_nexti, 136,
                MPI_COMM_WORLD, &MPIStatus);
    MPI_Sendrecv(e3, 1, slice_i1_start, rank_previ, 137,
                e3, 1, slice_i1_stopN, rank_nexti, 137,
                MPI_COMM_WORLD, &MPIStatus);

    MPI_Sendrecv(e1, 1, slice_h1_start, rank_prevh, 138,
                e1, 1, slice_h1_stopN, rank_nexth, 138,
                MPI_COMM_WORLD, &MPIStatus);
    MPI_Sendrecv(e2, 1, slice_h1_start, rank_prevh, 139,
                e2, 1, slice_h1_stopN, rank_nexth, 139,
                MPI_COMM_WORLD, &MPIStatus);
#endif

        ZSLOOP(0,0,0,1,0,0) {
                F1[h][i][j][B1] = 0. ;
                F1[h][i][j][B2] =  0.5*(e3[h][i][j] + e3[h][i][j+1]) ;
                F1[h][i][j][B3] = -0.5*(e2[h][i][j] + e2[h+1][i][j]) ;
        }
        ZSLOOP(0,0,0,0,0,1) {
                F2[h][i][j][B1] = -0.5*(e3[h][i][j] + e3[h][i+1][j]) ;
                F2[h][i][j][B2] = 0. ;
                F2[h][i][j][B3] =  0.5*(e1[h][i][j] + e1[h+1][i][j]) ;
        }
        ZSLOOP(0,0,0,0,0,0) {
                F3[h][i][j][B1] =  0.5*(e2[h][i][j] + e2[h][i+1][j]) ;
                F3[h][i][j][B2] = -0.5*(e1[h][i][j] + e1[h][i][j+1]) ;
                F3[h][i][j][B3] = 0. ;
        }
  
}

#endif


#if(DOKTOT)
void compute_ktot(double *_pi,double *_prh, double *_pr, int h, int i, int j, double Dt)
{
  double (*   pi)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_pi);
  double (*   prh)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_prh);
  double (*   pr)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_pr);
    
    double kappatot;

    kappatot =(gam-1.)*pr[h][i][j][UU]*pow(pr[h][i][j][RHO],-gam);
    
    //reset to actual internal energy/Entropy
    pr[h][i][j][KTOT] += (kappatot-pr[h][i][j][KTOT]);
    
    return;
}
/* this function sets the entropy
 * based on the initial conditions */
void init_entropy(){
    double (*   p)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_p);
    int i, j, h;
    ZSLOOP(0,0,-2,2,-2,2){
        p[h][i][j][KTOT] = (gam-1.)*p[h][i][j][UU]*pow(p[h][i][j][RHO],-gam);
    }
}
#endif
