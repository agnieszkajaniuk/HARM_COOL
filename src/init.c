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

/*
 *
 * generates initial conditions for a fishbone & moncrief disk
 * with exterior at minimum values for density & internal energy.
 *
 * cfg 8-10-01
 *
 */

#include "decs.h"
#include "random.h"
#include "assert.h"

#ifdef MPI_USED
#include "mpi.h"
#endif

#if( COOL_EOS )
#include "cooleos.h"
#endif

void coord_transform(double *pr,int h, int i, int j) ;
void debug2(double *pr, const char *where);


void init()
{
    int h,i,j ;
    double r,th,ph,sth,cth ;
    double ur,uh,up,u,rho ;
    double X[NDIM] ;
    struct of_geom geom ;

    /* for disk interior */
    double l,rin,lnh,expm2chi,up1 ;
    double DD,AA,SS,thin,sthin,cthin,DDin,AAin,SSin ;
    double kappa,hm1 ;

    /* for magnetic field */
//	double A[N3+2][N1+4][N2+4] ;
    double rho_av,rhomax,umax,beta,bsq_ij,bsq_max,norm,q,beta_act ;
    double rmax, lfish_calc(double rmax) ;


    /* some physics parameters */
    gam = 4./3. ;

    /* disk parameters (use fishbone.m to select new solutions) */
    //    a = 0.98;
    // rin =3.1 ; //for MBH=3 and Mdisk=0.1
    // rmax = 9.1 ;

    //new set of params

//    a = 0.98;
//    rin =3.1 ; //for MBH=3 and Mdisk=0.1
//    rmax = 9.1 ;
 

    //	a = 0.98;
    // rin =3.5 ; //for MBH=10 and Mdisk=1.0
    //	rmax = 9.8 ;

        a = 0.6; //for MBH=3, M=0.11
        rin = 4.0;
        rmax = 11.8;


    l = lfish_calc(rmax) ;

    kappa = 1.e-3 ;
    beta = 50.0 ;
//    beta = 100.0 ;

    /* some numerical parameters */
    lim = MC ;
    failed = 0 ;	/* start slow */
    cour = 0.9 ;
    dt = 1.e-5 ;
    R0 = 0.0 ;
    Rin = 0.98*(1. + sqrt(1. - a*a)) ;
//    Rout = 50. ;
 //   Rout = 1000. ;
    Rout = 100;

    t = 0. ;
    hslope = 0.3 ;

    set_arrays() ;
    set_grid() ;

    double (*   p)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_p);
    double (*   A)[SN1][SN2] = (double (*) [SN1][SN2])(_A);

    if(rank == 0)
    {
        coord(hstart,istart,jstart,CENT,X) ;
        bl_coord(X,&r,&th,&ph) ;
        fprintf(stderr,"rmin: " FMT_DBL_OUT "\n",r) ;
        fprintf(stderr,"rmin/rm: " FMT_DBL_OUT "\n",r/(1. + sqrt(1. - a*a))) ;

#if( COOL_EOS )
        printf("L_UNIT=%E, RHO_UNIT=%E, U_UNIT=%E, T_UNIT=%E\n", L_UNIT, RHO_UNIT, U_UNIT, T_UNIT);
#endif

    }

    /* output choices */
    tf = 20000.0 ; //for beta=50
   
    DTd = 10. ;	/* dumping frequency, in units of M */
    //DTd = 1.e-8 ;	/* dumping frequency, in units of M */
    // DTd = 2000. ;	/* dumping frequency, in units of M */
    DTl = 2. ;	/* logfile frequency, in units of M */
    DTi = 2.0E6 ; 	/* image file frequ., in units of M */
    //	DTi = 1.0 ; 	/* image file frequ., in units of M */
    DTr = 3600 ; 	/* restart file frequ., in real seconds */

    /* start diagnostic counters */
    dump_cnt = 0 ;
    image_cnt = 0 ;
    rdump_cnt = 0 ;
    defcon = 1. ;

    rhomax = 0. ;
    umax = 0. ;
/*
    ZSLOOPJ(0,0) {
      coord(0, 0, j, CENT, X);
      bl_coord(X, &r, &th, &ph);
	  printf("j=%3d th=%lf mpim=%E\n", j, th, M_PI-th);
	}
*/

//	ZSLOOP(0,N1-1,0,N2-1) {
    ZLOOP {
        coord(h,i,j,CENT,X) ;
        bl_coord(X,&r,&th,&ph) ;

        sth = sin(th) ;
        cth = cos(th) ;

        /* calculate lnh */
        DD = r*r - 2.*r + a*a ;
        AA = (r*r + a*a)*(r*r + a*a) - DD*a*a*sth*sth ;
        SS = r*r + a*a*cth*cth ;

        thin = M_PI/2. ;
        sthin = sin(thin) ;
        cthin = cos(thin) ;
        DDin = rin*rin - 2.*rin + a*a ;
        AAin = (rin*rin + a*a)*(rin*rin + a*a)
        - DDin*a*a*sthin*sthin ;
        SSin = rin*rin + a*a*cthin*cthin ;

        if(r >= rin) {
            lnh = 0.5*log((1. + sqrt(1. + 4.*(l*l*SS*SS)*DD/
            (AA*sth*AA*sth)))/(SS*DD/AA))
            - 0.5*sqrt(1. + 4.*(l*l*SS*SS)*DD/(AA*AA*sth*sth))
            - 2.*a*r*l/AA
            - (0.5*log((1. + sqrt(1. + 4.*(l*l*SSin*SSin)*DDin/
            (AAin*AAin*sthin*sthin)))/(SSin*DDin/AAin))
            - 0.5*sqrt(1. + 4.*(l*l*SSin*SSin)*DDin/
            (AAin*AAin*sthin*sthin))
            - 2.*a*rin*l/AAin ) ;
        }
        else
            lnh = 1. ;


        /* regions outside torus */
        if(lnh < SMALL || r < rin) {
            rho = 1.e-7*RHOMIN ;
            u = 1.e-7*UUMIN ;

            /* these values are demonstrably physical
               for all values of a and r */
            /*
                        ur = -1./(r*r) ;
                        uh = 0. ;
            up = 0. ;
            */

            ur = 0. ;
            uh = 0. ;
            up = 0. ;

            /*
            get_geometry(h,i,j,CENT,&geom) ;
                        ur = geom.gcon[0][1]/geom.gcon[0][0] ;
                        uh = geom.gcon[0][2]/geom.gcon[0][0] ;
                        up = geom.gcon[0][3]/geom.gcon[0][0] ;
            */

            p[h][i][j][RHO] = rho ;
            p[h][i][j][UU] = u ;
            p[h][i][j][U1] = ur ;
            p[h][i][j][U2] = uh ;
            p[h][i][j][U3] = up ;
        }
        /* region inside magnetized torus; u^i is calculated in
         * Boyer-Lindquist coordinates, as per Fishbone & Moncrief,
         * so it needs to be transformed at the end */
        else {
            hm1 = exp(lnh) - 1. ;
            rho = pow(hm1*(gam - 1.)/(kappa*gam),
            1./(gam - 1.)) ;
            u = kappa*pow(rho,gam)/(gam - 1.) ;
            ur = 0. ;
            uh = 0. ;

            /* calculate u^phi */
            expm2chi = SS*SS*DD/(AA*AA*sth*sth) ;
            up1 = sqrt((-1. + sqrt(1. + 4.*l*l*expm2chi))/2.) ;
            up = 2.*a*r*sqrt(1. + up1*up1)/sqrt(AA*SS*DD) +
            sqrt(SS/AA)*up1/sth ;


            p[h][i][j][RHO] = rho ;
            if(rho > rhomax) rhomax = rho ;
//			p[h][i][j][UU] = u*(1. + 4.e-2*(ranc(0)-0.5)) ;
            {
                unsigned int ii = i - (irank==0 ? 2 : 1) + di_rank_start;
                unsigned int jj = j - (jrank==0 ? 2 : 1) + dj_rank_start;
                unsigned int hh = N3 == 1 ? 0 : h - 1 + dh_rank_start;
//		hh = 0;

                RANDOMDEF rd = { 12345+hh*ii*jj+ii*jj+jj, 65435+hh*ii*jj+ii*jj+jj, 34221+hh*ii*jj+ii*jj+jj, 43543+hh*ii*jj+ii*jj+jj };
                p[h][i][j][UU] = u*(1. + 4.e-2*(UniformRandomRange(&rd, 0.0, 1.0)-0.5)) ;
            }
//			if(u > umax && r > rin) umax = u ;
            if(u > umax && r >= rin) umax = u ;
            p[h][i][j][U1] = ur ;
            p[h][i][j][U2] = uh ;

            p[h][i][j][U3] = up ;

            /* convert from 4-vel to 3-vel */
            coord_transform(p[h][i][j],h,i,j) ;
        }

        p[h][i][j][B1] = 0. ;
        p[h][i][j][B2] = 0. ;
        p[h][i][j][B3] = 0. ;


    }

#ifdef MPI_USED
    MPI_Allreduce(MPI_IN_PLACE,&rhomax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&umax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif

    /* Normalize the densities so that max(rho) = 1 */
    if(rank == 0)
        fprintf(stderr,"rhomax: " FMT_DBL_OUT "\n",rhomax) ;

//	ZSLOOP(0,N1-1,0,N2-1) {
    ZLOOP {
        p[h][i][j][RHO] /= rhomax ;
        p[h][i][j][UU]  /= rhomax ;
    }

    umax /= rhomax ;
    rhomax = 1. ;

    fixup((double*)p) ;

    bound_prim((double*)p) ;

#ifdef MPI_USED

    MPI_Sendrecv(p, 1, slice_iNPR1_stop, rank_nexti, 123,
                 p, 1, slice_iNPR1_startP, rank_previ, 123,
                 MPI_COMM_WORLD, &MPIStatus);

#endif

//	debug2((double *)p, "Init 1");

    /* first find corner-centered vector potential */
//	ZSLOOP(0,N1,0,N2) A[h][i][j] = 0. ;
    ZSLOOP(0,0,0,1,0,1) A[h][i][j] = 0. ;
//        ZSLOOP(0,N1,0,N2) {
    ZSLOOP(0,0,0,1,0,1) {
        /* vertical field version */
#if 0       
        coord(h,i,j,CENT,X) ;
        bl_coord(X,&r,&th,&ph) ;

		//DipoleField
        //A[h][i][j] = 0.5*sin(th)/(r*r) ;

		//Monopole
        A[h][i][j] = 1. - cos(th) ;
#else

        /* field-in-disk version */
        /* flux_ct */
        rho_av = 0.25*(
                     p[h][i][j][RHO] +
                     p[h][i-1][j][RHO] +
                     p[h][i][j-1][RHO] +
                     p[h][i-1][j-1][RHO]) ;

        q = rho_av/rhomax - 0.2 ;
        if(q > 0.) A[h][i][j] = q ;
#endif

    }

//    debug2((double *)p, "Init 1.01");


#ifdef MPI_USED
    MPI_Sendrecv(A, 1, slice_i1_start, rank_previ, 223,
                 A, 1, slice_i1_stopN, rank_nexti, 223,
                 MPI_COMM_WORLD, &MPIStatus);
#endif

    /* now differentiate to find cell-centered B,
       and begin normalization */
    bsq_max = 0. ;
    ZLOOP {
        get_geometry(h,i,j,CENT,&geom) ;

        /* flux-ct */
        p[h][i][j][B1] = -(A[h][i][j] - A[h][i][j+1]
        + A[h][i+1][j] - A[h][i+1][j+1])/(2.*dx[TH]*geom.g) ;
        p[h][i][j][B2] = (A[h][i][j] + A[h][i][j+1]
        - A[h][i+1][j] - A[h][i+1][j+1])/(2.*dx[RR]*geom.g) ;

        p[h][i][j][B3] = 0. ;

        bsq_ij = bsq_calc(p[h][i][j],&geom) ;
        if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
    }

#ifdef MPI_USED
    MPI_Allreduce(MPI_IN_PLACE,&bsq_max,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif


//    debug2((double *)p, "Init 1.0");

    if(rank == 0)
        fprintf(stderr,"initial bsq_max: " FMT_DBL_OUT "\n",bsq_max) ;

    /* finally, normalize to set field strength */
    beta_act = (gam - 1.)*umax/(0.5*bsq_max) ;
    if(rank == 0)
        fprintf(stderr,"initial beta: " FMT_DBL_OUT " (should be " FMT_DBL_OUT ")\n",beta_act,beta) ;
    norm = sqrt(beta_act/beta) ;
    bsq_max = 0. ;
    ZLOOP {
        p[h][i][j][B1] *= norm ;
        p[h][i][j][B2] *= norm ;

        get_geometry(h,i,j,CENT,&geom) ;
        bsq_ij = bsq_calc(p[h][i][j],&geom) ;
        if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
    }
#ifdef MPI_USED
    MPI_Allreduce(MPI_IN_PLACE,&bsq_max,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif

    beta_act = (gam - 1.)*umax/(0.5*bsq_max) ;

    if(rank == 0)
        fprintf(stderr,"final beta: " FMT_DBL_OUT " (should be " FMT_DBL_OUT ")\n",beta_act,beta) ;

//    debug2((double *)p, "Init 1.1");

    /* enforce boundary conditions */
    fixup((double*)p) ;

//    debug2((double *)p, "Init 1.2");

    bound_prim((double*)p) ;

//    debug2((double *)p, "Init 2");

#if( DO_FONT_FIX )
    set_Katm();
#endif

/*
	if(irank==0) ZSLOOPH(0,0) ZSLOOPJ(0,0) {
	    get_geometry(h,istart,j,CENT,&geom) ;
		double alpha = 1./sqrt(-geom.gcon[0][0]) ;
		fprintf(stdout, "j=%3d	alpha =%E\n", j, alpha);	
	}

	if(irank==isize-1) ZSLOOPH(0,0) ZSLOOPJ(0,0) {
	    get_geometry(h,istop,j,CENT,&geom) ;
		double alpha = 1./sqrt(-geom.gcon[0][0]) ;
		fprintf(stdout, "j=%3d	alpha =%E\n", j, alpha);	
	}
*/
    ZLOOP {
        if (isnan(p[h][i][j][UU])) {
            fprintf(stdout,"%d %d %d u=%g\n",h,i,j,p[h][i][j][UU]);
            exit(1);
        }
    }

    {
        double rhomin=1.e99, rhomax=0.0;
        double umin=1.e99, umax=0.0;

        ZLOOP {
            if(p[h][i][j][RHO] > rhomax) rhomax = p[h][i][j][RHO];
            if(p[h][i][j][RHO] < rhomin) rhomin = p[h][i][j][RHO];
            if(p[h][i][j][UU] > umax) umax = p[h][i][j][UU];
            if(p[h][i][j][UU] < umin) umin = p[h][i][j][UU];
        }

#ifdef MPI_USED
        MPI_Allreduce(MPI_IN_PLACE,&rhomin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,&umin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,&rhomax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,&umax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif
	    if(rank == 0)
	        printf("rhomin=%E, rhomax=%E, umin=%E, umax=%E\n", rhomin, rhomax, umin, umax);
    }

}


/* this version starts w/ BL 4-velocity and
 * converts to relative 4-velocities in modified
 * Kerr-Schild coordinates */

void coord_transform(double *pr,int hh, int ii, int jj)
{
  double X[NDIM],r,th,phi,ucon[NDIM],uconp[NDIM],trans[NDIM][NDIM],tmp[NDIM] ;
  double AA,BB,CC,discr ;
  double utconp[NDIM], dxdxp[NDIM][NDIM], dxpdx[NDIM][NDIM] ;
  struct of_geom geom ;
  struct of_state q ;
  int i,j,k,m ;

  coord(hh,ii,jj,CENT,X) ;
  bl_coord(X,&r,&th,&phi) ;
  blgset(hh,ii,jj,&geom) ;

  X[0] = 1;

  ucon[1] = pr[U1] ;
  ucon[2] = pr[U2] ;
  ucon[3] = pr[U3] ;

  AA =     geom.gcov[TT][TT] ;
  BB = 2.*(geom.gcov[TT][1]*ucon[1] +
           geom.gcov[TT][2]*ucon[2] +
           geom.gcov[TT][3]*ucon[3]) ;
  CC = 1. +
          geom.gcov[1][1]*ucon[1]*ucon[1] +
          geom.gcov[2][2]*ucon[2]*ucon[2] +
          geom.gcov[3][3]*ucon[3]*ucon[3] +
      2.*(geom.gcov[1][2]*ucon[1]*ucon[2] +
          geom.gcov[1][3]*ucon[1]*ucon[3] +
          geom.gcov[2][3]*ucon[2]*ucon[3]) ;

  discr = BB*BB - 4.*AA*CC ;
  ucon[TT] = (-BB - sqrt(discr))/(2.*AA) ;
  /* now we've got ucon in BL coords */

  /* transform to Kerr-Schild */
  /* make transform matrix */
  DLOOP trans[j][k] = 0. ;
  DLOOPA trans[j][j] = 1. ;
  trans[0][1] = 2.*r/(r*r - 2.*r + a*a) ;
  trans[3][1] = a/(r*r - 2.*r + a*a) ;

  /* transform ucon */
  DLOOPA tmp[j] = 0. ;
  DLOOP tmp[j] += trans[j][k]*ucon[k] ;
  DLOOPA ucon[j] = tmp[j] ;
  /* now we've got ucon in KS coords */

  /* transform to KS' coords */
  /* dr^\mu/dx^\nu jacobian, where x^\nu are internal coords */
  dxdxp_func(X, dxdxp);
  /* dx^\mu/dr^\nu jacobian */
  invert_matrix(dxdxp, dxpdx);
  
  for(i=0;i<NDIM;i++) {
    uconp[i] = 0;
    for(j=0;j<NDIM;j++){
      uconp[i] += dxpdx[i][j]*ucon[j];
    }
  }
  //old way of doing things for Gammie coords
  //ucon[1] *= (1./(r - R0)) ;
  //ucon[2] *= (1./(M_PI + (1. - hslope)*M_PI*cos(2.*M_PI*X[2]))) ;
  //ucon[3] *= 1.; //!!!ATCH: no need to transform since will use phi = X[3]

  get_geometry(hh, ii, jj, CENT, &geom);
  
  /* now solve for relative 4-velocity that is used internally in the code:
   * we can use the same u^t because it didn't change under KS -> KS' */
  ucon_to_utcon(uconp,&geom,utconp);
  
  pr[U1] = utconp[1] ;
  pr[U2] = utconp[2] ;
  pr[U3] = utconp[3] ;

  /* done! */
}

double lfish_calc(double r)
{
    return(
              ((pow(a,2) - 2.*a*sqrt(r) + pow(r,2))*
               ((-2.*a*r*(pow(a,2) - 2.*a*sqrt(r) + pow(r,2)))/
                sqrt(2.*a*sqrt(r) + (-3. + r)*r) +
                ((a + (-2. + r)*sqrt(r))*(pow(r,3) + pow(a,2)*(2. + r)))/
                sqrt(1 + (2.*a)/pow(r,1.5) - 3./r)))/
              (pow(r,3)*sqrt(2.*a*sqrt(r) + (-3. + r)*r)*(pow(a,2) + (-2. + r)*r))
          ) ;
}
