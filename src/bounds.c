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

void bound_x2dn_polefix( double *_prim );
void bound_x2up_polefix( double *_prim );


/* bound array containing entire set of primitive variables */

void bound_prim( double *_prim)
{
	double (*   prim)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_prim);
    char   (*pflag)[SN1][SN2] = (char (*) [SN1][SN2])(_pflag);

    int h,i,j,k ;
	void inflow_check(double *pr, int hh, int ii, int jj, int type );
#if( RESCALE )
        struct of_geom geom ;
#endif
        /* inner r boundary condition: u, gdet extrapolation */
//        for(j=0;j<N2;j++) {
		if(irank==0) ZSLOOPH(0,0) ZSLOOPJ(0,0) {
#if( RESCALE )
		get_geometry(h,istart,j,CENT,&geom) ;
		rescale(prim[istart][j],FORWARD, 1, h,istart,j,CENT,&geom) ;
#endif
                PLOOP prim[h][istart-1][j][k] = prim[h][istart][j][k] ;
                PLOOP prim[h][istart-2][j][k] = prim[h][istart][j][k] ;
                pflag[h][istart-1][j] = pflag[h][istart][j] ;
                pflag[h][istart-2][j] = pflag[h][istart][j] ;

#if( RESCALE )
		get_geometry(h,istart,j,CENT,&geom) ;
		rescale(prim[h][istart][j],REVERSE, 1, h,istart,j,CENT,&geom) ;
		get_geometry(h,istart-1,j,CENT,&geom) ;
		rescale(prim[h][istart-1][j],REVERSE, 1, h,istart-1,j,CENT,&geom) ;
		get_geometry(h,istart-2,j,CENT,&geom) ;
		rescale(prim[h][istart-2][j],REVERSE, 1, h,istart-2,j,CENT,&geom) ;
#endif

        }

        /* outer r BC: outflow */


//        for(j=0;j<N2;j++) {
		if(irank==isize-1) ZSLOOPH(0,0) ZSLOOPJ(0,0) {
#if( RESCALE )
		get_geometry(h,istop,j,CENT,&geom) ;
		rescale(prim[h][istop][j],FORWARD, 1, h,istop,j,CENT,&geom) ;
#endif

		PLOOP prim[h][istop+1  ][j][k] = prim[h][istop][j][k] ;
		PLOOP prim[h][istop+2][j][k] = prim[h][istop][j][k] ;
		pflag[h][istop+1  ][j] = pflag[h][istop][j] ;
		pflag[h][istop+2][j] = pflag[h][istop][j] ;


#if( RESCALE )
		get_geometry(h,istop,j,CENT,&geom) ;
		rescale(prim[h][istop][j],REVERSE, 1, h,istop,j,CENT,&geom) ;
		get_geometry(h,istop+1,j,CENT,&geom) ;
		rescale(prim[h][istop+1][j],REVERSE, 1, h,istop+1,j,CENT,&geom) ;
		get_geometry(h,istop+2,j,CENT,&geom) ;
		rescale(prim[h][istop+2][j],REVERSE, 1, h,istop+2,j,CENT,&geom) ;
#endif
        }

        /* make sure there is no inflow at the inner boundary */
//	for(i=-2;i<=-1;i++)  for(j=-2;j<N2+2;j++) { 
/*Original HARM, seems to be excessive check!!!
	if(irank==0) ZSLOOPH(0,0) for(i=-2;i<=-1;i++) ZSLOOPJ(-2,2) {
	  inflow_check(prim[h][istart-1][j],h,istart+i,j,0) ;
	  inflow_check(prim[h][istart-2][j],h,istart+i,j,0) ;
	}
*/
#if 0
	if(irank==0) ZSLOOPH(0,0) for(i=istart-2;i<=istart-1;i++) ZSLOOPJ(0,0) {
	  inflow_check(prim[h][i][j],h,i,j,0) ;
	}
#endif
        /* make sure there is no inflow at the outer boundary */
//        for(i=N1;i<=N1+1;i++)  for(j=-2;j<N2+2;j++) { 
/*Original HARM, seems to be excessive check!!!
	if(irank==isize-1) ZSLOOPH(0,0) for(i=istop+1;i<=istop+2;i++) ZSLOOPJ(-2,2) {
	  inflow_check(prim[h][istop+1  ][j],h,i,j,1) ;
	  inflow_check(prim[h][istop+2][j],h,i,j,1) ;
	}
*/
#if 0
	if(irank==isize-1) ZSLOOPH(0,0) for(i=istop+1;i<=istop+2;i++) ZSLOOPJ(0,0) {
	  inflow_check(prim[h][i][j],h,i,j,1) ;
	}
#endif


        /* polar BCs */
/*        for(i=-2;i<=N1+1;i++) { 
	  PLOOP {
                prim[h][i][-1  ][k] = prim[h][i][0   ][k] ;
                prim[h][i][-2  ][k] = prim[h][i][1   ][k] ;
                prim[h][i][N2  ][k] = prim[h][i][N2-1][k] ;
                prim[h][i][N2+1][k] = prim[h][i][N2-2][k] ;
	  }
	  pflag[h][i][-1  ] = pflag[h][i][0   ] ;
	  pflag[h][i][-2  ] = pflag[h][i][1   ] ;
	  pflag[h][i][N2  ] = pflag[h][i][N2-1] ;
	  pflag[h][i][N2+1] = pflag[h][i][N2-2] ;
        }
*/

#if(POLEFIX && POLEFIX < N2/2 && BL)
	bound_x2dn_polefix( _prim );
    bound_x2up_polefix( _prim );
#endif

#if N3 > 1
	  double (*   dq)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_dq);
	  static    MPI_Request request[4*N3];
	  static    MPI_Status status[4*N3];

      /* make sure b and u are antisymmetric at the poles */
	  if(jrank==0) ZSLOOPH(0,0) ZSLOOPI(-2,2) {
		  PLOOP {
                dq[h][i][jstart-1  ][k] = prim[h][i][jstart+0+THETAZEROPI   ][k] ;
                dq[h][i][jstart-2  ][k] = prim[h][i][jstart+1+THETAZEROPI   ][k] ;
		  }

		  pflag[h][i][jstart-1  ] = pflag[h][i][jstart+0+THETAZEROPI   ] ;
		  pflag[h][i][jstart-2  ] = pflag[h][i][jstart+1+THETAZEROPI   ] ;
      }

      /* make sure b and u are antisymmetric at the poles */
	  if(jrank==jsize-1) ZSLOOPH(0,0) ZSLOOPI(-2,2) {
		  PLOOP {
                dq[h][i][jstop+1  ][k] = prim[h][i][jstop-0-THETAZEROPI][k] ;
                dq[h][i][jstop+2][k] = prim[h][i][jstop-1-THETAZEROPI][k] ;
		  }

		  pflag[h][i][jstop+1  ] = pflag[h][i][jstop-0-THETAZEROPI] ;
		  pflag[h][i][jstop+2] = pflag[h][i][jstop-1-THETAZEROPI] ;
      }


	  int rnum = 0;
	  if(jrank==0) ZSLOOPH(0,0) {
		if(rank != slice_h_rankop[h])
		{
	      MPI_Isend(dq, //buffer that's being sent
    	            1,              //number of items sent
        	        slice_hopNPR2_m2[h],         //data type
            	    slice_h_rankop[h],//the rank of destination process
                	1000 + slice_h_hop[h],                //tag
	               	MPI_COMM_WORLD,     //communicator
    	           	&request[rnum++] //error
                	);


      		MPI_Irecv(prim, //buffer that's being received
        	        1,              //number of items received (same as those sent)
            	    slice_hopNPR2_m2[h],         //data type
                	slice_h_rankop[h],//the rank of source process
                	1000 + h,                //tag (should be same as in the send process)
                	MPI_COMM_WORLD,     //communicator
                	&request[rnum++] //error
                )	;
		} else {
			  ZSLOOPI(-2,2)  PLOOP {
                prim[h][i][jstart-2   ][k] = dq[slice_h_hop[h]][i][jstart-2  ][k];
                prim[h][i][jstart-1   ][k] = dq[slice_h_hop[h]][i][jstart-1  ][k];
		  }

		}
		
	  }

	  if(jrank==jsize-1) ZSLOOPH(0,0) {

		if(rank != slice_h_rankop[h])
		{
	      MPI_Isend(dq, //buffer that's being sent
    	            1,              //number of items sent
        	        slice_hopNPR2_n1[h],         //data type
            	    slice_h_rankop[h],//the rank of destination process
                	2000 + slice_h_hop[h],                //tag
	               	MPI_COMM_WORLD,     //communicator
    	           	&request[rnum++] //error
                	);

      		MPI_Irecv(prim, //buffer that's being received
        	        1,              //number of items received (same as those sent)
            	    slice_hopNPR2_n1[h],         //data type
                	slice_h_rankop[h],//the rank of source process
                	2000 + h,                //tag (should be same as in the send process)
                	MPI_COMM_WORLD,     //communicator
                	&request[rnum++] //error
                )	;
		} else {
			  ZSLOOPI(-2,2)  PLOOP {
                prim[h][i][jstop+1   ][k] = dq[slice_h_hop[h]][i][jstop+1  ][k];
                prim[h][i][jstop+2   ][k] = dq[slice_h_hop[h]][i][jstop+2  ][k];
		  }
		}
	  }


          MPI_Waitall(rnum, request, status);




#if !GHOSTZONESIGNCHANGEINMETRIC
          if(jrank==0) ZSLOOPH(0,0) ZSLOOPI(-2,2) for(j=jstart-2;j<jstart;j++) {
                        prim[h][i][j][U2] *= -1. ;
                        prim[h][i][j][B2] *= -1. ;
// Axisymmetric boundary, Same as reflective, except for the angular component of vphi or Bphi, which also changes sign
#if JBOUNDAXISYMMETRIC
                        prim[h][i][j][U3] *= -1. ;
                        prim[h][i][j][B3] *= -1. ;
#endif
                      }

          if(jrank==jsize-1) ZSLOOPH(0,0) ZSLOOPI(-2,2) for(j=jstop+1;j<=jstop+2;j++) {
                        prim[h][i][j][U2] *= -1. ;
                        prim[h][i][j][B2] *= -1. ;
#if JBOUNDAXISYMMETRIC
                        prim[h][i][j][U3] *= -1. ;
                        prim[h][i][j][B3] *= -1. ;
#endif
                }
#endif

#else //#if N3 > 1

	  if(jrank==0) ZSLOOPH(0,0) ZSLOOPI(-2,2) {
		  PLOOP {
                prim[h][i][jstart-1  ][k] = prim[h][i][jstart+0+THETAZEROPI][k] ;
                prim[h][i][jstart-2  ][k] = prim[h][i][jstart+1+THETAZEROPI][k] ;
		  }
		  pflag[h][i][jstart-1  ] = pflag[h][i][jstart+0+THETAZEROPI] ;
		  pflag[h][i][jstart-2  ] = pflag[h][i][jstart+1+THETAZEROPI] ;
      }

	  if(jrank==jsize-1) ZSLOOPH(0,0) ZSLOOPI(-2,2) {
		  PLOOP {
                prim[h][i][jstop+1  ][k] = prim[h][i][jstop-0-THETAZEROPI][k] ;
                prim[h][i][jstop+2][k] = prim[h][i][jstop-1-THETAZEROPI][k] ;
		  }
		  pflag[h][i][jstop+1  ] = pflag[h][i][jstop-0-THETAZEROPI] ;
		  pflag[h][i][jstop+2] = pflag[h][i][jstop-1-THETAZEROPI] ;
      }

        /* make sure b and u are antisymmetric at the poles */

#if !GHOSTZONESIGNCHANGEINMETRIC
          if(jrank==0) ZSLOOPH(0,0) ZSLOOPI(-2,2) for(j=jstart-2;j<jstart;j++) {
                        prim[h][i][j][U2] *= -1. ;
                        prim[h][i][j][B2] *= -1. ;

#if JBOUNDAXISYMMETRIC
                        prim[h][i][j][U3] *= -1. ;
                        prim[h][i][j][B3] *= -1. ;
#endif
                      }

          if(jrank==jsize-1) ZSLOOPH(0,0) ZSLOOPI(-2,2) for(j=jstop+1;j<=jstop+2;j++) {
                        prim[h][i][j][U2] *= -1. ;
                        prim[h][i][j][B2] *= -1. ;

#if JBOUNDAXISYMMETRIC
                        prim[h][i][j][U3] *= -1. ;
                        prim[h][i][j][B3] *= -1. ;
#endif
                }
#endif

#endif

	if(irank==0) ZSLOOPH(0,0) for(i=istart-2;i<=istart-1;i++) ZSLOOPJ(-2,2) {
	  inflow_check(prim[h][i][j],h,i,j,0) ;
	}

	if(irank==isize-1) ZSLOOPH(0,0) for(i=istop+1;i<=istop+2;i++) ZSLOOPJ(-2,2) {
	  inflow_check(prim[h][i][j],h,i,j,1) ;
	}

}

void inflow_check(double *pr, int hh, int ii, int jj, int type )
{
        struct of_geom geom ;
        double ucon[NDIM] ;
        int j,k ;
        double alpha,beta1,gamma,vsq ;

        get_geometry(hh,ii,jj,CENT,&geom) ;
        ucon_calc(pr, &geom, ucon) ;

        if( ((ucon[1] > 0.) && (type==0)) || ((ucon[1] < 0.) && (type==1)) ) { 
                /* find gamma and remove it from primitives */
	  if( gamma_calc(pr,&geom,&gamma) ) { 
	    fflush(stderr);
	    fprintf(stderr,"\ninflow_check(): gamma failure \n");
	    fflush(stderr);
	    fail(FAIL_GAMMA);
	  }
	    pr[U1] /= gamma ;
	    pr[U2] /= gamma ;
	    pr[U3] /= gamma ;
	    alpha = 1./sqrt(-geom.gcon[0][0]) ;
	    beta1 = geom.gcon[0][1]*alpha*alpha ;

	    /* reset radial velocity so radial 4-velocity
	     * is zero */
	    pr[U1] = beta1/alpha ;

	    /* now find new gamma and put it back in */
	    vsq = 0. ;
	    SLOOP vsq += geom.gcov[j][k]*pr[U1+j-1]*pr[U1+k-1] ;
	    if( fabs(vsq) < 1.e-13 )  vsq = 1.e-13;
	    if( vsq >= 1. ) { 
	      vsq = 1. - 1./(GAMMAMAX*GAMMAMAX) ;
	    }
	    gamma = 1./sqrt(1. - vsq) ;
	    pr[U1] *= gamma ;
	    pr[U2] *= gamma ;
	    pr[U3] *= gamma ;

	    /* done */
	  }
	  else
	    return ;

	}


void bound_x2dn_polefix( double *_prim )
{
  double (*   prim)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_prim);
  int h,i,j,k,jref ;
  
#if(POLEFIX && POLEFIX < N2/2 && BL)
  //copy all densities and B^phi in; interpolate linearly transverse velocity
  jref = POLEFIX;
  if(jrank==0) ZSLOOPH(0,0) ZSLOOPI(-2,2) {
      for(j=0;j<jref;j++) {
        PLOOP {
          if(k==B1 || k==B2 || (N3>1 && k==B3))
            //don't touch magnetic fields
            continue;
          else if(k==U2) {
            //linear interpolation of transverse velocity (both poles)
            prim[h][i][jstart+j][k] = (j+0.5)/(jref+0.5) * prim[h][i][jstart+jref][k];
          }
          else {
            //everything else copy (both poles)
            prim[h][i][jstart+j][k] = prim[h][i][jstart+jref][k];
          }
        }
      }
  }
#endif
}

void bound_x2up_polefix( double *_prim )
{
  double (*   prim)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_prim);
  int h,i,j,k,jref ;

#if(POLEFIX && POLEFIX < N2/2 && BL)
  //copy all densities and B^phi in; interpolate linearly transverse velocity
  jref = POLEFIX;
  if(jrank==jsize-1) ZSLOOPH(0,0) ZSLOOPI(-2,2) { 
      for(j=0;j<jref;j++) {
        PLOOP {
          if(k==B1 || k==B2 || (N3>1 && k==B3))
            //don't touch magnetic fields
            continue;
          else if(k==U2) {
            //linear interpolation of transverse velocity (both poles)
            prim[h][i][jstop-j][k] = (j+0.5)/(jref+0.5) * prim[h][i][jstop-jref][k];
          }
          else {
            //everything else copy (both poles)
            prim[h][i][jstop-j][k] = prim[h][i][jstop-jref][k];
          }
        }
      }
  }
#endif
  
}
