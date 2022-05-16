/* 
-------------------------------------------------------------------------------
    Copyright 2005 Scott C. Noble, Charles F. Gammie, 
                   Jonathan C. McKinney, and Luca Del Zanna


    This file is part of PVS-GRMHD.

    PVS-GRMHD is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    PVS-GRMHD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with PVS-GRMHD; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

-------------------------------------------------------------------------------
*/

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************

utoprim_5d.c: 
---------------

    Uses the 5D method: 
       -- solves for the independent variables rho,u,\tilde{u}^i using a 5D
          Newton-Raphson method 
       -- can be used (in principle) with a general equation of state, but this 
          implementation is not designed with the ability to easily switch EOS's;

  -- Currently returns with an error state (>0) if a negative rest-mass
      density or internal energy density is calculated.  You may want 
      to change this aspect of the code so that it still calculates the 
      velocity and so that you can floor the densities.  If you want to 
      change this aspect of the code please comment out the "return(retval)"
      statement after "retval = 5;" statement in Utoprim_new_body();


   -- Note: many of these routines (especially dudp_calc() routines) were 
            written and used for HARM by Jon McKinney and Charles Gammie.
            -- changes were made to make it more compatible and consistent 
               with the newer Utoprim*() routines;

******************************************************************************/


#include "decs.h"
#include "u2p_util.h"

#if( COOL_EOS )
#include "cooleos.h"
#endif

#define NEWT_DIM 5

#define NORMMETHOD (1)

// Tolerance for residual function to prevent trunc. error problems with 
// the backward substitution procedure:
#define TOL_RESID (1.e-4)    

#if 0
#define TT 0

/* loop over Primitive variables */
#define PLOOP for(k=0;k<NPR;k++)
/* loop over all Dimensions; second rank loop */
#define DLOOP for(j=0;j<NDIM;j++)for(k=0;k<NDIM;k++)
/* loop over all Dimensions; first rank loop */
#define DLOOPA for(j=0;j<NDIM;j++)
/* loop over all Space dimensions; second rank loop */
#define SLOOP for(j=1;j<NDIM;j++)for(k=1;k<NDIM;k++)
/* loop over all Space dimensions; first rank loop */
#define SLOOPA for(j=1;j<NDIM;j++)
#endif

/* loop over all for j and Space for k; second rank loop */
#define DSLOOP for(j=0;j<NDIM;j++)for(k=1;k<NDIM;k++)
/* loop over all for k and Space for j; second rank loop */
#define SDLOOP for(j=1;j<NDIM;j++)for(k=0;k<NDIM;k++)

FTYPE U_target[NPR];
FTYPE glob_gcov[NDIM][NDIM];
FTYPE glob_gcon[NDIM][NDIM];
FTYPE glob_gdet;


  static int general_newton_raphson( FTYPE x[], int n, 
				     void (*funcd) (FTYPE [], FTYPE [], FTYPE [], 
						    FTYPE [][NEWT_DIM], FTYPE *, 
						    FTYPE *, int) );

  
  static void func_5d(FTYPE prguess[], FTYPE dx[], FTYPE resid[],
		      FTYPE jac[][NEWT_DIM],
		      FTYPE *f, FTYPE *df, int n);

  
  static int LU_decompose( FTYPE A[][NEWT_DIM], int permute[] );
  static void LU_substitution( FTYPE A[][NEWT_DIM], FTYPE B[], int permute[] );
  static int dudp_calc_g(FTYPE *pr, FTYPE Am[][NEWT_DIM] );

  static void dgdvi_calc_g(FTYPE *pr,FTYPE *dgdvi);
  static void duidvj_calc_g(FTYPE *dgdvi,FTYPE duidvj[][NDIM]);
  static void dutdui_calc_g(FTYPE *ucov,FTYPE *dutdui) ;

  static void dbtdui_calc_g(FTYPE *dutdui, FTYPE *pr, FTYPE *dbtdui) ;
  static void dbiduj_calc_g(FTYPE *dbtdui,FTYPE *dutdui,FTYPE *ucon, FTYPE *b, 
			    FTYPE dbiduj[][NDIM]) ;
  static void dbsqdui_calc_g(FTYPE dbiduj[][NDIM],FTYPE *bcov, FTYPE *dbsqdui) ;
  static void dudduu_calc_g(FTYPE*dutdui, FTYPE dudduu[][NDIM]) ;
  static void dbdiduj_calc_g(FTYPE dbiduj[][NDIM],FTYPE dbdiduj[][NDIM]);

  static void calc_errx_5d( FTYPE *W, FTYPE *vsq, FTYPE x[] );



/**********************************************************************/ 
/******************************************************************

  Utoprim_5d():

  -- solves for the primitive variables, "prim[]", associated with the 
     conserved variables U[];

     It assumes that on input/output:


              /  rho u^t           \
         U =  |  T^t_t   + rho u^t |  sqrt(-det(g_{\mu\nu}))
              |  T^t_\mu           | 
              \   B^i              / 

             /    rho        \
	 P = |    uu         |
             | \tilde{u}^i   |
             \   B^i         /


     ala HARM. 

   Arguments:
       U[NPR]    = conserved variables (current values on input/output);
       gcov[NDIM][NDIM] = covariant form of the metric ;
       gcon[NDIM][NDIM] = contravariant form of the metric ;
       gdet             = sqrt( - determinant of the metric) ;
       prim[NPR] = primitive variables (guess on input, calculated values on 
                                        output if there are no problems);
  
   -- NOTE: for those using this routine for special relativistic MHD and are
            unfamiliar with metrics, merely set 
              gcov = gcon = diag(-1,1,1,1)  and gdet = 1.  ;

return:  (i*100 + j)  where 
         i = 0 ->  Newton-Raphson solver either was not called (yet or not used) 
                   or returned successfully;
             1 ->  Newton-Raphson solver did not converge to a solution with the 
                   given tolerances;
             2 ->  Newton-Raphson procedure encountered a numerical divergence 
                   (occurrence of "nan" or "+/-inf" ;
	     
         j = 0 -> success 
             1 -> failure: some sort of failure in Newton-Raphson; 
             2 -> failure: utsq<0 w/ initial p[] guess;
	     3 -> failure: W<0 or W>W_TOO_BIG
             4 -> failure: v^2 > 1 
             5 -> failure: rho,uu <= 0 ;

******************************************************************/

int Utoprim_5d(FTYPE U[NPR], FTYPE gcov[NDIM][NDIM], FTYPE gcon[NDIM][NDIM], 
	       FTYPE gdet, FTYPE prim[NPR])
{

  FTYPE prim_tmp[NPR], gamma_tmp;
  int ret, i,j,k, ret_gam; 


  // set global geometric functions:
  glob_gdet = gdet;
  for( i = 0; i < NDIM ; i++ ) { 
    for( j = 0; j < NDIM ; j++ ) { 
      glob_gcov[i][j] = gcov[i][j];
      glob_gcon[i][j] = gcon[i][j];
    }
  }


  // Set global conserved variable array:
  PLOOP U_target[k] = U[k];

  /* Set the primitive B-fields (same as conserved B-fields): */
  for (k = BCON1; k <= BCON3; k++)
    prim[k] = U[k] / gdet;  


  PLOOP prim_tmp[k] = prim[k];

  // Solve the 5D linear system with Newton's method:
  ret = general_newton_raphson(prim_tmp, NPR - 3, func_5d);

  if( ret ){
    ret = 100*ret + 1;
    return(ret); 

  }// otherwise got good solution
  else{
    // check densities for positivity
    if((prim_tmp[RHO]<0.0)||(prim_tmp[UU]<0.0)){
      ret = 5;
      // mostly we hit pr[UU]<0, or maybe never pr[RHO]<0 and pr[UU]<0 alot
      return(ret);
    }

    //Calculate gamma, Lorentz factor, to see if we have a superluminal velocity:
    ret_gam = gamma_calc_g(prim_tmp, gcov, &gamma_tmp );
    
    if( ret_gam == 1 ) { 
      ret = 4;
      return(ret);
    }
  }

  // If we got here, then "everything" is ok : 
  ret = 0;

  PLOOP prim[k] = prim_tmp[k];

  return( ret );

}

/*******************************************************************/ 
/*******************************************************************
 
   func_5d(): 

    --  Auxiliary function required by general_newton_raphson()  
        to calculate the residual, Newton step, etc. ; 

*******************************************************************/
static void func_5d(FTYPE prguess[], FTYPE dxm[], FTYPE resid[], 
		    FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n)
{
  FTYPE prstart[NPR];
  FTYPE pr[NPR];
  static FTYPE U_curr[NPR];
  
  int i = 0, j = 0, k = 0;
  int failreturn=0;
  FTYPE normtmp;
  FTYPE norm;
  FTYPE residnorm;
  FTYPE d;

  int numnormterms;

  FTYPE  max_subs_change, subs_tmp;
  static int permute[NEWT_DIM];
  static FTYPE Am[NEWT_DIM][NEWT_DIM];
  static FTYPE beta[NEWT_DIM];


  // Initialized LU_deocmpose() quantities:
  for( i = 0 ; i < NEWT_DIM ; i++ ) { 
    permute[i] = 0 ;
  }


  for( i = 0;  i <  BCON1; i++)   pr[i]  =  prguess[i];
  for( i = BCON1; i <= BCON3; i++)   pr[i]  =  U_target[i] / glob_gdet;

  // store old pr
  PLOOP prstart[k]=pr[k];

  // Find the conserved variables associated with the current guess for the 
  //  primitive variables:
  primtoU_g(pr, glob_gcov, glob_gcon, glob_gdet, U_curr);

  // Calculate the jacobian of the forward transformation:
  failreturn=dudp_calc_g(pr, Am);

  // normalize error = beta to "average" of jacobian elements or to \rho u^t 
#if(NORMMETHOD==0)
  norm=1.0/U_target[RHO];
#elif(NORMMETHOD==1)
  norm=0.0;
  numnormterms=0;
  for (j = 0; j < NPR - 3; j++){
    for (k = 0; k < NPR - 3; k++){
      if(Am[j][k] > NUMEPSILON){
	  norm+=fabs( Am[j][k] );
	  numnormterms++;
      }
    }
  }
  norm=((FTYPE)(numnormterms))/(norm); // (i.e. inverse of average)
#elif(NORMMETHOD==2)
  norm = 1.0;
#endif

  //Normalize Jacobian:
  for (j = 0; j < NPR - 3; j++) {
    for (k = 0; k < NPR - 3; k++) {
      jac[j][k] = Am[j][k];
      Am[j][k] *= (norm);
    }
  }
  
  // determine normalized error
  *df = residnorm = 0.;
  for (k = 0; k < NPR - 3; k++) {
    resid[k] =  U_curr[k] - U_target[k] ;
    beta[k] = -resid[k] *(norm);
    residnorm = 0.5*(fabs(U_curr[k]) + fabs(U_target[k]));
    residnorm = (residnorm == 0.) ? 1. : 1./residnorm;
    *df -= residnorm * resid[k]*resid[k];
  }

  *f = -0.5*(*df);

  // Need to return with normalization factor that normalizes U[] to 1 :
  normtmp = 0.;
  for (k = 0; k < NPR - 3; k++)
    normtmp += fabs(U_target[k]) ;

  normtmp *= norm;
  norm = normtmp;


  /*******************************************/
  /* Solve the linear system:                */
  /*******************************************/

  /*********************************************************/
  /**** With simple LU-decomp. and backward/forward subst. */

  for( i = 0 ; i < n ; i++ ) { 
    dxm[i] = beta[i];
  }

  // Get the LU matrix:
  if( LU_decompose( Am,  permute ) != 0  ) { 
    fprintf(stderr, "func_5d(): singular matrix encountered! \n");
    
  }

  // Solve the linear system: 
  LU_substitution( Am,  dxm, permute );

}


/***************************************************************/
/***************************************************************
   
   calc_errx_5d(): 
  
     -- finds W and v^2 to be used by general_newton_raphson()
        to calculate the error function, which determines when we 
        have converged to a solution;

***************************************************************/
static void calc_errx_5d( FTYPE *W, FTYPE *vsq, FTYPE x[] ) 
{
  int id, jd;
  FTYPE gsqtmp;


  // First find v^2:
  *vsq = 0. ;
  for(id=1;id<4;id++)
    for(jd=1;jd<4;jd++) 
      *vsq += glob_gcov[id][jd]*x[UTCON1+id-1]*x[UTCON1+jd-1] ;      
  
  gsqtmp = 1. + *vsq;
  *vsq /= gsqtmp;


#if( COOL_EOS )
  *W = gsqtmp * ( x[UU] + x[RHO] + pressure_rho0_u_CoolEOS(x[RHO],x[UU]) ) ;
#else
  *W = gsqtmp * ( GAMMA*x[UU] + x[RHO] ) ;
#endif  
}



/***************************************************************/
/***************************************************************

  dudp_calc_g():

       --  calculate du/dp analytically.  
       --  pr is the full (NPR) element primitive variables
       --  Am is be a 5x5 matrix containing dU(1-5)/dp(1-5).
            used only by Utoprim_5d();

***************************************************************/

static int dudp_calc_g(FTYPE *pr, FTYPE Am[][NEWT_DIM] )
{
  static FTYPE dutdui[NDIM] ;
  static FTYPE dbtdui[NDIM] ;
  static FTYPE dbsqdui[NDIM] ;
  static FTYPE dbiduj[NDIM][NDIM] ;
  static FTYPE dgdvi[NDIM];
  static FTYPE duidvj[NDIM][NDIM];
  static FTYPE dudduu[NDIM][NDIM];
  static FTYPE dbdiduj[NDIM][NDIM];
  static FTYPE tmp1[NDIM],tmp2[NDIM] ;	
  FTYPE eta,bsq ;
  int i,j,k,l ;
  static FTYPE tempA[NPR][NPR];
  static FTYPE bcon[NDIM], bcov[NDIM], ucov[NDIM], ucon[NDIM], ncov[NDIM];


  for( j=0; j < NEWT_DIM; j++)
    for( k=0; k < NEWT_DIM; k++){ tempA[j][k] = Am[j][k] = 0. ;}


  // Calculate: bcon, bcov, ucon, ucov
  ucon_calc_g(pr,glob_gcov,glob_gcon,ucon) ;
  lower_g(ucon,glob_gcov,ucov) ;
  ncov_calc(glob_gcon,ncov) ;
  bcon_calc_g(pr,ucon,ucov,ncov,bcon) ;
  lower_g(bcon,glob_gcov,bcov) ;

  // Calculate other auxiliary arrrays and parameters:

  // b^\mu b_\mu : 
  bsq = dot(bcon, bcov);
  // d u^t / d u^i
  dutdui_calc_g(ucov,dutdui) ;
  // d b^t / d u^i
  dbtdui_calc_g(dutdui,pr,dbtdui) ;
  // d b^{mu} / d u^j
  dbiduj_calc_g(dbtdui,dutdui,ucon,bcon,dbiduj) ;
  // d b^2 / d u^i
  dbsqdui_calc_g(dbiduj,bcov,dbsqdui) ;
  // d u^{mu} / d u^j
  dudduu_calc_g(dutdui,dudduu);
  // d b_{mu} / d u^j
  dbdiduj_calc_g(dbiduj,dbdiduj);

#if( COOL_EOS )
  eta = pr[RHO] + pr[UU] + pressure_rho0_u_CoolEOS(pr[RHO], pr[UU]) + bsq ;
#else
  eta = pr[RHO] + GAMMA*pr[UU] + bsq ;
#endif

  ///////////////////
  // now define Am
  //
  // Am is dU^i/dp^j = Am[i][j]

  // rho on rho
  Am[RHO][RHO] = ucon[TT] ;

  // u+momentums on rho
#if( COOL_EOS )
  double dpdrho = dpdrho_calc_CoolEOS(pr[RHO], pr[UU]);
  DLOOPA Am[UU+j][RHO] = (1.0 + dpdrho) * ucon[TT] * ucov[j] + delta(TT,j) * dpdrho;
#else
  DLOOPA Am[UU+j][RHO] = ucon[TT] * ucov[j];
#endif

  // U[rho] on u
  Am[RHO][UU] = 0.0;
  // u+momentums on u
#if( COOL_EOS )
  double dpdu = dpdu_calc_CoolEOS(pr[RHO], pr[UU]);
  DLOOPA Am[UU+j][UU] = (1.0 + dpdu) * ucon[TT] * ucov[j]  +  delta(TT,j) * dpdu;
#else
  DLOOPA Am[UU+j][UU] = GAMMA * ucon[TT] * ucov[j]  +  delta(TT,j) * (GAMMA-1.0);
#endif
  // rho on momentums
  SLOOPA Am[RHO][UTCON1+j-1]=pr[RHO]*dutdui[j];

  DSLOOP{//D=j S=k // u+momentum's on velocities
    Am[UU+j][UTCON1+k-1] = 
      dbsqdui[k] * ucon[TT] * ucov[j]
      +eta * (dutdui[k] * ucov[j] + ucon[TT] * dudduu[j][k])
      +delta(j,0) * 0.5 * dbsqdui[k]
      -(dbtdui[k] * bcov[j] + bcon[TT] * dbdiduj[j][k]);
  }

  /* this bit of legacy code can be uncommented if
     the rest-mass flux is added to the energy 
     flux, which may be numerically convenient  */
  for(k=0; k < NEWT_DIM; k++) Am[UU][k] += Am[RHO][k] ;



  // change of variables
  // d\gamma^2 / d V^i
  dgdvi_calc_g(pr,dgdvi);
  // d u^{mu} / d V^j (only  mu=1,2,3 used)
  duidvj_calc_g(dgdvi,duidvj);

  // convert to relative velocity (over all Am[U=RHO->UTCON3][p=v1,v2,v3])
  // 4vel -> rel4vel
  for(i=0; i < NEWT_DIM; i++) for(j=RHO; j<=UU; j++) tempA[i][j] = Am[i][j];

  for(i=0; i < NEWT_DIM; i++) SLOOPA for(l=1;l<=3;l++) 
    tempA[RHO+i][UTCON1+j-1] += Am[RHO+i][UTCON1+l-1]*duidvj[l][j];

  for(i=0; i < NEWT_DIM; i++) for(j=0; j < NEWT_DIM; j++) Am[i][j] = tempA[i][j];
	
	  
  // Am is dU^i/dp^j = Am[i][j]

  // N.B.: all the conserved variables contain a factor 
  // of \sqrt{det(g_{\mu\nu})} 
  for( j=0; j < NEWT_DIM; j++)
    for( k = 0; k < NEWT_DIM; k++) 
      Am[j][k] *= glob_gdet ;


  return(0) ;
}

static void dgdvi_calc_g(FTYPE *pr, FTYPE *dgdvi)
{
  int j,k;
  FTYPE gamma;

  gamma_calc_g(pr, glob_gcov, &gamma) ;

  SLOOPA dgdvi[j]=0.0;

  SLOOP dgdvi[j] += 1.0 / gamma * glob_gcov[j][k] * pr[UTCON1+k-1];

  // no such dgdvi[TT]

}

static void duidvj_calc_g(FTYPE *dgdvi, FTYPE duidvj[][NDIM])
{
  int j,k;
  FTYPE alpha,betacon[NDIM];
  
  alpha = 1.0 / sqrt(-glob_gcon[TT][TT]);
  SLOOPA betacon[j] = glob_gcon[TT][j]*alpha*alpha;

  SLOOPA duidvj[TT][j] = dgdvi[j]/alpha;
 
  SLOOP duidvj[j][k] = delta(j,k) - betacon[j]/alpha*dgdvi[k];

  // duidvj[j][TT] doesn't exist since no v0 primitive variable, 
  //  so should never be referenced

}

static void dutdui_calc_g(FTYPE *ucov, FTYPE *dutdui) 
{
	int j ;

	SLOOPA dutdui[j] = -ucov[j]/ucov[0] ;

	return ;
}

static void dbtdui_calc_g(FTYPE *dutdui, FTYPE *pr, FTYPE *dbtdui) 
{
	int j ;
	static FTYPE B[NDIM],Bcov[NDIM] ;

	B[0] = 0. ;
	SLOOPA B[j] = pr[j+BCON1-1] ;

	lower_g(B,glob_gcov,Bcov) ;
	SLOOPA dbtdui[j] = Bcov[j] + Bcov[0]*dutdui[j] ;

	return ;
}

// valid for all i, j=1,2,3
static void dbiduj_calc_g(FTYPE *dbtdui,FTYPE *dutdui,FTYPE *ucon, 
			  FTYPE *b, FTYPE dbiduj[][NDIM]) 
{
	int j,k ;

	DLOOP dbiduj[j][k] = 0. ;

	SLOOP dbiduj[j][k] = -b[j]*dutdui[k]/ucon[TT] 
			+ ucon[j]*dbtdui[k]/ucon[TT] ;

	SLOOPA dbiduj[j][j] += b[TT]/ucon[TT] ;

	SLOOPA dbiduj[TT][j] = dbtdui[j] ;

	return ;
}

static void dbsqdui_calc_g(FTYPE dbiduj[][NDIM],FTYPE *bcov, FTYPE *dbsqdui) 
{
	int j,k ;

	DLOOPA dbsqdui[j] = 0. ;
	DLOOP dbsqdui[j] += 2.*bcov[k]*dbiduj[k][j] ;

	return ;
}

// valid for all i, j=1,2,3
static void dudduu_calc_g(FTYPE*dutdui, FTYPE dudduu[][NDIM]) 
{
  int j,k;

  DSLOOP dudduu[j][k] = glob_gcov[j][k] + glob_gcov[j][TT]*dutdui[k];

}


// valid for all i, j=1,2,3
static void dbdiduj_calc_g(FTYPE dbiduj[][NDIM],FTYPE dbdiduj[][NDIM])
{
  int j,k;
  int l;

  DSLOOP dbdiduj[j][k] = 0.0;

  // just a lowered dbiduj (lower the 1st index)
  DSLOOP for(l=0;l<NDIM;l++) dbdiduj[j][k] += glob_gcov[j][l]*dbiduj[l][k];

}

/*************************************************************************/
/*************************************************************************

   LU_decompose():

       Performs a LU decomposition of the matrix A using Crout's method      
       with partial implicit pivoting.  The exact LU decomposition of the    
       matrix can be reconstructed from the resultant row-permuted form via  
       the integer array permute[]                                            
                                                                             
       The algorithm closely follows ludcmp.c of "Numerical Recipes  
       in C" by Press et al. 1992.                                           
                                                                             
       This will be used to solve the linear system  A.x = B                 
                                                                             
       Returns (1) if a singular matrix is found,  (0) otherwise.            

*************************************************************************/


static int LU_decompose( FTYPE A[][NEWT_DIM], int permute[] )
{

  const  FTYPE absmin = 1e-30; /* Value used instead of 0 for singular matrices */

  static FTYPE row_norm[NEWT_DIM];
  FTYPE  absmax, maxtemp, mintemp;

  int i, j, k, max_row;
  int n = NEWT_DIM;


  max_row = 0;

  /* Find the maximum elements per row so that we can pretend later
     we have unit-normalized each equation: */

  for( i = 0; i < n; i++ ) { 
    absmax = 0.;
    
    for( j = 0; j < n ; j++ ) { 
      
      maxtemp = fabs( A[i][j] ); 

      if( maxtemp > absmax ) { 
	absmax = maxtemp; 
      }
    }

    /* Make sure that there is at least one non-zero element in this row: */
    if( absmax == 0. ) { 
     fprintf(stderr, "LU_decompose(): row-wise singular matrix!\n");
      return(1);
    }

    row_norm[i] = 1. / absmax ;   /* Set the row's normalization factor. */
  }


  /* The following the calculates the matrix composed of the sum 
     of the lower (L) tridagonal matrix and the upper (U) tridagonal
     matrix that, when multiplied, form the original maxtrix.  
     This is what we call the LU decomposition of the maxtrix. 
     It does this by a recursive procedure, starting from the 
     upper-left, proceding down the column, and then to the next
     column to the right.  The decomposition can be done in place 
     since element {i,j} require only those elements with {<=i,<=j} 
     which have already been computed.
     See pg. 43-46 of "Num. Rec." for a more thorough description. 
  */

  /* For each of the columns, starting from the left ... */
  for( j = 0; j < n; j++ ) {

    /* For each of the rows starting from the top.... */

    /* Calculate the Upper part of the matrix:  i < j :   */
    for( i = 0; i < j; i++ ) {
      for( k = 0; k < i; k++ ) { 
	A[i][j] -= A[i][k] * A[k][j];
      }
    }

    absmax = 0.0;

    /* Calculate the Lower part of the matrix:  i <= j :   */

    for( i = j; i < n; i++ ) {

      for (k = 0; k < j; k++) { 
	A[i][j] -= A[i][k] * A[k][j];
      }

      /* Find the maximum element in the column given the implicit 
	 unit-normalization (represented by row_norm[i]) of each row: 
      */
      maxtemp = fabs(A[i][j]) * row_norm[i] ;

      if( maxtemp >= absmax ) {
	absmax = maxtemp;
	max_row = i;
      }

    }

    /* Swap the row with the largest element (of column j) with row_j.  absmax
       This is the partial pivoting procedure that ensures we don't divide
       by 0 (or a small number) when we solve the linear system.  
       Also, since the procedure starts from left-right/top-bottom, 
       the pivot values are chosen from a pool involving all the elements 
       of column_j  in rows beneath row_j.  This ensures that 
       a row  is not permuted twice, which would mess things up. 
    */
    if( max_row != j ) {

      /* Don't swap if it will send a 0 to the last diagonal position. 
	 Note that the last column cannot pivot with any other row, 
	 so this is the last chance to ensure that the last two 
	 columns have non-zero diagonal elements.
       */

      if( (j == (n-2)) && (A[j][j+1] == 0.) ) {
	max_row = j;
      }
      else { 
	for( k = 0; k < n; k++ ) { 

	  maxtemp       = A[   j   ][k] ; 
	  A[   j   ][k] = A[max_row][k] ;
	  A[max_row][k] = maxtemp; 

	}

	/* Don't forget to swap the normalization factors, too... 
	   but we don't need the jth element any longer since we 
	   only look at rows beneath j from here on out. 
	*/
	row_norm[max_row] = row_norm[j] ; 
      }
    }

    /* Set the permutation record s.t. the j^th element equals the 
       index of the row swapped with the j^th row.  Note that since 
       this is being done in successive columns, the permutation
       vector records the successive permutations and therefore
       index of permute[] also indexes the chronology of the 
       permutations.  E.g. permute[2] = {2,1} is an identity 
       permutation, which cannot happen here though. 
    */

    permute[j] = max_row;

    if( A[j][j] == 0. ) { 
      A[j][j] = absmin;
    }


  /* Normalize the columns of the Lower tridiagonal part by their respective 
     diagonal element.  This is not done in the Upper part because the 
     Lower part's diagonal elements were set to 1, which can be done w/o 
     any loss of generality.
  */
    if( j != (n-1) ) { 
      maxtemp = 1. / A[j][j]  ;
      
      for( i = (j+1) ; i < n; i++ ) {
	A[i][j] *= maxtemp;
      }
    }

  }

  return(0);

  /* End of LU_decompose() */

}


/************************************************************************/
/************************************************************************

   LU_substitution():

       Performs the forward (w/ the Lower) and backward (w/ the Upper)   
       substitutions using the LU-decomposed matrix A[][] of the original
       matrix A' of the linear equation:  A'.x = B.  Upon entry, A[][]   
       is the LU matrix, B[] is the source vector, and permute[] is the  
       array containing order of permutations taken to the rows of the LU
       matrix.  See LU_decompose() for further details. 		     
    								     
       Upon exit, B[] contains the solution x[], A[][] is left unchanged.
								     
************************************************************************/

static void LU_substitution( FTYPE A[][NEWT_DIM], FTYPE B[], int permute[] )
{
  int i, j ;
  int n = NEWT_DIM;
  FTYPE tmpvar,tmpvar2;

  
  /* Perform the forward substitution using the LU matrix. 
   */
  for(i = 0; i < n; i++) {

    /* Before doing the substitution, we must first permute the 
       B vector to match the permutation of the LU matrix. 
       Since only the rows above the currrent one matter for 
       this row, we can permute one at a time. 
    */
    tmpvar        = B[permute[i]];
    B[permute[i]] = B[    i     ];
    for( j = (i-1); j >= 0 ; j-- ) { 
      tmpvar -=  A[i][j] * B[j];
    }
    B[i] = tmpvar; 
  }
	   

  /* Perform the backward substitution using the LU matrix. 
   */
  for( i = (n-1); i >= 0; i-- ) { 
    for( j = (i+1); j < n ; j++ ) { 
      B[i] -=  A[i][j] * B[j];
    }
    B[i] /= A[i][i] ; 
  }

  /* End of LU_substitution() */

}


/**********************************************************************/ 
/************************************************************

  general_newton_raphson(): 

    -- performs Newton-Rapshon method on an arbitrary system.

    -- inspired in part by Num. Rec.'s routine newt();

    Arguements: 

       -- x[]   = set of independent variables to solve for;
       -- n     = number of independent variables and residuals;
       -- funcd = name of function that calculates residuals, etc.;

*****************************************************************/
static int general_newton_raphson( FTYPE x[], int n, 
				   void (*funcd) (FTYPE [], FTYPE [], FTYPE [], 
						  FTYPE [][NEWT_DIM], FTYPE *, 
						  FTYPE *, int))
{
  FTYPE f, df, dx[NEWT_DIM], x_old[NEWT_DIM], 
    resid[NEWT_DIM], jac[NEWT_DIM][NEWT_DIM];
  FTYPE errx, x_orig[NEWT_DIM];
  int    n_iter, id, jd, i_extra, doing_extra;
  FTYPE dW,dvsq,vsq_old,vsq,W,W_old;

  int   keep_iterating;


  // Initialize various parameters and variables:
  errx = 1. ; 
  df = f = 1.;
  i_extra = doing_extra = 0;
  for( id = 0; id < n ; id++)  x_old[id] = x_orig[id] = x[id] ;

  vsq_old = vsq = W = W_old = 0.;

  n_iter = 0;


  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 

    (*funcd) (x, dx, resid, jac, &f, &df, n);  /* returns with new dx, f, df */


    /* Save old values before calculating the new: */
    errx = 0.;
    for( id = 0; id < n ; id++) {
      x_old[id] = x[id] ;
    }

    /* Make the newton step: */
    for( id = 0; id < n ; id++) {
      x[id] += dx[id]  ;
    }


    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/
    /* For the new criterion, always look at error in "W" : */
    // METHOD specific:
    W_old = W;
    vsq_old = vsq;
    calc_errx_5d(&W, &vsq, x);
    errx  = (W==0.) ?  fabs(W-W_old) : fabs((W-W_old)/W);


    /*****************************************************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*  before stopping                                                          */
    /*****************************************************************************/
    
    if( (fabs(errx) <= NEWT_TOL) && (f <= TOL_RESID) 
	&& (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( ((fabs(errx) <= NEWT_TOL)&&(doing_extra == 0)) 
	|| (i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }   // END of while(keep_iterating)


  /*  Check for bad untrapped divergences : */
  if( (finite(f)==0) || (finite(df)==0) ) {
    return(2);
  }


  if( (fabs(errx) > MIN_NEWT_TOL) || (f > TOL_RESID)  ){
    return(1);
  } 
  if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ){
    return(0);
  }
  if( fabs(errx) <= NEWT_TOL ){
    return(0);
  }

  return(0);

}

