
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
#ifdef MPI_USED
#include "mpi.h"
#endif

#include <unistd.h>
#include <signal.h>
#include <time.h>
#include <fcntl.h>
#include <sys/time.h>

#include "decs.h"
#include "u2p_util.h"
#include "cache4d.h"
#include "random.h"
#include "mnewt.h"

#include "cooleos.h"


#define max(a,b)  (((a) > (b)) ? (a) : (b))
#define min(a,b)  (((a) < (b)) ? (a) : (b))

/* your choice of floating-point data type */
#define FTYPE double

#define USE_ROLLING_CACHE 0

#define NCACHEVAL 5		//e, p, Qnu, tau, Yee
#define CACHE_E   0
#define CACHE_P   1
#define CACHE_Qnu 2
#define CACHE_Tau 3
#define CACHE_Yee 4
static int dimmask_ALL[NCACHEVAL] = { 1, 1, 1, 1, 1 };
static int dimmask_E[NCACHEVAL] = { 1, 0, 0, 0, 0 };
static int dimmask_P[NCACHEVAL] = { 0, 1, 0, 0, 0 };


#define MAXLOS 1
#define DEBUGLOG 0
#define CACHEDEBUG 0

#ifdef MPI_USED
#define CACHETHREADS 4		//Number of cache threads
#else
#define CACHETHREADS 6		//Number of cache threads
#endif

void coord (int h, int i, int j, int loc, double *X);
void bl_coord (double *X, double *r, double *th, double *ph);



extern "C" void abundance_ (double *rho, double *T, double *HdFac, double *e,
                            double *p, double *Qnu, double *tau, double *xnee, double *xnep);

static int calc_cache_val (double rho, double T, double HdFac, double *val);


static Cache4d *cache = NULL;
static int cacheIx = -999;
static int ZIdx = 0;
double *LastT = NULL;
//static double _aLastT[N3+2][N1+4][N2+4];      /* Last T for each of cell*/
static double *_aLastT;		/* Last T for each of cell */
static Cache4d **cacheT;

static const double hh = 1.0E-6;

double rhoMin, rhoMax;
double TMin, TMax;

double init_rhoMax;


RANDOMDEF rd = { 12345, 65435, 34221, 43543 };


const double GNEWT = 6.6742867e-8;
const double MSUN = 1.9892e33;
const double YEAR = 365.25 * 24. * 3600.;
const double KBOL = 1.3e-16;
const double CL = 2.99792458e10;
const double MP = 1.67262163783e-24;
const double ME = 9.10956e-28;
const double MN = 1.67492729e-24;
const double MHE = 6.6465e-24;	//2p+2e+binding energy
const double EE = 4.80320680e-10;
const double SIGMATH = 6.65248e-25;
const double SIGMASBOL = 5.67051e-5;
const double MAV = (MP + MN + 2. * ME + MHE) / 5.;	// avrage particle mass

/* set mass of the black hole- sets the lenght scale, and the initial mass of the disk*/
const double M_BLH = 3. * MSUN;

//const double M_BLH = 62. * MSUN;
//      const double M_unit=1.5e-5*MSUN;
//const double M_BLH=10.*MSUN;
const double bhs = 0.6;
//const double bhs = 0.98;	//black hole spin, has to match the init.c val.



/* fundamental units */
const double L_UNIT = GNEWT * M_BLH / (CL * CL);
//const double M_UNIT = 5.e-5*MSUN;
//const double M_UNIT = 1.e-5 * MSUN;
//const double M_UNIT = 1.e-5 * MSUN;
//const double M_UNIT = 5.e-4 * MSUN;

//const double RHO_SCALE = 0.16;
const double RHO_SCALE = 1.5e-5;

const double M_UNIT = RHO_SCALE * MSUN;
const double T_UNIT = L_UNIT / CL;

/* characteristic fraction of radius over
   which absorption occurs */
//      const double L = 0.1 * r * L_UNIT;

//const double USCALE = 2.891985047341179325e+02;  //u is devided by rhomax 
//const double USCALE = 1.0e4;

//const double USCALE = 2.659700688070531192e-03;
const double USCALE = 1.0;


/* scaling ratio */
const double RHO_UNIT = M_UNIT / (L_UNIT * L_UNIT * L_UNIT);
const double U_UNIT = RHO_UNIT * CL * CL * USCALE;
const double P_UNIT = U_UNIT;



typedef struct
{
    double rho;
    double e;
} EVOLPAR_RHOE;


inline bool
checkisnan (double var)
{
    volatile double d = var;
    return d != d;
}

#if CACHEDEBUG
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
#endif

/************************************************************************/
int
calc_cache_val (double ro, double te, double HdFac, double *val)
/************************************************************************/
{
#if CACHEDEBUG
	extern unsigned long long getTicks();
	unsigned long long lastClock = getTicks();
#endif

    double e, p, Qnu, tau, xnee, xnep;

    abundance_ (&ro, &te, &HdFac, &e, &p, &Qnu, &tau, &xnee, &xnep);

    val[CACHE_E] = e;
    val[CACHE_P] = p;
    val[CACHE_Qnu] = Qnu;
    val[CACHE_Tau] = tau;
    val[CACHE_Yee] = (xnee-xnep)/(ro/1.66e15);

#if CACHEDEBUG
	printf("\ncalc_cache_val ro=%28.18e, T=%28.18e, e=%28.18e, p=%28.18e, Yee=%28.18e, time=%3lf\n", ro, te, e, p, val[CACHE_Yee], ((double)(getTicks()-lastClock))/1000);
#endif
    return 1;
}

/************************************************************************/
inline void
mnewt_mini_rho_e (double x[], void *usrfunpar, double y[])
/************************************************************************/
{
    double val[NCACHEVAL];

    EVOLPAR_RHOE *par = (EVOLPAR_RHOE *) usrfunpar;


    double rho = par->rho;
    double e = par->e;


    double currT = x[1];

    int vv = cache->getvalues (rho, currT, 0.0, SPLINE, val, NULL, ZIdx, dimmask_E);

    if (!vv)
    {
        printf ("Poza zakresem mnewt_mini_rho_e: rho=%E T=%E\n", rho, currT);
    }

    double ecur = val[CACHE_E];

//      printf("currT=%E ecur=%E e=%E\n", currT, ecur, e);

    double F = (e - ecur) / e;

    y[1] = F;
}


/************************************************************************/
inline void
mnewt_mini_rho_e_der (double x[], void *usrfunpar, double y[], double d[])
/************************************************************************/
{
    double val[NCACHEVAL];
    double der[NCACHEVAL];

    EVOLPAR_RHOE *par = (EVOLPAR_RHOE *) usrfunpar;


    double rho = par->rho;
    double e = par->e;


    double currT = x[1];

    int vv = cache->getvalues (rho, currT, 0.0, SPLINE, val, der, ZIdx, dimmask_E);

    if (!vv)
    {
        printf ("Poza zakresem mnewt_mini_rho_e: rho=%E T=%E\n", rho, currT);
    }

    double ecur = val[CACHE_E];
    double dercur = der[CACHE_E];

//      printf("currT=%E ecur=%E e=%E\n", currT, ecur, e);

    double F = (e - ecur) / e;

    y[1] = F;
    d[1] = -dercur / e;
}

/************************************************************************/
void
usrfun_mnewt_rho_e (double x[2], double lalpha[2][2], double lbetha[2],
                    void *usrfunpar)
/*Function called by mnewt in new version with eos_neu_trapping*/
/*Vertically averaged equations to calculate the steady state disk*/
/************************************************************************/
{
#if 0
    const double epsx1 = x[1] * hh;

    double xn[2];
    double F1P[2];
    double F1M[2];
    double F[2];

    mnewt_mini_rho_e (x, usrfunpar, F);

    xn[1] = x[1] + epsx1;
    mnewt_mini_rho_e (xn, usrfunpar, F1P);

    xn[1] = x[1] - epsx1;
    mnewt_mini_rho_e (xn, usrfunpar, F1M);
    /*
       {
       static int ntr;
       printf("usrfun_mnewt_dysk ntr=%6d [%E] -> %E\n", ++ntr, x[1], F[1]);
       }
     */

    lalpha[1][1] = (F1P[1] - F1M[1]) / (2.0 * epsx1);

    lbetha[1] = -F[1];
#else

    double F[2];
    double D[2];

    mnewt_mini_rho_e_der (x, usrfunpar, F, D);

    lalpha[1][1] = D[1];

    lbetha[1] = -F[1];
#endif

}



/************************************************************************/
double
pressure_rho0_u_CoolEOS (double rho0, double u)
/************************************************************************/
{
    double x[2];
    double minx[2];
    double maxx[2];
    double berr;
    double val[NCACHEVAL];


    double rho = rho0 * RHO_UNIT;
    double e = u * U_UNIT;
    /*
      if (checkisnan (rho0))
        {
          printf ("rho0=NAN pressure_rho0_e_CoolEOS: rho=%E e=%E\n", rho0, e);
          exit (0);
        }
    */

    if (cache == NULL)
    {
        printf ("Invalid pressure_rho0_u_CoolEOS call!\n");
        exit (0);
    }

    if (*LastT == 0.0)
    {
//        const double corrfact = 1.5e-1;
        const double corrfact = 1.0;
        *LastT = (gam - 1.) * u / rho0 / KBOL * MAV * CL * CL * corrfact;	//T [K] //Initial guess
		if(*LastT < TMin * (1.0 + hh))
			*LastT = TMin * (1.0 + hh);
    }

    if (rho < rhoMin)
    {
#if DEBUGLOG
        printf
        ("rho=%E below limit %E in pressure_rho0_u_CoolEOS. Using Gamma Low!\n",
         rho, rhoMin);
#endif
        return ((GAMMA - 1.) * u);
    }


    if (rho > rhoMax)
    {
//        printf ("Invalid rho=%E in pressure_rho0_u_CoolEOS call!\n", rho);
//        exit (0);
#if DEBUGLOG
        printf
        ("rho=%E above limit %E in pressure_rho0_u_CoolEOS. Using Gamma Low!\n",
         rho, rhoMax);
#endif
        return ((GAMMA - 1.) * u);
    }

    if (e < 0.0)
    {
//        printf("e is below ZERO!\n");
        return 0.0;
    }
    
    minx[1] = TMin * (1.0 + hh);
    maxx[1] = TMax * (1.0 - hh);


    EVOLPAR_RHOE usrfunpar;
    usrfunpar.rho = rho;
    usrfunpar.e = e;

//                      printf("rho=%E startT=%E e=%E\n", rho, LastT, e);


    x[1] = *LastT;

    int nlos = MAXLOS;
    while (nlos--)
    {
        mnewtc < 1 > (300, x, 1.0E-6, minx, maxx, usrfun_mnewt_rho_e,
                      &usrfunpar, &berr, 1, &rd);

        if (berr < 1.0E-6 || nlos == 0)
            break;

        x[1] =
            UniformRandomRangeLog (&rd, max (minx[1], *LastT * (1.0 - 1.0e-3)),
                                   min (maxx[1], *LastT * (1.0 + 1.0e-3)));
    }

    //printf("0. Zbieglo dla rho=%E T=%E berr=%E\n", rho, x[1], berr);

    if (berr > 1.0E-6)
    {
//        const double corrfact = 1.5e-1;
        const double corrfact = 1.0;
        double estT = (gam - 1.) * u / rho0 / KBOL * MAV * CL * CL * corrfact;	//T [K] //Initial guess


#if DEBUGLOG
        printf
        ("mnewtc pressure_rho0_u failed berr=%E rho=%E e=%E, estT=%E, Using min or Gamma Low!\n",
         berr, rho, e, estT);
#endif

//        if(estT < TMax)
//			x[1] = TMin;
//		else
	        return ((GAMMA - 1.) * u);
/*			
//	     return ((GAMMA - 1.) * u);
        return 0.0;
*/
        //getc(stdin);
    }

    int vv = cache->getvalues (rho, x[1], 0.0, SPLINE, val, NULL, ZIdx, dimmask_P);

    if (!vv)
    {
        printf ("Poza zakresem pressure_rho0_u_CoolEOS: rho=%E T=%E\n", rho, x[1]);
    }


    *LastT = x[1];

    return val[CACHE_P] / P_UNIT;

}

/************************************************************************/
void 
test_CoolEOS(double rho0)
/************************************************************************/
{
    double rho = rho0 * RHO_UNIT;
    double val[NCACHEVAL];
	int vv;

	printf("test_CoolEOS rho0=%E, rho=%E\n", rho0, rho);

//    vv = cache->getvalues (rho, TMin, 0.0, SPLINE, val, ZIdx, dimmask_ALL);
	vv = calc_cache_val (rho, TMin, cache->cdesc.zmin, val);
	printf("rho=%E, T=%E, E=%E, P=%E, ugl=%E\n", rho, TMin, val[CACHE_E], val[CACHE_P], TMin/(gam - 1.) *  rho0 * KBOL / MAV / CL / CL * U_UNIT);


        double Tdelta = (log10(TMax) - log10(TMin))/(cache->cdesc.y - 1);
        double T = pow(10.0, log10(TMin) + Tdelta*cache->cdesc.y/2);

        vv = calc_cache_val (rho, T, cache->cdesc.zmin, val);
        printf("rho=%E, T=%E, E=%E, P=%E, ugl=%E\n", rho, T, val[CACHE_E], val[CACHE_P], T/(gam - 1.) *  rho0 * KBOL / MAV / CL / CL * U_UNIT);



//    vv = cache->getvalues (rho, TMax, 0.0, SPLINE, val, ZIdx, dimmask_ALL);
	vv = calc_cache_val (rho, TMax, cache->cdesc.zmin, val);
	printf("rho=%E, T=%E, E=%E, P=%E, ugl=%E\n", rho, TMax, val[CACHE_E], val[CACHE_P],  TMax/(gam - 1.) *  rho0 * KBOL / MAV / CL / CL * U_UNIT);


}


/************************************************************************/
void
CoolVal_rho0_u_CoolEOS (double rho0, double u, int h, int i, int j, double coolvals[NCOOLVAL])
/************************************************************************/
{
    double x[2];
    double minx[2];
    double maxx[2];
    double berr;
    double val[NCACHEVAL];


    double rho = rho0 * RHO_UNIT;
    double e = u * U_UNIT;

    if (cache == NULL)
    {
        printf ("Invalid CoolVal_rho0_u_CoolEOS call!\n");
        exit (0);
    }

    if (rho < rhoMin)
    {
#if DEBUGLOG
        printf
        ("rho=%E below limit %E in CoolVal_rho0_u_CoolEOS. Using Gamma Low!\n",
         rho, rhoMin);
#endif
        coolvals[COOL_P] = (GAMMA - 1.) * u;	//p
        coolvals[COOL_LAMBDA_SIM] = SMALL;
        coolvals[COOL_Qnu] = SMALL;

        double HdFac = cache->VectorZ[ZIdx];
		double Hd = HdFac * CL * sqrt (coolvals[COOL_P] * P_UNIT / (rho * CL * CL + e));
		coolvals[COOL_Hd] = Hd / L_UNIT;
		coolvals[COOL_T] = (gam - 1.) * u / rho0 / KBOL * MAV * CL * CL;
		coolvals[COOL_Tau] = 0.0;
		coolvals[COOL_Yee] = 0.0;
        return;
    }


    if (rho > rhoMax)
    {
//        printf ("Invalid rho=%E in CoolVal_rho0_u_CoolEOS call!\n", rho);
//        exit (0);
#if DEBUGLOG
        printf
        ("rho=%E above limit %E in CoolVal_rho0_u_CoolEOS. Using Gamma Low!\n",
         rho, rhoMin);
#endif
        coolvals[COOL_P] = (GAMMA - 1.) * u;	//p
        coolvals[COOL_LAMBDA_SIM] = SMALL;
        coolvals[COOL_Qnu] = SMALL;

        double HdFac = cache->VectorZ[ZIdx];
		double Hd = HdFac * CL * sqrt (coolvals[COOL_P] * P_UNIT / (rho * CL * CL + e));
		coolvals[COOL_Hd] = Hd / L_UNIT;
		coolvals[COOL_T] = (gam - 1.) * u / rho0 / KBOL * MAV * CL * CL;
		coolvals[COOL_Tau] = 0.0;
		coolvals[COOL_Yee] = 0.0;
        return;
    }


    if (*LastT == 0.0)
    {
//        const double corrfact = 1.5e-1;
        const double corrfact = 1.0;
        *LastT = (gam - 1.) * u / rho0 / KBOL * MAV * CL * CL * corrfact;	//T [K] //Initial guess
		if(*LastT < TMin * (1.0 + hh))
			*LastT = TMin * (1.0 + hh);
    }

    minx[1] = TMin * (1.0 + hh);
    maxx[1] = TMax * (1.0 - hh);


    EVOLPAR_RHOE usrfunpar;
    usrfunpar.rho = rho;
    usrfunpar.e = e;

//                      printf("rho=%E startT=%E e=%E\n", rho, LastT, e);


    x[1] = *LastT;

    int nlos = MAXLOS;
    while (nlos--)
    {
        mnewtc < 1 > (300, x, 1.0E-6, minx, maxx, usrfun_mnewt_rho_e,
                      &usrfunpar, &berr, 1, &rd);

        if (berr < 1.0E-6 || nlos == 0)
            break;

        x[1] =
            UniformRandomRangeLog (&rd, max (minx[1], *LastT * (1.0 - 1.0e-3)),
                                   min (maxx[1], *LastT * (1.0 + 1.0e-3)));
    }

    //printf("0. Zbieglo dla rho=%E T=%E berr=%E\n", rho, x[1], berr);

    if (berr > 1.0E-6)
    {
//      printf ("mnewtc CoolVal_rho0_u failed berr=%E Using Gamma Low!\n",
//	      berr);
//        const double corrfact = 1.5e-1;
        const double corrfact = 1.0;
        double estT = (gam - 1.) * u / rho0 / KBOL * MAV * CL * CL * corrfact;	//T [K] //Initial guess

#if DEBUGLOG
        printf
        ("mnewtc pressure_rho0_u failed berr=%E rho=%E e=%E, estT=%E, Using min or Gamma Low!\n",
         berr, rho, e, estT);
#endif
//        if(estT < TMax)
//			x[1] = TMin;
//		else
		{

    	    coolvals[COOL_P] = (GAMMA - 1.) * u;	//p
	        coolvals[COOL_LAMBDA_SIM] = SMALL;
	        coolvals[COOL_Qnu] = SMALL;

	        double HdFac = cache->VectorZ[ZIdx];
			double Hd = HdFac * CL * sqrt (coolvals[COOL_P] * P_UNIT / (rho * CL * CL + e));
			coolvals[COOL_Hd] = Hd / L_UNIT;
			coolvals[COOL_T] = (gam - 1.) * u / rho0 / KBOL * MAV * CL * CL;
			coolvals[COOL_Tau] = 0.0;
			coolvals[COOL_Yee] = 0.0;
        	return;
		}
        //getc(stdin);
    }

    int vv = cache->getvalues (rho, x[1], 0.0, SPLINE, val, NULL, ZIdx, dimmask_ALL);

    if (!vv)
    {
        printf ("Poza zakresem CoolVal_rho0_u_CoolEOS: rho=%E T=%E\n", rho, x[1]);
    }

    *LastT = x[1];

    //Qnu
    {
        double HdFac = cache->VectorZ[ZIdx];
      
		double Hd = HdFac * CL * sqrt (val[CACHE_P] / (rho * CL * CL + val[CACHE_E]));

		double Hdisk;
		{
	    	double X[NDIM];
    	    double r, th, ph;
            double H;

	        coord (h, i, j, CENT, X);
    	    bl_coord (X, &r, &th, &ph);
		
			Hdisk = 0.5*r*L_UNIT; //disk thickness - assumption
		}

        /*convert lambda to emissivity in erg/cm^2/s */
        double lambda_cgs = val[CACHE_Qnu] / Hdisk;

        /*convert cooling rate from cgs units to the code units */
		coolvals[COOL_LAMBDA_SIM] = lambda_cgs / (RHO_UNIT * CL * CL * CL / L_UNIT);

		coolvals[COOL_Hd] = Hd / L_UNIT;

		coolvals[COOL_Qnu] = val[CACHE_Qnu] / (RHO_UNIT * CL * CL * CL / L_UNIT);
    }

    coolvals[COOL_P] = val[CACHE_P] / P_UNIT;	//p
    coolvals[COOL_T] = x[1];	//T
    
	if(val[CACHE_Tau] > 0) //Tau   
	    coolvals[COOL_Tau] = val[CACHE_Tau];
	else
	    coolvals[COOL_Tau] = 0.0;

	//Yee
     coolvals[COOL_Yee] = val[CACHE_Yee];
}



/************************************************************************/
double
dpdu_calc_CoolEOS (double rho0, double u)
/************************************************************************/
{
    double h = u * hh/U_UNIT;
    double uP = u + h;
    double uM = u - h;

    double ret;

    if (uM > 0)
    {
        double pP = pressure_rho0_u_CoolEOS (rho0, uP);

        double pM = pressure_rho0_u_CoolEOS (rho0, uM);

        ret = (pP - pM) / (2.0*h);
    }
    else
    {
        uM = u;
        double pP = pressure_rho0_u_CoolEOS (rho0, uP);

        double pM = pressure_rho0_u_CoolEOS (rho0, uM);

        ret = (pP - pM) / (h);
    }

    /*	if(checkisnan(ret))
    	{
            printf("ret=NAN dpdvsq_calc_CoolEOS 3: rho0=%E w=%E\n", rho0, w);
    		exit(0);
    	}*/

    return ret;

}

/************************************************************************/
double
dpdrho_calc_CoolEOS (double rho0, double u)
/************************************************************************/
{
    double h = rho0 * hh/RHO_UNIT;
    double rhoP = rho0 + h;
    double rhoM = rho0 - h;
    double ret;

    if (rhoM > 0)
    {
        double pP = pressure_rho0_u_CoolEOS (rhoP, u);

        double pM = pressure_rho0_u_CoolEOS (rhoM, u);

        ret = (pP - pM) / (2.0 * h);
    }
    else
    {
        rhoM = rho0;
        double pP = pressure_rho0_u_CoolEOS (rhoP, u);

        double pM = pressure_rho0_u_CoolEOS (rhoM, u);

        ret = (pP - pM) / (h);
    }

    /*	if(checkisnan(ret))
    	{
            printf("ret=NAN dpdvsq_calc_CoolEOS 3: rho0=%E w=%E\n", rho0, w);
    		exit(0);
    	}*/

    return ret;

}

/**************************************************************/
void
ProgramStop (int code)
/**************************************************************/
{
    printf("STOP BEGIN:\n");

	if(_aLastT)
        free(_aLastT);

	if(cacheT)
	{
	    for (int i = 0; i < SN1; i++)
			delete cacheT[i];
    	free(cacheT);
	}

#ifdef MPI_USED
    MPI_Finalize ();
#endif
    printf("STOP END\n");

    exit (0);
}

/************************************************************************/
int
InitCache_CoolEOS (int hh, int ii, int jj)
/************************************************************************/
{
    double (*aLastT)[SN1][SN2] = (double (*)[SN1][SN2]) (_aLastT);

    if (ii != cacheIx)
    {
//          printf("Calculations start for i=%d j=%d\n", ii, jj);
    }

    LastT = &aLastT[hh][ii][jj];

//          printf("Calculations start for i=%d j=%d, LastT=%E\n", ii, jj, *LastT);

    if (*LastT == 0.0)
    {	
	 if(jj > jstart
            && (*LastT = aLastT[hh][ii][jj - 1]) != 0.0)
    	{
    	}
/*    	else if (*LastT == 0.0 && ii > istart
             && (*LastT = aLastT[hh][ii - 1][jj]) != 0.0)
    	{
    	}
    	else if (*LastT == 0.0 && hh > hstart
             && (*LastT = aLastT[hh - 1][ii][jj]) != 0.0)
    	{
	}
*/
    }


    cache = cacheT[ii];
    cacheIx = ii;

    ZIdx = 0;
    return 0;
}

/************************************************************************/
int
SetupThreads_CoolEOS ()
/************************************************************************/
{
    if (CacheWorker::getInstance ().getCount () < CACHETHREADS)
    {
        CacheWorker::getInstance ().startThreads (CACHETHREADS, calc_cache_val,
                NCACHEVAL);
		//Renice main thread to get more cpu for eos threads
		nice(10);
    }
    return 0;
}
/************************************************************************/
int
SetupCache_CoolEOS ()
/************************************************************************/
{

    if(rank == 0)
    {
        printf("SetupCache_CoolEOS, L_UNIT=%E, RHO_UNIT=%E, U_UNIT=%E, SIM_UNIT=%E\n", L_UNIT, RHO_UNIT, U_UNIT, RHO_UNIT * CL * CL * CL / L_UNIT);
    }



    _aLastT = (double *) calloc (SN3 * SN1 * SN2, sizeof (double));

    cacheT = (Cache4d **) calloc (SN1, sizeof (Cache4d *));


    if (CacheWorker::getInstance ().getCount () < CACHETHREADS)
    {
        CacheWorker::getInstance ().startThreads (CACHETHREADS, calc_cache_val,
                NCACHEVAL);
    }


#if (USE_ROLLING_CACHE)		//Rolling cache
    for (int i = 0; i < SN1; i++)
    {
        CacheDesc cdesc;
        char fname[256];
        sprintf (fname, "cache/cache_%04d.dat", ix);

        rhoMin = 1.e2;
        rhoMax = 1.e13;

        cdesc.mode = CONSTGRID;
        cdesc.x = 256;		//ro
        cdesc.xmin = rhoMin;
        cdesc.xmax = rhoMax;
        cdesc.xrat = 8.0E-3;
        cdesc.xtype = ROLLING;

        TMin = 1.e2;
        TMax = 1.e14;

        cdesc.y = 256;		//T
        cdesc.ymin = TMin;
        cdesc.ymax = TMax;
        cdesc.yrat = 8.0E-3;
        cdesc.ytype = ROLLING;


        double HdFac;

        //Please note, Hd inside EOS is calculted as Hd=HdFac*sqrt(p/rho)
        {
            double X[NDIM];
            double r, th, ph;
            double cs, omega;

            coord (hstart, i, jstart, CENT, X);
            bl_coord (X, &r, &th, &ph);
            //  cs = sqrt((gam-1.0)*gam*u_eq/ro_eq)*CL;
            //  omega = CL*CL*CL/GNEWT/M_BLH/(a+pow(r, 1.5));

            //printf("r=%E\n", r);

            cs = 1.0;
            omega = CL * CL * CL / GNEWT / M_BLH / (bhs + pow (r, 1.5));

            HdFac = cs / omega;
        }


        cdesc.z = 1;		//HdFac
        cdesc.zmin = cdesc.zmax = HdFac;
        cdesc.ztype = CUSTOM;

        cdesc.d = NCACHEVAL;	//tablice: e, p
        cdesc.autosave = 100;  //Cache save after number of points

        cacheT[i] = new Cache4d (&cdesc, calc_cache_val, fname);
    }
#else

    for (int i = 0; i < SN1; i++)
    {
        char apdx[32] = "";
        int idx = i - (irank == 0 ? 2 : 1) + N1 / isize * irank;
        if (i < istart || i > istop)
        {
            sprintf (apdx, "_%d_%d", isize, irank);
        }

        //Please note, Hd inside EOS is calculted as Hd=HdFac*sqrt(p/rho)

        CacheDesc cdesc;
        char fname[256];
        sprintf (fname, "cache/cache_blh_%E_%d%s.dat", M_BLH * bhs, idx, apdx);

        cdesc.mode = CONSTGRID;
        cdesc.x = 72;		//ro

        rhoMin = 1.e2;
        cdesc.xmin = rhoMin;

        rhoMax = 1.e15;
        cdesc.xmax = rhoMax;
        cdesc.xtype = LOGARITHMIC;

        cdesc.y = 72;		//T
//        TMin = 1.e6;
        TMin = 1.e4;
        cdesc.ymin = TMin;

        TMax = 1.e18;
        cdesc.ymax = TMax;
        cdesc.ytype = LOGARITHMIC;

        cdesc.autosave = 100;  //Cache save after number of points
        //Please note, Hd inside EOS is calculted as Hd=HdFac*sqrt(p/rho)
        double HdFac;

        {
            double X[NDIM];
            double r, th, ph;
            double cs, omega;

            coord (hstart, i, jstart, CENT, X);
            bl_coord (X, &r, &th, &ph);
            //  cs = sqrt((gam-1.0)*gam*u_eq/ro_eq)*CL;
            //  omega = CL*CL*CL/GNEWT/M_BLH/(a+pow(r, 1.5));

//              printf("i=%04d r=%E\n", i, r);

            cs = 1.0;
            omega = CL * CL * CL / GNEWT / M_BLH / (bhs + pow (r, 1.5));

            HdFac = cs / omega;
              printf("i=%04d r=%E HdFac=%E\n", i, r, HdFac);
        }

        cdesc.z = 1;		//HdFac
        cdesc.zmin = HdFac;
        cdesc.zmax = HdFac;
        cdesc.ztype = CUSTOM;

        cdesc.d = NCACHEVAL;	//tablice: e, p, Qnu

        cacheT[i] = new Cache4d (&cdesc, calc_cache_val, fname);
    }
#endif

    //Maintance
    {
#ifndef WIN32

        signal (SIGINT, ProgramStop);	//Ctrl-C
        signal (SIGQUIT, SIG_IGN);
        signal (SIGHUP, SIG_IGN);

        signal (SIGTERM, ProgramStop);	//kill
        signal (SIGCHLD, SIG_IGN);
#endif



    }
    return 0;
}

/************************************************************************/
int
restart_write_CoolEOS (FILE * fp)
/************************************************************************/
{
    double (*aLastT)[SN1][SN2] = (double (*)[SN1][SN2]) (_aLastT);
    int h, i, j;

    for (h = 0; h < SN3; ++h)
        for (i = 0; i < SN1; ++i)
            for (j = 0; j < SN2; ++j)
                fprintf (fp, FMT_DBL_OUT, aLastT[h][i][j]);

	if(cacheT)
	{
	    for (int i = 0; i < SN1; i++)
			cacheT[i]->save();
	}

    return 0;
}

/************************************************************************/
int
restart_read_CoolEOS (FILE * fp)
/************************************************************************/
{
    double (*aLastT)[SN1][SN2] = (double (*)[SN1][SN2]) (_aLastT);
    int h, i, j;

    for (h = 0; h < SN3; ++h)
        for (i = 0; i < SN1; ++i)
            for (j = 0; j < SN2; ++j)
                fscanf (fp, "%lf", &(aLastT[h][i][j]));

    return 0;
}
