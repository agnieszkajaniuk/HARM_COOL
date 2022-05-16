/***********************************************************************************

***********************************************************************************/
#include <unistd.h>
#include <sys/stat.h>

#include "decs.h"
#include "u2p_util.h"

#ifdef MPI_USED
#include "mpi.h"
#endif

#include <vector>

#include "cooleos.h"
#include "tracers.h"

struct of_tracer {
    int h,i,j ;
    double X[NDIM] ;
};

std::vector<of_tracer> vtr;

static double tracers_Rin;
static double tracers_Rout;
//static double tracers_Rmax;
static double tracers_ThetaMin;
static double tracers_ThetaMax;

static double tracers_Rank_Rin;
static double tracers_Rank_Rout;
static double tracers_Rank_Phin;
static double tracers_Rank_Phout;

static void dump_one_tracer(struct of_tracer *tracer, int hh, int ii, int jj);

int search_tracer(double *XX, int *hh, int *ii, int *jj)
{
    double Xl[NDIM], Xr[NDIM];
    int i,j,h;
    int ok;

    coord(hstart,irank==0?istart-1:istart,jrank==0?jstart-1:jstart,CENT,Xl) ;
    if(XX[RR] < Xl[RR] || XX[TH] < Xl[TH])
        return 0;

    *hh = 0;

    ok = 0;
    ZSLOOPI(-1,0) {
        coord(hstart,i,jstart,CENT,Xl);
        coord(hstart,i+1,jstart,CENT,Xr);

        if(XX[RR] >= Xl[RR] && XX[RR] < Xr[RR])
        {
            *ii = i;
            ok = 1;
            break;
        }
    }

    if(!ok) return 0;

    ok = 0;
    ZSLOOPJ(-1,0) {
        coord(hstart,istart,j,CENT,Xl);
        coord(hstart,istart,j+1,CENT,Xr);

        if(XX[TH] >= Xl[TH] && XX[TH] < Xr[TH])
        {
            *jj = j;
            ok = 1;
            break;
        }
    }
    if(!ok) return 0;

#if N3>2
    ok = 0;
    ZSLOOPH(0,0) {
        coord(h,istart,jstart,FACE3,Xl);
        coord(h+1,istart,jstart,FACE3,Xr);

        if(XX[PH] >= Xl[PH] && XX[PH] < Xr[PH])
        {
            *hh = h;
            ok = 1;
            break;
        }
    }
    if(!ok) return 0;

#endif

    return 1;
}

char *tracer_name(const char *loc, int h, int i, int j, char *name, size_t size) {
#if	N3>2
    snprintf(name, size, "%s/tracer_%03d_%03d_%03d.dat", loc, i, j, h);
#else
    snprintf(name, size, "%s/tracer_%03d_%03d.dat", loc, i, j);
#endif
    return name;
}

void init_tracers(void)
{
    double (*   p)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_p);

    struct of_tracer tracer;

    int i,j,h;
    int ii,jj,hh;
    int cnt = 0;
    double X[NDIM];

    //Lets predefine some values
    tracers_Rin = exp(startx[RR] - 0.5*dx[RR]) + R0;
    tracers_Rout = 0.8*(exp(startx[RR] + N1*dx[RR]) + R0);
//	tracers_Rmax = 11.8; //???
    tracers_ThetaMin = 0.02*M_PI;
    tracers_ThetaMax = 0.98*M_PI;

    coord(hstart,istart,jstart,CENT,X);
    tracers_Rank_Rin = X[RR];
    coord(hstart,istop+1,jstart,CENT,X);
    tracers_Rank_Rout = X[RR];

#if	N3>2
    coord(hstart,istart,jstart,CENT,X);
    tracers_Rank_Phin = X[PH];
    coord(hstop+1,istart,jstart,CENT,X);
    tracers_Rank_Phout = X[PH];
#endif

    printf("Tracers Rin=%E, Rout=%E\n", tracers_Rin, tracers_Rout);

#ifdef MPI_USED
    MPI_Sendrecv(_p, 1, slice_iNPR1_start, rank_previ, 1123,
                 _p, 1, slice_iNPR1_stopN, rank_nexti, 1123,
                 MPI_COMM_WORLD, &MPIStatus);

#if	N3>2
    {
        MPI_Sendrecv(_p, 1, slice_hNPR1_start, rank_prevh, 1125,
                     _p, 1, slice_hNPR1_stopN, rank_nexth, 1125,
                     MPI_COMM_WORLD, &MPIStatus);
    }
#endif

#endif

    for(h=0; h<N3; h++)for(i=0; i<N1; i++)for(j=0; j<N2; j++) {
//printf("init_tracers i=%d, j=%d, h=%d\n", i,j,h);
                tracer.h = h;
                tracer.i = i;
                tracer.j = j;
                tracer.X[TT] = 0.;

                tracer.X[RR] = startx[RR] + (i + 0.5)*dx[RR] ;
                tracer.X[TH] = startx[TH] + (j + 0.5)*dx[TH] ;
                tracer.X[PH] = startx[PH] + h*dx[PH] ;

                if(search_tracer(tracer.X, &hh, &ii, &jj)) {
//printf("init_tracers ii=%d, jj=%d, hh=%d, rho=%E, rhomin=%E\n", ii,jj,hh, p[hh][ii][jj][RHO], RHOMIN);
//				    double X[NDIM], r, th, ph;

//		        	coord(hh,ii,jj,CENT,X);
//        			bl_coord(X,&r,&th,&ph) ;

//                    if(p[hh][ii][jj][RHO] > RHOMIN || r <= tracers_Rmax) {
                    if(p[hh][ii][jj][RHO] > RHOMIN) {
                        vtr.push_back(tracer);
                        cnt++;
                    }
                }
            }

#ifdef MPI_USED
    MPI_Allreduce(MPI_IN_PLACE,&cnt,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
    if(rank == 0)
        printf("init_tracers: Number of active tracers=%d\n", cnt);

}

//
//  Evaluate the piecewise linear function in 3D.
//
void linear_interpolation(double tracerX[NDIM], const int hh, const int ii, const int jj, const int VALDIM, double *val,
                          void (*val_fun)(int h, int i, int j, double *val))
{

    const int WRR=2;
    const int WTH=2;
    const int WPH=(N3 > 2 ? 2 : 1);

    double X[NDIM];
    double XX[WPH][WRR][WTH][NDIM];

    double valt[WPH][WRR][WTH][VALDIM];
    double tmp1[WPH][WRR][VALDIM];
    double tmp2[WPH][VALDIM];

    for(int l3=0; l3<WPH; l3++)
        for(int l1=0; l1<WRR; l1++)
            for(int l2=0; l2<WTH; l2++) {
                coord(hh+l3,ii+l1,jj+l2,CENT,XX[l3][l1][l2]) ;

                val_fun(hh+l3,ii+l1,jj+l2, valt[l3][l1][l2]);
            }

    for(int l3=0; l3<WPH; l3++)
        for(int l1=0; l1<WRR; l1++) {

            double T1 = 1.0/(XX[l3][l1][1][TH] - XX[l3][l1][0][TH]);
            double T2 = tracerX[TH] - XX[l3][l1][0][TH];

            for (int  k = 0; k < VALDIM; k++ )
            {
                double ypval = ( valt[l3][l1][1][k] - valt[l3][l1][0][k]) *T1;
                tmp1[l3][l1][k] = valt[l3][l1][0][k] + T2 * ypval;
            }
        }

    for(int l3=0; l3<WPH; l3++) {
        double T1 = 1.0/(XX[l3][1][1][RR] - XX[l3][0][0][RR]);
        double T2 = tracerX[RR] - XX[l3][0][0][RR];

        for (int  k = 0; k < VALDIM; k++ )
        {
            double ypval = ( tmp1[l3][1][k] - tmp1[l3][0][k]) *T1;
            tmp2[l3][k] = tmp1[l3][0][k] + T2 * ypval;
        }
    }
#if N3 <= 2
    for(int k=0; k<VALDIM; k++) {
        val[k] = tmp2[0][k];
    }
#else
    {
        double T1 = 1.0/(XX[1][1][1][PH] - XX[0][0][0][PH]);
        double T2 = tracerX[PH] - XX[0][0][0][PH];

        for (int  k = 0; k < VALDIM; k++ )
        {
            double ypval = ( tmp2[1][k] - tmp2[0][k]) *T1;
            val[k] = tmp2[0][k] + T2 * ypval;
        }
    }
#endif

}


void ucon_fun(int h, int i, int j, double *ucon)
{
    double (*   p)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_p);
    struct of_geom geom ;

    get_geometry(h,i,j,CENT,&geom) ;
    ucon_calc(p[h][i][j], &geom, ucon) ;

    ucon[RR] /= ucon[TT];
    ucon[TH] /= ucon[TT];
    ucon[PH] /= ucon[TT];
}

#ifdef MPI_USED
void Tracers_Sendrecv(std::vector<struct of_tracer> &vectr, int rank_to, int rank_from)
{
    int trsend = vectr.size();
    int trrecv = 0;

    MPI_Sendrecv(&trsend, 1, MPI_INT, rank_to, 1128,
                 &trrecv, 1, MPI_INT, rank_from, 1128,
                 MPI_COMM_WORLD, &MPIStatus);
//	printf("\nTracers_Sendrecv send %d to %d, recv %d from %d\n", trsend, rank_to, trrecv, rank_from);

    if(trsend > 0)
        MPI_Send(&vectr[0], trsend*sizeof(struct of_tracer), MPI_BYTE, rank_to, 1129, MPI_COMM_WORLD);

    if(trrecv > 0) {
        int cursize = vtr.size();
        vtr.resize(cursize + trrecv);

        MPI_Recv(&vtr[cursize], trrecv*sizeof(struct of_tracer), MPI_BYTE, rank_from, 1129, MPI_COMM_WORLD, &MPIStatus);
    }
#ifdef MPI_USED
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}
#endif

void update_tracers(double Dt)
{
    std::vector<struct of_tracer>::iterator iter;
    struct of_tracer *tracer;

    std::vector<of_tracer> vtr_Rin;
    std::vector<of_tracer> vtr_Rout;

    int k;
    int ii,jj,hh;
    int cnt = 0;
    double ucon[NDIM];

#ifdef MPI_USED
    MPI_Sendrecv(_p, 1, slice_iNPR1_start, rank_previ, 1123,
                 _p, 1, slice_iNPR1_stopN, rank_nexti, 1123,
                 MPI_COMM_WORLD, &MPIStatus);
#if	N3>2
    {
        MPI_Sendrecv(_p, 1, slice_hNPR1_start, rank_prevh, 1125,
                     _p, 1, slice_hNPR1_stopN, rank_nexth, 1125,
                     MPI_COMM_WORLD, &MPIStatus);
    }
#endif

#endif


    for (iter = vtr.begin(); iter != vtr.end(); ) {
        tracer = &(*iter);

        if(!search_tracer(tracer->X, &hh, &ii, &jj))
        {
            printf("\nupdate_tracers: invalid tracer i,j,h=%d,%d,%d r=%E, th=%E, in rank %d, removed\n",
                   tracer->i,tracer->j,tracer->h, tracer->X[RR], tracer->X[TH], rank);
            iter = vtr.erase(iter);
            continue;
        }

        linear_interpolation(tracer->X, hh, ii, jj, NDIM, ucon, ucon_fun);

        for(int k=1; k<NDIM; k++) {
            tracer->X[k] += ucon[k]*Dt;
        }

        double r, th, ph;
        bl_coord(tracer->X,&r,&th,&ph) ;

        if(r < tracers_Rin || th < 0 || th > M_PI) {
            dump_one_tracer(tracer, hh, ii, jj);

            char tc[256], tn[256];
            printf("\nupdate_tracers: Tracer In i=%d, j=%d, h=%d\n", tracer->i,tracer->j,tracer->h);

            tracer_name("tracers", tracer->h, tracer->i, tracer->j, tc, sizeof(tc));
            tracer_name("tracers/in", tracer->h, tracer->i, tracer->j, tn, sizeof(tn));
            rename (tc, tn);
            iter = vtr.erase(iter);
        } else if(r >= tracers_Rout) {
            dump_one_tracer(tracer, hh, ii, jj);

            char tc[256], tn[256];
            tracer_name("tracers", tracer->h, tracer->i, tracer->j, tc, sizeof(tc));

//            printf("\nupdate_tracers: i=%d, j=%d, h=%d, r=%E, th=%E, tracers_Rout=%E, tracers_ThetaMin=%E, tracers_ThetaMax=%E x[RR]=%E, X[TH]=%E, x=%E, y=%E\n",
//                   tracer->i,tracer->j,tracer->h, r,th, tracers_Rout, tracers_ThetaMin, tracers_ThetaMax, tracer->X[RR], tracer->X[TH], r*sin(th), r*cos(th));

            if(th >= tracers_ThetaMin && th <= tracers_ThetaMax)
            {
                printf("\nupdate_tracers: Tracer Out i=%d, j=%d, h=%d\n", tracer->i,tracer->j,tracer->h);
                tracer_name("tracers/out", tracer->h, tracer->i, tracer->j, tn, sizeof(tn));
            } else {
                printf("\nupdate_tracers: Tracer Jet i=%d, j=%d, h=%d\n", tracer->i,tracer->j,tracer->h);
                tracer_name("tracers/jet", tracer->h, tracer->i, tracer->j, tn, sizeof(tn));
            }
            rename (tc, tn);
            iter = vtr.erase(iter);
        } else if(irank > 0 && tracer->X[RR] < tracers_Rank_Rin) {
            vtr_Rin.push_back(*tracer);
            iter = vtr.erase(iter);
        } else if(irank < isize-1 && tracer->X[RR] >= tracers_Rank_Rout) {
            vtr_Rout.push_back(*tracer);
            iter = vtr.erase(iter);
        } else {
            ++cnt;
            ++iter;
        }
    }



#ifdef MPI_USED
//    MPI_Allreduce(MPI_IN_PLACE,&cnt,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
//    if(rank == 0)
//       printf("update_tarcers: Number of active tracers=%d\n", cnt);

#ifdef MPI_USED
    Tracers_Sendrecv(vtr_Rin, rank_previ, rank_nexti);
    Tracers_Sendrecv(vtr_Rout, rank_nexti, rank_previ);
#endif

#if	N3>2
    std::vector<of_tracer> vtr_Phin;
    std::vector<of_tracer> vtr_Phout;

    for (iter = vtr.begin(); iter != vtr.end(); ) {
        tracer = &(*iter);

        if(hrank == 0 && tracer->X[PH] < 0.) {
            tracer->X[PH] += 2.*M_PI;
            vtr_Phin.push_back(*tracer);
            iter = vtr.erase(iter);
        } else if(hrank == 0 && tracer->X[PH] > 2.*M_PI) {
            tracer->X[PH] -= 2.*M_PI;
            vtr_Phout.push_back(*tracer);
            iter = vtr.erase(iter);
        } else if(tracer->X[PH] < tracers_Rank_Phin) {
            vtr_Phin.push_back(*tracer);
            iter = vtr.erase(iter);
        } else if(tracer->X[PH] >= tracers_Rank_Phout) {
            vtr_Phout.push_back(*tracer);
            iter = vtr.erase(iter);
        } else {
            ++iter;
        }
    }
    Tracers_Sendrecv(vtr_Phin, rank_prevh, rank_nexth);
    Tracers_Sendrecv(vtr_Phout, rank_nexth, rank_prevh);
#endif

}


void rhoTYe_fun(int h, int i, int j, double *val)
{
    double (*   p)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_p);

    double rho = p[h][i][j][RHO] * RHO_UNIT;

#if( COOL_EOS )

    double (*   CoolVal)[SN1][SN2][NCOOLVAL] = (double (*) [SN1][SN2][NCOOLVAL])(_CoolVal);
    InitCache_CoolEOS(h, i, j);
    CoolVal_rho0_u_CoolEOS(p[h][i][j][RHO], p[h][i][j][UU], h, i, j, CoolVal[h][i][j]);

    double T = CoolVal[h][i][j][COOL_T];
    double ye = CoolVal[h][i][j][COOL_Yee];
#else
    const double CL = 2.99792458e10;
    const double KBOL = 1.38e-16;
    const double mn = 1.67492729e-24;
    const double me = 9.109565e-28;
    const double mp = 1.673e-24;
    const double mhe = 6.6465e-24; // 2p+2e+binding energy
    const double mav = (mp+mn+2.*me+mhe)/5.; //avrage particle mass

    double T = (gam-1.)*p[h][i][j][UU]/p[h][i][j][RHO]/KBOL*mav*CL*CL;
    double ye = -1.;

//   if (isnan(T)) {
//		printf("rhoTYe_fun T isnan! u=%E, rho=%E, h=%d, i=%d, j=%d\n", p[h][i][j][UU], p[h][i][j][RHO], h, i, j);
//   }

#endif
    val[0] = rho;
    val[1] = T;
    val[2] = ye;
}

void dump_one_tracer(struct of_tracer *tracer, int hh, int ii, int jj)
{
    char tc[256];

    tracer_name("tracers", tracer->h, tracer->i, tracer->j, tc, sizeof(tc));

    FILE *fp = fopen(tc, "a+");
    if(fp == NULL) {
        printf("Can't open file %s\n", tc);
        return;
    }

    double tim = t * T_UNIT;
    double val[3];

    linear_interpolation(tracer->X, hh, ii, jj, 3, val, rhoTYe_fun);
    double rho = val[0];
    double T   = val[1];
    double ye  = val[2];

//   if (isnan(T)) {
//		printf("dump_one_tracer T isnan! rho=%E, h=%d, i=%d, j=%d, hh=%d, ii=%d, jj=%d, istop=%d\n", rho, h, i, j, hh, ii, jj, istop);
//		exit(0);
//   }

    double r, th, ph;

    bl_coord(tracer->X,&r,&th,&ph) ;

#if N3>2
    if(ph > 2.*M_PI)
        ph -= 2.*M_PI;
    else if (ph < 0.)
        ph += 2.*M_PI;
#else
    ph = 0.0;
#endif

    fprintf(fp, FMT_DBL_OUT FMT_DBL_OUT FMT_DBL_OUT FMT_DBL_OUT FMT_DBL_OUT FMT_DBL_OUT FMT_DBL_OUT "\n", tim, rho, T, ye, r, th, ph);

    fclose(fp);

}

void dump_tracers(void) {
    std::vector<struct of_tracer>::iterator iter;
    struct of_tracer *tracer;

    int k;
    int ii,jj,hh;
    int cnt = 0;


#ifdef MPI_USED
    MPI_Sendrecv(_p, 1, slice_iNPR1_start, rank_previ, 1126,
                 _p, 1, slice_iNPR1_stopN, rank_nexti, 1126,
                 MPI_COMM_WORLD, &MPIStatus);
#if	N3>2
    {
        MPI_Sendrecv(_p, 1, slice_hNPR1_start, rank_prevh, 1127,
                     _p, 1, slice_hNPR1_stopN, rank_nexth, 1127,
                     MPI_COMM_WORLD, &MPIStatus);
    }
#endif

#endif

    for (iter = vtr.begin(); iter != vtr.end(); ) {
        tracer = &(*iter);
        if(!search_tracer(tracer->X, &hh, &ii, &jj))
        {
            printf("\ndump_tracers: invalid tracer in rank %d, removed\n", rank);
            iter = vtr.erase(iter);
            continue;
        }

        dump_one_tracer(tracer, hh, ii, jj);

        ++cnt;
        ++iter;
    }

#ifdef MPI_USED
//    MPI_Allreduce(MPI_IN_PLACE,&cnt,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
//    if(rank == 0)
//       printf("Number of active tracers=%d\n", cnt);

#ifdef MPI_USED
//    MPI_Barrier(MPI_COMM_WORLD);
#endif

}


/************************************************************************/
int
restart_write_tracers (FILE * fp)
/************************************************************************/
{
    std::vector<struct of_tracer>::iterator iter;
    int k;

    fprintf(fp, FMT_DBL_OUT, tracers_Rin);
    fprintf(fp, FMT_DBL_OUT, tracers_Rout);
    fprintf(fp, FMT_DBL_OUT, tracers_ThetaMin);
    fprintf(fp, FMT_DBL_OUT, tracers_ThetaMax);
    fprintf(fp, FMT_DBL_OUT, tracers_Rank_Rin);
    fprintf(fp, FMT_DBL_OUT, tracers_Rank_Rout);
    fprintf(fp, FMT_DBL_OUT, tracers_Rank_Phin);
    fprintf(fp, FMT_DBL_OUT, tracers_Rank_Phout);

    fprintf(fp, FMT_INT_OUT, (int)vtr.size());

    for (iter = vtr.begin(); iter != vtr.end(); ) {
        struct of_tracer *tracer = &(*iter);
        fprintf(fp, FMT_INT_OUT, tracer->i);
        fprintf(fp, FMT_INT_OUT, tracer->j);
        fprintf(fp, FMT_INT_OUT, tracer->h);
        for(k=0; k<NDIM; k++) fprintf(fp, FMT_DBL_OUT, tracer->X[k]);

		//Store tracer's file size
		{ 
	        char tname[256];
    	    tracer_name("tracers", tracer->h, tracer->i, tracer->j, tname, sizeof(tname));

			struct stat st;
			stat(tname, &st);
			int size = st.st_size;

	        fprintf(fp, FMT_INT_OUT, size);
		}
		++iter;
    }

    return 0;
}

/************************************************************************/
int
restart_read_tracers (FILE * fp)
/************************************************************************/
{
    std::vector<struct of_tracer>::iterator iter;
    int i,k;
    int cnt;

    fscanf(fp, "%lf", &tracers_Rin);
    fscanf(fp, "%lf", &tracers_Rout);
    fscanf(fp, "%lf", &tracers_ThetaMin);
    fscanf(fp, "%lf", &tracers_ThetaMax);
    fscanf(fp, "%lf", &tracers_Rank_Rin);
    fscanf(fp, "%lf", &tracers_Rank_Rout);
    fscanf(fp, "%lf", &tracers_Rank_Phin);
    fscanf(fp, "%lf", &tracers_Rank_Phout);

    fscanf(fp, "%d",  &cnt);

    for(i=0; i<cnt; i++) {
        struct of_tracer tracer;
        fscanf(fp, "%d",  &tracer.i);
        fscanf(fp, "%d",  &tracer.j);
        fscanf(fp, "%d",  &tracer.h);
        for(k=0; k<NDIM; k++) fscanf(fp, "%lf", &(tracer.X[k]));

	    //Restore tracer's file size
	    { 
		int size;
	        fscanf(fp, "%d", &size);

	        char tname[256];
    	    tracer_name("tracers", tracer.h, tracer.i, tracer.j, tname, sizeof(tname));
            if(access(tname, 0))
			{
				//Tracer file does not exists, skip this tracer
                printf("\nrestart_read_tracers: invalid tracer %s in rank %d (probably moved to out or in), tracer removed\n", tname, rank);
				continue;
			}

			truncate(tname, size);
		}

        vtr.push_back(tracer);
    }

    return 0;
}
