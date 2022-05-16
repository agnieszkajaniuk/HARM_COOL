#if !defined(MNEWT_H__INCLUDED_)
#define MNEWT_H__INCLUDED_
#include <math.h>
#include <stdio.h>

#include "random.h"

template <int N> inline
void lubksb(double a[N+1][N+1],int indx[N+1], double b[N+1])
{
    int i,ii=0,ip,j;
    double sum;

    for (i=1;i<=N;i++) {
        ip=indx[i];
        sum=b[ip];
        b[ip]=b[i];
        if (ii)
            for (j=ii;j<=i-1;j++)
                sum -= a[i][j]*b[j];
        else if (sum)
            ii=i;
        b[i]=sum;
    }
    for (i=N;i>=1;i--) {
        sum=b[i];
        for (j=i+1;j<=N;j++)
            sum -= a[i][j]*b[j];
        b[i]=sum/a[i][i];
    }
}


#define MNEWTTINY 1.0e-20

template <int N> inline
int ludcmp(double a[N+1][N+1],int indx[N+1],double *d)
{
    int i,imax=0,j,k;
    double big,dum,sum,temp;
    double vv[N+1];

    *d=1.0;
    for (i=1;i<=N;i++) {
        big=0.0;
        for (j=1;j<=N;j++)
            if ((temp=fabs(a[i][j])) > big)
                big=temp;

        //                if (big == 0.0) nrerror("Singular matrix in routine LUDCMP");
        if (big == 0.0) {
            printf("Singular matrix in routine LUDCMP\n");
	        for (j=1;j<=N;j++)
			{
    	        printf("a[%d][%d]=%E\n", i, j, a[i][j]);
			}
            return 0;
//			exit(0);
        }
        vv[i]=1.0/big;
    }
    for (j=1;j<=N;j++) {
        for (i=1;i<j;i++) {
            sum=a[i][j];
            for (k=1;k<i;k++)
                sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
        }
        big=0.0;
        for (i=j;i<=N;i++) {
            sum=a[i][j];
            for (k=1;k<j;k++)
                sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
            if ( (dum=vv[i]*fabs(sum)) >= big) {
                big=dum;
                imax=i;
            }
        }
        if (j != imax) {
            for (k=1;k<=N;k++) {
                dum=a[imax][k];
                a[imax][k]=a[j][k];
                a[j][k]=dum;
            }
            *d = -(*d);
            vv[imax]=vv[j];
        }
        indx[j]=imax;
        if (a[j][j] == 0.0)
            a[j][j]=MNEWTTINY;
        if (j != N) {
            dum=1.0/(a[j][j]);
            for (i=j+1;i<=N;i++)
                a[i][j] *= dum;
        }
    }

    return 1;
}

#undef MNEWTTINY


template <int N> inline
int mnewtc(int ntrial, double x[N+1], double tolf, double minx[N+1], double maxx[N+1],
           void (*usrfun)(double [N+1], double [N+1][N+1], double [N+1], void *), void *usrfunpar, double *berr, int nlos, RANDOMDEF *rd)
{
    int k, i, l;
	int indx[N+1];
    double errf;
    double d;
	double bet[N+1];
    double alpha[N+1][N+1];
    double bestval[N+1];
    double besterr = 1.0E100;


    int retval = 0;

    for (i=1;i<=N;i++)
   	    bestval[i] = x[i];

//    indx=ivector(1,n);
//    bet=vector(1,n);
//    alpha=matrix(1,n,1,n);

    for(l = 0; l < nlos; ++l) {
//        printf("New los %d\n", l+1);
        for (k=1;k<=ntrial;k++) {

            usrfun(x,alpha, bet, usrfunpar);
            errf=0.0;
            for (i=1;i<=N;i++)
                errf += fabs(bet[i]);

//            printf("mnet x1=%E x2=%E errf=%E\n", x[1], x[2], errf);
            //getc(stdin);


            if (errf <= tolf) {
                *berr = errf;
                retval = 1;
                goto koniec;
            }

/*            
            if(errf > 1.0e6)
            {
//	    		printf("Errf =%E is too large\n", errf);
    	    	break;
    	    }
*/
            if(errf < besterr) {
                if(k > 1)
                    retval = 2;

                besterr = errf;
                for (i=1;i<=N;i++)
                    bestval[i] = x[i];
            }

            if(ludcmp<N>(alpha,indx,&d))
	            lubksb<N>(alpha,indx,bet);

			int newvalcnt = 0;
            for (i=1;i<=N;i++) 
			{
                double newval = x[i] + bet[i];
//				printf("i = %d x[i]=%E bet[i]=%E -> %E\n", i, x[i], bet[i], newval);

                if(newval >= minx[i] && newval <= maxx[i] && newval != x[i]) {
                    x[i] = newval;
					newvalcnt++;
            }
			if(!newvalcnt)
				break;
            }
        }

        if(rd && l<nlos-1) 
        for (i=1;i<=N;i++) {
            //            double newval = bestval[i]*pow(10.0, (0.5-(double)rand()/RAND_MAX)*1.5);
            //            if(newval >= minx[i] && newval <= maxx[i])
            //                x[i] = newval;

            x[i] = UniformRandomRangeLog(rd, minx[i], maxx[i]);

        }

    }


    for (i=1;i<=N;i++)
   	    x[i] = bestval[i];

    *berr = besterr;


koniec:
//    free_matrix(alpha,1,n,1,n);
//    free_vector(bet,1,n);
//    free_ivector(indx,1,n);
    return retval;
}


#endif
