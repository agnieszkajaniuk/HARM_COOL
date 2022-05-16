#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*
A new method of interpolation and smooth curve fitting based on local
procedures. Hiroshi Akima, J. ACM, October 1970, 17(4), 589-602.
*/

int interpolate_akima(
    int si,            /* size of input arrays */
    char *xi, int dxi, /* x coordinates and stride */
    char *yi, int dyi, /* y coordinates and stride */
    int so,            /* size of output arrays */
    char *xo, int dxo, /* x coordinates of output and stride */
    char *yo, int dyo, /* y output coordinates and stride */
    double *p           /* buffer for polynomial coefficients of size 4*si+4 */
)
{
    int i;
    size_t s;
    double x0, x1, x2, x3;       /* extrapolated x values */
    double y0, y1, y2, y3;       /* extrapolated y values */
    double t0, t1, t2, t3;       /* temporary values */
    double d0, d1;
    double g0, g1;               /* gradients at extrapolated values */
    double *p0, *p1, *p2, *p3;   /* buffer pointers. p3 holds slopes */
    char *pxi, *pyi, *pxo, *pyo; /* data pointers */

    p0 = p;
    p1 = p + si + 1;
    p2 = p + si*2 + 2;
    p3 = p + si*3 + 3;

    /* slopes of input data */
    pxi = xi;
    pyi = yi;
    t0 = *((double *)pxi);
    t1 = *((double *)pyi);
    i = si-1;
    while (i--) {
        pxi += dxi;
        pyi += dyi;
        t0 = *((double *)pxi) - t0;
        if (t0 < 1e-12)
		{
			printf("t0 =%E < 1e-12\n", t0);
            return -1;
		}
        *p3++ = (*((double *)pyi) - t1) / t0;
        t0 = *((double *)pxi);
        t1 = *((double *)pyi);
    }
    p3 = p + si*3 + 3;

    /* extrapolate 2 points on left side */
    t0 = *((double *)(xi));
    t1 = *((double *)(xi + dxi));
    t2 = *((double *)(yi));
    t3 = *((double *)(yi + dyi));

    x1 = t0 + t1 - *((double *)(xi+dxi+dxi));
    x0 = x1 + t0 - t1;
    y1 = (t0 - x1) * (p3[1] - 2.0*p3[0]) + t2;
    g1 = (t2 - y1) / (t0 - x1);
    y0 = (x1 - x0) * (p3[0] - 2.0*g1) + y1;
    g0 = (y1 - y0) / (x1 - x0);

    /* extrapolate 2 points on right side */
    s = (size_t)xi + dxi*(si-1);
    t0 = *((double *)(s - dxi));
    t1 = *((double *)(s));
    x2 = t1 + t0 - *((double *)(s - dxi - dxi));
    x3 = x2 + t1 - t0;

    s = (size_t)yi + dyi*(si-1);
    t2 = *((double *)(s - dyi));
    t3 = *((double *)(s));
    y2 = (2.0*p3[si-2] - p3[si-3]) * (x2 - t1) + t3;
    p3[si-1] = (y2 - t3) / (x2 - t1);
    y3 = (2.0*p3[si-1] - p3[si-2]) * (x3 - x2) + y2;
    p3[si] = (y3 - y2) / (x3 - x2);

    /* slopes */
    t1 = g0;
    t2 = g1;
    t3 = *p3;
    i = si;
    while (i--) {
        t0 = t1;
        t1 = t2;
        t2 = t3;
        t3 = *(p3 + 1);
        d0 = t3 - t2;
        if (d0 < 0.0)
            d0 *= -1.0;
        d1 = t1 - t0;
        if (d1 < 0.0)
            d1 *= -1.0;
        if ((d0 + d1) < 1e-9) {
            *p3++ = 0.5 * (t1 + t2);
        } else {
            *p3++ = (d0*t1 + d1*t2) / (d0 + d1);
        }
    }
    /* polynomial coefficients */
    pxi = xi;
    pyi = yi;
    t0 = *((double *)pxi);
    t1 = *((double *)pyi);
    p3 = p + si*3 + 3;
    g1 = *p3++;
    i = si;
    while (i--) {
        pxi += dxi;
        pyi += dyi;
        d0 = (*((double *)pxi) - t0);
        d1 = (*((double *)pyi) - t1);
        t2 = d1 / d0;
        g0 = g1;
        g1 = *p3++;
        *p0++ = t1;
        *p1++ = (3.0*t2 - 2.0*g0 - g1) / d0;
        *p2++ = (g0 + g1 - 2.0*t2) / (d0*d0);
        t0 = *((double *)pxi);
        t1 = *((double *)pyi);
    }
    /* interpolate data */
    p0 = p;
    p1 = p + si + 1;
    p2 = p + si*2 + 2;
    p3 = p + si*3 + 3;
    pxi = xi;
    pyi = yi;
    pxo = xo;
    pyo = yo;
    si -= 2;
    i = -1;
    s = so;
    while (s--) {
        t0 = *((double *)pxo);
//		printf("akima 1: i=%d t0=%E *pxi=%E\n", i, t0, *((double *)pxi)); 
        while ((t0 > *((double *)pxi)) && (i < si)) {
            pxi += dxi;
            i++;
//		printf("akima 2: i=%d t0=%E *pxi=%E\n", i, t0, *((double *)pxi)); 
        }
//		printf("akima 3: i=%d t0=%E *pxi=%E\n", i, t0, *((double *)pxi)); 
        if (i < 0) {
            i = 0;
            pxi = xi + dxi;
        }
        t1 = t0 - *((double *)(pxi - dxi));
//		printf("akima 3: i=%d t0=%E t1=%E - %E %E %E %E\n", i, t0, t1, p0[i], p3[i]*t1, p1[i]*t1*t1, p2[i]*t1*t1*t1); 
	

        *((double *)pyo) = p0[i] + p3[i]*t1 + p1[i]*t1*t1 + p2[i]*t1*t1*t1;
        pyo += dyo;
        pxo += dxo;
    }
    return 0;
}

void spline_akima ( int ndim, int dimmask[], int ndata, double tdata[], double ydata[],
                    double tval, double yval[] )
{
	double p[4*ndata+4];
/*
    for(int i=0; i<ndata; i++)
	{
		printf("i=%d, data=%E, tval=%E, xnuc=%E\n", i, tdata[i], tval, ydata[4+i*ndim]);
	}
*/
    for(int i=0; i<ndim; i++)
    {
//		memset(&p[0], 0, (4*ndata+4)*sizeof(double));

//		printf("Calc for i=%d\n", i);
		if(!dimmask || dimmask[i])
		{
	        interpolate_akima(
    	        ndata,            /* size of input arrays */
        	    (char *)tdata, sizeof(double), /* x coordinates and stride */
            	(char *)(&ydata[i]), ndim*sizeof(double), /* y coordinates and stride */
	            1,            /* size of output arrays */
    	        (char *)&tval, sizeof(double), /* x coordinates of output and stride */
            	(char *)(&yval[i]), sizeof(double), /* y output coordinates and stride */
        	    p           /* buffer for polynomial coefficients of size 4*si+4 */
	        );
		}

//		if(i == 4)
//	        printf("Ret =%d for i=%d val=%E\n", ret, i, yval[i]);
    }
}

void spline_akima_der ( int ndim, int dimmask[], int ndata, double tdata[], double ydata[],
                    double tval, double yval[], double dval[] )
{
	double p[4*ndata+4];
	double tval3[3];
	double yval3[3];

    const double hh = 1.0E-6;
    const double h = tval * hh;


	tval3[0] = tval - h;
	tval3[1] = tval;
	tval3[2] = tval + h;

    for(int i=0; i<ndim; i++)
    {
//		printf("Calc for i=%d\n", i);
		if(!dimmask || dimmask[i])
		{
	        interpolate_akima(
    	        ndata,            /* size of input arrays */
        	    (char *)tdata, sizeof(double), /* x coordinates and stride */
            	(char *)(&ydata[i]), ndim*sizeof(double), /* y coordinates and stride */
	            3,            /* size of output arrays */
    	        (char *)tval3, sizeof(double), /* x coordinates of output and stride */
            	(char *)(yval3), sizeof(double), /* y output coordinates and stride */
        	    p           /* buffer for polynomial coefficients of size 4*si+4 */
	        );
			yval[i]  = yval3[1];
			dval[i]  = (yval3[2]-yval3[0]) / (2.0*h);
		}
    }
}


//********************************************************************

void dvec_bracket ( int n, double x[], double xval, int *left,
                    int *right )

//********************************************************************
//
//  Purpose:
//
//    DVEC_BRACKET searches a sorted array for successive brackets of a value.
//
//  Discussion:
//
//    If the values in the vector are thought of as defining intervals
//    on the real line, then this routine searches for the interval
//    nearest to or containing the given value.
//
//    It is always true that RIGHT = LEFT+1.
//
//    If XVAL < X[0], then LEFT = 1, RIGHT = 2, and
//      XVAL   < X[0] < X[1];
//    If X(1) <= XVAL < X[N-1], then
//      X[LEFT-1] <= XVAL < X[RIGHT-1];
//    If X[N-1] <= XVAL, then LEFT = N-1, RIGHT = N, and
//      X[LEFT-1] <= X[RIGHT-1] <= XVAL.
//
//    For consistency, this routine computes indices RIGHT and LEFT
//    that are 1-based, although it would be more natural in C and
//    C++ to use 0-based values.
//
//  Modified:
//
//    24 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, length of input array.
//
//    Input, double X[N], an array that has been sorted into ascending order.
//
//    Input, double XVAL, a value to be bracketed.
//
//    Output, int *LEFT, *RIGHT, the results of the search.
//
{
    int i;

    for ( i = 2; i <= n - 1; i++ ) {
        if ( xval < x[i-1] ) {
            *left = i - 1;
            *right = i;
            return;
        }

    }

    *left = n - 1;
    *right = n;

    return;
}

//**********************************************************************

void spline_linear_val2 (int ndim,  int ndata, double tdata[], double ydata[],
                         double tval, double *yval)
//**********************************************************************
//
//  Purpose:
//
//    SPLINE_LINEAR_VAL evaluates a piecewise linear spline at a point.
//
//  Discussion:
//
//    Because of the extremely simple form of the linear spline,
//    the raw data points ( TDATA(I), YDATA(I)) can be used directly to
//    evaluate the spline at any point.  No processing of the data
//    is required.
//
//  Modified:
//
//    24 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NDATA, the number of data points defining the spline.
//
//    Input, double TDATA[NDATA], YDATA[NDATA], the values of the independent
//    and dependent variables at the data points.  The values of TDATA should
//    be distinct and increasing.
//
//    Input, double TVAL, the point at which the spline is to be evaluated.
//
//    Output, double *YVAL, *YPVAL, the value of the spline and its first
//    derivative dYdT at TVAL.  YPVAL is not reliable if TVAL is exactly
//    equal to TDATA(I) for some I.
//
{
    int left;
    int right;
    double ypval;
    //
    //  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
    //  nearest to, TVAL.
    //
    dvec_bracket ( ndata, tdata, tval, &left, &right );
    //
    //  Now evaluate the piecewise linear function.
    //
    
/*
    for (int  i = 0; i < ndim; i++ )
    {     
        ypval = ( ydata[i+(right-1)*ndim] - ydata[i+(left-1)*ndim] )
             / ( tdata[right-1] - tdata[left-1] );

		yval[i] = ydata[i+(left-1)*ndim] +  ( tval - tdata[left-1] ) * ypval;
    }
*/    
    double T1 = 1.0/(tdata[right-1] - tdata[left-1]);
    double T2 = tval - tdata[left-1];
    
    for (int  i = 0; i < ndim; i++ )
    {     
        ypval = ( ydata[i+(right-1)*ndim] - ydata[i+(left-1)*ndim] ) * T1;

		yval[i] = ydata[i+(left-1)*ndim] + T2 * ypval;
    }
    
    return;
}
