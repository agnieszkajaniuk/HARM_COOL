#if !defined(AKIMA_H__INCLUDED_)
#define AKIMA_H__INCLUDED_


void spline_akima ( int ndim, int dimmask[], int ndata, double tdata[], double ydata[],
                    double tval, double yval[] );
                    
void spline_akima_der ( int ndim, int dimmask[], int ndata, double tdata[], double ydata[],
                    double tval, double yval[],  double dval[]);
                    
int interpolate_akima(
    int si,            /* size of input arrays */
    char *xi, int dxi, /* x coordinates and stride */
    char *yi, int dyi, /* y coordinates and stride */
    int so,            /* size of output arrays */
    char *xo, int dxo, /* x coordinates of output and stride */
    char *yo, int dyo, /* y output coordinates and stride */
    double *p           /* buffer for polynomial coefficients of size 4*si+4 */
);
                    
void dvec_bracket ( int n, double x[], double xval, int *left,
                    int *right );
void spline_linear_val2 (int ndim,  int ndata, double tdata[], double ydata[],
                         double tval, double *yval);
#endif