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


#include "u2p_defs.h"

extern void primtoU_g( FTYPE prim[], FTYPE gcov[NDIM][NDIM], FTYPE gcon[NDIM][NDIM], FTYPE gdet,  FTYPE U[] );
extern void ucon_calc_g(FTYPE prim[],FTYPE gcov[NDIM][NDIM],FTYPE gcon[NDIM][NDIM],FTYPE ucon[]);
extern void raise_g(FTYPE vcov[], FTYPE gcon[NDIM][NDIM], FTYPE vcon[]);
extern void lower_g(FTYPE vcon[], FTYPE gcov[NDIM][NDIM], FTYPE vcov[]);
extern void ncov_calc(FTYPE gcon[NDIM][NDIM],FTYPE ncov[]) ;
extern void bcon_calc_g(FTYPE prim[],FTYPE ucon[],FTYPE ucov[],FTYPE ncov[],FTYPE bcon[]); 
extern FTYPE pressure_rho0_u(FTYPE rho0, FTYPE u);
extern FTYPE pressure_rho0_w(FTYPE rho0, FTYPE w);
extern FTYPE pressure_rho0_w(FTYPE rho0, FTYPE w);

extern int gamma_calc_g(FTYPE *pr, FTYPE gcov[NDIM][NDIM], FTYPE *gamma);
