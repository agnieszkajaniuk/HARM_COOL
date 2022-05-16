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

#include <unistd.h>
#include <sys/stat.h>


/* restart functions; restart_init and restart_dump */

#include "decs.h"

#if( COOL_EOS )
#include "cooleos.h"
#endif

#if TRACERSDUMP
#include "tracers.h"
#endif

extern double tdump,timage,tlog ;
extern int nfailed;

/***********************************************************************/
/***********************************************************************
  restart_write():
     -- writes current state of primitive variables to the 
        checkpointing/restart file. 
     -- uses ASCII text format ;
     -- when changing this routine, be sure to make analogous changes 
        in restart_read();
************************************************************************/
void restart_write()
{
  double (*   p)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_p);
  FILE *fp ;
  int h,i,j,k ;
  char filename[256];
  char filename2[256];

/*
  if(rdump_cnt%2 == 0) {
    fp = fopen("dumps/rdump0","wt") ;
    fprintf(stderr,"RESTART  file=dumps/rdump0\n") ;
  }
  else {
    fp = fopen("dumps/rdump1","wt") ;
    fprintf(stderr,"RESTART file=dumps/rdump1\n") ;
  }
*/
  sprintf(filename, "dumps/rdump0.%d_%dx%dx%d", rank, hsize,isize,jsize);
  fp = fopen(filename,"wt") ;
   fprintf(stderr,"RESTART  file=%s\n",filename) ;

  if(fp == NULL) {
    fprintf(stderr,"Cannot open restart file\n") ;
    exit(2) ;
  }

  /*************************************************************
	  Write the header of the restart file: 
  *************************************************************/
  fprintf(fp, FMT_INT_OUT, N1       );
  fprintf(fp, FMT_INT_OUT, N2       );
  fprintf(fp, FMT_INT_OUT, N3       );
  fprintf(fp, FMT_DBL_OUT, t        );
  fprintf(fp, FMT_DBL_OUT, tf       );
  fprintf(fp, FMT_INT_OUT, nstep    );
  fprintf(fp, FMT_DBL_OUT, a        );
  fprintf(fp, FMT_DBL_OUT, gam      );
  fprintf(fp, FMT_DBL_OUT, cour     );
  fprintf(fp, FMT_DBL_OUT, DTd      );
  fprintf(fp, FMT_DBL_OUT, DTl      );
  fprintf(fp, FMT_DBL_OUT, DTi      );
  fprintf(fp, FMT_INT_OUT, DTr      );
  fprintf(fp, FMT_INT_OUT, dump_cnt );
  fprintf(fp, FMT_INT_OUT, image_cnt);
  fprintf(fp, FMT_INT_OUT, rdump_cnt);
  fprintf(fp, FMT_DBL_OUT, dt       );
  fprintf(fp, FMT_INT_OUT, lim      );
  fprintf(fp, FMT_INT_OUT, failed   );
  fprintf(fp, FMT_DBL_OUT, Rin      );
  fprintf(fp, FMT_DBL_OUT, Rout     );
  fprintf(fp, FMT_DBL_OUT, hslope   );
  fprintf(fp, FMT_DBL_OUT, R0       );
  fprintf(fp, FMT_DBL_OUT, mass_start);
  fprintf(fp, FMT_DBL_OUT, mass_in  );
  fprintf(fp, FMT_DBL_OUT, mass_out );
  fprintf(fp, FMT_DBL_OUT, tdump    );
  fprintf(fp, FMT_DBL_OUT, timage   );
  fprintf(fp, FMT_DBL_OUT, tlog     );
  fprintf(fp, FMT_DBL_OUT, defcon   );
  fprintf(fp, FMT_INT_OUT, nfailed  );
  fprintf(fp, FMT_DBL_OUT, Rbr      );
  fprintf(fp, FMT_DBL_OUT, cyl_R0   );
  fprintf(fp, FMT_DBL_OUT, cyl_Th0  );

  fprintf(fp,"\n");

  /*************************************************************
	  Write the body of the restart file: 
  *************************************************************/
//  ZSLOOP(-2,N1+1,-2,N2+1) {
  for(h=0;h<SN3;++h)for(i=0;i<SN1;++i)for(j=0;j<SN2;++j) {
    PLOOP fprintf(fp, FMT_DBL_OUT, p[h][i][j][k]); 
    fprintf(fp, "\n"); 
  }

  //Store ener.out file size
  if(rank==0)
  { 
	struct stat st;
	stat("ener.out", &st);
	int size = st.st_size;

    fprintf(fp, FMT_INT_OUT, size);
  }

#if TRACERSDUMP
	restart_write_tracers(fp);
#endif


#if( COOL_EOS )
	restart_write_CoolEOS(fp);
#endif


  fclose(fp) ;

#ifdef MPI_USED
		MPI_Barrier(MPI_COMM_WORLD);
#endif

  rdump_cnt++ ;

  sprintf(filename2, "dumps/rdump1.%d_%dx%dx%d", rank, hsize,isize,jsize);
  remove(filename2);

#ifdef MPI_USED
		MPI_Barrier(MPI_COMM_WORLD);
#endif

  rename(filename,filename2);
  return;
}

/***********************************************************************/
/***********************************************************************
  restart_init():
     -- main driver for setting initial conditions from a checkpoint 
        or restart file. 
     -- determines if there are any restart files to use and then 
         lets the  user choose if there are more than one file. 
     -- then calls initializes run with restart data;
************************************************************************/
int restart_init()
{
  FILE *fp, *fp1, *fp0 ;
//  char ans[100] ;
  //int strncmp(char *s1, char *s2, int n) ;
  int h,i,j,k ;
  char filename[256];

  /********************************************************************
   Check to see which restart files exist. 
   Use the only one that exists, else prompt user to decide 
     which one to use if we have a choice : 
  ********************************************************************/
  sprintf(filename, "dumps/rdump0.%d_%dx%dx%d", rank, hsize,isize,jsize);
  fp0 = fopen(filename,"r") ;
  sprintf(filename, "dumps/rdump1.%d_%dx%dx%d", rank, hsize,isize,jsize);
  fp1 = fopen(filename,"r") ;

  if( (fp0 == NULL) && (fp1 == NULL) ) {

    fprintf(stderr,"No restart file\n") ;
    return(0) ;
  }
  else { 
    fprintf(stderr,"\nRestart file exists! \n") ;
    if( fp0 == NULL ) { 
      fprintf(stderr,"Using dumps/rdump1 ... \n");
      fp = fp1;
    }
    else if( fp1 == NULL ) { 
      fprintf(stderr,"Using dumps/rdump0 ... \n");
      fp = fp0;
    }
    else { 
/*
      if(rank==0) fprintf(stderr,"Use dumps/rdump0 (0) or dumps/rdump1 (1)?   [0|1]  \n");
      fscanf(stdin,"%s",ans) ;
      if(strncmp(ans,"0",1) == 0) { 
	fp = fp0;
      }
      else{ 
	fp = fp1;
      }
*/
      fprintf(stdout,"Invalid dump files (two rdump0 and aand rdump1) for rank %d, please fix... \n", rank);
      fprintf(stdout,"For e.g remove all rdump0 files in case when all of rdump1 exist\n");
      fprintf(stdout,"Otherwise remove all rdump1 files\n");

	  exit(0);
    }
  }

  /* set up global arrays */
  set_arrays() ;


  /********************************************************************
   Now that we know we are restarting from a checkpoint file, then 
   we need to read in data, assign grid functions and define the grid: 
  ********************************************************************/
  restart_read(fp) ;

#if TRACERSDUMP
	restart_read_tracers(fp);
#endif

  double (*   p)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_p);
  double (*  ph)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_ph);

  /* set half-step primitives everywhere */
//  ZSLOOP(-2,N1+1,-2,N2+1)
//  ZSLOOP(-1,1,-2,2,-2,2) PLOOP ph[h][i][j][k] = p[h][i][j][k] ;
  for(h=0;h<SN3;++h)for(i=0;i<SN1;++i)for(j=0;j<SN2;++j) PLOOP ph[h][i][j][k] = p[h][i][j][k] ;

  /* set metric functions */
  set_grid() ;

#if( COOL_EOS )
  restart_read_CoolEOS(fp);
#endif

  fclose(fp) ;

  /* bound */
  bound_prim(_p) ;
  bound_prim(_ph) ;


#if( DO_FONT_FIX ) 
  set_Katm();
#endif 

  /***********************************************************************
    Make any changes to parameters in restart file  here: 
      e.g., cour = 0.4 , change in limiter...
  ************************************************************************/
  //lim = MC ;
  //cour = 0.9 ;
  //lim = VANL ;
  //tf = 4000. ;

  fprintf(stderr,"done with restart init.\n") ;

  /* done! */
  return(1) ;

}

/***********************************************************************/
/***********************************************************************
  restart_read():
     -- reads in data from the restart file, which is specified in 
         restart_init() but is usually named "dumps/rdump[0,1]" 
************************************************************************/
void restart_read(FILE *fp)
{
  double (*   p)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_p);
  int idum,h,i,j,k ;

  /*************************************************************
	  READ the header of the restart file: 
  *************************************************************/
  fscanf(fp, "%d", &idum  );
  if(idum != N1) {
    fprintf(stderr,"error reading restart file; N1 differs\n") ;
    exit(3) ;
  }
  fscanf(fp, "%d", &idum  );
  if(idum != N2) {
    fprintf(stderr,"error reading restart file; N2 differs\n") ;
    exit(4) ;
  }
  fscanf(fp, "%d", &idum  );
  if(idum != N3) {
    fprintf(stderr,"error reading restart file; N3 differs\n") ;
    exit(4) ;
  }

  fscanf(fp, "%lf", &t        );
  fscanf(fp, "%lf", &tf       );
  fscanf(fp, "%d",  &nstep    );
  fscanf(fp, "%lf", &a        );
  fscanf(fp, "%lf", &gam      );
  fscanf(fp, "%lf", &cour     );
  fscanf(fp, "%lf", &DTd      );
  fscanf(fp, "%lf", &DTl      );
  fscanf(fp, "%lf", &DTi      );
  fscanf(fp, "%d",  &DTr      );
  fscanf(fp, "%d",  &dump_cnt );
  fscanf(fp, "%d",  &image_cnt);
  fscanf(fp, "%d",  &rdump_cnt);
  fscanf(fp, "%lf", &dt       );
  fscanf(fp, "%d",  &lim      );
  fscanf(fp, "%d",  &failed   );
  fscanf(fp, "%lf", &Rin      );
  fscanf(fp, "%lf", &Rout     );
  fscanf(fp, "%lf", &hslope   );
  fscanf(fp, "%lf", &R0       );
  fscanf(fp, "%lf", &mass_start);
  fscanf(fp, "%lf", &mass_in  );
  fscanf(fp, "%lf", &mass_out );
  fscanf(fp, "%lf", &tdump    );
  fscanf(fp, "%lf", &timage   );
  fscanf(fp, "%lf", &tlog     );
  fscanf(fp, "%lf", &defcon   );
  fscanf(fp, "%d",  &nfailed  );
  fscanf(fp, "%lf",  &Rbr      );
  fscanf(fp, "%lf",  &cyl_R0  );
  fscanf(fp, "%lf",  &cyl_Th0  );

  /*************************************************************
	  READ the body of the restart file: 
  *************************************************************/
//  ZSLOOP(-2,N1+1,-2,N2+1)
//  ZSLOOP(-1,1,-2,2,-2,2)  PLOOP fscanf(fp, "%lf", &(p[h][i][j][k])); 
  for(h=0;h<SN3;++h)for(i=0;i<SN1;++i)for(j=0;j<SN2;++j) PLOOP fscanf(fp, "%lf", &(p[h][i][j][k]));

  //Restore ener.out file size
  if(rank==0)
  { 
	int size;
    fscanf(fp, "%d", &size);
    truncate("ener.out", size);
  }

  return ;
}

#undef FMT_DBL_OUT
#undef FMT_INT_OUT
