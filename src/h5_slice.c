#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "hdf5.h"

#define DATASETNAME "IntArray"
#define NX_SUB  3           /* hyperslab dimensions */
#define NY_SUB  4
#define NX 7           /* output buffer dimensions */
#define NY 7
#define NZ  3
#define RANK         2
#define RANK_OUT     3



#define FMT_DBL_OUT "%28.18e"
#define FMT_INT_OUT "%10d"

#define NDIM       (4)        /* number of total dimensions.  Never changes */



const int RR = 0; 
const int TH = 1; 
const int PH = 2; 
int          N1, N2, N3;
double startx[3], dx[3];

double t;
double tf;
int nstep;
double  a;
double  gam;
double  cour;
double  DTd;
double  DTl;
double  DTi;
int DTr;
int dump_cnt;
int image_cnt;
int rdump_cnt;
double dt;
int lim;
int failed;
double Rin;
double Rout;
double hslope;
double R0;

double  *r_1D;
double  *th_1D;
double  *ph_1D;

double  *X1;
double  *X2;
double  *X3;

double  *RHO, *UU, *GEOMG, *DIVB;
double  *QNU, *P;
double  *U1U2U3, *B1B2B3;

double  *UCON, *UCOV, *BCON, *BCOV;

int
main (int argc, char *argv[])
{
    hid_t       file, dataset;         /* handles */
    hid_t       dataspace;
    hid_t       memspace;
    herr_t		status;

    hsize_t     dimsm[4];              /* memory space dimensions */

    hsize_t      count[4];              /* size of the hyperslab in the file */
    hsize_t      offset[4];             /* hyperslab offset in the file */
    hsize_t      count_out[3];          /* size of the hyperslab in memory */
    hsize_t      offset_out[3];         /* hyperslab offset in memory */

	char         h5_file_name[256];
	char         h5_file_dir[256];
	char         coords_file_name[256];
	char         output_file_name[256];
	int          slicePHI = 1;
	int          sliceNum = 0;
	int			 coolEnabled = 1;


	if(argc < 2)
	{
		printf("Program usage: h5_slice file_name [PHI,THETA] number [COOL,NOCOOL]\n");
		printf("example: h5_slice dump0000.h5 PHI 0 COOL\n\n");
		exit(1);
	}

	snprintf(h5_file_name, sizeof(h5_file_name), "%s", argv[1]);
	if(access(h5_file_name, 0))
	{
		printf("File %s does not exist!\n", h5_file_name);
		exit(2);
	}

	{
		strcpy(h5_file_dir, h5_file_name);
		char *p = strrchr(h5_file_dir, '/');
		if(p != NULL)
			*(p+1) = '\0';
		else
			h5_file_dir[0] = '\0';
	}

	sprintf(coords_file_name, "%scoords.h5", h5_file_dir);

	if(access(coords_file_name, 0))
	{
		printf("File %s does not exist!\n", coords_file_name);
		exit(2);
	}


	if(argc > 2 && strcmp(argv[2], "THETA") == 0)
		slicePHI = 0;
	
	if(argc > 3)
		sliceNum = atoi(argv[3]);

	if(argc > 4 && strcmp(argv[4], "NOCOOL") == 0)
		coolEnabled = 0;

	{
		char *p = strrchr(h5_file_name, '/');
		if(p != NULL)
			strcpy(output_file_name, p + 1);
		else
			strcpy(output_file_name, h5_file_name);

		p = strrchr(output_file_name, '.');
		if(p != NULL)
			*(p+1) = '\0';
		
		sprintf(output_file_name+strlen(output_file_name), "%s#%d", slicePHI ? "PHI" : "THETA", sliceNum);
	}


    /*
     * Open the file and the dataset.
     */
    file = H5Fopen(coords_file_name, H5F_ACC_RDONLY, H5P_DEFAULT);

    {
        const char *datasetname[] = {"/N1", "/N2", "/N3", "/lim", "/DTr"};
        int   *datasetval[]  = {&N1, &N2, &N3, &lim, &DTr};

        for(unsigned int i = 0; i < sizeof(datasetname)/sizeof(datasetname[0]); i++)
        {
		    dataset = H5Dopen2(file, datasetname[i], H5P_DEFAULT);
			status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, datasetval[i]);
		    H5Dclose(dataset);
        }
	}

	printf("DadaGrid: N1xN2xN3=%dx%dx%d\n", N1, N2, N3);

	if(slicePHI && (sliceNum < 0 || sliceNum >= N3))
	{
		printf("Invalid slice num %d for slicePHI\n", sliceNum);
		exit(4);
	}
	if(!slicePHI && (sliceNum < 0 || sliceNum >= N2))
	{
		printf("Invalid slice num %d for sliceTHETA\n", sliceNum);
		exit(4);
	}
	


    {
        const char *datasetname[] = {"/startx1", "/startx2", "/startx3", "/dx1", "/dx2", "/dx3", "/a", "/Rin", "/Rout", "/hslope",
									"/gam", "/cour", "/tf", "/DTd", "/DTl", "/DTi", "/R0"};

        double *datasetval[]  = {&startx[RR], &startx[TH], &startx[PH], &dx[RR], &dx[TH], &dx[PH], &a, &Rin, &Rout, &hslope,
										&gam, &cour, &tf, &DTd, &DTl, &DTi, &R0};

        for(unsigned int i = 0; i < sizeof(datasetname)/sizeof(datasetname[0]); i++)
        {
		    dataset = H5Dopen2(file, datasetname[i], H5P_DEFAULT);
			status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, datasetval[i]);
		    H5Dclose(dataset);
        }
	}

	r_1D = (double*)malloc(N1*sizeof(double));
	th_1D = (double*)malloc(N2*sizeof(double));
	ph_1D = (double*)malloc(N3*sizeof(double));

    {
        const char *datasetname[] = {"/r_1D", "/th_1D", "/ph_1D"};
        double *datasetval[]  = {r_1D, th_1D, ph_1D};

        for(unsigned int i = 0; i < sizeof(datasetname)/sizeof(datasetname[0]); i++)
        {
		    dataset = H5Dopen2(file, datasetname[i], H5P_DEFAULT);
			status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, datasetval[i]);
		    H5Dclose(dataset);
        }
	}

	X1 = (double*)malloc(N1*sizeof(double));
	X2 = (double*)malloc(N2*sizeof(double));
	X3 = (double*)malloc(N3*sizeof(double));

    {
        const char *datasetname[] = {"/X1", "/X2", "/X3"};
        double *datasetval[]  = {X1, X2, X3};

        for(unsigned int i = 0; i < sizeof(datasetname)/sizeof(datasetname[0]); i++)
        {
		    dataset = H5Dopen2(file, datasetname[i], H5P_DEFAULT);
			status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, datasetval[i]);
		    H5Dclose(dataset);
        }
	}



    /*
     * Define the memory dataspace.
     */
    dimsm[0] = N1;
    dimsm[1] = slicePHI ? N2 : N3;
    memspace = H5Screate_simple(2,dimsm,NULL);

    /*
     * Define memory hyperslab.
     */
    offset_out[0] = 0;
    offset_out[1] = 0;
    count_out[0]  = N1;
    count_out[1]  = slicePHI ? N2 : N3;
    status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);


    /*
     * Define file hyperslab.
     */
    offset[0] = 0;
    offset[1] = !slicePHI ? sliceNum : 0;
    offset[2] = slicePHI ? sliceNum : 0;
    count[0]  = N1;
    count[1]  = slicePHI ? N2 : 1;
    count[2]  = !slicePHI ? N3 : 1;

	GEOMG = (double*)malloc(N1*(slicePHI ? N2 : N3)*sizeof(double));

	{

   		dataset = H5Dopen2(file, "/gdet", H5P_DEFAULT);

   		dataspace = H5Dget_space(dataset);    /* dataspace handle */

   		status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

    	/*
   	 	* Read data from hyperslab in the file into the hyperslab in
   		* memory and display.
   		*/
   		status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, GEOMG);

   		H5Sclose(dataspace);
   		H5Dclose(dataset);
	}

    H5Fclose(file);







    /*
     * Open the file and the dataset.
     */
    file = H5Fopen(h5_file_name, H5F_ACC_RDONLY, H5P_DEFAULT);

    {
        const char *datasetname[] = {"/nstep", "/failed", "/dump_cnt", "/image_cnt", "/rdump_cnt"};
        int   *datasetval[]  = {&nstep, &failed, &dump_cnt, &image_cnt, &rdump_cnt};

        for(unsigned int i = 0; i < sizeof(datasetname)/sizeof(datasetname[0]); i++)
        {
		    dataset = H5Dopen2(file, datasetname[i], H5P_DEFAULT);
			status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, datasetval[i]);
		    H5Dclose(dataset);
        }
	}

    {
        const char *datasetname[] = {"/t", "/dt"};
        double *datasetval[]  = {&t, &dt};

        for(unsigned int i = 0; i < sizeof(datasetname)/sizeof(datasetname[0]); i++)
        {
		    dataset = H5Dopen2(file, datasetname[i], H5P_DEFAULT);
			status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, datasetval[i]);
		    H5Dclose(dataset);
        }
	}




	{
        const char *datasetname[] = {"/Rho", "/Energy", "/DivB", "/Qnu", "/P"};
        double **datasetval[]  = {&RHO, &UU, &DIVB, &QNU, &P};

        for(unsigned int i = 0; i < sizeof(datasetname)/sizeof(datasetname[0]); i++)
		{

			if(H5Oexists_by_name(file, datasetname[i], H5P_DEFAULT) < 0)
			{
				printf("dataset %s not present\n", datasetname[i]);
				continue;
			}

    		dataset = H5Dopen2(file, datasetname[i], H5P_DEFAULT);
			if(dataset < 0)
			{
				printf("dataset %s not present\n", datasetname[i]);
				continue;
			}

    		dataspace = H5Dget_space(dataset);    /* dataspace handle */

    		status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

			*datasetval[i] = (double*)malloc(N1*(slicePHI ? N2 : N3)*sizeof(double));

	    	/*
    	 	* Read data from hyperslab in the file into the hyperslab in
     		* memory and display.
     		*/
    		status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, *datasetval[i]);

    		H5Sclose(dataspace);
    		H5Dclose(dataset);
		}
    }

    H5Sclose(memspace);


    /*
     * Define the memory dataspace.
     */
    dimsm[0] = N1;
    dimsm[1] = slicePHI ? N2 : N3;
    dimsm[2] = 3;
    memspace = H5Screate_simple(3,dimsm,NULL);

    /*
     * Define memory hyperslab.
     */
    offset_out[0] = 0;
    offset_out[1] = 0;
    offset_out[2] = 0;
    count_out[0]  = N1;
    count_out[1]  = slicePHI ? N2 : N3;
    count_out[2]  = 3;
    status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);


    /*
     * Define file hyperslab.
     */
    offset[0] = 0;
    offset[1] = !slicePHI ? sliceNum : 0;
    offset[2] = slicePHI ? sliceNum : 0;
    offset[3] = 0;
    count[0]  = N1;
    count[1]  = slicePHI ? N2 : 1;
    count[2]  = !slicePHI ? N3 : 1;
    count[3]  = 3;

	{
        const char *datasetname[] = {"/U1U2U3", "/B1B2B3"};
        double **datasetval[]  = {&U1U2U3, &B1B2B3};

        for(unsigned int i = 0; i < sizeof(datasetname)/sizeof(datasetname[0]); i++)
		{
    		dataset = H5Dopen2(file, datasetname[i], H5P_DEFAULT);

    		dataspace = H5Dget_space(dataset);    /* dataspace handle */

    		status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

			*datasetval[i] = (double*)malloc(3*N1*(slicePHI ? N2 : N3)*sizeof(double));

	    	/*
    	 	* Read data from hyperslab in the file into the hyperslab in
     		* memory and display.
     		*/
    		status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, *datasetval[i]);

    		H5Sclose(dataspace);
    		H5Dclose(dataset);
		}
    }

    H5Sclose(memspace);


    /*
     * Define the memory dataspace.
     */
    dimsm[0] = N1;
    dimsm[1] = slicePHI ? N2 : N3;
    dimsm[2] = NDIM;
    memspace = H5Screate_simple(3,dimsm,NULL);

    /*
     * Define memory hyperslab.
     */
    offset_out[0] = 0;
    offset_out[1] = 0;
    offset_out[2] = 0;
    count_out[0]  = N1;
    count_out[1]  = slicePHI ? N2 : N3;
    count_out[2]  = NDIM;
    status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);


    /*
     * Define file hyperslab.
     */
    offset[0] = 0;
    offset[1] = !slicePHI ? sliceNum : 0;
    offset[2] = slicePHI ? sliceNum : 0;
    offset[3] = 0;
    count[0]  = N1;
    count[1]  = slicePHI ? N2 : 1;
    count[2]  = !slicePHI ? N3 : 1;
    count[3]  = NDIM;


	{
        const char *datasetname[] = {"/ucon", "/ucov", "/bcon", "/bcov"};
        double **datasetval[]  = {&UCON, &UCOV, &BCON, &BCOV};

        for(unsigned int i = 0; i < sizeof(datasetname)/sizeof(datasetname[0]); i++)
		{
    		dataset = H5Dopen2(file, datasetname[i], H5P_DEFAULT);

    		dataspace = H5Dget_space(dataset);    /* dataspace handle */

    		status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

			*datasetval[i] = (double*)malloc(NDIM*N1*(slicePHI ? N2 : N3)*sizeof(double));

	    	/*
    	 	* Read data from hyperslab in the file into the hyperslab in
     		* memory and display.
     		*/
    		status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, *datasetval[i]);

    		H5Sclose(dataspace);
    		H5Dclose(dataset);
		}
    }

    H5Fclose(file);


	//Write output file	
	{
		remove(output_file_name);
		FILE *fp = fopen(output_file_name, "w+");
		if(fp==NULL) 
		{
			printf("error opening dump file %s\n", output_file_name) ;
			exit(3) ;
		}
		

        /***************************************************************
          Write header information :
        ***************************************************************/
        fprintf(fp, FMT_DBL_OUT, t        );
        fprintf(fp, FMT_INT_OUT, N1       );
		if(slicePHI)
        	fprintf(fp, FMT_INT_OUT, N2       );
		if(!slicePHI)
        	fprintf(fp, FMT_INT_OUT, N3       );
        fprintf(fp, FMT_DBL_OUT, startx[RR]);

		if(slicePHI)
	        fprintf(fp, FMT_DBL_OUT, startx[TH]);
		if(!slicePHI)
        	fprintf(fp, FMT_DBL_OUT, startx[PH]);
        fprintf(fp, FMT_DBL_OUT, dx[RR]    );
		if(slicePHI)
	        fprintf(fp, FMT_DBL_OUT, dx[TH]    );
		if(!slicePHI)
        	fprintf(fp, FMT_DBL_OUT, dx[PH]    );

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

        fprintf(fp,"\n") ;

        /***************************************************************
          Write header information :
        ***************************************************************/
    	int jjmax = slicePHI ? N2 : N3;
	
		for(int ii=0; ii < N1; ii++)
		for(int jj=0; jj < jjmax; jj++)
		{
            fprintf(fp, FMT_DBL_OUT, X1[ii]       );

			if(slicePHI)
    	        fprintf(fp, FMT_DBL_OUT, X2[jj]       );

			if(!slicePHI)
	            fprintf(fp, FMT_DBL_OUT, X3[jj]       );

            fprintf(fp, FMT_DBL_OUT, r_1D[ii]          );

			if(slicePHI)
	            fprintf(fp, FMT_DBL_OUT, th_1D[jj]         );

			if(!slicePHI)
	            fprintf(fp, FMT_DBL_OUT, ph_1D[jj]         );

            fprintf(fp, FMT_DBL_OUT, RHO[ii*jjmax + jj]          );
            fprintf(fp, FMT_DBL_OUT,  UU[ii*jjmax + jj]          );

			for(int k=0; k<3; k++) 
	            fprintf(fp, FMT_DBL_OUT,  U1U2U3[3*(ii*jjmax + jj) + k]          );

			for(int k=0; k<3; k++) 
	            fprintf(fp, FMT_DBL_OUT,  B1B2B3[3*(ii*jjmax + jj) + k]          );

			if(coolEnabled)
			{
				if(QNU != NULL)
		            fprintf(fp, FMT_DBL_OUT,  QNU[ii*jjmax + jj]          );
				else
		            fprintf(fp, FMT_DBL_OUT,  0.0          );

				if(P != NULL)
		            fprintf(fp, FMT_DBL_OUT,  P[ii*jjmax + jj]          );
				else
		            fprintf(fp, FMT_DBL_OUT,  0.0          );
	            fprintf(fp, FMT_DBL_OUT,  0.0          );
	            fprintf(fp, FMT_DBL_OUT,  0.0          );
			}

            fprintf(fp, FMT_DBL_OUT,  DIVB[ii*jjmax + jj]          );

            if(!failed) 
			{
				for(int k=0; k<NDIM; k++) 
		            fprintf(fp, FMT_DBL_OUT,  UCON[NDIM*(ii*jjmax + jj) + k]          );

				for(int k=0; k<NDIM; k++) 
		            fprintf(fp, FMT_DBL_OUT,  UCOV[NDIM*(ii*jjmax + jj) + k]          );
				for(int k=0; k<NDIM; k++) 
		            fprintf(fp, FMT_DBL_OUT,  BCON[NDIM*(ii*jjmax + jj) + k]          );
				for(int k=0; k<NDIM; k++) 
		            fprintf(fp, FMT_DBL_OUT,  BCOV[NDIM*(ii*jjmax + jj) + k]          );

				//vmin, vmax
            	fprintf(fp, FMT_DBL_OUT, 0.0);
            	fprintf(fp, FMT_DBL_OUT, 0.0);
        	    fprintf(fp, FMT_DBL_OUT, 0.0 );
    	        fprintf(fp, FMT_DBL_OUT, 0.0 );

	            fprintf(fp, FMT_DBL_OUT, GEOMG[ii*jjmax + jj] );

			}

            fprintf(fp,"\n") ;

		}

		fclose(fp);
	}

    status++;

    return 0;
}
