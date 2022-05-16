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

#include "decs.h"

#if( COOL_EOS )
#include "cooleos.h"
#endif

#if HDFDUMP
#include "hdf5.h"
#endif


#if HDFDUMP
int WriteCoordsHDF(const char *FileName);
int WriteBH_HorizonHDF(const char *fullFileName);
void write_xdmf_xml(double ttime, int hdfNo);
void write_BH_Horizon_xdmf_xml(void);

#if N3<=2
#define DUMPZLOOP for(i=istart;i<=istop;i++)for(j=jstart;j<=jstop;j++)for(h=hstart;h<=hstop;h++)
#else
#define DUMPZLOOP for(i=istart;i<=istop;i++)for(j=jstart;j<=jstop;j++)for(h=hstart;h<=(hrank==hsize-1?hstop+1:hstop);h++)
#endif


void dumph5(int dump_cnt)
{
    double (*   p)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_p);
    double (*   gdet)[SN2][NPG] = (double (*) [SN2][NPG])(_gdet);
#if( COOL_EOS )
    double (*   CoolVal)[SN1][SN2][NCOOLVAL] = (double (*) [SN1][SN2][NCOOLVAL])(_CoolVal);
#endif

    /*
     * HDF5 APIs definitions
     */
    hid_t       file_id, dset_id;         /* file and dataset identifiers */
    hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
    hsize_t     dimsf[4];                 /* dataset dimensions */
    double     *data;                    /* pointer to data buffer to write */
    hsize_t		count[4];	          /* hyperslab selection parameters */
    hsize_t		offset[4];
    hid_t		plist_id;                 /* property list identifier */
    int         h,i,j,k;
    herr_t		status;
    char 	    fullFileName[256];
//	double 		r,th,ph,X[NDIM] ;


    if(dump_cnt == 0)
    {
        WriteCoordsHDF("dumps/coords.h5");
		WriteBH_HorizonHDF("dumps/bh_horizon.h5");
    }

    sprintf(fullFileName,"dumps/dump%03d.h5",dump_cnt) ;

    if(rank == 0)
        fprintf(stderr,"DUMP     file=%s\n",fullFileName) ;

    if(rank == 0)
        write_xdmf_xml(t, dump_cnt);


#if	N3>2
	{
		MPI_Sendrecv(p, 1, slice_hNPR1_start, rank_prevh, 125,
                p, 1, slice_hNPR1_stopN, rank_nexth, 125,
                MPI_COMM_WORLD, &MPIStatus);
	}
#endif

    /*
     * Set up file access property list with parallel I/O access
     */
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    if(plist_id < 0) printf("HDF Error 1\n");

    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    /*
     * Create a new file collectively and release property list identifier.
    */

    file_id = H5Fcreate(fullFileName, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    if(file_id < 0) printf("HDF Error 2\n");

    H5Pclose(plist_id);

    /*
     * Create the dataspace for the dataset.
     */

    dimsf[0] = N1;
    dimsf[1] = N2;
#if N3<=2
    dimsf[2] = N3;
#else
    dimsf[2] = N3+1;
#endif
    dimsf[3] = 4;

    /*
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
    count[0] = (istop - istart + 1);
    count[1] =  (jstop - jstart + 1);
#if N3<=2
    count[2] =  (hstop - hstart + 1);
#else
    count[2] =  (hstop - hstart + 1)+(hrank==hsize-1?1:0);
#endif
    count[3] =  4;

    offset[0] = istart - (irank==0 ? 2 : 1) + di_rank_start;
    offset[1] = jstart - (jrank==0 ? 2 : 1) + dj_rank_start;
    offset[2] = (N3 == 1 ? 0 : hstart - 1 + dh_rank_start);
    offset[3] = 0;


#ifdef MPI_USED
	if(N3>2)
	{
		MPI_Sendrecv(p, 1, slice_hNPR1_stop, rank_nexth, 124,
                p, 1, slice_hNPR1_startP, rank_prevh, 124,
                MPI_COMM_WORLD, &MPIStatus);
	}
#endif

    /*
    	 * Initialize data buffer
    	 */
    data = (double *) malloc(count[0] * count[1] * count[2] * count[3] * sizeof(double));
	if(data == NULL)
	{
		fprintf(stderr,"can not allocate %lld bytes rank=%d\n", count[0] * count[1] * count[2] * count[3] * sizeof(double), rank);
		exit(-1);
	}

    //Scalar primitives
    for(k=0; k<U1; k++)
    {
        const char *primitives[2] = { "/Rho", "/Energy" };

#if !HDFDUMP_Rho
		if(k==0) continue;
#endif

#if !HDFDUMP_Energy
		if(k==1) continue;
#endif


        filespace = H5Screate_simple(3, dimsf, NULL);

        /*
        * Create the dataset with default properties and close filespace.
        */
        dset_id = H5Dcreate(file_id, primitives[k], H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);

        memspace = H5Screate_simple(3, count, NULL);

        /*
         * Select hyperslab in the file.
         */
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

        int ii = 0;
        DUMPZLOOP data[ii++] = p[h][i][j][k];

        /*
         * Create property list for collective dataset write.
         */
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace,
                          plist_id, data);
        H5Pclose(plist_id);

        /*
         * Close/release resources.
         */
        H5Dclose(dset_id);
        H5Sclose(filespace);
        H5Sclose(memspace);
    }

#if HDFDUMP_Temp
    {
        const char *primitives[1] = { "/Temp" };

        filespace = H5Screate_simple(3, dimsf, NULL);

        /*
        * Create the dataset with default properties and close filespace.
        */
        dset_id = H5Dcreate(file_id, primitives[0], H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);

        memspace = H5Screate_simple(3, count, NULL);

        /*
         * Select hyperslab in the file.
         */
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

	    const double CL = 2.99792458e10;
    	const double KBOL = 1.38e-16;
    	const double mn = 1.67492729e-24;
    	const double me = 9.109565e-28;
    	const double mp = 1.673e-24;
    	const double mhe = 6.6465e-24; // 2p+2e+binding energy
    	const double mav = (mp+mn+2.*me+mhe)/5.; //avrage particle mass

        int ii = 0;
        DUMPZLOOP {
			double T = (gam-1.)*p[h][i][j][UU]/p[h][i][j][RHO]/KBOL*mav*CL*CL;
			data[ii++] = T;
		}
        /*
         * Create property list for collective dataset write.
         */
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace,
                          plist_id, data);
        H5Pclose(plist_id);

        /*
         * Close/release resources.
         */
        H5Dclose(dset_id);
        H5Sclose(filespace);
        H5Sclose(memspace);
    }
#endif


#if HDFDUMP_DivB
	//Div
	{

#ifdef MPI_USED
	MPI_Sendrecv(p, 1, slice_iNPR1_stop, rank_nexti, 124,
                p, 1, slice_iNPR1_startP, rank_previ, 124,
                MPI_COMM_WORLD, &MPIStatus);

	if(N3>1)
	{
		MPI_Sendrecv(p, 1, slice_hNPR1_stop, rank_nexth, 124,
                p, 1, slice_hNPR1_startP, rank_prevh, 124,
                MPI_COMM_WORLD, &MPIStatus);
	}
#endif

        int ii = 0;
        DUMPZLOOP
        {
#if N3 == 1
             data[ii++] = fabs( 0.5*(
                                 p[h][i][j][B1]*gdet[i][j][CENT]
                                 + p[h][i][j-1][B1]*gdet[i][j-1][CENT]
                                 - p[h][i-1][j][B1]*gdet[i-1][j][CENT]
                                 - p[h][i-1][j-1][B1]*gdet[i-1][j-1][CENT]
                             )/dx[RR] +
                             0.5*(
                                 p[h][i][j][B2]*gdet[i][j][CENT]
                                 + p[h][i-1][j][B2]*gdet[i-1][j][CENT]
                                 - p[h][i][j-1][B2]*gdet[i][j-1][CENT]
                                 - p[h][i-1][j-1][B2]*gdet[i-1][j-1][CENT]
                             )/dx[TH]) ;
#else
//			coord(h,i,j,CENT,X) ;
//			bl_coord(X,&r,&th,&ph) ;

			data[ii++] = fabs( 0.25*(
								+ p[h-1][i][j][B1]*gdet[i][j][CENT]
    	                        + p[h-1][i][j-1][B1]*gdet[i][j-1][CENT]
        	                    - p[h-1][i-1][j][B1]*gdet[i-1][j][CENT]
            	                - p[h-1][i-1][j-1][B1]*gdet[i-1][j-1][CENT]

                                 + p[h][i][j][B1]*gdet[i][j][CENT]
                                 + p[h][i][j-1][B1]*gdet[i][j-1][CENT]
                                 - p[h][i-1][j][B1]*gdet[i-1][j][CENT]
                                 - p[h][i-1][j-1][B1]*gdet[i-1][j-1][CENT]
                             )/dx[RR] +
                             0.25*(
                                 + p[h-1][i][j][B2]*gdet[i][j][CENT]
                                 + p[h-1][i-1][j][B2]*gdet[i-1][j][CENT]
                                 - p[h-1][i][j-1][B2]*gdet[i][j-1][CENT]
                                 - p[h-1][i-1][j-1][B2]*gdet[i-1][j-1][CENT]

                                 + p[h][i][j][B2]*gdet[i][j][CENT]
                                 + p[h][i-1][j][B2]*gdet[i-1][j][CENT]
                                 - p[h][i][j-1][B2]*gdet[i][j-1][CENT]
                                 - p[h][i-1][j-1][B2]*gdet[i-1][j-1][CENT]
                             )/dx[TH] +
                             0.25*(
                                 + p[h][i][j][B3]*gdet[i][j][CENT]
                                 + p[h][i-1][j][B3]*gdet[i-1][j][CENT]
                                 - p[h-1][i][j][B3]*gdet[i][j][CENT]
                                 - p[h-1][i-1][j][B3]*gdet[i-1][j][CENT]

                                 + p[h][i][j-1][B3]*gdet[i][j][CENT]
                                 + p[h][i-1][j-1][B3]*gdet[i-1][j-1][CENT]
                                 - p[h-1][i][j-1][B3]*gdet[i][j-1][CENT]
                                 - p[h-1][i-1][j-1][B3]*gdet[i-1][j-1][CENT]
                             )/dx[PH]) ;
#endif

#if( COOL_EOS )
	        InitCache_CoolEOS(h, i, j);

			CoolVal_rho0_u_CoolEOS(p[h][i][j][RHO], p[h][i][j][UU], h, i, j, CoolVal[h][i][j]);
#endif
		}

        filespace = H5Screate_simple(3, dimsf, NULL);

        /*
        * Create the dataset with default properties and close filespace.
        */
        dset_id = H5Dcreate(file_id, "/DivB", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);

        memspace = H5Screate_simple(3, count, NULL);

        /*
         * Select hyperslab in the file.
         */
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

        /*
         * Create property list for collective dataset write.
         */
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace,
                          plist_id, data);
        H5Pclose(plist_id);

        /*
         * Close/release resources.
         */
        H5Dclose(dset_id);
        H5Sclose(filespace);
        H5Sclose(memspace);

	}
#endif

    //Scalar primitives
#if( COOL_EOS )

    for(k=0; k<7; k++)
    {
        const char *primitives[7] = { "/Lambda_sim", "/P", "/Hd", "/Qnu", "/Tau", "/T", "/Yee"};

        filespace = H5Screate_simple(3, dimsf, NULL);

        /*
        * Create the dataset with default properties and close filespace.
        */
        dset_id = H5Dcreate(file_id, primitives[k], H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);

        memspace = H5Screate_simple(3, count, NULL);

        /*
         * Select hyperslab in the file.
         */
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

        int ii = 0;
        DUMPZLOOP data[ii++] = CoolVal[h][i][j][k];

        /*
         * Create property list for collective dataset write.
         */
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace,
                          plist_id, data);
        H5Pclose(plist_id);

        /*
         * Close/release resources.
         */
        H5Dclose(dset_id);
        H5Sclose(filespace);
        H5Sclose(memspace);
    }
#endif

    //vector primitives
    for(k=0; k<2; k++)
    {
        const char *primitives[2] = { "/U1U2U3gdet", "/B1B2B3gdet" };

#if !HDFDUMP_U1U2U3gdet
		if(k==0) continue;
#endif

#if !HDFDUMP_B1B2B3gdet
		if(k==1) continue;
#endif


        dimsf[3] = 3;
	    count[3] = 3;

        filespace = H5Screate_simple(4, dimsf, NULL);

        /*
        * Create the dataset with default properties and close filespace.
        */
        dset_id = H5Dcreate(file_id, primitives[k], H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);

        memspace = H5Screate_simple(4, count, NULL);

        /*
         * Select hyperslab in the file.
         */
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

        int ii = 0;
        DUMPZLOOP
        {
	        struct of_geom geom ;
            get_geometry(h,i,j,CENT,&geom) ;

            double *vv = &p[h][i][j][U1+k*3];
            data[ii++] = *(vv++)*geom.g;
            data[ii++] = *(vv++)*geom.g;
            data[ii++] = *(vv++)*geom.g;
        }

        /*
         * Create property list for collective dataset write.
         */
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace,
                          plist_id, data);
        H5Pclose(plist_id);

        /*
         * Close/release resources.
         */
        H5Dclose(dset_id);
        H5Sclose(filespace);
        H5Sclose(memspace);
    }


    if(!failed)
    {
        struct of_state *qstate;

        qstate = (of_state *) malloc(count[0] * count[1] * count[2] * sizeof(struct of_state));
		
		if(qstate == NULL)
		{
			fprintf(stderr,"can not allocate %lld bytes rank=%d\n", count[0] * count[1] * count[2] * sizeof(struct of_state), rank);
			exit(-2);
		}

        int ii = 0;
        DUMPZLOOP
        {
	        struct of_geom geom ;
            get_geometry(h,i,j,CENT,&geom) ;
            get_state(p[h][i][j],&geom,&qstate[ii]) ;
			++ii;
        }

        for(k=0; k<4; k++)
        {
            const char *primitives[4] = { "/ucon", "/ucov", "/bcon", "/bcov"};

#if !HDFDUMP_Ucon
		if(k==0) continue;
#endif

#if !HDFDUMP_Ucov
		if(k==1) continue;
#endif
#if !HDFDUMP_Bcon
		if(k==2) continue;
#endif

#if !HDFDUMP_Bcov
		if(k==3) continue;
#endif

            dimsf[3] = NDIM;
			count[3] = NDIM;

            filespace = H5Screate_simple(4, dimsf, NULL);

            /*
            * Create the dataset with default properties and close filespace.
            */
            dset_id = H5Dcreate(file_id, primitives[k], H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Sclose(filespace);

            memspace = H5Screate_simple(4, count, NULL);

            /*
             * Select hyperslab in the file.
             */
            filespace = H5Dget_space(dset_id);
            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

            int ii = 0;
            DUMPZLOOP
            {
                double *vv = (double *)&qstate[ii] + k*NDIM;
                data[ii*4] = *(vv++);
                data[ii*4+1] = *(vv++);
                data[ii*4+2] = *(vv++);
                data[ii*4+3] = *(vv++);
				++ii;
            }

            /*
             * Create property list for collective dataset write.
             */
            plist_id = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

            status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace,
                              plist_id, data);
            H5Pclose(plist_id);

            /*
             * Close/release resources.
             */
            H5Dclose(dset_id);
            H5Sclose(filespace);
            H5Sclose(memspace);
        }

        free(qstate);
    }


    free(data);

    H5Fclose(file_id);

#ifdef MPI_USED
    MPI_Barrier(MPI_COMM_WORLD);
#endif


    if(rank == 0)
    {

	    file_id = H5Fopen(fullFileName, H5F_ACC_RDWR, H5P_DEFAULT);

        /* Create the data space for the attribute. */
        dimsf[0] = 1;
        {
            const char *datasetname[] = {"/t", "/dt"};
            const double datasetval[]  = {t, dt};

            for(unsigned int i = 0; i < sizeof(datasetname)/sizeof(datasetname[0]); i++)
            {
		        filespace = H5Screate_simple(1, dimsf, NULL);
		        dset_id = H5Dcreate(file_id, datasetname[i], H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &datasetval[i]);
		        H5Sclose(filespace);
                H5Dclose(dset_id);
            }
        }

       {
            const char *datasetname[] = {"/nstep", "/failed", "/dump_cnt", "/image_cnt", "/rdump_cnt"};
            const int   datasetval[]  = {nstep, failed, dump_cnt, image_cnt, rdump_cnt};

            for(unsigned int i = 0; i < sizeof(datasetname)/sizeof(datasetname[0]); i++)
            {
		        filespace = H5Screate_simple(1, dimsf, NULL);
		        dset_id = H5Dcreate(file_id, datasetname[i], H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                status = H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &datasetval[i]);
		        H5Sclose(filespace);
                H5Dclose(dset_id);
            }
        }

	    H5Fclose(file_id);

    }

    status++;

    return;
}


/**************************************************************/
int WriteCoordsHDF(const char *fullFileName)
/**************************************************************/
{
    /*
     * HDF5 APIs definitions
     */
    hid_t       file_id, dset_id;         /* file and dataset identifiers */
    hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
    hsize_t     dimsf[5];                 /* dataset dimensions */
    double     *data;                    /* pointer to data buffer to write */
    hsize_t		count[5];	          /* hyperslab selection parameters */
    hsize_t		offset[5];
    hid_t		plist_id;                 /* property list identifier */
    int         h,i,j,k;
    herr_t		status;


    if(rank == 0)
        fprintf(stderr,"DUMP     file=%s\n", fullFileName) ;

    /*
     * Set up file access property list with parallel I/O access
     */
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    if(plist_id < 0) printf("HDF Error 1\n");

    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    /*
     * Create a new file collectively and release property list identifier.
    */

    file_id = H5Fcreate(fullFileName, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    if(file_id < 0) printf("HDF Error 2\n");

    H5Pclose(plist_id);

    /*
     * Create the dataspace for the dataset.
     */

    dimsf[0] = N1;
    dimsf[1] = N2;
#if N3<=2
    dimsf[2] = N3;
#else
    dimsf[2] = N3+1;
#endif
    dimsf[3] = NDIM;
    dimsf[4] = NDIM;

    /*
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
    count[0] = (istop - istart + 1);
    count[1] =  (jstop - jstart + 1);
#if N3<=2
    count[2] =  (hstop - hstart + 1);
#else
    count[2] =  (hstop - hstart + 1)+(hrank==hsize-1?1:0);
#endif
    count[3] = NDIM;
    count[4] = NDIM;

    offset[0] = istart - (irank==0 ? 2 : 1) + di_rank_start;
    offset[1] = jstart - (jrank==0 ? 2 : 1) + dj_rank_start;
    offset[2] = (N3 == 1 ? 0 : hstart - 1 + dh_rank_start);
    offset[3] =0;
    offset[4] =0;


    /*
    	 * Initialize data buffer
    	 */
    data = (double *) malloc(count[0]*count[1]*count[2] * sizeof(double));
	if(data == NULL)
	{
		fprintf(stderr,"can not allocate %lld bytes rank=%d\n", count[0] * count[1] * count[2] * sizeof(double), rank);
		exit(-1);
	}

    //Create 1D vectors of coordinates
    for(k = 0; k<3; k++)
    {
        const char *datasetname[3] = {"/r_1D", "/th_1D", "/ph_1D"};


        filespace = H5Screate_simple(1, &dimsf[k], NULL);

        /*
        * Create the dataset with default properties and close filespace.
        */
        dset_id = H5Dcreate(file_id, datasetname[k], H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);

        memspace = H5Screate_simple(1, &count[k], NULL);

        /*
         * Select hyperslab in the file.
         */
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset[k], NULL, &count[k], NULL);


        int ii = 0;
        double X[NDIM];
        double r, th, ph;


        if(k == 0)
        {
            ZSLOOPI(0,0)
            {
                coord(hstart, i, jstart, CENT, X);
                bl_coord(X, &r, &th, &ph);
                data[ii++] = r;
            }
        }
        else if(k == 1)
        {
            ZSLOOPJ(0,0)
            {
                coord(hstart, istart, j, CENT, X);
                bl_coord(X, &r, &th, &ph);
                data[ii++] = th;
            }
        }
        else if(k == 2)
        {
#if N3>2
            ZSLOOPH(0,1)
#else
            ZSLOOPH(0,0)
#endif
            {
                coord(h, istart, jstart, CENT, X);
                bl_coord(X, &r, &th, &ph);
                data[ii++] = ph;
            }
        }


        /*
         * Create property list for collective dataset write.
         */
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace,
                          plist_id, data);
        H5Pclose(plist_id);

        /*
         * Close/release resources.
         */
        H5Dclose(dset_id);
        H5Sclose(filespace);
        H5Sclose(memspace);
    }

    for(k = 0; k<3; k++)
    {
        const char *datasetname[3] = {"/X1", "/X2", "/X3", };


        filespace = H5Screate_simple(1, &dimsf[k], NULL);

        /*
        * Create the dataset with default properties and close filespace.
        */
        dset_id = H5Dcreate(file_id, datasetname[k], H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);

        memspace = H5Screate_simple(1, &count[k], NULL);

        /*
         * Select hyperslab in the file.
         */
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset[k], NULL, &count[k], NULL);


        int ii = 0;
        double X[NDIM];


        if(k == 0)
        {
            ZSLOOPI(0,0)
            {
                coord(hstart, i, jstart, CENT, X);
                data[ii++] = X[RR];
            }
        }
        else if(k == 1)
        {
            ZSLOOPJ(0,0)
            {
                coord(hstart, istart, j, CENT, X);
                data[ii++] = X[TH];
            }
        }
        else if(k == 2)
        {
#if N3>2
            ZSLOOPH(0,1)
#else
            ZSLOOPH(0,0)
#endif
            {
                coord(h, istart, jstart, CENT, X);
                data[ii++] = X[PH];
            }
        }


        /*
         * Create property list for collective dataset write.
         */
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace,
                          plist_id, data);
        H5Pclose(plist_id);

        /*
         * Close/release resources.
         */
        H5Dclose(dset_id);
        H5Sclose(filespace);
        H5Sclose(memspace);
    }

    //Create 3D vectors of coordinates and gdet
    for(k = 0; k<4; k++)
    {
        const char *datasetname[4] = {"/r", "/th", "/ph", "/gdet"};


        filespace = H5Screate_simple(3, dimsf, NULL);

        /*
        * Create the dataset with default properties and close filespace.
        */
        dset_id = H5Dcreate(file_id, datasetname[k], H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);

        memspace = H5Screate_simple(3, count, NULL);

        /*
         * Select hyperslab in the file.
         */
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);


        int ii = 0;
        double X[NDIM];
        double r, th, ph;


        if(k == 0)
        {
            DUMPZLOOP
            {
                coord(hstart, i, jstart, CENT, X);
                bl_coord(X, &r, &th, &ph);
                data[ii++] = r;
            }
        }
        else if(k == 1)
        {
            DUMPZLOOP
            {
                coord(hstart, istart, j, CENT, X);
                bl_coord(X, &r, &th, &ph);
                data[ii++] = th;
            }
        }
        else if(k == 2)
        {
            DUMPZLOOP
            {
                coord(h, istart, jstart, CENT, X);
                bl_coord(X, &r, &th, &ph);
                data[ii++] = ph;
            }
        }
        else if(k == 3)
        {
            DUMPZLOOP
            {
                struct of_geom geom;
                get_geometry(h,i,j,CENT,&geom) ;
                data[ii++] = geom.g;
            }
        }

        /*
         * Create property list for collective dataset write.
         */
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace,
                          plist_id, data);
        H5Pclose(plist_id);

        /*
         * Close/release resources.
         */
        H5Dclose(dset_id);
        H5Sclose(filespace);
        H5Sclose(memspace);
    }

    //Create 3D vectors of coordinates
    for(k = 0; k<3; k++)
    {
        const char *datasetname[4] = {"/x", "/y", "/z"};


        filespace = H5Screate_simple(3, dimsf, NULL);

        /*
        * Create the dataset with default properties and close filespace.
        */
        dset_id = H5Dcreate(file_id, datasetname[k], H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);

        memspace = H5Screate_simple(3, count, NULL);

        /*
         * Select hyperslab in the file.
         */
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);


        int ii = 0;
        double X[NDIM];
        double r, th, ph;


        if(k == 0)
        {
            DUMPZLOOP
            {
                coord(h, i, j, CENT, X);
                bl_coord(X, &r, &th, &ph);
                data[ii++] = r*sin(th)*cos(ph);
            }
        }
        else if(k == 1)
        {
            DUMPZLOOP
            {
                coord(h, i, j, CENT, X);
                bl_coord(X, &r, &th, &ph);
                data[ii++] = r*sin(th)*sin(ph);
            }
        }
        else if(k == 2)
        {
            DUMPZLOOP
            {
                coord(h, i, j, CENT, X);
                bl_coord(X, &r, &th, &ph);
                data[ii++] = r*cos(th);
            }
        }

        /*
         * Create property list for collective dataset write.
         */
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace,
                          plist_id, data);
        H5Pclose(plist_id);

        /*
         * Close/release resources.
         */
        H5Dclose(dset_id);
        H5Sclose(filespace);
        H5Sclose(memspace);
    }

	/* Costas Supplement 17/1/2017
         * for saving also the gcov gcon
   	 */
#if HDFDUMP_Metric 
//Create 3D vectors of coordinates and gdet

    data = (double *) realloc(data, count[0]*count[1]*count[2]*count[3]*count[4] * sizeof(double));
    if(data == NULL)
    {
    	fprintf(stderr,"can not reallocate %lld bytes rank=%d\n", count[0] * count[1] * count[2] * count[3] * count[4]  * sizeof(double), rank);
    	exit(-1);
    }

#if HDFDUMP_gcov
    #if HDFDUMP_gcon
         const char *datasetname[HDFDUMP_gcov+HDFDUMP_gcon] = {"/gcov", "/gcon"};
    #else
         const char *datasetname[HDFDUMP_gcov+HDFDUMP_gcon] = {"/gcov"};
    #endif
#else 
         const char *datasetname[HDFDUMP_gcov+HDFDUMP_gcon] = {"/gcon"};
#endif
int k_lim=HDFDUMP_gcov+HDFDUMP_gcon;
       if (k_lim==0){fprintf(stderr,"\n\n Warning (HDF5_DUMP): You choose to save metric, but you have not chosen a covariant or contravariant form to be dumped!\n\n");}

    for(k = 0; k<k_lim; k++)
    {
        filespace = H5Screate_simple(5, dimsf, NULL);

        /*
        * Create the dataset with default properties and close filespace.
        */
        dset_id = H5Dcreate(file_id, datasetname[k], H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);

        memspace = H5Screate_simple(5, count, NULL);

        /*
         * Select hyperslab in the file.
         */
	//fprintf(stderr,"We reached here A rank=%d\n",rank);
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
	
        int ii = 0,mm = 0,ll = 0, ff1=0,ff2=0;


	DUMPZLOOP
        {
        	struct of_geom geom;
                get_geometry(h,i,j,CENT,&geom) ;
		for(mm=0;mm<NDIM;mm++){
        		for(ll=0;ll<NDIM;ll++){
                                if ((k==0)&&(HDFDUMP_gcov==1)){
                			data[ii++] = geom.gcov[mm][ll];
					//fprintf(stderr,"this is the geom.gcon[%d][%d]=%e\n",ll,mm,geom.gcov[mm][ll]);
					ff1++;
                                }
                                else{
					data[ii++] = geom.gcon[mm][ll];
					//fprintf(stderr,"this is the geom.gcon[%d][%d]=%e\n",ll,mm,geom.gcon[mm][ll]);
					ff2++;
				}
			}
		}    
        }

        fprintf(stderr,"rank=%d this is the ff1=%d and ff2=%d\n",rank,ff1,ff2);

        /*
         * Create property list for collective dataset write.
         */
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

       status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace,plist_id, data);
        H5Pclose(plist_id);

        /*
         * Close/release resources.
         */
        H5Dclose(dset_id);
        H5Sclose(filespace);
        H5Sclose(memspace);
    }
#endif

    
#if (HDFDUMP_mtrc_trns)

if (HDFDUMP_mtrc_trns==1){
	data = (double *) realloc(data, count[0]*count[1]*count[2]*count[3]*count[4] * sizeof(double));
	//const char *datasetname2[] = {"/dxpdx"};
	int ii=0,mm = 0,ll = 0;
        double X[NDIM];
        double dxpdx_res[NDIM][NDIM];
	
	filespace = H5Screate_simple(5, dimsf, NULL);

        /*
        * Create the dataset with default properties and close filespace.
        */
        dset_id = H5Dcreate(file_id, "/dxpdx", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//fprintf(stderr,"We reached here B rank=%d\n",rank);
        H5Sclose(filespace);
        memspace = H5Screate_simple(5, count, NULL);

        /*
         * Select hyperslab in the file.
         */
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

	DUMPZLOOP
	        {
			coord(h, i, j, CENT, X);
			dxdxp_func(X,dxpdx_res);
			for(mm=0;mm<NDIM;mm++){
	        		for(ll=0;ll<NDIM;ll++){
                			data[ii++] =dxpdx_res[mm][ll];
				}
			}
	        }

        /*
         * Create property list for collective dataset write.
         */
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

	status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace,plist_id, data);
	H5Pclose(plist_id);

        /*
         * Close/release resources.
         */
        H5Dclose(dset_id);
        H5Sclose(filespace);
        H5Sclose(memspace);
}
#endif

    H5Fclose(file_id);
    free(data);


#ifdef MPI_USED
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if(rank == 0)
    {
	    file_id = H5Fopen(fullFileName, H5F_ACC_RDWR, H5P_DEFAULT);

        /* Create the data space for the attribute. */
        dimsf[0] = 1;
        {
            const char *datasetname[] = {"/N1", "/N2", "/N3", "/lim", "/DTr"};
            const int   datasetval[]  = {N1, N2, N3, lim, DTr};

            for(unsigned int i = 0; i < sizeof(datasetname)/sizeof(datasetname[0]); i++)
            {
		        filespace = H5Screate_simple(1, dimsf, NULL);
		        dset_id = H5Dcreate(file_id, datasetname[i], H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                status = H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &datasetval[i]);
                H5Dclose(dset_id);
		        H5Sclose(filespace);
            }
        }

        {
            const char *datasetname[] = {"/startx1", "/startx2", "/startx3", "/dx1", "/dx2", "/dx3", "/a", "/Rin", "/Rout", "/hslope", 
										"/gam", "/cour", "/tf", "/DTd", "/DTl", "/DTi", "/R0"};
            const double datasetval[]  = {startx[RR], startx[TH], startx[PH], dx[RR], dx[TH], dx[PH], a, Rin, Rout, hslope, 
										gam, cour, tf, DTd, DTl, DTi, R0};

            for(unsigned int i = 0; i < sizeof(datasetname)/sizeof(datasetname[0]); i++)
            {
		        filespace = H5Screate_simple(1, dimsf, NULL);
		        dset_id = H5Dcreate(file_id, datasetname[i], H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &datasetval[i]);
                H5Dclose(dset_id);
		        H5Sclose(filespace);
            }
        }

	    H5Fclose(file_id);

    }
    status++;

    return 0;
}


/**************************************************************/
int WriteBH_HorizonHDF(const char *fullFileName)
/**************************************************************/
{
    /*
     * HDF5 APIs definitions
     */
    hid_t       file_id, dset_id;         /* file and dataset identifiers */
    hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
    hsize_t     dimsf[3];                 /* dataset dimensions */
    double     *data;                    /* pointer to data buffer to write */
    hsize_t		count[3];	          /* hyperslab selection parameters */
    hsize_t		offset[3];
    hid_t		plist_id;                 /* property list identifier */
    int         h,i,j,k;
    herr_t		status;
    
    double      bh_horizon = 1. + sqrt(1. - a*a);


    if(rank == 0)
        fprintf(stderr,"DUMP     file=%s\n", fullFileName) ;

    /*
     * Set up file access property list with parallel I/O access
     */
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    if(plist_id < 0) printf("HDF Error 1\n");

    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    /*
     * Create a new file collectively and release property list identifier.
    */

    file_id = H5Fcreate(fullFileName, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    if(file_id < 0) printf("HDF Error 2\n");

    H5Pclose(plist_id);

    /*
     * Create the dataspace for the dataset.
     */

    dimsf[0] = N1;
    dimsf[1] = N2;
#if N3<=2
    dimsf[2] = N3;
#else
    dimsf[2] = N3+1;
#endif

    /*
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
    count[0] = (istop - istart + 1);
    count[1] =  (jstop - jstart + 1);
#if N3<=2
    count[2] =  (hstop - hstart + 1);
#else
    count[2] =  (hstop - hstart + 1)+(hrank==hsize-1?1:0);
#endif

    offset[0] = istart - (irank==0 ? 2 : 1) + di_rank_start;
    offset[1] = jstart - (jrank==0 ? 2 : 1) + dj_rank_start;
    offset[2] = (N3 == 1 ? 0 : hstart - 1 + dh_rank_start);


    /*
    	 * Initialize data buffer
    	 */
    data = (double *) malloc(count[0]*count[1]*count[2] * sizeof(double));
	if(data == NULL)
	{
		fprintf(stderr,"can not allocate %lld bytes rank=%d\n", count[0] * count[1] * count[2] * sizeof(double), rank);
		exit(-1);
	}

    //Create 1D vectors of coordinates
    for(k = 0; k<3; k++)
    {
        const char *datasetname[3] = {"/r_1D", "/th_1D", "/ph_1D"};


        filespace = H5Screate_simple(1, &dimsf[k], NULL);

        /*
        * Create the dataset with default properties and close filespace.
        */
        dset_id = H5Dcreate(file_id, datasetname[k], H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);

        memspace = H5Screate_simple(1, &count[k], NULL);

        /*
         * Select hyperslab in the file.
         */
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset[k], NULL, &count[k], NULL);


        int ii = 0;
        double X[NDIM];
        double r, th, ph;


        if(k == 0)
        {
            ZSLOOPI(0,0)
            {
                coord(hstart, i, jstart, CENT, X);
                bl_coord(X, &r, &th, &ph);
                
                r = r/Rout*bh_horizon;
                
                data[ii++] = r;
            }
        }
        else if(k == 1)
        {
            ZSLOOPJ(0,0)
            {
                coord(hstart, istart, j, CENT, X);
                bl_coord(X, &r, &th, &ph);
                data[ii++] = th;
            }
        }
        else if(k == 2)
        {
#if N3>2
            ZSLOOPH(0,1)
#else
            ZSLOOPH(0,0)
#endif
            {
                coord(h, istart, jstart, CENT, X);
                bl_coord(X, &r, &th, &ph);
                data[ii++] = ph;
            }
        }


        /*
         * Create property list for collective dataset write.
         */
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace,
                          plist_id, data);
        H5Pclose(plist_id);

        /*
         * Close/release resources.
         */
        H5Dclose(dset_id);
        H5Sclose(filespace);
        H5Sclose(memspace);
    }


    //Create 3D vectors of coordinates and gdet
    for(k = 0; k<4; k++)
    {
        const char *datasetname[4] = {"/r", "/th", "/ph", "/bh_horizon"};


        filespace = H5Screate_simple(3, dimsf, NULL);

        /*
        * Create the dataset with default properties and close filespace.
        */
        dset_id = H5Dcreate(file_id, datasetname[k], H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);

        memspace = H5Screate_simple(3, count, NULL);

        /*
         * Select hyperslab in the file.
         */
        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);


        int ii = 0;
        double X[NDIM];
        double r, th, ph;


        if(k == 0)
        {
            DUMPZLOOP
            {
                coord(hstart, i, jstart, CENT, X);
                bl_coord(X, &r, &th, &ph);
                r = r/Rout*bh_horizon;
                data[ii++] = r;
            }
        }
        else if(k == 1)
        {
            DUMPZLOOP
            {
                coord(hstart, istart, j, CENT, X);
                bl_coord(X, &r, &th, &ph);
                data[ii++] = th;
            }
        }
        else if(k == 2)
        {
            DUMPZLOOP
            {
                coord(h, istart, jstart, CENT, X);
                bl_coord(X, &r, &th, &ph);
                data[ii++] = ph;
            }
        }
        else if(k == 3)
        {
            DUMPZLOOP
            {
                coord(hstart, i, jstart, CENT, X);
                bl_coord(X, &r, &th, &ph);
                r = r/Rout;
                
                if(r <= bh_horizon)
	                data[ii++] = 1.;
	            else
	                data[ii++] = 0.;
            }
        }

        /*
         * Create property list for collective dataset write.
         */
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace,
                          plist_id, data);
        H5Pclose(plist_id);

        /*
         * Close/release resources.
         */
        H5Dclose(dset_id);
        H5Sclose(filespace);
        H5Sclose(memspace);
    }

    free(data);
    H5Fclose(file_id);

#ifdef MPI_USED
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if(rank == 0)
    {
	    file_id = H5Fopen(fullFileName, H5F_ACC_RDWR, H5P_DEFAULT);

        /* Create the data space for the attribute. */
        dimsf[0] = 1;
        {
            const char *datasetname[] = {"/N1", "/N2", "/N3"};
            const int   datasetval[]  = {N1, N2, N3};

            for(unsigned int i = 0; i < sizeof(datasetname)/sizeof(datasetname[0]); i++)
            {
		        filespace = H5Screate_simple(1, dimsf, NULL);
		        dset_id = H5Dcreate(file_id, datasetname[i], H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                status = H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &datasetval[i]);
                H5Dclose(dset_id);
		        H5Sclose(filespace);
            }
        }

	    H5Fclose(file_id);

    }
    status++;

	write_BH_Horizon_xdmf_xml();

    return 0;
}

void write_xdmf_xml(double ttime, int hdfNo)
{
    char hdf_file[128];
    FILE *xmf = 0;
    int  isxmf = (access("dumps/HARM_model.xmf", 0) == 0);


    sprintf(hdf_file, "dump%03d.h5",dump_cnt);

    /*
     * Open the file and write the XML description of the mesh..
     */
    xmf = fopen("dumps/HARM_model.xmf", hdfNo == 0 || !isxmf ? "w" : "r+");

    if(hdfNo == 0 || !isxmf)
    {
        fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
        fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
        fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
        fprintf(xmf, " <Domain>\n");
        fprintf(xmf, " <Grid GridType=\"Collection\" CollectionType=\"Temporal\">\n");
    }
    else
    {
        fseek(xmf, -strlen(" </Grid>\n" " </Domain>\n" "</Xdmf>\n"), SEEK_END);
    }

    fprintf(xmf, "   <Grid Name=\"Model for t=%.1lf\" GridType=\"Uniform\">\n", ttime);
    fprintf(xmf, "   <Time Value=\"%.1lf\"/>\n", ttime);
    fprintf(xmf, "     <Topology TopologyType=\"3DSMesh\" NumberOfElements=\"%d %d %d\"/>\n", N1, N2, N3>2?N3+1:N3);
    fprintf(xmf, "     <Geometry GeometryType=\"X_Y_Z\">\n");
    fprintf(xmf, "       <DataItem Name=\"x\" Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", N1, N2, N3>2?N3+1:N3);
    fprintf(xmf, "        coords.h5:/x\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Name=\"y\" Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", N1, N2, N3>2?N3+1:N3);
    fprintf(xmf, "        coords.h5:/y\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Name=\"z\" Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", N1, N2, N3>2?N3+1:N3);
    fprintf(xmf, "        coords.h5:/z\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Geometry>\n");
#if HDFDUMP_Rho
    fprintf(xmf, "     <Attribute Name=\"Density\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", N1, N2, N3>2?N3+1:N3);
    fprintf(xmf, "        %s:/Rho\n", hdf_file);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");

#if 0
    fprintf(xmf, "     <Attribute Name=\"RhoGdet\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem ItemType=\"Function\" Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\"\n", N1, N2, N3);
    fprintf(xmf, "        Function=\"($0 * $1)\">\n");
    fprintf(xmf, "         <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", N1, N2, N3);
    fprintf(xmf, "             %s:/Rho\n", hdf_file);
    fprintf(xmf, "         </DataItem>\n");
    fprintf(xmf, "         <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", N1, N2, N3);
    fprintf(xmf, "             coords.h5:/gdet\n");
    fprintf(xmf, "         </DataItem>\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
#endif

#endif
#if HDFDUMP_Energy
    fprintf(xmf, "     <Attribute Name=\"Energy\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", N1, N2, N3>2?N3+1:N3);
    fprintf(xmf, "        %s:/Energy\n", hdf_file);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
#endif
#if HDFDUMP_DivB
    fprintf(xmf, "     <Attribute Name=\"DivB\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", N1, N2, N3>2?N3+1:N3);
    fprintf(xmf, "        %s:/DivB\n", hdf_file);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
#endif
#if HDFDUMP_U1U2U3gdet
    fprintf(xmf, "     <Attribute Name=\"U1U2U3gdet\" AttributeType=\"Vector\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", N1, N2, N3>2?N3+1:N3);
    fprintf(xmf, "        %s:/U1U2U3gdet\n", hdf_file);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
#endif
#if HDFDUMP_B1B2B3gdet
    fprintf(xmf, "     <Attribute Name=\"B1B2B3gdet\" AttributeType=\"Vector\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", N1, N2, N3>2?N3+1:N3);
    fprintf(xmf, "        %s:/B1B2B3gdet\n", hdf_file);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
#endif
#if HDFDUMP_Ucon
    fprintf(xmf, "     <Attribute Name=\"ucon\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", N1, N2, N3>2?N3+1:N3);
    fprintf(xmf, "        %s:/ucon\n", hdf_file);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
#endif
#if HDFDUMP_Ucov
    fprintf(xmf, "     <Attribute Name=\"ucov\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", N1, N2, N3>2?N3+1:N3);
    fprintf(xmf, "        %s:/ucov\n", hdf_file);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
#endif
#if HDFDUMP_Bcon
    fprintf(xmf, "     <Attribute Name=\"bcon\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", N1, N2, N3>2?N3+1:N3);
    fprintf(xmf, "        %s:/bcon\n", hdf_file);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
#endif
#if HDFDUMP_Bcov
    fprintf(xmf, "     <Attribute Name=\"bcov\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", N1, N2, N3>2?N3+1:N3);
    fprintf(xmf, "        %s:/bcov\n", hdf_file);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
#endif

#if 0
    fprintf(xmf, "     <Attribute Name=\"gdet\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", N1, N2, N3);
    fprintf(xmf, "        coords.h5:/gdet\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
#endif

    fprintf(xmf, "   </Grid>\n");
    fprintf(xmf, " </Grid>\n");
    fprintf(xmf, " </Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);
}

void write_BH_Horizon_xdmf_xml(void)
{
    char hdf_file[128];
    FILE *xmf = 0;


    sprintf(hdf_file, "dump%03d.h5",dump_cnt);

    /*
     * Open the file and write the XML description of the mesh..
     */
    xmf = fopen("dumps/BH_Horizon.xmf", "w");

    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, " <Domain>\n");
    fprintf(xmf, "   <Grid Name=\"BH Horizon\" GridType=\"Uniform\">\n");
    fprintf(xmf, "     <Topology TopologyType=\"3DSMesh\" NumberOfElements=\"%d %d %d\"/>\n", N1, N2, N3>2?N3+1:N3);
    fprintf(xmf, "     <Geometry GeometryType=\"X_Y_Z\">\n");
    fprintf(xmf, "       <DataItem Name=\"r\" Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", N1, N2, N3>2?N3+1:N3);
    fprintf(xmf, "        bh_horizon.h5:/r\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Name=\"theta\" Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", N1, N2, N3>2?N3+1:N3);
    fprintf(xmf, "        bh_horizon.h5:/th\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Name=\"phi\" Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", N1, N2, N3>2?N3+1:N3);
    fprintf(xmf, "        bh_horizon.h5:/ph\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Geometry>\n");
    fprintf(xmf, "     <Attribute Name=\"BH_Horizon\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", N1, N2, N3>2?N3+1:N3);
    fprintf(xmf, "        bh_horizon.h5:/bh_horizon\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    fprintf(xmf, "   </Grid>\n");
    fprintf(xmf, " </Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);
}


#else
void dump(FILE *fp)
{
    double (*   p)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_p);
    double (*   gdet)[SN2][NPG] = (double (*) [SN2][NPG])(_gdet);
#if( COOL_EOS )
    double (*   CoolVal)[SN1][SN2][NCOOLVAL] = (double (*) [SN1][SN2][NCOOLVAL])(_CoolVal);
#endif

    int h,i,j,k ;
    double divb ;
    double X[NDIM] ;
    double r,th,ph,vmin,vmax ;
    struct of_geom geom ;
    struct of_state q ;

//    return;

#ifdef MPI_USED
	MPI_Sendrecv(p, 1, slice_iNPR1_stop, rank_nexti, 124,
                p, 1, slice_iNPR1_startP, rank_previ, 124,
                MPI_COMM_WORLD, &MPIStatus);

	if(N3>1)
	{
		MPI_Sendrecv(p, 1, slice_hNPR1_stop, rank_nexth, 124,
                p, 1, slice_hNPR1_startP, rank_prevh, 124,
                MPI_COMM_WORLD, &MPIStatus);
	}

#endif

    if (rank == 0)
    {
        /***************************************************************
          Write header information :
        ***************************************************************/
        fprintf(fp, FMT_DBL_OUT, t        );
        fprintf(fp, FMT_INT_OUT, N1       );
#if N3 == 1 || defined(PHSLICE)
        fprintf(fp, FMT_INT_OUT, N2       );
#endif
#if N3 > 1 && !defined(PHSLICE)
        fprintf(fp, FMT_INT_OUT, N3       );
#endif
        fprintf(fp, FMT_DBL_OUT, startx[RR]);
#if N3 == 1 || defined(PHSLICE)
        fprintf(fp, FMT_DBL_OUT, startx[TH]);
#endif
#if N3 > 1 && !defined(PHSLICE)
        fprintf(fp, FMT_DBL_OUT, startx[PH]);
#endif
        fprintf(fp, FMT_DBL_OUT, dx[RR]    );
#if N3 == 1 || defined(PHSLICE)
        fprintf(fp, FMT_DBL_OUT, dx[TH]    );
#endif
#if N3 > 1 && !defined(PHSLICE)
        fprintf(fp, FMT_DBL_OUT, dx[PH]    );
#endif
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
    }

#ifdef MPI_USED
    //Wait for previous task
	if (irank > 0 && hrank == 0)
	{
		int coords[] = {0, irank-1};
	    int wrank;
		extern MPI_Comm cart_grid;
		MPI_Cart_rank(cart_grid, coords, &wrank);
   	    MPI_Recv(&divb, 1, MPI_DOUBLE, wrank, 140, MPI_COMM_WORLD, &MPIStatus);
	}

#endif

    ZSLOOPI(0,0)
    {
#if defined(THSLICE)
	j = jstart + THSLICE;
#else

#ifdef MPI_USED
        //Wait for previous task
	if (hrank > 0)
	{
		int coords[] = {hrank-1, irank};
	    int wrank;
		extern MPI_Comm cart_grid;
		MPI_Cart_rank(cart_grid, coords, &wrank);
   	    MPI_Recv(&divb, 1, MPI_DOUBLE, wrank, 150, MPI_COMM_WORLD, &MPIStatus);
	}
#endif


	ZSLOOPJ(0,0)
#endif
	{

#if N3 > 1 && defined(PHSLICE)
        h = PHSLICE - dh_rank_start + 1;
        
        if(h >= hstart && h <= hstop )
#else
        ZSLOOPH(0,0)
#endif
        {
            coord(h,i,j,CENT,X) ;
            bl_coord(X,&r,&th,&ph) ;

            fprintf(fp, FMT_DBL_OUT, X[1]       );
#if N3 == 1 || defined(PHSLICE)
            fprintf(fp, FMT_DBL_OUT, X[2]       );
#endif
#if N3 > 1 && !defined(PHSLICE)
            fprintf(fp, FMT_DBL_OUT, X[3]       );
#endif
            fprintf(fp, FMT_DBL_OUT, r          );
#if N3 == 1 || defined(PHSLICE)
            fprintf(fp, FMT_DBL_OUT, th         );
#endif
#if N3 > 1 && !defined(PHSLICE)
            fprintf(fp, FMT_DBL_OUT, ph         );
#endif
            PLOOP if(k!=KTOT) fprintf(fp, FMT_DBL_OUT, p[h][i][j][k] );

#if( COOL_EOS )

	        InitCache_CoolEOS(h, i, j);

			CoolVal_rho0_u_CoolEOS(p[h][i][j][RHO], p[h][i][j][UU], h, i, j, CoolVal[h][i][j]);

//        CLOOP fprintf(fp, FMT_DBL_OUT, CoolVal[h][i][j][k] );
            //Lambda_sim
            fprintf(fp, FMT_DBL_OUT, CoolVal[h][i][j][COOL_LAMBDA_SIM] );
            //P
            fprintf(fp, FMT_DBL_OUT, CoolVal[h][i][j][COOL_P] );
			//T
            fprintf(fp, FMT_DBL_OUT, CoolVal[h][i][j][COOL_T] );
			//Hd
            fprintf(fp, FMT_DBL_OUT, CoolVal[h][i][j][COOL_Hd] );
			//Qnu
            fprintf(fp, FMT_DBL_OUT, CoolVal[h][i][j][COOL_Qnu] );
			//Tau
            fprintf(fp, FMT_DBL_OUT, CoolVal[h][i][j][COOL_Tau] );
			//Yee
            fprintf(fp, FMT_DBL_OUT, CoolVal[h][i][j][COOL_Yee] );
#else
            fprintf(fp, FMT_DBL_OUT, 0.0 );
            fprintf(fp, FMT_DBL_OUT, 0.0 );
            fprintf(fp, FMT_DBL_OUT, 0.0 );
            fprintf(fp, FMT_DBL_OUT, 0.0 );
            fprintf(fp, FMT_DBL_OUT, 0.0 );
            fprintf(fp, FMT_DBL_OUT, 0.0 );
            fprintf(fp, FMT_DBL_OUT, 0.0 );
#endif

            /* divb flux-ct defn; corner-centered.  Use
            only interior corners */
#if N3 == 1
            if(i > 0 && j > 0 && i < N1 && j < N2) {
                divb = fabs( 0.5*(
                                 p[h][i][j][B1]*gdet[i][j][CENT]
                                 + p[h][i][j-1][B1]*gdet[i][j-1][CENT]
                                 - p[h][i-1][j][B1]*gdet[i-1][j][CENT]
                                 - p[h][i-1][j-1][B1]*gdet[i-1][j-1][CENT]
                             )/dx[RR] +
                             0.5*(
                                 p[h][i][j][B2]*gdet[i][j][CENT]
                                 + p[h][i-1][j][B2]*gdet[i-1][j][CENT]
                                 - p[h][i][j-1][B2]*gdet[i][j-1][CENT]
                                 - p[h][i-1][j-1][B2]*gdet[i-1][j-1][CENT]
                             )/dx[TH]) ;
            }
            else divb = 0. ;
#else
			{
                divb = fabs( 0.25*(
								+ p[h-1][i][j][B1]*gdet[i][j][CENT]
    	                        + p[h-1][i][j-1][B1]*gdet[i][j-1][CENT]
        	                    - p[h-1][i-1][j][B1]*gdet[i-1][j][CENT]
            	                - p[h-1][i-1][j-1][B1]*gdet[i-1][j-1][CENT]

                                 + p[h][i][j][B1]*gdet[i][j][CENT]
                                 + p[h][i][j-1][B1]*gdet[i][j-1][CENT]
                                 - p[h][i-1][j][B1]*gdet[i-1][j][CENT]
                                 - p[h][i-1][j-1][B1]*gdet[i-1][j-1][CENT]
                             )/dx[RR] +
                             0.25*(
                                 + p[h-1][i][j][B2]*gdet[i][j][CENT]
                                 + p[h-1][i-1][j][B2]*gdet[i-1][j][CENT]
                                 - p[h-1][i][j-1][B2]*gdet[i][j-1][CENT]
                                 - p[h-1][i-1][j-1][B2]*gdet[i-1][j-1][CENT]

                                 + p[h][i][j][B2]*gdet[i][j][CENT]
                                 + p[h][i-1][j][B2]*gdet[i-1][j][CENT]
                                 - p[h][i][j-1][B2]*gdet[i][j-1][CENT]
                                 - p[h][i-1][j-1][B2]*gdet[i-1][j-1][CENT]
                             )/dx[TH] +
                             0.25*(
                                 + p[h][i][j][B3]*gdet[i][j][CENT]
                                 + p[h][i-1][j][B3]*gdet[i-1][j][CENT]
                                 - p[h-1][i][j][B3]*gdet[i][j][CENT]
                                 - p[h-1][i-1][j][B3]*gdet[i-1][j][CENT]

                                 + p[h][i][j-1][B3]*gdet[i][j][CENT]
                                 + p[h][i-1][j-1][B3]*gdet[i-1][j-1][CENT]
                                 - p[h-1][i][j-1][B3]*gdet[i][j-1][CENT]
                                 - p[h-1][i-1][j-1][B3]*gdet[i-1][j-1][CENT]
                             )/dx[PH]) ;
			}
#endif
            fprintf(fp, FMT_DBL_OUT, divb     );

            if(!failed) {
                get_geometry(h,i,j,CENT,&geom) ;
                get_state(p[h][i][j],&geom,&q) ;

                for(k=0; k<NDIM; k++) fprintf(fp,FMT_DBL_OUT,q.ucon[k]) ;
                for(k=0; k<NDIM; k++) fprintf(fp,FMT_DBL_OUT,q.ucov[k]) ;
                for(k=0; k<NDIM; k++) fprintf(fp,FMT_DBL_OUT,q.bcon[k]) ;
                for(k=0; k<NDIM; k++) fprintf(fp,FMT_DBL_OUT,q.bcov[k]) ;

                vchar(p[h][i][j],&q,&geom,1,&vmax,&vmin) ;
                fprintf(fp, FMT_DBL_OUT, vmin );
                fprintf(fp, FMT_DBL_OUT, vmax );

                vchar(p[h][i][j],&q,&geom,2,&vmax,&vmin) ;
                fprintf(fp, FMT_DBL_OUT, vmin );
                fprintf(fp, FMT_DBL_OUT, vmax );

                fprintf(fp, FMT_DBL_OUT, geom.g );

				{//Now find v^2:
  					double vsq = 0. ;
  					for(int id=1;id<4;id++)
    					for(int jd=1;jd<4;jd++) 
      						vsq += geom.gcov[id][jd]*q.ucon[id]*q.ucon[jd];      
  
  					double gsqtmp = 1. + vsq;
  					vsq /= gsqtmp;

	                fprintf(fp, FMT_DBL_OUT, sqrt(vsq) );
				}
				
/*				
				{ //ZAMO velocities, http://arxiv.org/abs/0906.1950v1
			        double vzamo[NDIM] ;
					double gtttilde = geom.gcov[TT][TT] - geom.gcov[TT][PH]*geom.gcov[TT][PH]/geom.gcov[PH][PH];
					if(gtttilde < 0.)
						vzamo[TT] = sqrt(-1.*gtttilde)*q.ucon[TT];
					else
						vzamo[TT] = 0.;

					vzamo[RR] = sqrt(geom.gcov[RR][RR])*q.ucon[RR];
					vzamo[TH] = sqrt(geom.gcov[TH][TH])*q.ucon[TH];
					vzamo[PH] = geom.gcov[TT][PH]/sqrt(geom.gcov[PH][PH])*q.ucon[TT] + 
    	                            sqrt(geom.gcov[PH][PH])*q.ucon[PH];

//					for(k=0;k<NDIM;k++) fprintf(fp,FMT_DBL_OUT,vzamo[k]) ;

					double vrr = 0., vth = 0., vph = 0.;
					if(vzamo[TT] != 0.)
					{
						vrr = vzamo[RR]/vzamo[TT];
						vth = vzamo[TH]/vzamo[TT];
						vph = vzamo[PH]/vzamo[TT];
					}
					fprintf(fp, FMT_DBL_OUT, vzamo[TT] );
					fprintf(fp, FMT_DBL_OUT, vrr );
					fprintf(fp, FMT_DBL_OUT, vth );
					fprintf(fp, FMT_DBL_OUT, vph );

				}
*/			
				{
					void kstobl(double r, double th, double *ucon);

					double uur,PS_cs,PS_pressure, PS_Mach;
			
					PS_pressure = p[h][i][j][UU] * (gam - 1.0);   // pressure of gas
			
					PS_cs = sqrt(4.0/3.0*PS_pressure / ( p[h][i][j][RHO] + 4.0 * PS_pressure)); //sound_speed
			

/*
					dxdxp_func(X, dxdxp);
				
					uur = dxdxp[1][1]*q.ucon[RR];    //radial velocity
*/
					//Transform velocieties from KS to BL
//					kstobl(r, th, q.ucon);
					uur = sqrt(geom.gcov[RR][RR])*q.ucon[RR]/q.ucon[TT];    //radial velocity
			
					PS_Mach = -uur/PS_cs;
			
					fprintf(fp,FMT_DBL_OUT,PS_pressure) ;
					fprintf(fp,FMT_DBL_OUT,PS_cs) ;
					fprintf(fp,FMT_DBL_OUT,uur) ;
					fprintf(fp,FMT_DBL_OUT,PS_Mach) ;
			    }
            }

            fprintf(fp,"\n") ;
        }

	}

	fflush(fp);

#ifdef MPI_USED
	//Signal the next task
	if (hrank < hsize - 1)
	{	
		int coords[] = {hrank+1, irank};
	    int wrank;
		extern MPI_Comm cart_grid;
		MPI_Cart_rank(cart_grid, coords, &wrank);
        MPI_Send(&divb, 1, MPI_DOUBLE, wrank, 150, MPI_COMM_WORLD);
	}
#endif

    }

#ifdef MPI_USED
    //Signal the next task
    if (irank < isize - 1 && hrank == 0)
	{	
		int coords[] = {0, irank+1};
	    int wrank;
		extern MPI_Comm cart_grid;
		MPI_Cart_rank(cart_grid, coords, &wrank);
        MPI_Send(&divb, 1, MPI_DOUBLE, wrank, 140, MPI_COMM_WORLD);
	}

    MPI_Barrier(MPI_COMM_WORLD);
#endif

#if 0
	h=0;
	j=jstart+(jstop-jstart)/2;
	for(i=istart;i<istop;i++)
	{
       coord(h,i,j,CENT,X) ;
       bl_coord(X,&r,&th,&ph) ;
       get_geometry(h,i,j,CENT,&geom) ;

       get_state(p[h][i][j],&geom,&q) ;
/*
        {
			double vobserver[NDIM] ;
			double u_tilde[NDIM] ;
			u_tilde[0] = 0. ;
			u_tilde[1] = p[h][i][j][U1];
			u_tilde[2] = p[h][i][j][U2];
			u_tilde[3] = p[h][i][j][U3];

  			double u_tilde_sq = 0. ;
  			for(int ii=0;ii<NDIM;ii++)
    		for(int jj=0;jj<NDIM;jj++)
		      u_tilde_sq += geom.gcov[ii][jj]*u_tilde[ii]*u_tilde[jj] ;
			u_tilde_sq = fabs(u_tilde_sq) ;

			double lapse = 1./sqrt(-geom.gcon[TT][TT]) ;
			double shift[NDIM] ;

			shift[0] = 0.0;
			for(int jj=1;jj<NDIM;jj++) shift[jj] = geom.gcon[TT][jj]*lapse*lapse ;

  			for(int ii=0;ii<NDIM;ii++) vobserver[ii] = u_tilde[ii]*lapse - shift[ii] ;

//  			for(int ii=0;ii<NDIM;ii++) vobserver[ii] = q.ucon ;
		    printf("r=%g vr/c=%g, vth/c=%g, vphi/c=%g rho=%g\n", r, vobserver[RR], vobserver[TH], vobserver[PH], p[h][i][j][RHO]);
		}
*/
/*
	{
  		// Now find v^2:

  		double vsq = 0. ;
  		for(int id=1;id<4;id++)
    		for(int jd=1;jd<4;jd++) 
      			vsq += geom.gcov[id][jd]*q.ucon[id]*q.ucon[jd];      
  
  		double gsqtmp = 1. + vsq;
  		vsq /= gsqtmp;

	    printf("r=%g v/c=%g rho=%g\n", r, sqrt(vsq), p[h][i][j][RHO]);
	}
*/

	{
  		// Now find v^2:

  		double vsq = 0. ;
  		for(int id=1;id<4;id++)
    		for(int jd=1;jd<4;jd++) 
      			vsq += geom.gcov[id][jd]*q.ucon[id]*q.ucon[jd];      
  
//  		double gsqtmp = 1. + vsq;
//  		vsq /= gsqtmp;

	    printf("r=%g v/c=%g rho=%g\n", r, sqrt(vsq), p[h][i][j][RHO]);
	}

	}
#endif
}
#endif
                         