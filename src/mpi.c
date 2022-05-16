#ifdef MPI_USED
#include "mpi.h"
#endif

#include "decs.h"
#include "mpi_ranking.c"

#ifdef MPI_USED
MPI_Comm cart_grid;
void set_MPI(int argc,char *argv[]);
void set_MPI_Slices(void);
#endif


#ifdef MPI_USED
/*****************************************************************/
void set_MPI(int argc,char *argv[])
{
	int coords[2];
	int dims[2] = { 0, 0}, period[2] = { 1, 0};   //boundry in Phi is periodic
    int kos_mpi_status=0,dh,di,help_dim;

	if(N3 == 1)
		dims[0] = 1; 
/*
	MPI_Dims_create is not aware of the application topology. E.g., if the application uses a grid of 1200 x 500 elements and 60 MPI 
	processes then the expected Cartesian process topology should be 12 x 5 and not 10 x 6 as returned by existing MPI_Dims_create.

	In case of nonequal grids it is better to provide custom process topology:
	mpirun -np 128 ./harm 16x1x8
*/

  /* Cartesian grid creation */
    if(argc >= 2)
	{
		int js;
		sscanf(argv[1], "%dx%dx%d", &isize, &js, &hsize);
		dims[0] = hsize;
		dims[1] = isize;

		if(js != 1)
		{
			printf("invalid slab jsize div != 1\n");
			exit(0);
		}

	} 
	else
	{
        kos_mpi_status=kos_mpi_ranking(nprocs,N1,N2,N3,dims);
        if (kos_mpi_status==-1){
                exit(123);
        }
        printf("Kostas mesagge - Kostas routine ranks division: dims[0]=%d dims[1]=%d\n",dims[0],dims[1]);

		hsize = dims[0];
		isize = dims[1];
	}

	if(rank == 0)
		printf("dims : %d x %d x %d\n", isize, jsize, hsize);

	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, 1, &cart_grid);
	MPI_Cart_coords(cart_grid,rank,2,coords);



   	hrank = coords[0];
	irank = coords[1];

	if(N3>1)
	{
		hstart = 1;
        dh_rank_start = N3/hsize*hrank;
        int AA=N3-hsize*(N3/hsize),BB=AA*(N3/hsize+1);
        int dh_rank_start_kostas=(hrank<=AA)?hrank*(N3/hsize+1):BB+(hrank-AA)*(N3/hsize);
           
		hstop = hstart + (hrank < hsize-1 ? N3/hsize-1: N3-dh_rank_start-1);//This is the problematic command
        if (N3% hsize ==0){
		    dh=N3/hsize;
		}
		else{
           if(hrank<(N3-hsize*(N3/hsize))){
               dh=N3/hsize+1;
           }
           else{
               dh=N3/hsize;
           }
        }
        dh_rank_start = dh_rank_start_kostas;
        hstop=hstart+dh-1;
        printf("\nKostas mesagge - MPI hstop: rank=%d hstop=%d hsize=%d quand=%d\n",rank,hstop,hsize,(N3-(hsize+1)*(N3/(hsize+1))));
        printf("Kostas mesagge - MPI hstop: rank=%d irank=%d dh=%d\n",rank,irank,dh);
        printf("Kostas mesagge - MPI hstop: dh_rank_start = %d dh_rank_start_kostas=%d\n",dh_rank_start,dh_rank_start_kostas);
                
	}
	else
	{
		hstart = hstop = 0;
	}


    int AAi=isize*(N1/isize+1)-N1,BBi=AAi*(N1/isize);
    int di_rank_start_kostas=(irank<=AAi)?irank*(N1/isize):BBi+(irank-AAi)*(N1/isize+1);
    printf("Kostas mesagge - MPI hstop: AAi = %d BBi=%d di_rank_start_kostas=%d\n",AAi,BBi,di_rank_start_kostas);
	istart = 1 + (irank==0 ? 1 : 0);
    di_rank_start = N1/isize*irank;
	istop = istart + (irank < isize-1 ? N1/isize-1: N1-N1/isize*irank-1);
    if (N1% isize ==0){di=N1/isize;}
    else{
        if(irank<=(isize*(1+N1/isize)-N1-1)){ //the first rank has two ghost so load it with less points
            di=N1/isize;
        }
        else{
            di=N1/isize+1;
        }
    }
    di_rank_start = di_rank_start_kostas;
    istop=istart+di-1;
    printf("Kostas mesagge - MPI istop: rank=%d istart=%d istop=%d isize=%d quand_i=%d\n",rank,istart,istop,isize,isize*(1+N1/isize)-N1);
    printf("Kostas mesagge - MPI istop: rank=%d hrank=%d di=%d\n",rank,hrank,di);
         

	if(istop - istart + 1 < 1)
	{
		printf("invalid slab div!\n");
		exit(0);
	}

	MPI_Cart_shift(cart_grid, 0, 1, &rank_prevh, &rank_nexth);
	MPI_Cart_shift(cart_grid, 1, 1, &rank_previ, &rank_nexti);

	SN1 = istop-istart+1 + (irank==0 ? 2 : 1) + (irank==isize-1 ? 2 : 1); 
	SN2 = jstop-jstart+1 + (jrank==0 ? 2 : 1) + (jrank==jsize-1 ? 2 : 1); 
	SN3 = hstop-hstart+1 + (N3>1 ? 1 : 0) + (N3>1 ? 1 : 0); 


	printf("rank %d SN1=%d SN2=%d SN3=%d\n",rank, SN1, SN2, SN3);

//	printf("rank %d PH=[%d, %d]\n", rank, hstart, hstop);
//	printf("rank %d RHO=[%d, %d]\n", rank, istart, istop);
//	printf("rank %d TH=[%d, %d]\n", rank, jstart, jstop);


   //////////////////////////////////////////////////////////////////////////
   // Create MPI default types for 2D slices

	{
//		MPI_Datatype slice_iNPR1_start;
//		MPI_Datatype slice_iNPR1_startP;

//		MPI_Datatype slice_iNPR1_stop;
//		MPI_Datatype slice_iNPR1_stopN;

    	int bigsizes[4]  = {SN3,SN1,SN2,NPR};
	    int subsizes[4]  = {N3>1?hstop-hstart+1+2:1,1,SN2,NPR};
    	int starts[4] = {N3>1?hstart-1:0,istart,0,0};

        starts[1] = istart;
	    MPI_Type_create_subarray(4, bigsizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &slice_iNPR1_start);
    	MPI_Type_commit(&slice_iNPR1_start);

        starts[1] = istart-1;
	    MPI_Type_create_subarray(4, bigsizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &slice_iNPR1_startP);
    	MPI_Type_commit(&slice_iNPR1_startP);

        starts[1] = istop;
	    MPI_Type_create_subarray(4, bigsizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &slice_iNPR1_stop);
    	MPI_Type_commit(&slice_iNPR1_stop);

        starts[1] = istop+1;
	    MPI_Type_create_subarray(4, bigsizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &slice_iNPR1_stopN);
    	MPI_Type_commit(&slice_iNPR1_stopN);

	}
/*
	{
//		MPI_Datatype slice_iNPR2_start;
//		MPI_Datatype slice_iNPR2_startP;

//		MPI_Datatype slice_iNPR2_stop;
//		MPI_Datatype slice_iNPR2_stopN;

    	int bigsizes[4]  = {N3+2,N1+4,N2+4,NPR};
	    int subsizes[4]  = {hstop-hstart+1,2,N2+4,NPR};
    	int starts[4] = {hstart,istart,0,0};

        starts[1] = istart;
	    MPI_Type_create_subarray(4, bigsizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &slice_iNPR2_start);
    	MPI_Type_commit(&slice_iNPR2_start);

        starts[1] = istart-2;
	    MPI_Type_create_subarray(4, bigsizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &slice_iNPR2_startP);
    	MPI_Type_commit(&slice_iNPR2_startP);

        starts[1] = istop-1;
	    MPI_Type_create_subarray(4, bigsizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &slice_iNPR2_stop);
    	MPI_Type_commit(&slice_iNPR2_stop);

        starts[1] = istop+1;
	    MPI_Type_create_subarray(4, bigsizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &slice_iNPR2_stopN);
    	MPI_Type_commit(&slice_iNPR2_stopN);

	}
*/
	{
//		MPI_Datatype slice_i1_start;
//		MPI_Datatype slice_i1_startP;

//		MPI_Datatype slice_i1_stop;
//		MPI_Datatype slice_i1_stopN;

    	int bigsizes[3]  = {SN3,SN1,SN2};
	    int subsizes[3]  = {N3>1?hstop-hstart+1+2:1,1,SN2};
    	int starts[3] = {N3>1?hstart-1:0,istart,0};

        starts[1] = istart;
	    MPI_Type_create_subarray(3, bigsizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &slice_i1_start);
    	MPI_Type_commit(&slice_i1_start);

        starts[1] = istart-1;
	    MPI_Type_create_subarray(3, bigsizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &slice_i1_startP);
    	MPI_Type_commit(&slice_i1_startP);

        starts[1] = istop;
	    MPI_Type_create_subarray(3, bigsizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &slice_i1_stop);
    	MPI_Type_commit(&slice_i1_stop);

        starts[1] = istop+1;
	    MPI_Type_create_subarray(3, bigsizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &slice_i1_stopN);
    	MPI_Type_commit(&slice_i1_stopN);

	}
	{
//		MPI_Datatype slice_i1c_start;
//		MPI_Datatype slice_i1c_startP;

//		MPI_Datatype slice_i1c_stop;
//		MPI_Datatype slice_i1c_stopN;

    	int bigsizes[3]  = {SN3,SN1,SN2};
	    int subsizes[3]  = {N3>1?hstop-hstart+1+2:1,1,SN2};
    	int starts[3] = {N3>1?hstart-1:0,istart,0};

        starts[1] = istart;
	    MPI_Type_create_subarray(3, bigsizes, subsizes, starts, MPI_ORDER_C, MPI_CHAR, &slice_i1c_start);
    	MPI_Type_commit(&slice_i1c_start);

        starts[1] = istart-1;
	    MPI_Type_create_subarray(3, bigsizes, subsizes, starts, MPI_ORDER_C, MPI_CHAR, &slice_i1c_startP);
    	MPI_Type_commit(&slice_i1c_startP);

        starts[1] = istop;
	    MPI_Type_create_subarray(3, bigsizes, subsizes, starts, MPI_ORDER_C, MPI_CHAR, &slice_i1c_stop);
    	MPI_Type_commit(&slice_i1c_stop);

        starts[1] = istop+1;
	    MPI_Type_create_subarray(3, bigsizes, subsizes, starts, MPI_ORDER_C, MPI_CHAR, &slice_i1c_stopN);
    	MPI_Type_commit(&slice_i1c_stopN);

	}

	//N3 direction
#if N3 > 1
	{
		{
//			MPI_Datatype slice_hNPR1_start;
//			MPI_Datatype slice_hNPR1_startP;

//			MPI_Datatype slice_hNPR1_stop;
//			MPI_Datatype slice_hNPR1_stopN;

        	int bigsizes[4]  = {SN3,SN1,SN2,NPR};
		    int subsizes[4]  = {1,SN1,SN2,NPR};
    		int starts[4] = {hstart,0,0,0};

	        starts[0] = hstart;
		    MPI_Type_create_subarray(4, bigsizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &slice_hNPR1_start);
    		MPI_Type_commit(&slice_hNPR1_start);

        	starts[0] = hstart-1;
		    MPI_Type_create_subarray(4, bigsizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &slice_hNPR1_startP);
    		MPI_Type_commit(&slice_hNPR1_startP);

	        starts[0] = hstop;
	 	   	MPI_Type_create_subarray(4, bigsizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &slice_hNPR1_stop);
    		MPI_Type_commit(&slice_hNPR1_stop);

	        starts[0] = hstop+1;
		    MPI_Type_create_subarray(4, bigsizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &slice_hNPR1_stopN);
	    	MPI_Type_commit(&slice_hNPR1_stopN);
		}

		{
//			MPI_Datatype slice_h1_start;
//			MPI_Datatype slice_h1_startP;

//			MPI_Datatype slice_h1_stop;
//			MPI_Datatype slice_h1_stopN;

	    	int bigsizes[3]  = {SN3,SN1,SN2};
		    int subsizes[3]  = {1,SN1,SN2};
    		int starts[3] = {hstart,0,0};

	        starts[0] = hstart;
		    MPI_Type_create_subarray(3, bigsizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &slice_h1_start);
    		MPI_Type_commit(&slice_h1_start);

	        starts[0] = hstart-1;
		    MPI_Type_create_subarray(3, bigsizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &slice_h1_startP);
    		MPI_Type_commit(&slice_h1_startP);

	        starts[0] = hstop;
		    MPI_Type_create_subarray(3, bigsizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &slice_h1_stop);
    		MPI_Type_commit(&slice_h1_stop);

	        starts[0] = hstop+1;
		    MPI_Type_create_subarray(3, bigsizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &slice_h1_stopN);
    		MPI_Type_commit(&slice_h1_stopN);
    	}
	}
#endif

}


/*****************************************************************/
void set_MPI_Slices(void)
{

   //////////////////////////////////////////////////////////////////////////
   // Create MPI types for 2D slices
	//N3 direction
#if N3 > 1
	{

		//Find opposite hrank and h 
		for(int h=hstart; h<=hstop; h++)
		{
		    double r,th,ph, opphi;
		    double X[NDIM] ;

			int coords[] = {hrank, irank};
		    int oph, ophrank, wrank;
			extern MPI_Comm cart_grid;

			//Get actual phi
	        coord(h,istart,jstart,CENT,X) ;
    	    bl_coord(X,&r,&th,&ph) ;

			opphi = ph + M_PI;
			if(opphi >= M_PI*2.) opphi -= M_PI*2.;

			oph = ophrank = -1;

			int rs = N3/hsize;

			for(int hh=0; hh<N3; hh++)
			{
#if N3 > 2
				double cph = startx[PH] + (hh + 0.5)*dx[PH];
#else				
				double cph = startx[PH] + hh*dx[PH];
#endif
				if(fabs(cph-opphi) < dx[PH]/2.)
				{
					ophrank = hh/rs; 
					if(ophrank > hsize-1) ophrank = hsize-1;

					oph = hh - ophrank*rs + 1;
					break;
				}          
			}


           int AA=N3-hsize*(N3/hsize),BB=AA*(N3/hsize+1);
           int ophrank_kostas=-1,oph_kostas=-1,wrank_kostas;
           int xx,oxx;

           if (N3%2!=0){
               printf("\n\nError, please put an even resolution in phi direction!");
               exit(123);
           }
           xx=dh_rank_start+h;
           //oxx=(xx+N3/2)%(N3);
           oxx=(xx+N3/2==N3)?N3:((xx+N3/2)%(N3));
           if (oxx<=BB){
                ophrank_kostas=(oxx-1)/((N3/hsize)+1);
                oph_kostas=(oxx-1)%((N3/hsize)+1)+1;
           }
           else{
                ophrank_kostas=AA+(oxx-1-BB)/(N3/hsize);
                oph_kostas=((oxx-1-BB)%(N3/hsize))+1;
           }

			if(oph == -1 || ophrank < 0)
			{
				printf("Can't find opposite h for ph=%E\n", ph);
				exit(0);
			}			

            //ophrank = ophrank_kostas;
            coords[0] = ophrank_kostas;
            MPI_Cart_rank(cart_grid, coords, &wrank_kostas);

			coords[0] = ophrank;
			MPI_Cart_rank(cart_grid, coords, &wrank);
                           
            ophrank = ophrank_kostas;
            oph = oph_kostas;
            wrank= wrank_kostas;

			//printf("h = %d phi=%E opphi =%E, ophrank=%d, wrank=%d, oph=%d\n", h, ph, opphi, ophrank, wrank, oph);
            printf("KD: h = %d phi=%E rank=%d opphi =%E, ophrank_kos=%d, wrank=%d, oph_kos=%d xx=%d oxx=%d\n", h, ph,rank, opphi, ophrank_kostas, wrank_kostas, oph_kostas,xx,oxx); // DON'T touch string format if python debugging tool is to work.

#if 1
			{
				slice_h_rankop[h] = wrank;
				slice_h_hop[h] = oph;

	    		int bigsizes[4]  = {SN3,SN1,SN2,NPR};
			    int subsizes[4]  = {1,SN1,2,NPR};
	    		int starts[4] = {h,0,jstart-2,0};


	        	starts[2] = jstart-2;
		    	MPI_Type_create_subarray(4, bigsizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &slice_hopNPR2_m2[h]);
    			MPI_Type_commit(&slice_hopNPR2_m2[h]);

	        	starts[2] = jstop+1;
		    	MPI_Type_create_subarray(4, bigsizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &slice_hopNPR2_n1[h]);
    			MPI_Type_commit(&slice_hopNPR2_n1[h]);
			}

#endif
		}


	}
#endif

}

#endif

void debugcalc(double *_pr, double val[NPR])
{
    double (*   pr)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_pr);
    int h,i,j,k;

//	fprintf(stdout, "debugcalc rank=%d irank=%d, hrank=%d\n", rank, irank, hrank);

#ifdef MPI_USED
/*
    //Wait for previous task
	if (irank > 0 && hrank == 0)
	{
		int coords[] = {hsize-1, irank-1};
	    int wrank;
		extern MPI_Comm cart_grid;
		MPI_Cart_rank(cart_grid, coords, &wrank);
   	    MPI_Recv(val, NPR, MPI_DOUBLE, wrank, 140, MPI_COMM_WORLD, &MPIStatus);
	}
*/
#endif

    ZSLOOPI(-2,2)
    {

#ifdef MPI_USED


        //Wait for previous task
	if (hrank > 0)
	{
		int coords[] = {hrank-1, irank};
	    int wrank;
		extern MPI_Comm cart_grid;
		MPI_Cart_rank(cart_grid, coords, &wrank);
   	    MPI_Recv(val, NPR, MPI_DOUBLE, wrank, 150, MPI_COMM_WORLD, &MPIStatus);
	}

		if (hrank == 0)
		{	
			if(i == istart && irank > 0)
			{
				int coords[] = {hsize-1, irank-1};
	    		int wrank;
				extern MPI_Comm cart_grid;
				MPI_Cart_rank(cart_grid, coords, &wrank);
   	    		MPI_Recv(val, NPR, MPI_DOUBLE, wrank, 120, MPI_COMM_WORLD, &MPIStatus);
			}
			else if(i > istart -2)
			{
				int coords[] = {hsize-1, irank};
	    		int wrank;
				extern MPI_Comm cart_grid;
				MPI_Cart_rank(cart_grid, coords, &wrank);
   	    		MPI_Recv(val, NPR, MPI_DOUBLE, wrank, 110, MPI_COMM_WORLD, &MPIStatus);
			}
		}

/*
    //Wait for previous task
    	if (irank > 0 && hrank == 0)
    	{
        	int coords[] = {hsize-1, irank-1};
        	int wrank;
        	extern MPI_Comm cart_grid;
        	MPI_Cart_rank(cart_grid, coords, &wrank);
			fprintf(stdout, "1. debugcalc waiting for wrank=%d\n", wrank);
        	MPI_Recv(val, NPR, MPI_DOUBLE, wrank, 140, MPI_COMM_WORLD, &MPIStatus);
			fprintf(stdout, "1. debugcalc waiting for wrank=%d done\n", wrank);
	    }

        //Wait for previous task
        if (hrank > 0)
        {
            int coords[] = {hrank-1, irank};
            int wrank;
            extern MPI_Comm cart_grid;
            MPI_Cart_rank(cart_grid, coords, &wrank);
			fprintf(stdout, "2. debugcalc waiting for wrank=%d\n", wrank);
            MPI_Recv(val, NPR, MPI_DOUBLE, wrank, 150, MPI_COMM_WORLD, &MPIStatus);
			fprintf(stdout, "2. debugcalc waiting for wrank=%d done\n", wrank);
        }
*/

/*
		if (hsize > 1 && hrank > 0)
        {
            int coords[] = {hsize-1, isize-1};
            int wrank;
            extern MPI_Comm cart_grid;
            MPI_Cart_rank(cart_grid, coords, &wrank);
			fprintf(stdout, "3. debugcalc waiting for wrank=%d\n", wrank);
            MPI_Recv(val, NPR, MPI_DOUBLE, wrank, 150, MPI_COMM_WORLD, &MPIStatus);
			fprintf(stdout, "3. debugcalc waiting for wrank=%d done\n", wrank);
        }
*/
#endif


        ZSLOOPH(0,0)
        {

	        ZSLOOPJ(-2,2)
            {

                //val[RHO] += pr[h][i][j][RHO];
                //val[UU] += pr[h][i][j][UU];
                PLOOP val[k] += pr[h][i][j][k];

//		for(k=0;k<3;k++) val[k] += pr[h][i][j][k];

//		val[0] += pr[h][i][j][0];
//		val[1] += pr[h][i][j][1];
//		val[2] += pr[h][i][j][2];
                /*
                	int arraidx =  h*SN1*SN2*NPR + i*SN2*NPR+j*NPR;
                	if(&pr[h][i][j][0] != _pr + arraidx)
                	{
                		fprintf(stderr, "0 addr addr not match for rank=%d, i=%d,j=%d, h=%d : %p vs %p placemt=%ld, _pr=%p\n", rank, i, j, h, &pr[h][i][j][0], _pr + arraidx, arraidx, _pr);
                		assert(&pr[h][i][j][0] == _pr + arraidx);
                	}

                	arraidx =  h*SN1*SN2*NPR + i*SN2*NPR+j*NPR + 1;
                	if(&pr[h][i][j][1] != _pr + arraidx)
                	{
                		fprintf(stderr, "1 addr addr not match for rank=%d, i=%d,j=%d, h=%d : %p vs %p placemt=%ld, _pr=%p\n", rank, i, j, h, &pr[h][i][j][1], _pr + arraidx, arraidx, _pr);
                		assert(&pr[h][i][j][1] == _pr + arraidx);
                	}

                	arraidx =  h*SN1*SN2*NPR + i*SN2*NPR+j*NPR + 2;
                	if(&pr[h][i][j][2] != _pr + arraidx)
                	{
                		fprintf(stderr, "2 addr addr not match for rank=%d, i=%d,j=%d, h=%d : %p vs %p placemt=%ld, _pr=%p\n", rank, i, j, h, &pr[h][i][j][2], _pr + arraidx, arraidx, _pr);
                		assert(&pr[h][i][j][2] == _pr + arraidx);
                	}
                */
            }
        }

#ifdef MPI_USED
/*        //Signal the next task
        if (hrank < hsize - 1)
        {
            int coords[] = {hrank+1, irank};
            int wrank;
            extern MPI_Comm cart_grid;
            MPI_Cart_rank(cart_grid, coords, &wrank);
			fprintf(stdout, "4. debugcalc sending to wrank=%d\n", wrank);
            MPI_Send(&val, NPR, MPI_DOUBLE, wrank, 150, MPI_COMM_WORLD);
			fprintf(stdout, "4. debugcalc sending to wrank=%d done\n", wrank);
        }

	    if (hrank == hsize - 1)
        {
            int coords[] = {0, irank+1};
            int wrank;
            extern MPI_Comm cart_grid;
            MPI_Cart_rank(cart_grid, coords, &wrank);
			fprintf(stdout, "5. debugcalc waiting for wrank=%d\n", wrank);
            MPI_Recv(val, NPR, MPI_DOUBLE, wrank, 150, MPI_COMM_WORLD, &MPIStatus);
			fprintf(stdout, "5. debugcalc waiting for wrank=%d done\n", wrank);
        }
*/

		//Signal the next task
		if (hrank < hsize - 1)
		{	
			int coords[] = {hrank+1, irank};
	    	int wrank;
			extern MPI_Comm cart_grid;
			MPI_Cart_rank(cart_grid, coords, &wrank);
    	    MPI_Send(val, NPR, MPI_DOUBLE, wrank, 150, MPI_COMM_WORLD);
		}
		if (hrank == hsize - 1) 
		{	
			if(i == istop && irank < isize - 1)
			{
				int coords[] = {0, irank+1};
		    	int wrank;
				extern MPI_Comm cart_grid;
				MPI_Cart_rank(cart_grid, coords, &wrank);
    	    	MPI_Send(val, NPR, MPI_DOUBLE, wrank, 120, MPI_COMM_WORLD);
			} 
			else if(i == istop + 2 && irank == isize - 1)
			{
			}
			else
			{
				int coords[] = {0, irank};
		    	int wrank;
				extern MPI_Comm cart_grid;
				MPI_Cart_rank(cart_grid, coords, &wrank);
    	    	MPI_Send(val, NPR, MPI_DOUBLE, wrank, 110, MPI_COMM_WORLD);
			}
		}


#endif
    }

#ifdef MPI_USED
/*
    //Signal the next task
    if (irank < isize - 1 && hrank == 0)
    {
        int coords[] = {0, irank+1};
        int wrank;
        extern MPI_Comm cart_grid;
        MPI_Cart_rank(cart_grid, coords, &wrank);
		fprintf(stdout, "6. debugcalc sending to wrank=%d\n", wrank);
        MPI_Send(val, NPR, MPI_DOUBLE, wrank, 140, MPI_COMM_WORLD);
		fprintf(stdout, "6. debugcalc sending to wrank=%d done\n", wrank);
    }


*/
/*
    //Signal the next task
    if (irank < isize - 1 && hrank == hsize-1)
	{	
		int coords[] = {0, irank+1};
	    int wrank;
		extern MPI_Comm cart_grid;
		MPI_Cart_rank(cart_grid, coords, &wrank);
        MPI_Send(val, NPR, MPI_DOUBLE, wrank, 140, MPI_COMM_WORLD);
	}
*/
    if(nprocs > 1)
    {
        if (irank == isize - 1 && hrank == hsize - 1)
            MPI_Send(val, NPR, MPI_DOUBLE, 0, 160, MPI_COMM_WORLD);

        if (rank == 0)
        {
            int coords[] = {hsize-1, isize-1};
            int wrank;
            extern MPI_Comm cart_grid;
            MPI_Cart_rank(cart_grid, coords, &wrank);
            MPI_Recv(val, NPR, MPI_DOUBLE, wrank, 160, MPI_COMM_WORLD, &MPIStatus);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
#endif

}


void debug2(double *_pr, const char *where)
{
    double (*   pr)[SN1][SN2][NPR] = (double (*) [SN1][SN2][NPR])(_pr);
/*
    int h,i,j,k;
    {
        double rho = 0.0;
        double u = 0.0;

        ZSLOOP(-1,1,-2,2,-2,2) {
            rho += pr[h][i][j][RHO];
            u += pr[h][i][j][UU];
        }
//	ZLOOP{rho += pr[h][i][j][RHO]; u += pr[h][i][j][UU];}

#ifdef MPI_USED
        MPI_Allreduce(MPI_IN_PLACE,&rho,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,&u,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif

        if(rank == 0)
            fprintf(stdout,"1. %s - rho="FMT_DBL_OUT" u="FMT_DBL_OUT"\n", where, rho, u) ;
    }
*/	
    double val[NPR];
    memset(val, 0, sizeof(val));
    /*
        ZSLOOP(-1,1,-2,2,-2,2)
    	{
    		//val[RHO] += pr[h][i][j][RHO];
    		//val[UU] += pr[h][i][j][UU];
    //    	PLOOP val[k] += pr[h][i][j][k];

    		for(k=0;k<3;k++) val[k] += pr[h][i][j][k];

    //		val[0] += pr[h][i][j][0];
    //		val[1] += pr[h][i][j][1];
    //		val[2] += pr[h][i][j][2];

    	}
    */
    debugcalc(_pr, val);

//#ifdef MPI_USED
//    MPI_Allreduce(MPI_IN_PLACE,val,NPR,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
//#endif
    if(rank == 0)
    {
	int k;
        PLOOP fprintf(stdout,"%s - val[%i]=" FMT_DBL_OUT "\n", where, k, val[k]) ;
    }
}
