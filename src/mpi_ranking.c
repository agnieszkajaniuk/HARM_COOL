#ifdef MPI_USED
#include "mpi.h"
#endif

int divisors(int numb, int mat[]); //finds the divisors
int divs_criteria(int NN1,int NN2,int NN3,int mat_length, int mat[],double divs_crit[]); // calculate the two criteria, NN2 doesn't participate in the grid points scattering. The smallest the better
const float fac_load=1.5; // system dependend constant. How slower is mpi communication for one point comparing to the calculations done because of it 

int kos_mpi_ranking(int num_total_process,int NN1, int NN2, int NN3, int prop_dim[])
{
     int *divs_mat;
     int divs_length;
     double *divs_crit;
     int divs_length_crit;

     
     int i;
     double min_crit_value;
     int proper_division;

     divs_mat = (int *)malloc(2*num_total_process*sizeof(int));
     divs_length=divisors(num_total_process,divs_mat);
     divs_mat = (int *)realloc(divs_mat,2*divs_length*sizeof(int));

     divs_crit = (double *)malloc(5*divs_length*sizeof(double));
     divs_length_crit = divs_criteria(NN1,NN2,NN3,divs_length,divs_mat,divs_crit);
     divs_crit = (double *)realloc(divs_crit,5*divs_length_crit*sizeof(double));

     printf("\n MPI process analyzing\n");
     printf("------------------------\n\n");
     printf("\nTotal of process:     %d\n",num_total_process);
     printf("Grid Size [N1,N2,N3]: %d x %d x %d\n",NN1,NN2,NN3);
     printf("Division pairs found: %d \n",divs_length);
     printf("Division pairs elligible for analysis: %d \n",divs_length_crit);


     printf("\n AA) [   nx , ny   ]     -> criteria=[ r_load  , surf_dim , fin_crit]\n");
     printf("----------------------------------------------------------------------\n");
    

     for(i=0;i<divs_length_crit;i++){
          if (i==0){proper_division=i;min_crit_value=divs_crit[5*i+4];}
          else{
               if(divs_crit[5*i+4]<min_crit_value){proper_division=i;min_crit_value=divs_crit[5*i+4];}
               else if(divs_crit[5*i+4]==min_crit_value){printf("\nWARNING: another best ranking division is found! It is ignored...\n");}
          }
          printf("%3d) [ %4d , %-4d ]     ->          [%f , %f , %f]\n",i+1,(int)floor(divs_crit[5*i]),(int)floor(divs_crit[5*i+1]),divs_crit[5*i+2],divs_crit[5*i+3],divs_crit[5*i+4]);
     }

     prop_dim[1]=(int)floor(divs_crit[5*proper_division]);
     prop_dim[0]=(int)floor(divs_crit[5*proper_division+1]);

     printf("\n Best resolution\n");
     printf("-----------------\n");
     printf("aa=%d n1=%d n3=%d total_criteria=%lf\n\n\n",proper_division+1,(int)floor(divs_crit[5*proper_division]),(int)floor(divs_crit[5*proper_division+1]),divs_crit[5*proper_division+4]);

     
     free(divs_mat);
     free(divs_crit);
     if(divs_length_crit==0){ return -1;}
     else {
          return divs_length_crit;
     }
}


int divs_criteria(int NN1,int NN2,int NN3,int mat_length, int mat[],double divs_crit[])
{
     int i,n1,n3,mn1,mn3,mx1,mx3;
     double r_load; // the unbalance point criterio
     double surf_load; // the surface communication criterio
     double tot_load; // the combination of the above criteria
     int num_ellig_crit=0; // the total number of elligable configuration


     for(i=0;i<mat_length;i++){
          n1=mat[2*i];
          n3=mat[2*i+1];

          if ((n1>NN1)||(n3>NN3)){
               r_load=-1.;surf_load=-1.;tot_load=-1.;
          }
          else{
               mn1=floor(NN1/n1);
               mn3=floor(NN3/n3);
               mx1=(NN1%n1==0)?mn1:mn1+1;
               mx3=(NN3%n3==0)?mn3:mn3+1;

               r_load=fabs((float)(mx3*mx1)/(mn3*mn1)-1.);

               surf_load= (float) (n1*NN3+n3*NN1)/(NN1*NN3);

               tot_load=r_load+fac_load*surf_load;

               divs_crit[5*num_ellig_crit]=(float) n1;
               divs_crit[5*num_ellig_crit+1]=(float) n3;	
               divs_crit[5*num_ellig_crit+2]=r_load;
               divs_crit[5*num_ellig_crit+3]=surf_load;
               divs_crit[5*num_ellig_crit+4]=tot_load;

               num_ellig_crit++;
          }

          //divs_crit[3*i]=r_load;
          //divs_crit[3*i+1]=surf_load;
          //divs_crit[3*i+2]=tot_load;
          //printf("N1=%d, N3= %d n1=%d n3=%d mn1=%d mx1=%d mn3=%d mx3=%d rload=%lf sload=%f\n",NN1,NN3,n1,n3,mn1,mx1,mn3,mx3,r_load,surf_load);
     }
     
     return num_ellig_crit;
}


int divisors(int numb, int mat[])
{
    int i,met=0;
    double met_rec;

    if (numb<=0){
        printf("error, divisors_routine are not made for zero or negative numbers!\n");
        exit(12);
    }
    
    for(i=1;i<=sqrt(numb);i++){
        if (numb % i ==0){
            mat[met]=i;
	    mat[met+1]=numb/i;
            met=met+2;
	}
    }

    met_rec=0.5*met;
    for(i=0;i<met_rec;i++){
        if(mat[2*i+1]!=mat[2*i]){ //so we will not the sqrt two times
            mat[met+2*i]=mat[2*i+1];
	    mat[met+2*i+1]=mat[2*i];
        }
        else met=met-1;
    }

    return met;
}
