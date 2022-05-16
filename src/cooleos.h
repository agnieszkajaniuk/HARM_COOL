#ifndef COOLEOS_H
#define COOLEOS_H
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

int SetupThreads_CoolEOS();
int SetupCache_CoolEOS();
int InitCache_CoolEOS(int hh, int ii, int jj);


double pressure_rho0_u_CoolEOS(double rho0, double u);

double dpdu_calc_CoolEOS( double rho0, double u );
double dpdrho_calc_CoolEOS( double rho0, double u );

void CoolVal_rho0_u_CoolEOS(double rho0, double u, int h, int i, int j, double coolvals[NCOOLVAL]);

int restart_write_CoolEOS(FILE *fp);
int restart_read_CoolEOS(FILE *fp);

void test_CoolEOS(double rho0);

extern const double RHO_UNIT;
extern const double U_UNIT;
extern const double P_UNIT;
extern const double T_UNIT;
extern const double L_UNIT;


#endif