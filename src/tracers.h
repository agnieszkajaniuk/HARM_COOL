#if !defined(TRACERS_H__INCLUDED_)
#define TRACERS_H__INCLUDED_

void init_tracers(void);
void update_tracers(double Dt);
void dump_tracers(void);
int restart_write_tracers (FILE * fp);
int restart_read_tracers (FILE * fp);

#endif
