#ifndef Cache4d_H

#define Cache4d_H

#define INI_MAX_THREAD 128

#pragma pack(1)

typedef enum {CONSTGRID=0, PREPARED} CACHEMODE;

typedef enum {LIN, SPLINE} INTERPMODE;

typedef enum {LOGARITHMIC=0, LINEAR, CUSTOM, NONUNIFORM, ROLLING} AXISTYPE;

typedef int (*CALCFUN)(double, double, double,  double *);

typedef struct {

    CACHEMODE mode;
    char  ver[13];          //"AGACACHE_2.0"

    int 	x;                //rozmiar X
    double	xmin;
    double  xmax;
    double  xrat;
    AXISTYPE xtype;

    int 	y;                //rozmiar Y
    double	ymin;
    double  ymax;
    double  yrat;
    AXISTYPE ytype;

    int 	z;                //rozmiar Z
    double	zmin;
    double  zmax;
    double  zrat;
    AXISTYPE ztype;

    int     d;                //liczba wartosci w punkcie xyz

    int     unused;
    int     autosave;         //Czestotliowsc zapisywania cache, liczba policzonych punktow
    int     wholecachecalculated;   //policzono caly cache

    double  *xcust;
    double  *ycust;
    double  *zcust;

    char   reserved[76];


}
CacheDesc;


struct proc_message_struct {

    long long idx;
    double 	  x;
    double 	  y;
    double    z;
};

typedef enum {READY=0, BUSY, DONE, INVALID} PROCSTATE;

class CacheWorker {

//to bedzie w pamieci dzielonej dlatego bez obiektow, dziedziczenia, ...
struct proc_type_struct {
    int proc_nr;
    PROCSTATE procstate;
    int pipe_read_child;
    int pipe_read_main;
    int pipe_write_child;
    int pipe_write_main;
};



private:
    /* No instantiation. */
    CacheWorker() {
    	shm_id_ = -1;
    	sem_id_ = -1;
		threadstarted = 0;
	}
    virtual ~CacheWorker();

public:
    static CacheWorker& getInstance()
    {
        /* Have a static local variable representing the unique instance.  Since
         * it's static, there is only one instance of this variable.  It's also only
         * initialized when getInstance is called.
         */
        static CacheWorker theInstance;
        return theInstance;
    }
public:

    int startThreads(int threads, CALCFUN calcfun_, int D_);
	void waitForOne();
	void waitForAll();

	const int getCount() {
		return threadstarted;
	}

	PROCSTATE getProcessState(int number)
	{
		if(number >= 0 && number < threadstarted)
			return proc_tab[number].procstate;
		return INVALID;
	}

	void setProcessState(int number, PROCSTATE state)
	{
		if(number >= 0 && number < threadstarted)
			proc_tab[number].procstate = state;
	}


	int getWorkerPipe(int number)
	{
		if(number >= 0 && number < threadstarted)
			return proc_tab[number].pipe_write_child;
		return -1;
	}

	int getPipe(int number)
	{
		if(number >= 0 && number < threadstarted)
			return proc_tab[number].pipe_read_main;
		return -1;
	}

private:
    void childMain(proc_type_struct *_proc_tab);


private:
    int   threadstarted;
    proc_type_struct *proc_tab;            //Tabela opisujaca procesy

	CALCFUN calcfun;
	int D;

    int   shm_id_;
    int   sem_id_;
};



class   Cache4d  {

public:
    int X;
    int Y;
    int Z;
    int D;
    CacheDesc cdesc;

    //Obiekty w pamieci dzielonej
    double *VectorX;
    double *VectorY;
    double *VectorZ;
    char   *MCounted;            //Trojwymiarowa mapa policzonych wartosci
//    long long MVStart;
    double *MValues;               //Trojwymiarowa macierz wartosci

    char filename[256];
    int  fd;
    CALCFUN calcfun;

private:

    int  needsave;

    void setindices(double *vect, int size, double min, double max, AXISTYPE type, double rat, double *vcust);

    void normalize(double *x, double *y, double *z);
    void calcbox(int size, double *vect, double val, INTERPMODE mode, int *left, int *right);


    int  initSharedMem();

    int checkIsOnGrid(double val, double *vector, int npoints);

    int  calculed;
	time_t last_save;

    int     ncached;
    double *cachedmem;

public:


    Cache4d(CacheDesc *desc, CALCFUN lcalcfun, const char *filename);  //opis, funkcja liczaca punkty, nazwa pliku z cachem
    virtual ~Cache4d();
    void test(void);

	int save();

    void calculateCache();  //Oblicza caly cache

    int  getvalues(double x, double y, double z, INTERPMODE mode, double *retval, double *derval = NULL, int zindic = -1, int *dimmask = NULL); //zwraca wartosci zinterpolowane

    int getivalues(int ix, int iy, int iz, double *retval, int onlycalc);

    int ImportFile(char *filename);  //iportuje dane z innego cache'a
    void ScanImport();  //scanuje katalog import w poszukiwania cachy do import
    void exportCache();

private:
    int load(char *filename);
    int save(char *filename);

    void calcCachePoints(struct proc_message_struct *messages, int count, int waitForFinish); //oblicza wybrane punkty

};

void Cache4d_childStop(int);
void PutOnConsole(const char *t1, const char *t2, const char *t3);

#pragma pack()

#endif

