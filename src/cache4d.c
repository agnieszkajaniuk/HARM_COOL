#include <stdio.h>  // Standard Input/Output routines
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>        // for floating point mathematics
#include <signal.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fcntl.h>
#include <dirent.h>
#include <sys/mman.h>
#if !defined __CYGWIN__
#include <sys/prctl.h>
#endif
#include <fcntl.h>
#include <time.h>

#ifndef WIN32
#include <sys/ipc.h>
#include <sys/msg.h>
#include <sys/sem.h>
#include <pwd.h>
#include <sys/shm.h>
#include <sys/wait.h>
#endif

#include "akima.h"
#include "cache4d.h"

#define  USE_SHARED_MEMORY 0

const int MaxCached = 128;

//#define TEST

void PutOnConsole(const char *t1, const char *t2, const char *t3) {

    printf("%s %s %s\n", t1, t2, t3);
}


Cache4d :: Cache4d(CacheDesc *desc, CALCFUN lcalcfun, const char *lfilename) {

    if(desc != NULL && access(lfilename, 0)) {
        memcpy(&cdesc, desc, sizeof(CacheDesc));

        strcpy(cdesc.ver, "AGACACHE_2.0");
        cdesc.wholecachecalculated = 0;

    } else {

        FILE *fcache;
        if((fcache = fopen(lfilename, "r+b"))==NULL) {

            printf("Cannot open file '%s' (errno %d %s)\n", lfilename, errno, strerror(errno));
            return;
        }

        if(fread(&cdesc, sizeof(CacheDesc), 1, fcache) != 1) {
            printf("Nie mozna odczytac 1. czesci cache'a\n");
            fclose(fcache);
            return;
        }

        fclose(fcache);

		cdesc.xcust = desc->xcust;
		cdesc.ycust = desc->ycust;
		cdesc.zcust = desc->zcust;
    }

    X = cdesc.x;
    Y = cdesc.y;
    Z = cdesc.z;
    D = cdesc.d;

    calcfun = lcalcfun;

    strcpy(filename, lfilename);

    //Inicjalizujemy pamiec dzielona
    fd = -1;
//    fid = -1;

    calculed = 0;

    needsave = 0;
    MCounted = NULL;

    VectorX = (double *)malloc(X * sizeof(double));
    VectorY = (double *)malloc(Y * sizeof(double));
    VectorZ = (double *)malloc(Z * sizeof(double));

    //Inicjalizacja

    setindices(VectorX, cdesc.x, cdesc.xmin, cdesc.xmax, cdesc.xtype, cdesc.xrat, cdesc.xcust);
    setindices(VectorY, cdesc.y, cdesc.ymin, cdesc.ymax, cdesc.ytype, cdesc.yrat, cdesc.ycust);
    setindices(VectorZ, cdesc.z, cdesc.zmin, cdesc.zmax, cdesc.ztype, cdesc.zrat, cdesc.zcust);


    if(initSharedMem())
    {
        exit(0);
    }




    //ladujemy
#if !USE_SHARED_MEMORY
        if(!access(filename, 0))
            if(!load(filename))  {
                ;//Cache4d_parentStop(0);
            }
#endif
}


Cache4d :: ~Cache4d() {

    if(cdesc.mode == CONSTGRID)
    {
#if USE_SHARED_MEMORY
        if (munmap(MValues,(long long)D*X*Y*Z*sizeof(double)) == -1)
        {
            int myerr = errno;
            printf("ERROR (main): munmap failed (errno %d %s)\n", myerr, strerror(myerr));
        }

        if (munmap(MCounted,X*Y*Z*sizeof(char)) == -1)
        {
            int myerr = errno;
            printf("ERROR (main): munmap failed (errno %d %s)\n", myerr, strerror(myerr));
        }
        close(fd);
#else
        save();
		free(MValues);
    	free(MCounted);
#endif
    }

    free(VectorX);
    free(VectorY);
    free(VectorZ);

}


int Cache4d :: load(char *filename) {

#if !USE_SHARED_MEMORY

    FILE *fcache;

    CacheDesc  fcachedesc;

    printf("Odczytujemy dane z cache'a...\n");

    if((fcache = fopen(filename, "r+b"))==NULL) {

        printf("Cannot open file %s\n", filename);
        return 0;
    }
    if(fread(&fcachedesc, sizeof(CacheDesc), 1, fcache) != 1) {
        printf("Nie mozna odczytac 1. czesci cache'a\n");
        fclose(fcache);
        return 0;
    }

    //Sprawdzamy czy naglowek z pliku pasuje do danych z pamieci

    if(cdesc.x != fcachedesc.x || cdesc.xmin != fcachedesc.xmin || cdesc.xmax != fcachedesc.xmax || cdesc.xtype != fcachedesc.xtype ||
            cdesc.y != fcachedesc.y || cdesc.ymin != fcachedesc.ymin || cdesc.ymax != fcachedesc.ymax || cdesc.ytype != fcachedesc.ytype ||
            cdesc.z != fcachedesc.z || cdesc.zmin != fcachedesc.zmin || cdesc.zmax != fcachedesc.zmax || cdesc.ztype != fcachedesc.ztype  ) {
        printf("Definicja cache'a nie pasuje do pliku!\n");
        fclose(fcache);
        exit(0);
        return 0;
    }


    if(fread(VectorX, X*sizeof(double), 1, fcache) != 1) {
        printf("Nie mozna odczytac 2. czesci cache'a\n");
        fclose(fcache);
        return 0;
    }

    if(fread(VectorY, Y*sizeof(double), 1, fcache) != 1) {
        printf("Nie mozna odczytac 3. czesci cache'a\n");
        fclose(fcache);
        return 0;
    }

    if(fread(VectorZ, Z*sizeof(double), 1, fcache) != 1) {
        printf("Nie mozna odczytac 4. czesci cache'a\n");
        fclose(fcache);
        return 0;
    }

    if(fread(MValues, X*Y*Z*D*sizeof(double), 1, fcache) != 1) {
        printf("Nie mozna odczytac 5. czesci cache'a\n");
        fclose(fcache);
        return 0;
    }

    if(strcmp(fcachedesc.ver, "AGACACHE_1.0") == 0) {

        int *tmp = new int[X*Y*Z];

        if(fread(tmp, X*Y*Z*sizeof(int), 1, fcache) != 1) {
            printf("Nie mozna odczytac 6. czesci cache'a\n");
            fclose(fcache);
            return 0;
        }

        for(int i = 0; i<X*Y*Z; ++i)
            MCounted[i] = tmp[i];

        delete []tmp;

    } else {
        if(fread(MCounted, X*Y*Z*sizeof(char), 1, fcache) != 1) {
            printf("Nie mozna odczytac 6. czesci cache'a\n");
            fclose(fcache);
            return 0;
        }
    }

    fclose(fcache);
    int total = 0;
    for(int iz = 0; iz < Z; ++iz)
        for(int ix = 0; ix < X; ++ix)
            for(int iy = 0; iy < Y; ++iy) {
                int idx = iz*(X*Y) + ix*Y + iy;
                if(MCounted[idx] == 1) {
                    ++total;
                }
            }



    if(total == X*Y*Z)
        cdesc.wholecachecalculated = 1;

    printf("Grid: %dx%dx%d\n", X, Y, Z);
    printf("X: %E - %E\n", cdesc.xmin, cdesc.xmax);
    printf("Y: %E - %E\n", cdesc.ymin, cdesc.ymax);
    printf("Z: %E - %E\n", cdesc.zmin, cdesc.zmax);

    printf("Total %d points of %d, whl=%d\n", total, X*Y*Z, cdesc.wholecachecalculated);

    printf("Koniec czytania cache'a.\n");
#endif
	last_save = time(NULL);

    return 1;
}


int Cache4d :: save(char *filename) {

#if !USE_SHARED_MEMORY
    FILE *fcache;

    printf("Zapisujemy dane do cache'a '%s'\n", filename);

    if((fcache = fopen(filename, "w+b"))==NULL) {

        printf("Cannot open file %s\n", filename);
        return 0;
    }

    if(fwrite(&cdesc, sizeof(CacheDesc), 1, fcache) != 1) {
        printf("Nie mozna zapisac 1. czesci cache'a\n");
        fclose(fcache);
        return 0;
    }

    if(fwrite(VectorX, X*sizeof(double), 1, fcache) != 1) {
        printf("Nie mozna zapisac 2. czesci cache'a\n");
        fclose(fcache);
        return 0;
    }

    if(fwrite(VectorY, Y*sizeof(double), 1, fcache) != 1) {
        printf("Nie mozna zapisac 3. czesci cache'a\n");
        fclose(fcache);
        return 0;
    }

    if(fwrite(VectorZ, Z*sizeof(double), 1, fcache) != 1) {
        printf("Nie mozna zapisac 4. czesci cache'a\n");
        fclose(fcache);
        return 0;
    }

    if(fwrite(MValues, X*Y*Z*D*sizeof(double), 1, fcache) != 1) {
        printf("Nie mozna zapisac 5. czesci cache'a\n");
        fclose(fcache);
        return 0;
    }

    if(fwrite(MCounted, X*Y*Z*sizeof(char), 1, fcache) != 1) {
        printf("Nie mozna zapisac 6. czesci cache'a\n");
        fclose(fcache);
        return 0;
    }

    fclose(fcache);

    printf("Koniec pisania cache'a.\n");
#endif
    return 1;
}

int Cache4d :: save() {
	int ret = 1;
    if(needsave) {
       ret = save(filename);
       needsave = 0;
    }
	last_save = time(NULL);
	return ret;
}

void Cache4d :: setindices(double *vect, int size, double min, double max, AXISTYPE type, double rat, double *vcust) {

    double delta;

    int r;

    switch(type) {
    case LOGARITHMIC:
        delta = (log10(max) - log10(min))/(size - 1);

        for(r = 0 ; r < size; r++)
            vect[r] = pow(10.0, log10(min) + delta*r);

        break;
    case LINEAR:
        delta = (max - min)/(size - 1);

        for(r = 0 ; r < size; r++)
            vect[r] = min + delta*r;

        break;
    case NONUNIFORM:
        delta = (max - min)*(rat - 1.0)/(pow(rat, size) - 1.0)*rat;

        vect[0] = min;

        for(r = 1 ; r < size; r++) {
            vect[r] = vect[r-1] + delta;
            delta *= rat;
        }

        break;
    case CUSTOM:
        if(size == 1)
            vect[0] = min;
        else
            memcpy(vect, vcust, size*sizeof(double));
        break;

    case ROLLING:
        vect[size/2] = (min + max)/2.0;
        for(r=size/2+1; r<size; r++)
            vect[r] = vect[r-1]*(1.0 + rat);
        for(r=size/2-1; r>=0; r--)
            vect[r] = vect[r+1]*(1.0 - rat);
        break;

    default:
        break;
    }

    /*
    	printf("Vector size = %d type=%d\n", size, type);

       	for(r = 1 ; r < size; r++)
    	{
    		printf("r=%3d %E, %E -> %E\n", r, vect[r-1], vect[r], (vect[r]/vect[r-1] - 1.0)*100.0);
    	}
    */
}


void Cache4d :: normalize(double *x, double *y, double *z) {
    //Wywalamy co wylazi poza zakres

    if(*x < VectorX[0])
        *x = VectorX[0];
    if(*x > VectorX[X-1])
        *x = VectorX[X-1];
    if(*y < VectorY[0])
        *y = VectorY[0];
    if(*y > VectorY[Y-1])
        *y = VectorY[Y-1];
    if(*z < VectorZ[0])
        *z = VectorZ[0];
    if(*z > VectorZ[Z-1])
        *z = VectorZ[Z-1];
}

void Cache4d :: calcbox(int size, double *vect, double val, INTERPMODE mode, int *left, int *right) {

    int t_left;
    int t_right;

    if(mode == LIN)
    {
        dvec_bracket (size, vect, val, &t_left, &t_right);
        *left = t_left - 1;
        *right = t_right;
    }
    else
    {
        dvec_bracket (size, vect, val, &t_left, &t_right);
        t_left--;
        t_right--;

        t_left -= 1;

//        t_left -= 2;

        if(t_left < 0) t_left = 0;

        t_right += 1;

//        t_right += 2;
        if(t_right >= size) t_right = size-1;

        *left = t_left;
        *right = t_right;
    }
}




int Cache4d :: initSharedMem()
{
#if USE_SHARED_MEMORY

    size_t page_size = (size_t) sysconf (_SC_PAGESIZE);

    long long shmsize = sizeof(CacheDesc);

    if(cdesc.mode == CONSTGRID)
    {
        shmsize += (page_size - sizeof(CacheDesc)%page_size) +
                   X*Y*Z*sizeof(char) + (page_size - (X*Y*Z*sizeof(char))%page_size) +
                   (long long)X*Y*Z*D*sizeof(double);
    }

    int fileExists = (access(filename, 0) == 0);


    /* * Create file we are using for mmap. The file must be * size of memory we wish to map. */
    fd = open(filename, O_RDWR| O_CREAT, S_IRUSR| S_IWUSR );

    if (fd == -1)
    {
        int myerr = errno;
        printf("ERROR: open %s failed (errno %d %s)\n", filename, myerr, strerror(myerr));
        return EXIT_FAILURE;
    }

    if(!fileExists)
    {
        if(write(fd, &cdesc, sizeof(CacheDesc)) != sizeof(CacheDesc)) {
            printf("Nie mozna zapisac 1. czesci cache'a\n");
            close(fd);
            return EXIT_FAILURE;
        }

        if (lseek(fd, shmsize, SEEK_SET) == -1)
        {
            int myerr = errno;
            printf("ERROR: lseek failed (errno %d %s)\n", myerr, strerror(myerr));
            return EXIT_FAILURE;
        }
        write(fd, "", 1);
    }

    if(cdesc.mode == CONSTGRID)
    {
        //Wykrawamy z duzego kawalka pamieci dzielonej
        long long offest = sizeof(CacheDesc) + (page_size - sizeof(CacheDesc)%page_size);
        int prot = (PROT_READ| PROT_WRITE);
        int flags = MAP_SHARED;
        //int flags = MAP_PRIVATE|MAP_POPULATE;

//        printf("offest=%lld allocation granularity=%ld\n", offest, page_size);

        MCounted = (char *)mmap(NULL, X*Y*Z*sizeof(char), prot, flags, fd, offest);
        if (MCounted == MAP_FAILED)
        {
            int myerr = errno;
            printf("ERROR (main): mmap MCounted failed (errno %d %s)\n", myerr, strerror(myerr));

            exit(0);
        }

        madvise(MCounted, X*Y*Z*sizeof(char), MADV_RANDOM|MADV_WILLNEED);

        long long MVStart = sizeof(CacheDesc) + (page_size - sizeof(CacheDesc)%page_size) +
                            X*Y*Z*sizeof(char) + (page_size - (X*Y*Z*sizeof(char))%page_size);

        MValues = (double *)mmap(NULL, (long long)X*Y*Z*D*sizeof(double), prot, flags, fd, MVStart);
        if (MValues == MAP_FAILED)
        {
            int myerr = errno;
            printf("ERROR (main): mmap MValues failed (errno %d %s)\n", myerr, strerror(myerr));
            exit(0);
        }
        madvise(MValues, (long long)X*Y*Z*D*sizeof(double), MADV_RANDOM|MADV_WILLNEED);
    }
#else
    MValues = (double *)malloc(X*Y*Z*D*sizeof(double));
    MCounted = (char *)malloc(X*Y*Z*sizeof(char));
	memset(MValues, 0, X*Y*Z*D*sizeof(double));
	memset(MCounted, 0, X*Y*Z*sizeof(char));

//    long long idx = 0*(X*Y) + 24*Y + 32;
//    printf("initSharedMem %d %d %d ? '0x%x' idx=%lld, MCounted=%p MCounted[idx]=%p\n", X, Y, Z, MCounted[idx], idx, MCounted, &MCounted[idx]);

#endif

    if(cdesc.mode == PREPARED)
    {
        cachedmem = (double *)malloc((3+D)*MaxCached*sizeof(double));
		memset(cachedmem, 0, (3+D)*MaxCached*sizeof(double));
        ncached = 0;
    }


    return 0;
}


void Cache4d :: calculateCache() {

#define MAXCALCPOINTS 200
    struct proc_message_struct messages[MAXCALCPOINTS];


    needsave = 0;

    int counter = 0;


    printf("Calculating cache...\n");

    for(int iz = 0; iz < Z; ++iz)
        for(int ix = 0; ix < X; ++ix)
            for(int iy = 0; iy < Y; ++iy) {
                long long idx = iz*(X*Y) + ix*Y + iy;


                if(MCounted[idx] == 0)  //Doliczamy brakujacy w cache'u punkt
                {

//                    printf("calculateCache: %4d %4d %4d\n", ix, iy, iz);
                    messages[counter].x = VectorX[ix];
                    messages[counter].y = VectorY[iy];
                    messages[counter].z = VectorZ[iz];
                    messages[counter].idx = idx;

                    // printf("after calculateCache: %4d %4d %4d\n", messages[counter].ix, messages[counter].iy, messages[counter].iz);

                    ++counter;

                    if(counter >= MAXCALCPOINTS) {

                        needsave = 1;
                        calcCachePoints(messages, counter, 0);
                        counter = 0;
                    }
                } else {
                    //                                           printf("Mam policzone     %d %d %d\n", ix, iy, iz);


                }
            }


    calcCachePoints(messages, counter, 1);


    save();


    printf("Caly cache policzony \n");

}

void Cache4d :: calcCachePoints(struct proc_message_struct *messages, int count, int waitForFinish) {


    ncached = 0;
    //Obliczamy caly brakujacy cache
//    printf("Doliczam punkty w cache'u: %d\n", count);

    for(int curmsg = 0; curmsg < count;) {

//        printf("zleceam policzenie %d punktu, threads %d\n", curmsg+1, threadstarted);

#ifndef WIN32

        //obsluga semaforow, czekamy az zwolni sie ktorys z procesow
        if (CacheWorker::getInstance().getCount() > 0) {

            int watek;

//            printf("ktorys watek jest wolny/skonczyl\n");

            for (watek = 0; watek < CacheWorker::getInstance().getCount(); watek++) {

                if (CacheWorker::getInstance().getProcessState(watek) == DONE)
                {
//                    printf("watek %d jest DONE\n", watek);

                    if(cdesc.mode == CONSTGRID)
                    {
                        long long idx;
                        if (read(CacheWorker::getInstance().getPipe(watek), &idx, sizeof(long long)) != sizeof(long long))
                        {
                            PutOnConsole("Blad odczytu pipe 1", "", "");
                            CacheWorker::getInstance().setProcessState(watek, INVALID);
                            continue;
                        }

//                        printf("czytam dane punku %lld\n", idx);
                        int bytesToRead = D*sizeof(double);
                        if (read(CacheWorker::getInstance().getPipe(watek), &MValues[idx*D], bytesToRead) != bytesToRead)
                        {
                            printf("Blad odczytu pipe 2 %d - %s, bytesToRead=%d, idx=%lld\n", errno, strerror(errno), bytesToRead, idx);
                            CacheWorker::getInstance().setProcessState(watek, INVALID);
                            continue;
                        }

                        MCounted[idx] = 1;
				        ++calculed;
						needsave = 1;
                    }
                    else if(cdesc.mode == PREPARED)
                    {

                        if(ncached + 1 < MaxCached)
                        {
                            long long idx;
                            if (read(CacheWorker::getInstance().getPipe(watek), &idx, sizeof(long long)) != sizeof(long long))
                            {
                                PutOnConsole("Blad odczytu pipe 3", "", "");
                                CacheWorker::getInstance().setProcessState(watek, INVALID);
                                continue;
                            }

                            int bytesToRead = (3+D)*sizeof(double);
                            if (read(CacheWorker::getInstance().getPipe(watek), &cachedmem[ncached*(3+D)], bytesToRead) != bytesToRead)
                            {
                                PutOnConsole("Blad odczytu pipe 4", "", "");
                                CacheWorker::getInstance().setProcessState(watek, INVALID);
                                continue;
                            }
                            ncached++;
                        }
                    }
                    CacheWorker::getInstance().setProcessState(watek, READY);
                }


                if (CacheWorker::getInstance().getProcessState(watek) == READY)
                    break;
            }

            if (watek >= CacheWorker::getInstance().getCount()) {

				CacheWorker::getInstance().waitForOne();

//				PutOnConsole("PANIC!!! brak wolnego watku, semafor podniesiony", "", "");
                //Cache4d_parentStop(0);

                continue;
            }

            CacheWorker::getInstance().setProcessState(watek, BUSY); //jeszcze sie nie skonczyl
//            printf("punkt %lld zlecony watkowi %d\n", messages[curmsg].idx, watek);

            if (write(CacheWorker::getInstance().getWorkerPipe(watek), &messages[curmsg], sizeof(struct proc_message_struct)) != sizeof(struct proc_message_struct)) {
                PutOnConsole("Blad zapisu do pipe", "", "");
                CacheWorker::getInstance().setProcessState(watek, INVALID);
                continue;
            }
        } else {  //Liczymy w jednym (glownym) watku
#endif

            //printf("licze w jednym watku %d %d %d\n", messages[curmsg].ix, messages[curmsg].iy, messages[curmsg].iz);

            //            printf("%E, %E, %E\n",    VectorX[(messages[curmsg]).ix], VectorY[(messages[curmsg]).iy], VectorZ[(messages[curmsg]).iz]);
            /*
            int prot = (PROT_READ| PROT_WRITE);
            int flags = MAP_SHARED;
            long long offest = MVStart + (messages[curmsg]).idx*D;

            double *MValues = (double *)mmap(NULL, D*sizeof(double), prot, flags, fd, offest);
            if (MValues == MAP_FAILED)
            {
            	int myerr = errno;
                printf("ERROR (parent): mmap failed (errno %d %s)\n", myerr, strerror(myerr));
                exit(0);
            }
            */
            if(cdesc.mode == CONSTGRID)
            {
                calcfun(messages[curmsg].x, messages[curmsg].y, messages[curmsg].z, &MValues[messages[curmsg].idx*D]);
                MCounted[messages[curmsg].idx] = 1;
		        ++calculed;
				needsave = 1;
            }
            else if(cdesc.mode == PREPARED)
            {

                double mvalues[D];

                calcfun(messages[curmsg].x, messages[curmsg].y, messages[curmsg].z, mvalues);

                if(ncached + 1 < MaxCached)
                {
                    int bytesToRead = (3+D)*sizeof(double);
                    cachedmem[ncached*(3+D)] = messages[curmsg].x;
                    cachedmem[ncached*(3+D)+1] = messages[curmsg].y;
                    cachedmem[ncached*(3+D)+2] = messages[curmsg].z;
                    memcpy(&cachedmem[ncached*(3+D)+3], mvalues, D*sizeof(double));
                    ncached++;
                }
            }
            /*
            if (munmap(MValues,D*sizeof(double)) == -1)
            {
            	int myerr = errno;
            	printf("ERROR (main): munmap failed (errno %d %s)\n", myerr, strerror(myerr));
            }
            */

#ifndef WIN32
        }
#endif
		++curmsg;

//        printf(".");
    }


    //Czekamy az wszystkie sie policza
#ifndef WIN32

    if (CacheWorker::getInstance().getCount() > 0 && waitForFinish) {

        //Czekamy az sie policza wszystkie watki
//        CacheWorker::getInstance().waitForAll();

        for (int watek = 0; watek < CacheWorker::getInstance().getCount();)
        {
            if (CacheWorker::getInstance().getProcessState(watek) == BUSY)
			{
				CacheWorker::getInstance().waitForOne();
				continue;
			}
			

            if (CacheWorker::getInstance().getProcessState(watek) == DONE)
            {
                if(cdesc.mode == CONSTGRID)
                {
                    long long idx;
                    if (read(CacheWorker::getInstance().getPipe(watek), &idx, sizeof(long long)) != sizeof(long long))
                    {
                        PutOnConsole("Blad odczytu pipe 5", "", "");
                        CacheWorker::getInstance().setProcessState(watek, INVALID);
                        continue;
                    }

                    int bytesToRead = D*sizeof(double);
                    if (read(CacheWorker::getInstance().getPipe(watek), &MValues[idx*D], bytesToRead) != bytesToRead)
                    {
                        printf("Blad odczytu pipe 6 %d - %s, bytesToRead=%d\n", errno, strerror(errno), bytesToRead);
                        CacheWorker::getInstance().setProcessState(watek, INVALID);
                        continue;
                    }

                    MCounted[idx] = 1;
			        ++calculed;
					needsave = 1;
                }
                else if(cdesc.mode == PREPARED)
                {

                    if(ncached + 1 < MaxCached)
                    {
                        long long idx;
                        if (read(CacheWorker::getInstance().getPipe(watek), &idx, sizeof(long long)) != sizeof(long long))
                        {
                            PutOnConsole("Blad odczytu pipe 7", "", "");
                            CacheWorker::getInstance().setProcessState(watek, INVALID);
                            continue;
                        }

                        int bytesToRead = (3+D)*sizeof(double);
                        if (read(CacheWorker::getInstance().getPipe(watek), &cachedmem[ncached*(3+D)], bytesToRead) != bytesToRead)
                        {
                            PutOnConsole("Blad odczytu pipe 8", "", "");
                            CacheWorker::getInstance().setProcessState(watek, INVALID);
                            continue;
                        }
                        ncached++;
                    }
                }
                CacheWorker::getInstance().setProcessState(watek, READY);
            }
			++watek;
        }

    }
#endif

    //        printf("Calculated point no: %d\n", calculed);
    if(calculed >= cdesc.autosave) {
         save();
         calculed = 0;
     }

//    printf("\n");

    return;
}


void Cache4d_childStop(int code) {

    exit(code);
}

int Cache4d :: checkIsOnGrid(double val, double *vector, int npoints) {

    for (int i = 0; i  < npoints; ++i) {

        //printf("val = %ld, i = %3d, g= %ld\n", val, i, vector[i]);
        if(vector[i] == val)
            return i;
        if(vector[i] > val)
            return -1;
    }

    return -1;
}


int Cache4d :: getvalues(double x, double y, double z, INTERPMODE mode, double *retval, double *dervaly, int zindic, int *dimmask) {

//	printf("getvalues %E %E %E\n", x, y, z);
//	printf("Cache4d getvalues %p\n", this);

    if(cdesc.mode == PREPARED)
    {
        while(1)
        {
            for(int i = 0; i <ncached; i++)
            {
                double *cacheitem = &cachedmem[i*(3+D)];

                if(cacheitem[0] == x && cacheitem[1] == y && cacheitem[2] == z)
                {
                    memcpy(retval, &cacheitem[3], D*sizeof(double));
                    return 1;
                }
            }
            printf("Nie ma wartosci w cachu dla x=%E, y=%E, z=%E, ncached=%d\n", x, y, z, ncached);
            struct proc_message_struct messages;
            messages.x = x;
            messages.y = y;
            messages.z = z;
            calcCachePoints(&messages, 1, 1);
        }
    }
    else if(cdesc.mode == CONSTGRID)
    {
        bool rolling = false;
        /*
            int prot = (PROT_READ| PROT_WRITE);
            int flags = MAP_SHARED;
        */
        int min_x, max_x, min_y, max_y, min_z, max_z;

        if(cdesc.xtype == ROLLING)
        {
			if(x < VectorX[1] || x > VectorX[X-2])
            	rolling = true;
        }
        else if(x < VectorX[0] || x > VectorX[X-1])
        {
                printf("getvalues(): poza zakresem\n");
                printf("%E %E %E\n", x, y, z);
                return 0;
        }


        if(cdesc.ytype == ROLLING)
        {
			if(y < VectorY[1] || y > VectorY[Y-2])
            	rolling = true;
        }
        else if(y < VectorY[0] || y > VectorY[Y-1])
        {
                printf("getvalues(): poza zakresem\n");
                printf("%E %E %E\n", x, y, z);
                return 0;
        }

        if(cdesc.ztype == ROLLING)
        {
			if(z < VectorZ[1] || z > VectorZ[Z-2])
            	rolling = true;
        }
		else if(zindic == -1 && (z < VectorZ[0] || z > VectorZ[Z-1]))
        {
            if(cdesc.ztype == ROLLING)
                rolling = true;
            else
            {
                printf("getvalues(): poza zakresem\n");
                printf("%E %E %E\n", x, y, z);
                return 0;
            }
        }


                                                  
        if(rolling)
        {
            memset(MCounted, 0, X*Y*Z*sizeof(char));

            printf("Rolling cache %s...\n", filename);
            if(cdesc.xtype == ROLLING)
            {
                cdesc.xmin = cdesc.xmax = x;
                setindices(VectorX, X, x, x, cdesc.xtype, cdesc.xrat, cdesc.xcust);
            }
            if(cdesc.ytype == ROLLING)
            {
                cdesc.ymin = cdesc.ymax = y;
                setindices(VectorY, Y, y, y, cdesc.ytype, cdesc.yrat, cdesc.ycust);
            }
            if(cdesc.ztype == ROLLING)
            {
                cdesc.zmin = cdesc.zmax = z;
                setindices(VectorZ, Z, z, z, cdesc.ztype, cdesc.zrat, cdesc.zcust);
            }

            lseek(fd, 0l, SEEK_SET);
            if(write(fd, &cdesc, sizeof(CacheDesc))<=0)
			{
				;
			}

        }

        //Ustalamy obszar potrzebny do obliczen
        calcbox(X, VectorX, x, mode, &min_x, &max_x);
        calcbox(Y, VectorY, y, mode, &min_y, &max_y);

        //Sprawdzamy czy punkt nie lezy dokladnie na na siatce z
        if(zindic != -1)
            max_z = min_z = zindic;
        else {
            //       printf("%lf nie jest na siatce wspolrzednych Z\n", z);
            calcbox(Z, VectorZ, z, mode, &min_z, &max_z);
        }

        //Sprawdzamy i ewentualnie doliczamy cache

        int newxsize = max_x - min_x + 1;
        int newysize = max_y - min_y + 1;
        int newzsize = max_z - min_z + 1;

//    	printf("x: %d %d, pts=%d\n", min_x, max_x, newxsize);
//	    printf("y: %d %d, pts=%d\n", min_y, max_y, newysize);
//	    printf("z: %d %d, pts=%d\n", min_z, max_z, newzsize);


        //printf("Cache: %s - %d\n", filename, cdesc.wholecachecalculated);

        if(!cdesc.wholecachecalculated) {
            int counter = 0;

            struct proc_message_struct messages[newxsize*newysize*newzsize];

            for(int iz = min_z; iz <= max_z; ++iz)
                for(int ix = min_x; ix <= max_x; ++ix)
                    for(int iy = min_y; iy <= max_y; ++iy) {
                        long long idx = iz*(X*Y) + ix*Y + iy;
						char w = MCounted[idx];

                        if(MCounted[idx] == 0)  //Doliczamy brakujacy w cache'u punkt
                        {
//                            printf("cache ok     value %d %d %d ? 0x%x (0x%x) sizeof=%d\n", ix, iy, iz, MCounted[idx], w, sizeof(MCounted[idx]));

                            messages[counter].x = VectorX[ix];
                            messages[counter].y = VectorY[iy];
                            messages[counter].z = VectorZ[iz];
                            messages[counter].idx = idx;

                            ++counter;
                        } else if(MCounted[idx] != 1) {    //nie ma wartosci
                            printf("cache broken value %d %d %d ? 0x%x (0x%x) sizeof=%lu\n", ix, iy, iz, MCounted[idx], w, sizeof(MCounted[idx]));
                            return 0;
                        }


                        /*else {

                        for(int i = 0;i <D;++i)
                        printf("chache val[%d]=%E\n", i, *(&MValues[idx*D]+i));
                        }*/
                    }

            //Doliczamy brakujace punkty
            if(counter > 0)  {
                //            printf("Doliczam brakujace %d punktow...\n", counter);
                calcCachePoints(messages, counter, 1);
            }
        }



        double ztmp[newzsize*D];
        double xtmp[newxsize*D];
        double xtmp_dery[newxsize*D];

        double *pout;

        if(zindic == -1)
            pout = ztmp;
        else
            pout = retval;

        if(mode == LIN)
        {
            for(int i = min_z, cnt1=0; i <= max_z; ++i, ++cnt1)
            {
                for(int j = min_x, cnt2=0; j <= max_x; ++j, ++cnt2)
                {
                    long long idx = (long long)D*(i*X*Y+j*Y + min_y);
                    /*
                    			    long long offest = MVStart + (long long)D*(i*X*Y+j*Y + min_y);
                    			    double *MValues = (double *)mmap(NULL, newysize*D*sizeof(double), prot, flags, fd, offest);
                    			    if (MValues == MAP_FAILED)
                    			    {
                        	    		int myerr = errno;
                    			        printf("ERROR (main): mmap failed (errno %d %s)\n", myerr, strerror(myerr));
                    			        exit(0);
                    			    }
                    */
                    spline_linear_val2(D, newysize, VectorY+min_y, &MValues[idx], y, &xtmp[cnt2*D]);
                    /*
                    				if (munmap(MValues,newysize*D*sizeof(double)) == -1)
                    				{
                    					int myerr = errno;
                    					printf("ERROR (main): munmap failed (errno %d %s)\n", myerr, strerror(myerr));
                    				}
                    */
                }
                spline_linear_val2(D, newxsize, VectorX+min_x, xtmp, x, &pout[cnt1*D]);
            }

            // Do the final interpolation
            if(zindic == -1)
                spline_linear_val2(D, newzsize, VectorZ+min_z, ztmp, z, retval);
        }
        else
        {
            for(int i = min_z, cnt1=0; i <= max_z; ++i, ++cnt1)
            {
                for(int j = min_x, cnt2=0; j <= max_x; ++j, ++cnt2)
                {
                    long long idx = (long long)D*(i*X*Y+j*Y + min_y);
                    /*
                    			    long long offest = MVStart + (long long)D*(i*X*Y+j*Y + min_y);

                    			    double *MValues = (double *)mmap(NULL, newysize*D*sizeof(double), prot, flags, fd, offest);
                    			    if (MValues == MAP_FAILED)
                    			    {
                        	    		int myerr = errno;
                    			        printf("ERROR (main): mmap failed (errno %d %s)\n", myerr, strerror(myerr));
                    			        exit(0);
                    			    }
                    */
//                    spline_overhauser_val3(D, newysize, VectorY+min_y, &MValues[idx], y, &xtmp[cnt2*D]);
					if(dervaly)
    	                spline_akima_der(D, dimmask, newysize, VectorY+min_y, &MValues[idx], y, &xtmp[cnt2*D], &xtmp_dery[cnt2*D]);
					else
    	                spline_akima(D, dimmask, newysize, VectorY+min_y, &MValues[idx], y, &xtmp[cnt2*D]);
                    /*
                    				if (munmap(MValues,newysize*D*sizeof(double)) == -1)
                    				{
                    					int myerr = errno;
                    					printf("ERROR (main): munmap failed (errno %d %s)\n", myerr, strerror(myerr));
                    				}
                    */
                }
//                spline_overhauser_val3(D, newxsize, VectorX+min_x, xtmp, x, &pout[cnt1*D]);
                spline_akima(D, dimmask, newxsize, VectorX+min_x, xtmp, x, &pout[cnt1*D]);
				if(dervaly)
	                spline_akima(D, dimmask, newxsize, VectorX+min_x, xtmp_dery, x, dervaly);
            }

            // Do the final interpolation
            if(zindic == -1)
            {
//                spline_overhauser_val3(D, newzsize, VectorZ+min_z, ztmp, z, retval);
                	spline_akima(D, dimmask, newzsize, VectorZ+min_z, ztmp, z, retval);
            }
        }
    }

	if(needsave && time(NULL) > last_save + 15*60)
		save();

    return 1;
}

int Cache4d :: getivalues(int ix, int iy, int iz, double *retval, int onlycalc) {


    long long idx = iz*(X*Y) + ix*Y + iy;


    if(!cdesc.wholecachecalculated && MCounted[idx] == 0)  //Doliczamy brakujacy w cache'u punkt
    {
        if(onlycalc)
            return 0;

        int counter = 0;

        struct proc_message_struct messages[1];

        messages[counter].x = VectorX[ix];
        messages[counter].y = VectorY[iy];
        messages[counter].z = VectorZ[iz];
        messages[counter].idx = idx;

        ++counter;
        calcCachePoints(messages, counter, 1);

    }
    /*
        memcpy(retval, &MValues[idx*D], D*sizeof(double));
    */
    return 1;

}


void Cache4d :: test(void) {

    /* int i;

     printf("----------------------------------------------------------\n");

     printf("Zakres x: %E - %E\n", cdesc.xmin, cdesc.xmax);

     for(i = 0; i < X; ++i)
         printf("%d  %E\n", i, VectorX[i]);

     printf("----------------------------------------------------------\n");
     */
    if(cdesc.mode == CONSTGRID)
    {

        int total = 0;
        int totalok = 0;

        for(int iz = 0; iz < Z; ++iz)
            for(int ix = 0; ix < X; ++ix)
                for(int iy = 0; iy < Y; ++iy) {
                    long long idx = iz*(X*Y) + ix*Y + iy;

                    if(MCounted[idx] != 0) {
                        ++total;
                    }

                    if(MCounted[idx] == 1) {
                        ++totalok;
                    }

                }
        if(totalok == X*Y*Z)
            cdesc.wholecachecalculated = 1;

        printf("total: %d\n", total);
    }
}

void Cache4d :: ScanImport() {
    char buf[512];
    DIR   *dirptr;
    struct dirent *direntptr;

    dirptr = opendir("import");
    if (dirptr==NULL) {
        return;
    }

    char *fname;

    fname = strrchr(filename, '/');

    if(fname != NULL)
        fname = fname + 1;
    else
        fname = filename;

    while((direntptr=readdir(dirptr))!=NULL) {

        if(!strcmp(direntptr->d_name, fname))  //Szukamy pliku o tej samej nazwie jak nasz cache
        {
            sprintf(buf, "import/%s", direntptr->d_name);

            if(ImportFile(buf)) {  //Poprawnie zaimportowane

                char buf2[256];

                sprintf(buf2, "%s_imp", buf);

                if(rename(buf, buf2) != 0) {
                    printf("Nie mozna zmienic nazwy pliku %s na %s\n", buf, buf2);
                    return;

                }
            }
        }
    }

    closedir(dirptr);

}


int Cache4d :: ImportFile(char *lfilename) {

#if 0
    FILE *fcache;

    CacheDesc  fcachedesc;



    printf("Odczytujemy dane z pliku %s...\n", lfilename);

    if((fcache = fopen(lfilename, "r+b"))==NULL) {

        printf("Cannot open file %s\n", lfilename);
        return 0;
    }

    if(fread(&fcachedesc, sizeof(CacheDesc), 1, fcache) != 1) {
        printf("Nie mozna odczytac 1. czesci cache'a\n");
        fclose(fcache);
        return 0;
    }

    if(strcmp(fcachedesc.ver, "AGACACHE_1.1") != 0)
    {
        printf("Nieprawidlowa versja cache'a: %s w pliku %s\n", fcachedesc.ver, lfilename);
        fclose(fcache);
        return 0;
    }


    //Sprawdzamy czy naglowek z pliku pasuje do danych z pamieci

    if(cdesc.x != fcachedesc.x || cdesc.xmin != fcachedesc.xmin || cdesc.xmax != fcachedesc.xmax || cdesc.xtype != fcachedesc.xtype ||
            cdesc.y != fcachedesc.y || cdesc.ymin != fcachedesc.ymin || cdesc.ymax != fcachedesc.ymax || cdesc.ytype != fcachedesc.ytype ||
            cdesc.z != fcachedesc.z || cdesc.zmin != fcachedesc.zmin || cdesc.zmax != fcachedesc.zmax || cdesc.ztype != fcachedesc.ztype  ) {
        printf("Definicja cache'a nie pasuje do pliku %s!\n", lfilename);
//        fclose(fcache);
//        return 0;
    }

    double *lVectorX = new double[X];
    double *lVectorY = new double[Y];
    double *lVectorZ = new double[Z];
    double *lMValues = new double[X*Y*Z*D];               //Trojwymiarowa macierz wartosci
    char   *lMCounted = new char[X*Y*Z];            //Trojwymiarowa macierz

    if(fread(lVectorX, X*sizeof(double), 1, fcache) != 1) {
        printf("Nie mozna odczytac 2. czesci cache'a\n");
        fclose(fcache);
        return 0;
    }

    if(fread(lVectorY, Y*sizeof(double), 1, fcache) != 1) {
        printf("Nie mozna odczytac 3. czesci cache'a\n");
        fclose(fcache);
        return 0;
    }

    if(fread(lVectorZ, Z*sizeof(double), 1, fcache) != 1) {
        printf("Nie mozna odczytac 4. czesci cache'a\n");
        fclose(fcache);
        return 0;
    }

    if(fread(lMValues, X*Y*Z*D*sizeof(double), 1, fcache) != 1) {
        printf("Nie mozna odczytac 5. czesci cache'a\n");
        fclose(fcache);
        return 0;
    }

    if(fread(lMCounted, X*Y*Z*sizeof(char), 1, fcache) != 1) {
        printf("Nie mozna odczytac 6. czesci cache'a\n");
        fclose(fcache);
        return 0;
    }

    fclose(fcache);

    int total;
    int count = 0;

    total = 0;

    for(int iz = 0; iz < Z; ++iz)
        for(int ix = 0; ix < X; ++ix)
            for(int iy = 0; iy < Y; ++iy) {
                long long idx = iz*(X*Y) + ix*Y + iy;
                if(lMCounted[idx] != 0 && MCounted[idx] == 0)  //Przepisujemy brakujacy w cache'u punkt
                {
                    memcpy(&MValues[idx*D], &lMValues[idx*D], D*sizeof(double));
                    MCounted[idx] = lMCounted[idx];

                    count++;
                }

                if(MCounted[idx] != 0) {
                    ++total;

                }
            }


    delete []lVectorX;
    delete []lVectorY;
    delete []lVectorZ;
    delete []lMValues;
    delete []lMCounted;

    if(total == X*Y*Z)
        cdesc.wholecachecalculated = 1;

    if(count > 0)
        save(filename);

    printf("Wczytano dane %d punktow z pliku %s, total: %d\n", count, lfilename, total);

    //    getc(stdin);
#endif
    return 1;
}

void Cache4d :: exportCache() {
#if 0
#ifdef WIN32
    mkdir("export");
#else
    mkdir("export", 0777);
#endif

    for(int d = 0; d < D; ++d) {
        for(int iz = 0; iz < Z; ++iz) {

            char buf[256];

            sprintf(buf, "export/cacheval_%E_%02d.dat", VectorZ[iz], d+1);
            FILE *fp = fopen(buf, "w+");

            if(fp == NULL) {
                printf("Nie moge otwiorzyc pliku %s\n", buf);
                return;
            }

            printf("Trwa eksport danych do pliku %s....\n", buf);

            for(int ix = 0; ix < X; ++ix) {
                for(int iy = 0; iy < Y; ++iy) {
                    int idx = iz*(X*Y) + ix*Y + iy;

                    if(MCounted[idx] == 1) {
                        fprintf(fp, "%E %E",  VectorX[ix], VectorY[iy]);
                        fprintf(fp, " %E",  *(&MValues[idx*D] + d));
                        fprintf(fp, "\n");
                    }
                }

//                fprintf(fp, "\n");
            }

            fclose(fp);
        }

    }

    int total = 0;
    for(int iz = 0; iz < Z; ++iz)
        for(int ix = 0; ix < X; ++ix)
            for(int iy = 0; iy < Y; ++iy) {
                int long long = iz*(X*Y) + ix*Y + iy;

                if(MCounted[idx] != 0) {
                    ++total;

                    //                    if(*(&MValues[idx*D]) == 0) {
                    //
                    //                        printf("kicha w cacheu %d %d %d\n", iz, iy, ix);
                    //                    }

                }
            }

    printf("Tworzenie pliku do wydruku cacheplot.d...\n");


    FILE *fp = fopen("export/cacheplot.d", "w+");

    if(fp == NULL) {
        printf("Nie moge otwiorzyc pliku %s\n", "export/cacheplot.d");
        return;
    }

    fprintf(fp,
            "set data style lines;\n"
            "set function style lines;\n"
            "set logscale xy;\n"
            "set terminal postscript landscape color;\n"
            "set size .5,.5;\n");


    for(int d = 0; d < D; ++d) {
        for(int iz = 0; iz < Z; iz += 4) {

            char buf[256];

            sprintf(buf, "%02d_%04d.ps", d+1, iz/4+1);

            fprintf(fp, "\nset output '%s';\n"
                    "set multiplot;\n", buf);

            fprintf(fp, "set origin 0.0,0.5;\n");
            fprintf(fp, "splot 'cacheval_%E_%02d.dat';\n", VectorZ[iz], d+1);

            if(iz+1 < Z) {
                fprintf(fp, "set origin 0.5,0.5;\n");
                fprintf(fp, "splot 'cacheval_%E_%02d.dat';\n", VectorZ[iz+1], d+1);
            }

            if(iz+2 < Z) {
                fprintf(fp, "set origin 0.0,0.0;\n");
                fprintf(fp, "splot 'cacheval_%E_%02d.dat';\n", VectorZ[iz+2], d+1);
            }

            if(iz+3 < Z) {
                fprintf(fp, "set origin 0.5,0.0;\n");
                fprintf(fp, "splot 'cacheval_%E_%02d.dat';\n", VectorZ[iz+3], d+1);
            }

            fprintf(fp, "unset multiplot;\n");
        }

    }

    fclose(fp);


    printf("Zawartosc cache'a %d punktow z %d\n", total, X*Y*Z);

#endif
}

int CacheWorker :: startThreads(int threads, CALCFUN calcfun_, int D_) {


    //Inicjalizacja
    if(threads > INI_MAX_THREAD)
        threads = INI_MAX_THREAD;

    if(threadstarted > 0)
        return 1;

    calcfun = calcfun_;
    D = D_;

#ifndef WIN32

    if(threads > 0)
    {
        if(proc_tab == NULL)
        {
            if ((shm_id_ = shmget(IPC_PRIVATE, threads*sizeof(proc_type_struct), IPC_CREAT | 0666)) == -1) {
                printf("Blad funkcji shmget size = %ld\n", threads*sizeof(proc_type_struct));
                return 1;
            }


            if ((proc_tab = (proc_type_struct *)shmat(shm_id_, 0, IPC_CREAT | 0666)) == NULL) {
                printf("Blad funkcji shmmat shm_id_ = %d\n", shm_id_);
                return 1;
            }

            for (int i=0; i < INI_MAX_THREAD; i++) {
                proc_tab[i].proc_nr = -1;
                proc_tab[i].procstate = INVALID;
            }


            if ((sem_id_ = semget(IPC_PRIVATE , 1, IPC_CREAT | IPC_EXCL | 0777)) == -1) {
                printf("Blad funkcji semget\n");
                return 1;
            }

            union semun {
                int val;
                struct semid_ds  * buf;
                ushort  * array;
            } semafor_arg;

            semafor_arg.val = 0;//;

            if (semctl(sem_id_, 0, SETVAL, semafor_arg) == -1) {
                printf("Blad ustawiania semaforow sem_id_ = %d\n", sem_id_);
                return 1;
            }
            
//            printf("shm_id=%d, sem_id=%d, proc_tab=%p\n", shm_id_, sem_id_, proc_tab);
        }
        
        //Uruchamiamy watki realizujace obliczenia
        if(threads > threadstarted)
        {
            int filedes[2];

            printf("Starting up %d threads\n", threads-threadstarted);

            for (; threadstarted < threads; threadstarted++) {

                //(*proc_tab)[watek].proc_nr = -1;

                if (pipe(filedes)) {
                    PutOnConsole("Blad funkcji pipe", "", "");
                    return 0;
                }

                proc_tab[threadstarted].pipe_read_child = filedes[0];
                proc_tab[threadstarted].pipe_write_child = filedes[1];


                if (pipe(filedes)) {
                    PutOnConsole("Blad funkcji pipe", "", "");
                    return 0;
                }
                proc_tab[threadstarted].pipe_read_main = filedes[0];
                proc_tab[threadstarted].pipe_write_main = filedes[1];


                int fork_nr = fork();
                if (fork_nr == -1) {
                    PutOnConsole("Blad funkcji fork", "", "");
                    return 0;
                }

                if (fork_nr != 0) {
                    //ojciec
//                    close(proc_tab[threadstarted].pipe_read_child);
                    fcntl(proc_tab[threadstarted].pipe_write_child, F_SETFL, O_NDELAY);
//                    close(proc_tab[threadstarted].pipe_write_main);

                } else {
                    //dziecko
//                    close(proc_tab[threadstarted].pipe_write_child);
//                    close(proc_tab[threadstarted].pipe_read_main);


/*
                    for(int i=0; i<threadstarted; i++)
                    {
                        close(proc_tab[i].pipe_read_child);
                        close(proc_tab[i].pipe_write_main);
                    }
*/
                    signal(SIGINT, SIG_IGN);
                    signal(SIGQUIT, SIG_IGN);
                    signal(SIGHUP, SIG_IGN);
#if !defined __CYGWIN__
					prctl(PR_SET_PDEATHSIG, SIGTERM);
					//signal(SIGTERM, SIG_IGN);
#endif
                    umask(0);

//	                printf("Thread %d started, pid=%d, proc_tab=%p\n", threadstarted, getpid(), proc_tab);
	                
                    proc_tab[threadstarted].proc_nr = getpid();

//                printf("Thread %d started\n", CHILD_NUMBER);

                    childMain(&proc_tab[threadstarted]);
                }
            }
        }
    }

#endif

    return 1;
}

void CacheWorker :: childMain(proc_type_struct *_proc_tab) {

#ifndef WIN32

    char mainPipe[64];
//    int  mainfd = -1;

    struct proc_message_struct message;
    //potomny proces


    signal(SIGTERM, Cache4d_childStop);


    _proc_tab->procstate = READY;


/*
    struct sembuf sem_ar;
    sem_ar.sem_num = 0 ;
    sem_ar.sem_op = 1 ;
    sem_ar.sem_flg = 0 ; //SEM_UNDO ;
/

    if (semop(sem_id_, &sem_ar , 1 ) == -1 ) {
        printf("PANIC!!! Blad semop DZIECKO errno %d\n", errno);
        int ppid = getppid();
        kill(ppid, SIGTERM);
    }
*/
    for (;;) {
        memset(&message, 0, sizeof(message));
        if (read(_proc_tab->pipe_read_child, &message, sizeof(message)) != sizeof(message)) {
            printf("Blad odczytu z pipe dziecka pid %d '%s'\n",  _proc_tab->proc_nr, strerror(errno));
            sleep(1);
            continue;
        }

/*
        if(strcmp(mainPipe, message.mainPipe) != 0)
        {
            if(mainfd != -1)
                close(mainfd);

            printf("About to open %s\n", message.mainPipe);
            mainfd = open(message.mainPipe, O_WRONLY);
//		    printf("Returned from open %s\n", name);
            if ( mainfd < 0)
            {
                printf("open err: %s\n", strerror(errno));
                exit(-1);
            }

            strcpy(mainPipe, message.mainPipe);
        }
*/

//        printf("------   Liczy watek %d pointidx=%lld\n", _proc_tab->proc_nr, message.idx);

        //printf("%E, %E, %E\n",           VectorX[message.ix], VectorY[message.iy], VectorZ[message.iz]);


        //            calcfun(VectorX[(messages[curmsg]).ix], VectorY[(messages[curmsg]).iy], VectorZ[(messages[curmsg]).iz], &MValues[(messages[curmsg]).idx*D]);

        double mvalues[D];

        calcfun(message.x, message.y, message.z, mvalues);

//        printf("------Policzyl watek %d pointidx=%lld\n", _proc_tab->proc_nr, message.idx);

        //printf("---policzyl watek %d\n", CHILD_NUMBER+1);
        if (write(_proc_tab->pipe_write_main, &message.idx, sizeof(long long)) != sizeof(long long))
        {
            printf("Blad zapisu pipe 1 fd = %d (errno %d %s)\n", _proc_tab->pipe_write_main, errno, strerror(errno));
            _proc_tab->procstate = INVALID;
            continue;
        }
        /*
                if (write(_proc_tab->pipe_write_main, &message.x, sizeof(double)) != sizeof(double) ||
                        write(_proc_tab->pipe_write_main, &message.y, sizeof(double)) != sizeof(double) ||
                        write(_proc_tab->pipe_write_main, &message.z, sizeof(double)) != sizeof(double))
                {
                    printf("Blad zapisu pipe 2.1 (errno %d %s)\n", errno, strerror(errno));
                    (*proc_tab)[CHILD_NUMBER].procstate = INVALID;
                    continue;
                }
        */
        int bytesToWrite = D*sizeof(double);
        if (write(_proc_tab->pipe_write_main, mvalues, bytesToWrite) != bytesToWrite)
        {
            printf("Blad zapisu pipe 2.2 (errno %d %s)\n", errno, strerror(errno));
            _proc_tab->procstate = INVALID;
            continue;
        }

        _proc_tab->procstate = DONE;

#if 0
        //struct sembuf sem_ar;
		struct sembuf sem_ar;
        sem_ar.sem_num = 0 ;
        sem_ar.sem_op = 1 ;
        sem_ar.sem_flg = 0 ; //SEM_UNDO ;

        if (semop(sem_id_, &sem_ar , 1 ) == -1 ) {
            printf("Blad semop DZIECKO (errno %d %s)\n", errno, strerror(errno));
            sleep(1);
            continue;
        }
#endif
    }

#endif

}

void CacheWorker :: waitForOne() {
#if 0
    struct sembuf sem_ar;
    sem_ar.sem_num = 0 ;
    sem_ar.sem_op = -1 ;
    sem_ar.sem_flg = 0 ;
    char Err_text[256];

    if (semop(sem_id_, &sem_ar , 1 ) == -1 ) {
        snprintf(Err_text, sizeof(Err_text), "%d", errno);
        PutOnConsole("Blad semop OJCIEC", strerror(errno), Err_text);
        return;
    }
#else
	usleep(100);
#endif
}

void CacheWorker :: waitForAll() {

    struct sembuf sem_ar;
    sem_ar.sem_num = 0 ;
    sem_ar.sem_op = -threadstarted ;
    sem_ar.sem_flg = 0 ;


    //Czekamy az sie policza wszystkie watki
    if (semop(sem_id_, &sem_ar , 1 ) == -1 ) {   //opuszczam semafory
        char Err_text[512];
        snprintf(Err_text, sizeof(Err_text), "%d", errno);
        PutOnConsole("Blad semop OJCIEC", strerror(errno), Err_text);
        return;
    }

    sem_ar.sem_op = threadstarted;

    //watki policzone, podnosimy semafory
    if (semop(sem_id_, &sem_ar , 1 ) == -1 ) {   //podnosze wszystkie opuszczone semafory
        char Err_text[512];
        snprintf(Err_text, sizeof(Err_text), "%d", errno);
        PutOnConsole("Blad semop OJCIEC", strerror(errno), Err_text);
        return;
    }


}


CacheWorker::~CacheWorker() {
#ifndef WIN32
    if(threadstarted > 0)
	    for(int i=0; i<threadstarted; i++)
    	{
        	kill(proc_tab[i].proc_nr, SIGKILL);
        }

    for(int i=0; i<threadstarted; i++)
    {
        close(proc_tab[i].pipe_read_child);
        close(proc_tab[i].pipe_read_main);
        close(proc_tab[i].pipe_write_child);
        close(proc_tab[i].pipe_write_main);
    }

    if(proc_tab != NULL) {
        struct shmid_ds shm_buf;

        if (shmdt(proc_tab) == -1) {
            PutOnConsole("Blad systemowy 'deatach shared memory'\n", "", "");
        }

        if (shmctl(shm_id_, IPC_RMID, &shm_buf)) {
            PutOnConsole("Blad systemowy 'remove shared memory'\n", "", "");
        }

        proc_tab = NULL;
    }

    if(sem_id_ != -1)
    {
        if (semctl(sem_id_, 0, IPC_RMID) == -1)
        {
            PutOnConsole("Blad systemowy 'remove sem'\n", "", "");
        }
        sem_id_ = -1;
    }
#endif
}

#ifdef TEST

void calc_cache_val_test(double x, double y, double z, double *val) {
    printf("Licze %E, %E, %E\n", x, y, z);

    sleep(2);

    val[0] = 100.0;
    val[1] = 200.0;
    val[2] = 300.0;
}

int main(void) {

    CacheDesc cdesc;

    cdesc.x = 128;
    cdesc.xmin      = 5E5;
    cdesc.xmax      = 1E13;
    cdesc.xtype = LOGARITHMIC;

    cdesc.y = 128;
    cdesc.ymin      = 5E7;
    cdesc.ymax      = 1E13;
    cdesc.ytype = LOGARITHMIC;

    cdesc.z = 32;
    cdesc.zmin      = 0.01;
    cdesc.zmax      = 20.0;
    cdesc.ztype = LOGARITHMIC; //LINEAR;

    cdesc.d  = 3;        //liczba tablic

    cdesc.autosave = 100;  //Liczba punktow po ktorych policzeniu zapisywany jest cache

    Cache4d cache(&cdesc, calc_cache_val_test, "c1.dat", 4);

    //	cache.calculateCache();

    //   cache.test();

    double val[3];
    double x, y, z;

    x = 1E9;
    y = 1E9;
    z = 0.02;

    cache.getvalues(x, y, z, val);
    printf("%E, %E, %E\n", val[0], val[1], val[2]);

    return 1;
}

#endif







