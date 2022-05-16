#include <stdio.h>
#include <math.h>


#include "random.h"


/* Use random seeds to reset z,w,jsr,jcong*/
void setRandomTable(RANDOMDEF *rd, UL i1,UL i2,UL i3,UL i4) {
    rd->z=i1;
    rd->w=i2,
          rd->jsr=i3;
    rd->jcong=i4;
}


double UniformRandom(RANDOMDEF *rd) {  //[-1, 1]


#define znew   (rd->z=36969*(rd->z&65535)+(rd->z>>16))
#define wnew   (rd->w=18000*(rd->w&65535)+(rd->w>>16))
#define MWC    ((znew<<16)+wnew )
#define SHR3  (rd->jsr^=(rd->jsr<<17), rd->jsr^=(rd->jsr>>13), rd->jsr^=(rd->jsr<<5))
#define CONG  (rd->jcong=69069*rd->jcong+1234567)
#define KISS  ((MWC^CONG)+SHR3)

    return ((double)((int) KISS))/0x7FFFFFFF;
}

double UniformRandomRange(RANDOMDEF *rd, double min, double max) {

#define znew   (rd->z=36969*(rd->z&65535)+(rd->z>>16))
#define wnew   (rd->w=18000*(rd->w&65535)+(rd->w>>16))
#define MWC    ((znew<<16)+wnew )
#define SHR3  (rd->jsr^=(rd->jsr<<17), rd->jsr^=(rd->jsr>>13), rd->jsr^=(rd->jsr<<5))
#define CONG  (rd->jcong=69069*rd->jcong+1234567)
#define KISS  ((MWC^CONG)+SHR3)

    double u = ((double)(KISS))/0xFFFFFFFF;  //[0, 1]

    return min + (max - min)*u;

}

double UniformRandomRangeLog(RANDOMDEF *rd, double min, double max) {

#define znew   (rd->z=36969*(rd->z&65535)+(rd->z>>16))
#define wnew   (rd->w=18000*(rd->w&65535)+(rd->w>>16))
#define MWC    ((znew<<16)+wnew )
#define SHR3  (rd->jsr^=(rd->jsr<<17), rd->jsr^=(rd->jsr>>13), rd->jsr^=(rd->jsr<<5))
#define CONG  (rd->jcong=69069*rd->jcong+1234567)
#define KISS  ((MWC^CONG)+SHR3)

    double u = ((double)(KISS))/0xFFFFFFFF;  //[0, 1]


    return pow(10.0, log10(min) + (log10(max) - log10(min))*u);

}
