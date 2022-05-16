
#if !defined(RANDOM_H__INCLUDED_)
#define RANDOM_H__INCLUDED_


typedef unsigned int UL;

typedef struct {
    UL z, w, jsr, jcong;
}
RANDOMDEF;

void setRandomTable(RANDOMDEF *rd, UL i1,UL i2,UL i3,UL i4);
double UniformRandom(RANDOMDEF *rd);  //[-1, 1]
double UniformRandomRange(RANDOMDEF *rd, double min, double max);
double UniformRandomRangeLog(RANDOMDEF *rd, double min, double max);



#endif
