#ifndef __GENERAL_DEFINITION
#define __GENERAL_DEFINITION

#include <assert.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>

#define BFHEAD_LENGTH (sizeof(int)*2)    //file header size

typedef char Block[];
//-------------------All-------------------------------------
#define MAXREAL         1e20f
#define FLOATZERO       1e-8f


#define MIN(a, b) (((a) < (b))? (a) : (b)  )
#define MAX(a, b) (((a) > (b))? (a) : (b)  )

#define TIMESFACTOR 5931.0f
#define FloatToInt(a) (int)((double)(sqrt(a))*TIMESFACTOR);

extern int blocksize;

//-------------------Class and other data types--------------
class BlockFile;  //for BlockFile definition
class Cache;
class Cacheable {  
	// inherit this class if you wish to use an external memory structure with a cache
public:
	BlockFile *file;
	Cache *cache;
};

enum SECTION {OVERLAP, INSIDE, S_NONE};
enum R_OVERFLOW {SPLIT, REINSERT, NONE};
enum R_DELETE {NOTFOUND,NORMAL,ERASED};
typedef float *floatptr;

struct SortMbr  {
    int dimension;
    float *mbr;
    float *center;
    int index;
};

//-----Global Functions--------------------------------------
void error(char *_errmsg, bool _terminate);
bool section(int dimension,float *mbr1,float *mbr2);
float area(int dimension, float *mbr);
float margin(int dimension, float *mbr);
float overlap(int dimension, float *r1, float *r2);
float* overlapRect(int dimension, float *r1, float *r2);
bool is_overlapRect(float *r1, float *r2, int dimension);
float MINDIST(float *bounces1, float *bounces2, int dim);
float MINDIST_nosqrt(float *bounces1, float *bounces2, int dim);
float MAXDIST(float *bounces1, float *bounces2, int dim);
float MAXDIST_nosqrt(float *bounces1, float *bounces2, int dim);
float MINMINDIST(float *bounces1, float *bounces2, int dim);
float MAX_P_DIST(float *bounces1, float *bounces2, int dim);

bool inside(float &p, float &lb, float &ub);
void enlarge(int dimension, float **mbr, float *r1, float *r2);
bool inside(float *p, float *mbr, int dimension);
bool insideRect(float *p, float *mbr, int dimension);

float calDiagonal(float *mbr, int dimension);

int sort_lower_mbr(const void *d1, const void *d2);
int sort_upper_mbr(const void *d1, const void *d2);
int sort_center_mbr(const void *d1, const void *d2);


#endif

