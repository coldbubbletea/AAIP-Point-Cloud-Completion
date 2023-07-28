#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gendef.h"
#include "float_cast.h"
#include<algorithm>

using namespace std;

void error(char *t, bool ex) {
    fprintf(stderr, t);
    if (ex) 
    {
        printf("unknown error\n");
        exit(0);
    }
}

bool section(int dimension,float *mbr1,float *mbr2) {
    for (int i=0;i<dimension;i++) 
		if (mbr1[2*i]>mbr2[2*i+1]||mbr1[2*i+1]<mbr2[2*i])
	    	return false;
    return true;
}


float area(int dimension, float *mbr) {
    float sum = 1.0;
    for (int i = 0; i < dimension; i++)
		sum *= mbr[2*i+1] - mbr[2*i];
    return sum;
}

float margin(int dimension, float *mbr) {
    float *ml, *mu, *m_last, sum;
    sum = 0.0;
    m_last = mbr + 2*dimension;
    ml = mbr;
    mu = ml + 1;
    while (mu < m_last) {
		sum += *mu - *ml;
		ml += 2;
		mu += 2;
    }
    return sum;
}

bool inside(float &p, float &lb, float &ub) {
    return (p >= lb && p <= ub);
}

bool inside(float *v, float *mbr, int dimension) {
    for (int i = 0; i < dimension; i++)
		if (!inside(v[i], mbr[2*i], mbr[2*i+1])) return false;
    return true;
}

bool insideRect(float *v, float *mbr, int dimension) {
    for (int i = 0; i < dimension; i++)
    {
		if (!inside(v[2*i], mbr[2*i], mbr[2*i+1])) return false;
        if (!inside(v[2*i+1], mbr[2*i], mbr[2*i+1])) return false;
    }
    return true;
}

// calcutales the overlapping rectangle between r1 and r2
// if rects do not overlap returns null
bool is_overlapRect(float *r1, float *r2, int dimension) {
	for (int i=0; i<dimension; i++) {
	    if ((r1[i*2]>r2[i*2+1]) || (r1[i*2+1]<r2[i*2])) { // non overlapping
            return false;
		}
	}
	return true;
}

// calcutales the overlapping rectangle between r1 and r2
// if rects do not overlap returns null
float* overlapRect(int dimension, float *r1, float *r2) {
	float *overlap = new float[2*dimension];
	for (int i=0; i<dimension; i++) {
	    if ((r1[i*2]>r2[i*2+1]) || (r1[i*2+1]<r2[i*2])) { // non overlapping
	        delete [] overlap;
			return NULL;
		}
		overlap[2*i] = MAX(r1[i*2], r2[i*2]);
	    overlap[2*i+1] = MIN(r1[i*2+1], r2[i*2+1]);
	}
	return overlap;
}

float overlap(int dimension, float *r1, float *r2) {
	// calcutales the overlapping area of r1 and r2
	// calculate overlap in every dimension and multiplicate the values
    float *r1pos, *r2pos, *r1last, r1_lb, r1_ub, r2_lb, r2_ub;
    float sum = 1.0;
    r1pos = r1; r2pos = r2;
    r1last = r1 + 2 * dimension;
    while (r1pos < r1last) {
		r1_lb = *(r1pos++);
		r1_ub = *(r1pos++);
		r2_lb = *(r2pos++);
		r2_ub = *(r2pos++);
        // calculate overlap in this dimension

        if (inside(r1_ub, r2_lb, r2_ub)) {
        	// upper bound of r1 is inside r2 
            if (inside(r1_lb, r2_lb, r2_ub))
            	// and lower bound of r1 is inside
                sum *= (r1_ub - r1_lb);
            else
                sum *= (r1_ub - r2_lb);
		} else {
            if (inside(r1_lb, r2_lb, r2_ub))
	    	// and lower bound of r1 is inside
				sum *= (r2_ub - r1_lb);
	    	else {
				if (inside(r2_lb, r1_lb, r1_ub)&&inside(r2_ub, r1_lb, r1_ub))
	        		sum *= (r2_ub - r2_lb);		// r1 contains r2
				else
					sum = 0.0;					// r1 and r2 do not overlap
	    	}
		}
    }
    return sum;
}

void enlarge(int dimension, float **mbr, float *r1, float *r2) {
	// enlarge r in a way that it contains s
    *mbr = new float[2*dimension];
    for (int i = 0; i < 2*dimension; i += 2) {
		(*mbr)[i]   = MIN(r1[i],   r2[i]);
		(*mbr)[i+1] = MAX(r1[i+1], r2[i+1]);
    }
}

int sort_lower_mbr(const void *d1, const void *d2) {
    SortMbr *s1, *s2;
    s1 = (SortMbr *) d1;
    s2 = (SortMbr *) d2;
    int dimension = s1->dimension;
    float erg = s1->mbr[2*dimension] - s2->mbr[2*dimension];
    if (erg < 0.0)
		return -1;
    else if (erg == 0.0)
		return 0;
    else 
		return 1;
}

int sort_upper_mbr(const void *d1, const void *d2) {
    SortMbr *s1, *s2;
    
    s1 = (SortMbr *) d1;
    s2 = (SortMbr *) d2;
    int dimension = s1->dimension;
    float erg = s1->mbr[2*dimension+1] - s2->mbr[2*dimension+1];
    if (erg < 0.0)
		return -1;
    else if (erg == 0.0)
		return 0;
    else 
		return 1;
}

int sort_center_mbr(const void *d1, const void *d2) {
    SortMbr *s1, *s2;
    float d, e1, e2;

    s1 = (SortMbr *) d1;
    s2 = (SortMbr *) d2;
    int dimension = s1->dimension;

    e1 = e2 = 0.0;
    for (int i = 0; i < dimension; i++) {
        d = ((s1->mbr[2*i] + s1->mbr[2*i+1]) / 2.0f) - s1->center[i];
        e1 += d*d;
        d = ((s2->mbr[2*i] + s2->mbr[2*i+1]) / 2.0f) - s2->center[i];
        e2 += d*d;
    }
    
    if (e1 < e2)
		return -1;
    else if (e1 == e2)
		return 0;
    else 
		return 1;
}

float MINDIST(float *bounces1, float *bounces2, int dim) {
    float r,summe=0.0;
    for (int i = 0; i < dim; i++) {
    	r=0;
		if (bounces1[2*i+1]<bounces2[2*i])
			r=bounces2[2*i]-bounces1[2*i+1];
		else if (bounces2[2*i+1]<bounces1[2*i])
			r=bounces1[2*i]-bounces2[2*i+1];
		summe += r*r;
    }
    return FloatToInt(summe);
}

float MINDIST_nosqrt(float *bounces1, float *bounces2, int dim) {
    float r,summe=0.0;
    for (int i = 0; i < dim; i++) {
    	r=0;
		if (bounces1[2*i+1]<bounces2[2*i])
			r=bounces2[2*i]-bounces1[2*i+1];
		else if (bounces2[2*i+1]<bounces1[2*i])
			r=bounces1[2*i]-bounces2[2*i+1];
		summe += r*r;
    }
    return summe;
}

float MAXDIST(float *bounces1, float *bounces2, int dim) {
    float r,summe=0.0;
    for (int i = 0; i < dim; i++) {
        if (abs(bounces2[2*i+1]-bounces1[2*i])>abs(bounces1[2*i+1]-bounces2[2*i]))
            r=bounces2[2*i+1]-bounces1[2*i];
        else
            r=bounces1[2*i+1]-bounces2[2*i];
		summe += r*r;
    }
    return sqrt(summe);
}

float MIN_DIST_NOSQRT(float *bounces1, float *bounces2, int dim) {
    float r,summe=0.0;
    for (int i = 0; i < dim; i++) {
    	r=0;
		if (bounces1[2*i+1]<bounces2[2*i])
			r=bounces2[2*i]-bounces1[2*i+1];
		else if (bounces2[2*i+1]<bounces1[2*i])
			r=bounces1[2*i]-bounces2[2*i+1];
		summe += r*r;
    }
    return (int)(summe);
}

float MAXDIST_nosqrt(float *bounces1, float *bounces2, int dim) {
    float r,summe=0.0;
    for (int i = 0; i < dim; i++) {
        if (abs(bounces2[2*i+1]-bounces1[2*i])>abs(bounces1[2*i+1]-bounces2[2*i]))
            r=bounces2[2*i+1]-bounces1[2*i];
        else
            r=bounces1[2*i+1]-bounces2[2*i];
		summe += r*r;
    }
    return (int)(summe);
}

float MAX_P_DIST(float *bounces1, float *bounces2, int dim) {
    float r,summe=0.0;
    for (int i = 0; i < dim; i++) {
    	r=0;
		if (bounces1[2*i+1]<bounces2[2*i])
			r=bounces2[2*i]-bounces1[2*i+1];
		else if (bounces2[2*i+1]<bounces1[2*i])
			r=bounces1[2*i]-bounces2[2*i+1];
        summe+=r;
    }
    return FloatToInt(summe);
}

float MINMINDIST(float *bounces1, float *bounces2, int dim) {
    float r,summe=0.0;
    for (int i = 0; i < dim; i++) {
    	r=0;
		if (bounces1[2*i+1]<bounces2[2*i])
			r=bounces2[2*i]-bounces1[2*i+1];
		else if (bounces2[2*i+1]<bounces1[2*i])
			r=bounces1[2*i]-bounces2[2*i+1];
		if (r<summe)
            summe=r;
    }
    return r;
}

float calDiagonal(float* mbr, int dim) {
    float sum=0.0f;
    for (int i=0; i<dim; ++i) {
        sum+=(mbr[2*i+1]-mbr[2*i])*(mbr[2*i+1]-mbr[2*i]);
    }
 
    return sqrt(sum);
}
