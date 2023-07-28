#pragma once
#include "OptAssign.h"

class COptConstraint :
    public COptAssign
{
public:
    COptConstraint(void);
    ~COptConstraint(void);

    void computeAssign(vector<int> &constraint);
};
