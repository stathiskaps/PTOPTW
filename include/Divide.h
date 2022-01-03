#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <map>
#include <climits>
#include <stdlib.h> 
//#include <unistd.h>
#include "List.h"
#include "OP.h"


#ifndef DIVIDE_H
#define DIVIDE_H


class Divide{
    private:
    std::vector<TA*> mAttractions;

    public:
    Divide();
    Divide(std::vector<TA*>);
    ~Divide();

    int calcMeanTimePoint();
    std::vector<OP> exec();
    static bool compareByProfit(const TA* a, const TA* b);



};

#endif