#pragma once

#include <iostream>
#include "Custom.h"
#include "OP.h"

struct Sol{
    List<TA> walk;
    std::vector<List<TA>> buckets;
}

class IteratedLocalSearch{
    private:

    struct Interval{
        double start_time;
        double end_time;
    };

    Main();
    LocalSearch();
    Perturbation();

    public:
    IteratedLocalSearch();
    ~IteratedLocalSearch();
    IteratedLocalSearch(OP&);

};