#include <iostream>
#include <vector>
#include "OP.h"

#ifndef PROCESSOR_H
#define PROCESSOR_H

class Processor {
private:
	std::vector<std::vector<Cluster>> m_clusters;
	std::vector<std::vector<double>> m_ttMatrix;
public:
	Processor();
	Processor(std::vector<std::vector<Cluster>>, std::vector<std::vector<double>>);
	~Processor();
	std::vector<std::vector<Solution>> exec();
};

#endif