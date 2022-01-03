#include <iostream>
#include <vector>
#include "OP.h"

#ifndef UNIFIER_H
#define UNIFIER_H

class Unifier {
private:
	std::vector<std::vector<Solution>> m_solutions;
	std::vector<std::vector<double>> m_ttMatrix;
	int m_tMax;
	TA* m_depot;

public:
	Unifier();
	Unifier(std::vector<std::vector<Solution>>, int, std::vector<std::vector<double>>, TA*);
	~Unifier();

	std::vector<Solution> exec();
};

#endif