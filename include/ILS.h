#include <iostream>
#include "OP.h"
#include "List.h"
#include "Solution.h"

#ifndef ILS_H
#define ILS_H

class ILS{
private:
	virtual std::tuple<int, double, int, int> getBestPos(TA*, ListTA&, std::vector<std::vector<double>>);
	virtual Solution updateTimes(Solution solution, int startIndex, bool smart, std::vector<std::vector<double>> ttMatrix);
	Solution LocalSearch(Solution solution, double avgPoint);
	Solution Shake(Solution solution, int S, int R, int numOfPois, std::vector<std::vector<double>> ttMatrix);
	std::tuple<double, double, double> calcTimeEventCut(ListTA&);
	ListTA setBucketActivityDurations(ListTA walk);
public:
    ILS();
    ~ILS();
    void exec();
    int Insert(ListTA&, ListTA&, std::vector<std::vector<double>>);
	virtual std::tuple<bool, std::string> validate(ListTA&, std::vector<std::vector<double>>);
	Solution Solve(ListTA&, TA*, TA*, std::vector<std::vector<double>>);
	

};

class ILS_OPTW : public ILS {
private:
	std::tuple<int, double, int, int> getBestPos(TA*, ListTA&, std::vector<std::vector<double>>) override;
	Solution updateTimes(Solution, int, bool, std::vector<std::vector<double>>) override;
public:
	using ILS::ILS; //inherit constructor
	std::tuple<bool, std::string> validate(ListTA&, std::vector<std::vector<double>>) override;
};

class ILS_TOPTW : public ILS {
private:
	std::tuple<int, double, int, int> getBestPos(TA*, ListTA&, std::vector<std::vector<double>>) override;
	Solution updateTimes(Solution, int, bool, std::vector<std::vector<double>>) override;
public:
	using ILS::ILS; //inherit constructor
	std::tuple<bool, std::string> validate(ListTA&, std::vector<std::vector<double>>) override;
};



#endif