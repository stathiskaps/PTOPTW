#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <stdexcept>
#include <tuple>
#include <climits>
#include <numeric>
#include <signal.h>
#include <fmt/ranges.h>
#include <plog/Log.h> // Step1: include the headers
#include "plog/Initializers/RollingFileInitializer.h"
#include "List.h"
#include "Solution.h"
#include "OP.h"

class ILS{
private:
	int mBucketsNum;

	virtual std::tuple<int, double, int, int> getBestPos(TA*, ListTA, std::vector<std::vector<double>>);
	virtual Solution updateTimes(Solution solution, int startIndex, bool smart, std::vector<std::vector<double>> ttMatrix);
	Solution LocalSearch(Solution solution, double avgPoint, OP&);
	Solution Shake(Solution solution, int S, int R, int numOfPois, std::vector<std::vector<double>> ttMatrix);
	std::tuple<double, double, double> calcTimeEventCut(ListTA&);
	virtual Walk updateMaxShifts(Walk, std::vector<std::vector<double>>);
	ListTA setBucketActivityDurations(ListTA&, double);
	Solution construct(Solution, std::vector<std::vector<double>>);
	std::vector<std::vector<TA*>> getBuckets(std::vector<TA*>, int);
	std::vector<double> getTimeCuts(std::vector<std::vector<TA*>>);
public:
    ILS();
	ILS(int);
    ~ILS();
	virtual std::tuple<bool, std::string> validate(ListTA&, std::vector<std::vector<double>>);
	Solution Solve(OP&);
	

};

class ILS_OPTW : public ILS {
private:
	std::tuple<int, double, int, int> getBestPos(TA*, ListTA, std::vector<std::vector<double>>) override;
	Solution updateTimes(Solution, int, bool, std::vector<std::vector<double>>) override;
	Walk updateMaxShifts(Walk, std::vector<std::vector<double>>) override;
public:
	using ILS::ILS; //inherit constructor
	std::tuple<bool, std::string> validate(ListTA&, std::vector<std::vector<double>>) override;
};

class ILS_TOPTW : public ILS {
private:
	std::tuple<int, double, int, int> getBestPos(TA*, ListTA, std::vector<std::vector<double>>) override;
	Solution updateTimes(Solution, int, bool, std::vector<std::vector<double>>) override;
	Walk updateMaxShifts(Walk, std::vector<std::vector<double>>) override;
public:
	using ILS::ILS; //inherit constructor
	std::tuple<bool, std::string> validate(ListTA&, std::vector<std::vector<double>>) override;
};

