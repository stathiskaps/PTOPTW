#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <stdexcept>
#include <tuple>
#include <chrono>
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

	virtual std::tuple<int, double, int, int> getBestPos(TA*, ListTA, std::vector<std::vector<double>>&);
	virtual void updateTimes(Solution&, int, bool, std::vector<std::vector<double>>&);
	void LocalSearch(std::vector<Solution>&, std::vector<double>, OP&);
	void Shake(std::vector<Solution>&, std::vector<ShakeParameters>, OP&);
	std::tuple<double, double, double> calcTimeEventCut(ListTA&);
	virtual void updateMaxShifts(Walk&, std::vector<std::vector<double>>&);
	ListTA setBucketActivityDurations(ListTA&, double);
	void construct(Solution&, std::vector<std::vector<double>>&);
	std::vector<std::vector<TA*>> getBuckets(std::vector<TA*>, int);
	int collectScore(std::vector<Solution>);
	std::vector<double> getTimeCuts(std::vector<std::vector<TA*>>);
	std::tuple<int, int> getMinMaxLength(std::vector<Solution> solutions);
	Solution connectSolutions(std::vector<Solution>);
public:
    ILS();
	ILS(int);
    ~ILS();
	virtual void validate(ListTA&, std::vector<std::vector<double>>);
	Solution Solve(OP&);
	

};

class ILS_OPTW : public ILS {
private:
	std::tuple<int, double, int, int> getBestPos(TA*, ListTA, std::vector<std::vector<double>>&) override;
	void updateTimes(Solution&, int, bool, std::vector<std::vector<double>>&) override;
	void updateMaxShifts(Walk&, std::vector<std::vector<double>>&) override;
public:
	using ILS::ILS; //inherit constructor
	void validate(ListTA&, std::vector<std::vector<double>>) override;
};

class ILS_TOPTW : public ILS {
private:
	std::tuple<int, double, int, int> getBestPos(TA*, ListTA, std::vector<std::vector<double>>&) override;
	void updateTimes(Solution&, int, bool, std::vector<std::vector<double>>&) override;
	void updateMaxShifts(Walk&, std::vector<std::vector<double>>&) override;
public:
	using ILS::ILS; //inherit constructor
	void validate(ListTA&, std::vector<std::vector<double>>) override;
};

