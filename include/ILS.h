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
#include "Divider.h"
#include "OP.h"
#include "Custom.h"
//#include "boost/geometry.hpp"
//#include "ygor/YgorClustering.hpp"

template <typename T>
using Vector2D = std::vector<std::vector<T>>;

class ILS{
private:
	int mBucketsNum;

	//virtual std::tuple<int, double, int, int> getBestPos(TA*, ListTA, std::vector<std::vector<double>>&);
	virtual std::tuple<CustomListTA::iterator, double, int, int> getBestPos(const TA&, const CustomList<TA>&, const Vector2D<double>&);
	//virtual void updateTimes(Solution&, int, bool, std::vector<std::vector<double>>&);
	virtual void updateTimes(CustomSolution&, CustomList<TA>::iterator&, const bool, const Vector2D<double>&);
	void LocalSearch(std::vector<Solution>&, std::vector<double>, OP&);
	void SplitSearch(Solution&, std::vector<double>, OP&);
	void Shake(std::vector<Solution>&, std::vector<ShakeParameters>, OP&);
	void Shake(Solution&, int&, int&, OP&);
	void Shake(CustomSolution&, int&, int&, OP&);
	std::tuple<double, double, double> calcTimeEventCut(ListTA&);
	//virtual void updateMaxShifts(Walk&, std::vector<std::vector<double>>&);
	virtual void updateMaxShifts(const CustomList<TA>&, const Vector2D<double>&);
	ListTA setBucketActivityDurations(ListTA&, double, std::vector<double>);
	void construct(CustomSolution&, const Vector2D<double>&);
	//void construct(Solution&, std::vector<std::vector<double>>&);
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
	inline virtual void validate(const CustomList<TA>&, const Vector2D<double>&);
	inline void print(const CustomList<TA>&);
	//Solution Solve(OP&);
	void SolveNew(OP&);
	Solution Preprocess(OP&);
	

};

class ILS_OPTW : public ILS {
private:
	std::tuple<CustomListTA::iterator, double, int, int> getBestPos(const TA&, const CustomList<TA>&, const Vector2D<double>&) override;
	//std::tuple<int, double, int, int> getBestPos(TA*, ListTA, std::vector<std::vector<double>>&) override;
	void updateTimes(CustomSolution&, CustomList<TA>::iterator&, const bool, const Vector2D<double>&) override;
	//void updateTimes(Solution&, int, bool, std::vector<std::vector<double>>&) override;
	void updateMaxShifts(const CustomList<TA>&, const Vector2D<double>&) override;
	//void updateMaxShifts(Walk&, std::vector<std::vector<double>>&) override;
public:
	using ILS::ILS; //inherit constructor
	void validate(ListTA&, std::vector<std::vector<double>>) override;
};

