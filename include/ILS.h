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

using Walks = std::vector<CustomList<TA>>;

class ILS{
private:
	struct Bin{
		CustomList<TA> unvisited;
		TimeWindow tw;
	};

	int mBucketsNum;

	bool compareTimeWindowCenter(const CustomList<TA>::iterator&, const CustomList<TA>::iterator&);
	virtual std::tuple<CustomList<TA>::iterator, double, int, int> getBestPos(const TA&, const CustomList<TA>&, const Vector2D<double>&);
	virtual std::tuple<Walks::iterator, CustomList<TA>::iterator, double, int, int> getBestPos(const TA&, Walks&, const Vector2D<double>&);
	virtual void updateTimes(CustomList<TA>&, const CustomList<TA>::iterator&, const bool, const Vector2D<double>&);
	void SplitSearch(std::vector<CustomSolution>&, const std::vector<double>&, OP&);
	std::vector<CustomSolution> splitSolution(CustomSolution&, const std::vector<double>&);
	inline Point getWeightedCentroid(const CustomList<TA>::iterator& first, const CustomList<TA>::iterator& last);
	void Shake(CustomSolution&, int&, int&, OP&, const int&);
	std::tuple<double, double, double> calcTimeEventCut(ListTA&);
	virtual void updateMaxShifts(const CustomList<TA>&, const Vector2D<double>&);
	ListTA setBucketActivityDurations(ListTA&, double, std::vector<double>);
	void construct(CustomSolution&, const Vector2D<double>&);
	std::vector<std::vector<TA*>> getBuckets(std::vector<TA*>, int);
	int collectScore(std::vector<Solution>);
	int collectScores(std::vector<CustomSolution>);
	std::vector<double> getTimeCuts(std::vector<std::vector<TA*>>);
	std::tuple<int, int> getMinMaxLength(std::vector<Solution> solutions);
	CustomSolution connectSolutions(std::vector<CustomSolution>&, const size_t);
	std::vector<double> Preprocessing(std::vector<TA*>, int, double);
	inline int collectProfit (const CustomList<TA>::iterator&, const CustomList<TA>::iterator&) const;


	
public:
    ILS();
	ILS(int);
    ~ILS();
	virtual void validate(const CustomList<TA>&, const Vector2D<double>&);
	virtual void validate(const Walks&, const Vector2D<double>&);
	inline void print(const CustomList<TA>&);
	void SolveNew(OP&);
	Solution Preprocess(OP&);
	//template <typename Container, typename ConstIterator>
	//typename Container::iterator remove_constness(Container& c, ConstIterator it) const
	//{
	//	return c.erase(it, it);
	//}

};

class ILS_TOPTW : public ILS {
private:
	std::tuple<CustomList<TA>::iterator, double, int, int> getBestPos(const TA&, const CustomList<TA>&, const Vector2D<double>&) override;
	std::tuple<Walks::iterator, CustomList<TA>::iterator, double, int, int> getBestPos(const TA&, Walks&, const Vector2D<double>&) override;
	void updateTimes(CustomList<TA>&, const CustomList<TA>::iterator&, const bool, const Vector2D<double>&) override;
	void updateMaxShifts(const CustomList<TA>&, const Vector2D<double>&) override;
public:
	using ILS::ILS; //inherit constructor
	void validate(const CustomList<CustomList<TA>>&, const Vector2D<double>&);
	void validate(const CustomList<TA>&, const Vector2D<double>&) override;
};

