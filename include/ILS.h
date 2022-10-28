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
#include <chrono>
#include <plog/Log.h> // Step1: include the headers
#include "plog/Initializers/RollingFileInitializer.h"
#include "List.h"
#include "Solution.h"
#include "Divider.h"
#include "OP.h"
#include "Custom.h"
//#include "boost/geometry.hpp"
#include "ygor/YgorClustering.hpp"

template <typename T>
using Vector2D = std::vector<std::vector<T>>;

using Walks = std::vector<List<TA>>;

using MapOfActivities = std::map<std::string, Activity>;

class ILS{
private:
	struct Bin{
		List<TA> unvisited;
		TimeWindow tw;
	};

	struct Interval {
		double start_time;
		double end_time;
	};

	//std::map<std::string, std::vector<ActivityInBucket>> registry;

	int mBucketsNum;
	void dbScan(OP& op);
	void AddStartDepots(std::vector<Solution>&, const int, const OP&);
	void AddEndDepots(std::vector<Solution>&, const std::vector<double>&, const int, OP&);
	bool compareTimeWindowCenter(const List<TA>::iterator&, const List<TA>::iterator&);
	virtual std::tuple<List<TA>::iterator, double, int, int> getBestPos(const TA&, const List<TA>&, const Vector2D<double>&);
	virtual std::tuple<Walks::iterator, List<TA>::iterator, double, int, int> getBestPos(const TA&, Walks&, const Vector2D<double>&);
	virtual void updateTimes(List<TA>&, const List<TA>::iterator&, const bool, const Vector2D<double>&);
	void SplitSearch(std::vector<Solution>&, const std::vector<double>&, OP&, std::map<std::string, Activity>&);
	void SplitSearch2(Solution&, const std::vector<double>&, OP&, std::map<std::string, Activity>&);
	std::vector<Bin> splitUnvisited(List<TA>&, std::map<std::string, Activity>&);
	bool hasWeightedCentroid(const Solution& sol, const int, const int);
	std::vector<Solution> splitSolution(Solution&, const std::vector<double>&, std::map<std::string, Activity>&);
	inline Point getWeightedCentroid(const List<TA>::iterator& first, const List<TA>::iterator& last);
	void Shake(Solution&, int&, int&, OP&, const int&);
	virtual void updateMaxShifts(const List<TA>&, const Vector2D<double>&);
	std::vector<std::string> construct(Solution&, const Vector2D<double>&);
	int collectScores(std::vector<Solution>);
	Solution connectSolutions(std::vector<Solution>&, const size_t);
	std::vector<double> Preprocessing(std::vector<TA*>, int, double);
	inline int collectProfit (const List<TA>::iterator&, const List<TA>::iterator&) const;
	inline std::map<std::string, Activity> initializeRegistry(const List<TA>&, const std::vector<double>&);


	
public:
    ILS();
	ILS(int);
    ~ILS();
	virtual void validate(const List<TA>&, const Vector2D<double>&);
	virtual void validate(const Walks&, const Vector2D<double>&);
	inline void print(const List<TA>&);
	void SolveNew(OP&);
	//template <typename Container, typename ConstIterator>
	//typename Container::iterator remove_constness(Container& c, ConstIterator it) const
	//{
	//	return c.erase(it, it);
	//}

};

class ILS_TOPTW : public ILS {
private:
	std::tuple<List<TA>::iterator, double, int, int> getBestPos(const TA&, const List<TA>&, const Vector2D<double>&) override;
	std::tuple<Walks::iterator, List<TA>::iterator, double, int, int> getBestPos(const TA&, Walks&, const Vector2D<double>&) override;
	void updateTimes(List<TA>&, const List<TA>::iterator&, const bool, const Vector2D<double>&) override;
	void updateMaxShifts(const List<TA>&, const Vector2D<double>&) override;
public:
	using ILS::ILS; //inherit constructor
	void validate(const List<List<TA>>&, const Vector2D<double>&);
	void validate(const List<TA>&, const Vector2D<double>&) override;
};

