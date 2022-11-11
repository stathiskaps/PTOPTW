#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
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

#define last_solution i==solutions.size()-1
#define first_solution i==0
#define middle_solution i>0 && i<solutions.size()-1

template <typename T>
using Vector2D = std::vector<std::vector<T>>;

template <typename T>
void printVector2D(Vector2D<T> v){
	for(auto &i : v){
		for(auto &j : i){
			std::cout << j << " ";
		}
		std::cout << std::endl;
	}
}

using Walks = std::vector<List<TA>>;
using Walk = List<TA>;

using MapOfActivities = std::map<std::string, Activity>;

class ILS{
private:
	Walk walk;
	struct Usage{
		u_int16_t imported;
		u_int16_t solved;
		u_int16_t improved;

		Usage() : imported(1), solved(1), improved(1) {}
	};

	struct Bin{
		List<TA> unvisited;
		TimeWindow tw;
	};

	struct Interval {
		double start_time;
		double end_time;
		List<TA> unvisited;
	};

	struct SR{
		int S = 1;
		int R = 1;
	};

	//std::map<std::string, std::vector<ActivityInBucket>> registry;

	int mBucketsNum;
	void dbScan(OP& op);
	void AddStartDepots(std::vector<Solution>&, const std::vector<ILS::Interval>&, const int, const OP&);
	void AddEndDepots(std::vector<Solution>&, const std::vector<Interval>&, const int, OP&);
	bool compareTimeWindowCenter(const List<TA>::iterator&, const List<TA>::iterator&);



	// virtual std::tuple<List<TA>::iterator, double, int, int> getBestPos(const TA&, const List<TA>&, const Vector2D<double>&);
	virtual std::tuple<Walks::iterator, List<TA>::iterator, double, int, int> getBestPos(const TA&, Walks&, const Vector2D<double>&);
	virtual void updateTimes(List<TA>&, const List<TA>::iterator&, const bool, const Vector2D<double>&);
	virtual std::tuple<bool, double> insertionIsValid(const TA&, const TA&, const TA&, const Vector2D<double>&);
	virtual std::tuple<bool, double> insertionIsValid(const TA& , const TA&, const Vector2D<double>&);
	std::map<std::string, std::vector<Usage>> initRegistry(List<TA>&, std::vector<ILS::Interval>);
	std::map<std::string, std::vector<double>> getActivities(List<TA>&, std::vector<ILS::Interval>);
	std::vector<Interval> getIntervals(std::vector<TA*>, int, double, double);
	void SplitSearch(std::vector<Solution>&, const std::vector<Interval>&, OP&, std::map<std::string, std::vector<Usage>>&);
	void SplitSearch2(Solution&, const std::vector<double>&, OP&, std::map<std::string, Activity>&);
	std::vector<Bin> splitUnvisited(List<TA>&, std::map<std::string, Activity>&);
	void gatherUnvisited(std::vector<Solution>&, List<TA>&);
	std::vector<List<TA>> splitUnvisitedList(std::vector<Solution>&, List<TA>&, int, std::map<std::string, std::vector<ILS::Usage>>&, std::map<std::string, std::vector<double>>);
	void LocalSearch(Solution&, ILS::Interval, OP&, std::map<std::string, std::vector<ILS::Usage>>&);
	void Controller(std::vector<Solution>&, std::vector<ILS::Interval>, OP&, std::map<std::string, std::vector<ILS::Usage>>&);
	bool hasWeightedCentroid(const Solution& sol, const int, const int);
	std::vector<Solution> splitSolution(Solution&, const std::vector<double>&, std::map<std::string, Activity>&);
	inline Point getWeightedCentroid(const List<TA>::iterator& first, const List<TA>::iterator& last);
	void SplitShake(std::vector<Solution>&, std::vector<ILS::SR>&, OP&, const int&);
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
	void validate(const std::vector<Solution>&, const Vector2D<double>&);
	inline void print(const List<TA>&);
	void SolveNew(OP&);

};

class ILS_TOPTW : public ILS {
private:
	// std::tuple<List<TA>::iterator, double, int, int> getBestPos(const TA&, const List<TA>&, const Vector2D<double>&) override;
	std::tuple<Walks::iterator, List<TA>::iterator, double, int, int> getBestPos(const TA&, Walks&, const Vector2D<double>&) override;
	void updateTimes(List<TA>&, const List<TA>::iterator&, const bool, const Vector2D<double>&) override;
	void updateMaxShifts(const List<TA>&, const Vector2D<double>&) override;
	std::tuple<bool, double> insertionIsValid(const TA&, const TA&, const TA&, const Vector2D<double>&) override;
	std::tuple<bool, double> insertionIsValid(const TA& , const TA&, const Vector2D<double>&) override;
public:
	using ILS::ILS; //inherit constructor
	void validate(const List<List<TA>>&, const Vector2D<double>&);
	void validate(const List<TA>&, const Vector2D<double>&) override;
};

