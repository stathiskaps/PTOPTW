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
#include <algorithm>
#include <plog/Log.h> // Step1: include the headers
#include "plog/Initializers/RollingFileInitializer.h"
#include "List.h"
#include "Solution.h"
#include "Divider.h"
#include "OP.h"
#include "Custom.h"
#include "Graphics.h"
//#include "boost/geometry.hpp"
#include "ygor/YgorClustering.hpp"

#ifndef ILS_H
#define ILS_H

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

class ILS{
private:

	struct Usage{
		u_int16_t imported;
		u_int16_t solved;
		u_int16_t improved;

		Usage() : imported(1), solved(1), improved(1) {}
	};

	struct TimeCuts{
		uint32_t left;
		uint32_t right;
	};

	struct Bin{
		size_t left_cut_index;
		size_t right_cut_index;
		size_t count;

		Bin(size_t left, size_t right): left_cut_index(left), right_cut_index(right), count(0) {}
	};

	struct SR{
		int S = 1;
		int R = 1;
	};

	struct Bounds{
		double minLat, minLon, maxLat, maxLon;
	};

	struct Metrics{
		std::vector<double> local_search;
		double split_unvisited;
		double shake;
		int final_pos, middle_pos;
	};


	int mIntervalsNum;
	void dbScan(OP& op);
	void AddStartDepots(std::vector<Solution>&, const std::vector<TimeWindow>&, const int, const OP&);
	void AddEndDepots(std::vector<Solution>&, const std::vector<TimeWindow>&, const int, OP&);
	bool compareTimeWindowCenter(const List<TA>::iterator&, const List<TA>::iterator&);
	static ILS* currentInstance;

	static void drawCallback(){
		currentInstance->draw();
	}

	// virtual std::tuple<List<TA>::iterator, double, int, int> getBestPos(const TA&, const List<TA>&, const Vector2D<double>&);
	virtual std::tuple<Walks::iterator, List<TA>::iterator, double, int, int> getBestPos(const TA&, Walks&, const Vector2D<double>&, 
		const std::vector<double>, const TimeWindow, const bool);
	virtual void updateTimes(List<TA>&, const List<TA>::iterator&, const bool, const Vector2D<double>&, const TimeWindow);
	virtual std::tuple<bool, double> insertionBetweenIsValid(const TA&, const TA&, const TA&, const Vector2D<double>&);
	virtual std::tuple<bool, double> insertionBeforeIsValid(const TA& , const TA&, const double, const Vector2D<double>&);
	virtual std::tuple<bool, double> insertionAfterIsValid(const TA&, const TA&, const double, const Vector2D<double>&);
	std::map<std::string, std::vector<Usage>> initRegistry(List<TA>&, std::vector<TimeWindow>);
	std::map<std::string, std::vector<double>> getActivities(List<TA>&, std::vector<TimeWindow>);
	void SplitSearch(std::vector<Solution>&, const std::vector<TimeWindow>&, OP&, std::map<std::string, std::vector<Usage>>&);
	void gatherUnvisited(std::vector<Solution>&, List<TA>&);
	std::vector<List<TA>> splitUnvisitedList(std::vector<Solution>&, List<TA>&, int, std::map<std::string, std::vector<ILS::Usage>>&, std::map<std::string, std::vector<double>>);
	bool hasWeightedCentroid(const Solution& sol, const int, const int);
	inline Point getWeightedCentroid(const List<TA>::iterator&, const List<TA>::iterator&, const int);
	int SplitShake(std::vector<Solution>&, std::vector<ILS::SR>&, OP&, const int&, const std::vector<TimeWindow>);
	int Shake(Solution&, int&, int&, OP&, const int&, const TimeWindow);
	virtual void updateMaxShifts(const List<TA>&, const Vector2D<double>&, const TimeWindow);
	std::vector<std::string> construct(Solution&, const Vector2D<double>&, const std::vector<Point>, const TimeWindow, const bool);
	int collectScores(std::vector<Solution>);
	Solution connectSolutions(std::vector<Solution>&, const size_t);
	inline int collectProfit (const List<TA>::iterator&, const List<TA>::iterator&) const;
	void printSolution(const std::string, const Solution& sol);
	void printSolutions(const std::string, const std::vector<Solution>& sols);
	void PrepareForShake(std::vector<Solution>&);
	void RemoveDummyNodes(std::vector<Solution>&);
	void InitSolutions(std::vector<Solution>&, const std::vector<TimeWindow> intervals, const OP& op);
	std::tuple<bool, double> CandidateEndDepotIsValid(const List<TA>&, const TA, TimeWindow);
	std::tuple<bool, double> CandidateStartDepotIsValid(const List<TA>&, const TA&, const double, const Vector2D<double>&);
	void drawSolutions(const std::vector<Solution>& solutions);
	std::vector<Point> getTargets(const std::vector<Solution>&, const int, const OP&);
	std::vector<TimeWindow> getIntervals(std::vector<TA>, int, double, double);
	void setupDrawCallback();

	
protected:
	Metrics metrics;
public:
	std::vector<Solution> best_solutions;
	Solution best_solution;
	void draw();
	void drawSolution(const Solution&);
    ILS();
	ILS(int);
    ~ILS();
	virtual void validate(const List<TA>&, const Vector2D<double>&, const bool);
	virtual void validate(const Walks&, const Vector2D<double>&, const bool);
	void validate(const std::vector<Solution>&, const Vector2D<double>&, const bool);
	void Solve(OP&);

};

class ILS_TOPTW : public ILS {
private:
	// std::tuple<List<TA>::iterator, double, int, int> getBestPos(const TA&, const List<TA>&, const Vector2D<double>&) override;
	std::tuple<Walks::iterator, List<TA>::iterator, double, int, int> getBestPos(const TA&, Walks&, const Vector2D<double>&, 
	const std::vector<double>, const TimeWindow, const bool) override;
	void updateTimes(List<TA>&, const List<TA>::iterator&, const bool, const Vector2D<double>&, const TimeWindow) override;
	void updateMaxShifts(const List<TA>&, const Vector2D<double>&, const TimeWindow) override;
	std::tuple<bool, double> insertionBetweenIsValid(const TA&, const TA&, const TA&, const Vector2D<double>&) override;
	std::tuple<bool, double> insertionBeforeIsValid(const TA& , const TA&, const double, const Vector2D<double>&) override;
	std::tuple<bool, double> insertionAfterIsValid(const TA&, const TA&, const double, const Vector2D<double>&) override;
public:
	using ILS::ILS; //inherit constructor
	void validate(const List<TA>&, const Vector2D<double>&, const bool) override;
};



#endif