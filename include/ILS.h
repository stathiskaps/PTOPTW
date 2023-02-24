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
#include <fstream>
#include <algorithm>
#include "List.h"
#include "Solution.h"
#include "OP.h"
#include "Custom.h"
#include "Graphics.h"

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
public:
	struct Configuration{
		bool time_limited_execution = false, write_results = false, write_solution = false, graphics = false;
		double execution_time_limit = 0;
	};

	void draw();
	void drawSolution(const Solution&);
    ILS();
	ILS(int, std::string, Configuration);
    ~ILS();
	std::pair<int, double> Solve(OP&);
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

	struct Metrics{
		std::vector<double> local_search;
		double split_unvisited;
		double shake;
		double validation_time;
		double total_execution_time;
		int final_pos, middle_pos;
		int second_phase_counter;
		int second_phase_window_sum;
		int second_phase_improved;
		int comparisons;
	};

	int mIntervalsNum;
	std::string mInstance;
	static ILS* currentInstance;
	Configuration mConf;
	std::vector<Solution> best_solutions;
	Solution best_solution;

	static void drawCallback(){
		currentInstance->draw();
	}

	void AddStartDepots(std::vector<Solution>&, const std::vector<TimeWindow>&, const int, const OP&);
	bool compareTimeWindowCenter(const List<TA>::iterator&, const List<TA>::iterator&);
	std::tuple<Walks::iterator, List<TA>::iterator, double, int, int> getBestPos(const TA&, Walks&, const Vector2D<double>&, 
		const std::vector<double>, const TimeWindow, const bool);
	void updateBasicTimes(List<TA>::iterator&, const Vector2D<double>&, const TimeWindow);
	bool updateTimes(List<TA>&, const List<TA>::iterator&, const bool, const Vector2D<double>&, const TimeWindow);
	void updateMaxShifts(const List<TA>&, const Vector2D<double>&, const TimeWindow);
	std::map<std::string, std::vector<Usage>> initRegistry(List<TA>&, std::vector<TimeWindow>);
	std::map<std::string, std::vector<double>> getActivities(List<TA>&, std::vector<TimeWindow>);
	void SplitSearch(std::vector<Solution>&, List<TA>& pool, const std::vector<TimeWindow>&, OP&, std::map<std::string, std::vector<Usage>>&);
	void SearchBetween(std::vector<Solution>&, const size_t, List<TA>& pool, const std::vector<TimeWindow>&, OP&, std::map<std::string, std::vector<Usage>>&);
	void gatherUnvisited(std::vector<Solution>&, List<TA>&);
	std::vector<List<TA>> splitUnvisitedList(std::vector<Solution>&, List<TA>&, int, std::map<std::string, std::vector<ILS::Usage>>&, std::map<std::string, std::vector<double>>);
	bool hasWeightedCentroid(const Solution& sol, const int, const int);
	inline Point getWeightedCentroid(const List<TA>::iterator&, const List<TA>::iterator&, const int);
	int SplitShake(std::vector<Solution>&, std::vector<ILS::SR>&, OP&, const std::vector<TimeWindow>);
	int Shake(Solution&, int&, int&, OP&, const TimeWindow);
	std::vector<std::string> construct(Solution&, const Vector2D<double>&, const std::vector<Point>, const TimeWindow, const bool);
	int collectScores(const std::vector<Solution>&) const;
	Solution connectSolutions(std::vector<Solution>&, const size_t);
	inline int collectProfit (const List<TA>::iterator&, const List<TA>::iterator&) const;
	void printSolution(const std::string, const Solution& sol);
	void printSolutions(const std::string, const std::vector<Solution>&);
	void PrepareForShake(std::vector<Solution>&);
	void RemoveDummyNodes(std::vector<Solution>&);
	void InitSolutions(std::vector<Solution>&, const std::vector<TimeWindow>, const OP&);
	std::tuple<bool, double> CandidateStartDepotIsValid(const List<TA>&, const TA&, const double, const Vector2D<double>&);
	std::vector<Point> getTargets(const std::vector<Solution>&, const int, const OP&);
	std::vector<TimeWindow> getIntervals(std::vector<TA>, int, double, double);
	TA getValidPreviousTA(std::vector<Solution>&, const int, const size_t);
	void RemoveUnfeasibleVisits(std::vector<Solution>&, const int, const size_t);
	void connectAndValidateSolutions(const std::vector<Solution>&, const OP&);
	size_t countNodes(const std::vector<Solution>&);
	void printMetrics();
	void setupDrawCallback();
	void checkSolutions(std::vector<Solution>&,  const std::vector<TimeWindow>& intervals, const OP&);
	std::vector<std::string> fixWalk(List<TA>&, const OP&, TimeWindow time_budget);
	void validateDirectedSolution(const Solution&, const OP& op, const bool);
	void validateTimes(const std::vector<Solution>&, const Vector2D<double>&, const bool);
	void validateTimes(const Solution&, const Vector2D<double>&, const bool);
	void validateTimes(const List<TA>&, const Vector2D<double>&, const bool);

protected:
	Metrics metrics;


};



#endif