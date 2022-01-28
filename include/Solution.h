#pragma once
#include <stdexcept>
#include "List.h"



class Solution {
	friend struct Cluster;
	friend class ILS;
	friend class ILS_OPTW;
	friend class ILS_TOPTW;

	ListTA mUnvisited;
	Walk mWalk;

public:
	Solution(ListTA unvisited) {
		mUnvisited = unvisited;
	};
	Solution(){
	};
	Solution(TA* depot, ListTA unvisited, double openTime, double closeTime) {
		mUnvisited = unvisited;
		mWalk.pushClone(depot);
		mWalk.pushClone(depot);

		mWalk.first()->arrTime = openTime;
		mWalk.first()->startOfVisitTime = openTime;
		mWalk.first()->depTime = openTime;
		mWalk.last()->timeWindow.closeTime = closeTime;
		mWalk.last()->maxShift = closeTime - mWalk.last()->depTime;
	};
	Solution(Walk walk, ListTA unvisited, double openTime, double closeTime) {
		mUnvisited = unvisited;
		mWalk = walk;

		mWalk.first()->arrTime = openTime;
		mWalk.first()->startOfVisitTime = openTime;
		mWalk.first()->depTime = openTime;
		mWalk.last()->timeWindow.closeTime = closeTime;
		mWalk.last()->maxShift = closeTime - mWalk.last()->depTime;
	};

	Solution(Walk walk, ListTA unvisited) {
		mUnvisited = unvisited;
		mWalk = walk;
	}

	Solution(Point p, double startTime, double endTime) {
		if (startTime > endTime) {
			throw std::invalid_argument("arguments are invalid");
		}
		TA* depot = new TA(p);
		depot->depTime = startTime;
		depot->timeWindow = TimeWindow{ startTime, endTime };
		mWalk.push(depot);
		depot = new TA(p);
		depot->timeWindow = TimeWindow{ startTime, endTime };
		depot->maxShift = endTime - startTime;
	}

	Solution(TA* start, TA* end, ListTA unvisited, double startTime, double endTime) {
		TA* startDepot = start->clone();
		TA* endDepot = end->clone();

		startDepot->depTime = startTime;
		startDepot->timeWindow = TimeWindow{ startTime, endTime };
		
		endDepot->timeWindow = TimeWindow{ startTime, endTime };
		endDepot->maxShift = endTime - startTime;

		mWalk.push(startDepot);
		mWalk.push(endDepot);

		mUnvisited = unvisited;
	}


	~Solution() {};
	void SetRoute(Walk route) { mWalk = route; }
	Walk GetRoute() { return mWalk; }
	
	void copy(Solution sol) {
		this->mUnvisited = sol.mUnvisited.copy();
		this->mWalk = sol.mWalk.copy();
	}

	void reset() {
		this->mUnvisited.empty();
		this->mWalk.empty();
	}

	void print() {
		std::cout << "Unvisited:\t"; mUnvisited.print();
		std::cout << "Route:\t" ; mWalk.print();
	}

	int getScore() {
		return mWalk.collectProfit();
	}

	void print(std::string msg) {
		std::cout << "||======================================================================||" << std::endl;
		std::cout << msg << std::endl;
		std::cout << "----------------------------------------------------------------------" << std::endl;
		mUnvisited.print("Unvisited:");
		mWalk.print("Route:");
		std::cout << "||======================================================================||" << std::endl;
	}

	TA* toTA() {
		mWalk.removeNodesWithId(DEFAULT_DEPOT_ID);
		std::vector<TA> nodes = mWalk.toVec();
		double duration = mWalk.getDuration();
		if (nodes.size() == 0) {
			return nullptr;
		} else if (nodes.size() == 1) {
			return new Sight(
				id::generate(),
				nodes.front().point,
				duration,
				getScore(),
				mWalk.first()->depTime,
				mWalk.last()->timeWindow.closeTime
			);
		}

		return new Route(
			id::generate(),
			Point(),
			nodes.front().point,
			nodes.back().point,
			duration,
			getScore(),
			mWalk.first()->depTime,
			mWalk.last()->timeWindow.closeTime);
		
	}

	

};