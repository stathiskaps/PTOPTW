#pragma once
#include <stdexcept>
#include "List.h"
#include "Custom.h"

class CustomSolution {
	friend class ILS;
	friend class ILS_OPTW;
	friend class ILS_TOPTW;
	CustomList<TA> m_unvisited, m_walk;
	CustomList<CustomList<TA>> m_walks;
public:
	CustomSolution() {}
	CustomSolution(CustomList<TA> unvisited, CustomList<TA> walk) : m_unvisited(unvisited), m_walk(walk) {}
	explicit CustomSolution(TA start, TA end, CustomList<TA> unvisited, double startTime, double endTime) : m_walk{ start, end }, m_unvisited(unvisited) {
		m_walk.front().depTime = startTime;
		m_walk.front().timeWindow = TimeWindow{ startTime, endTime };
		m_walk.back().timeWindow = TimeWindow{ startTime, endTime };
		m_walk.back().maxShift = endTime - startTime;
	}

	CustomSolution(TA start, TA end, CustomList<TA> unvisited, double startTime, double endTime, int walksNum) : m_unvisited(unvisited) {
		for (int i = 0; i < walksNum; ++i) {
			CustomList<TA> walk{ start, end };
			walk.front().depTime = startTime;
			walk.front().timeWindow = TimeWindow{ startTime, endTime };
			walk.back().timeWindow = TimeWindow{ startTime, endTime };
			walk.back().maxShift = endTime - startTime;
			m_walks.push_back(walk);
		}
	}

	inline int getScore() {
		int sum{}; for (auto& p : m_walk) sum += p.profit; return sum;
	}

	inline int getScores() {
		int sum{}; 
		for (auto& w : m_walks) {
			for (auto& p : w) {
				sum += p.profit;
			}
		}
		return sum;
	}

	inline int getVisits() {
		int visits{};
		for (auto& w : m_walks) visits += w.size() - 2;
		return visits;
	}

	inline int getMinWalkSize() const{
		int min_size = INT_MAX;
		for (auto& walk : m_walks) {
			if (walk.size() < min_size) min_size = walk.size();
		}
		return min_size;
	}

	CustomList<TA>& getUnvisited() { return m_unvisited; }
	CustomList<TA>& getWalk() { return m_walk; }

};

class Solution {
	friend struct Cluster;
	friend class ILS;
	friend class ILS_OPTW;
	friend class ILS_TOPTW;

	ListTA mUnvisited;
	ListTA mWalk;

public:
	Solution(ListTA unvisited) {
		mUnvisited = unvisited;
		mWalk = ListTA();
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
	Solution(ListTA walk, ListTA unvisited, double openTime, double closeTime) {
		mUnvisited = unvisited;
		mWalk = walk;

		mWalk.first()->arrTime = openTime;
		mWalk.first()->startOfVisitTime = openTime;
		mWalk.first()->depTime = openTime;
		mWalk.last()->timeWindow.closeTime = closeTime;
		mWalk.last()->maxShift = closeTime - mWalk.last()->depTime;
	};

	Solution(ListTA walk, ListTA unvisited) {
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
		mWalk.pushBack(depot);
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

		mWalk.pushBack(startDepot);
		mWalk.pushBack(endDepot);

		mUnvisited = unvisited;
	}


	~Solution() {};
	void SetRoute(ListTA route) { mWalk = route; }
	ListTA GetRoute() { return mWalk; }
	
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
		mWalk.removeById(DEFAULT_DEPOT_ID);
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