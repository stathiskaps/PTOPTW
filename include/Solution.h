#pragma once
#include <stdexcept>
#include <climits>
#include "List.h"
#include "Custom.h"

class Solution {
	friend class ILS;
	friend class ILS_OPTW;
	friend class ILS_TOPTW;
	List<TA> m_unvisited;
	std::vector<List<TA>> m_walks;
public:
	Solution() {}

	Solution(List<TA> unvisited): m_unvisited(unvisited) {}

	//todo: check what explicit is
	Solution(TA start, TA end, List<TA> unvisited, double startTime, double endTime, int walksNum) : m_unvisited(unvisited) {
		for (int i = 0; i < walksNum; ++i) {
			List<TA> walk{ start, end };
			walk.front().depTime = startTime;
			walk.front().timeWindow = TimeWindow{ startTime, endTime };
			walk.back().timeWindow = TimeWindow{ startTime, endTime };
			walk.back().maxShift = endTime - startTime;
			m_walks.push_back(walk);
		}
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

	List<TA>& getUnvisited() { return m_unvisited; }

};