#include <stdexcept>
#include "List.h"

#ifndef SOLUTION_H
#define SOLUTION_H

class Solution {
	friend class OP;
	friend class ILSOP;
	friend class OPTW;
	friend class TOPTW;
	friend struct Cluster;

	ListTA mUnvisited;
	Walk mWalk;
	int mScore, mOpenTime, mCloseTime;

public:
	Solution(ListTA unvisited) {
		mUnvisited = unvisited;
		mScore = 0;
		mOpenTime = 0;
		mCloseTime = 0;
	};
	Solution(){
		mScore = 0;
		mOpenTime = 0;
		mCloseTime = 0;
	};
	Solution(TA depot, ListTA unvisited, double openTime, double closeTime) {
		mScore = 0;
		mOpenTime = openTime;
		mCloseTime = closeTime;
		mUnvisited = unvisited;
		mWalk.pushNew(&depot);
		mWalk.pushNew(&depot);

		mWalk.first()->arrTime = openTime;
		mWalk.first()->startOfVisitTime = openTime;
		mWalk.first()->depTime = openTime;
		mWalk.last()->timeWindow.closeTime = closeTime;
		mWalk.last()->maxShift = closeTime - openTime;
	};
	~Solution() {};
	void SetScore(int score) { mScore = score; }
	int GetScore() { return mScore; }
	void SetRoute(Walk route) { mWalk = route; }
	Walk GetRoute() { return mWalk; }
	
	void copy(Solution sol) {
		this->mUnvisited = sol.mUnvisited.copy();
		this->mWalk = sol.mWalk.copy();
		this->mScore = sol.mScore;
		this->mOpenTime = sol.mOpenTime;
		this->mCloseTime = sol.mCloseTime;
	}

	void reset() {
		this->mUnvisited.empty();
		this->mWalk.empty();
		this->mScore = 0;
		this->mOpenTime = 0;
		this->mCloseTime = 0;
	}

	void print() {
		std::cout << "Unvisited:\t"; mUnvisited.print();
		std::cout << "Route:\t" ; mWalk.print();
		std::cout << "Score: " << mScore << std::endl;
		std::cout << "OpenTime:" << mOpenTime << std::endl;
		std::cout << "CloseTime:" << mCloseTime << std::endl;
	}

	void print(std::string msg) {
		std::cout << msg << std::endl;
		mUnvisited.print("Unvisited:");
		mWalk.print("Route:");
		std::cout << "Score: " << mScore << std::endl;
		std::cout << "OpenTime:" << mOpenTime << std::endl;
		std::cout << "CloseTime:" << mCloseTime << std::endl;
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
				mScore,
				mOpenTime,
				mCloseTime
			);
		}

		return new Route(
			id::generate(),
			Point(),
			nodes.front().point,
			nodes.back().point,
			duration,
			mScore,
			mOpenTime,
			mCloseTime);
		
	}

	

};

#endif