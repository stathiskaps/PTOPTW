#include <iostream>
#include <vector>
#include <tuple>
#include <climits>
#include <numeric>
#include <signal.h>
#include "json.hpp"
#include "Solution.h"
#include "TaskManager.h"

#ifndef OP_H
#define OP_H

//We can use OP to model a problem of Tourist Attractions or Clusters
//This OP is used to solve the TOP, so no time windows, no multiple routes
class OP {
private:
	std::tuple<double, double, double> calcTimeEventCut();
	virtual std::tuple<int, double, int, int> getBestPos(TA*);
	virtual bool updateTimes(int, bool);
	int Insert();
	int LocalSearch(double);
	void Local(double);
	void setBucketActivityDurations();
	Walk Construct();
	void Shake(int, int, int);
	void SaveSolution(Solution);
	void print(std::string);
	void SendSolution(Solution);
	static bool compareByProfit(const TA a, const TA b);
protected:
	TA* mDepot;
	Solution mProcessSolution, mBestSolution;
	std::vector<std::vector<double>> mTtMatrix;
	double mStartTime, mEndTime;
public:
	OP();
	OP(std::vector<TA*>, std::vector<std::vector<double>>, TA*, double, double);
	OP(std::vector<TA*>, Walk, std::vector<std::vector<double>>, double, double);
	~OP();
	Solution solve();
	virtual bool validate();

};

class OPTW : public OP {
private:
	std::tuple<int, double, int, int> getBestPos(TA*) override;
	bool updateTimes(int, bool) override;
public:
	using OP::OP; //inherit constructor	
	bool validate() override;
};

class TOPTW : public OP {
private:
	std::tuple<int, double, int, int> getBestPos(TA*) override;
	bool updateTimes(int, bool) override;
public:
	using OP::OP; //inherit constructor
	bool validate() override;
};

#endif
