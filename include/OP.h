#pragma once

#include <iostream>
#include <vector>
#include "List.h"

class OP {
	friend class ILS;
	friend class ILS_OPTW;
	friend class ILS_TOPTW;
private:
	std::vector<TA*> mAttractions;
	std::vector<Point> mPoints;
	std::vector<std::vector<double>> mTravelTimes;
	TA mStartDepot, mEndDepot;
	size_t m_walks_num;
	TimeWindow mTimeWindow;

	std::tuple<std::vector<std::vector<double>>, double> calcTravelTimes(std::vector<Point>& points);

public:
	OP();
	~OP();
	OP(std::vector<TA*>, std::vector<Point>, TA, TA, int, double, double);
	OP(std::vector<TA*>, std::vector<Point>, TA, TA, int, double, double, std::vector<std::vector<double>>);
	void AddPointToGraph(Point&);
	void PrintTravelTimes(std::string);

};