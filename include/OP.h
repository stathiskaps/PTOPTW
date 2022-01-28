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
	TA* mStartDepot;
	TA* mEndDepot;

	std::tuple<std::vector<std::vector<double>>, double> calcTravelTimesMatrix2(std::vector<Point>& points);
	double GetEuclideanDistance2(int x1, int y1, int x2, int y2);	//returns euclidean distance

public:
	OP();
	~OP();
	OP(std::vector<TA*>, std::vector<Point>, TA*, TA*);
	void AddPointToGraph(Point);

};