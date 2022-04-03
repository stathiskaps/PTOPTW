#include "OP.h"

OP::OP() {
	
}

OP::~OP() {
	
}

OP::OP(std::vector<TA*> attractions, std::vector<Point> points, TA* startDepot, TA* endDepot, int walksNum) : mAttractions(attractions), mPoints(points), m_walks_num(walksNum) {
	std::tuple<std::vector<std::vector<double>>, double> tuple = calcTravelTimes(points); //TODO: delete pointer
	mTravelTimes = std::get<0>(tuple);
	mStartDepot = startDepot->clone();
	mEndDepot = endDepot->clone();
}

double OP::GetEuclideanDistance(int x1, int y1, int x2, int y2) {
	double eu_dist;
	double nearest;
	eu_dist = pow(x2 - x1, 2) + pow(y2 - y1, 2);
	eu_dist = sqrt(eu_dist);
	nearest = (double)roundf(eu_dist * 100) / 100;
	return nearest;
}

std::tuple<std::vector<std::vector<double>>, double> OP::calcTravelTimes(std::vector<Point>& points)
{
	size_t pointsSize = points.size();
	std::vector<std::vector<double>> ttMatrix;
	double totalTravelTime = 0;
	double meanTravelTime;
	double val;
	std::cout << std::endl;
	for (int i = 0; i < pointsSize; ++i) {
		std::vector<double> vec;
		for (int j = 0; j < pointsSize; j++) {
			val = GetEuclideanDistance(points.at(i).pos.lat, points.at(i).pos.lon, points.at(j).pos.lat, points.at(j).pos.lon);
			vec.push_back(val);
			totalTravelTime += val;
		}
		ttMatrix.push_back(vec);
	}

	meanTravelTime = totalTravelTime / ((double)pointsSize * pointsSize);

	return std::make_tuple(ttMatrix, meanTravelTime);
}

void OP::AddPointToGraph(Point& p) {
	size_t pointsSize = mPoints.size();
	p.id = pointsSize;
	mPoints.push_back(p);
	std::vector<double> vec;

	for (size_t i = 0; i < pointsSize; ++i) {
		double dist = mPoints[i].euclidean_distance(p);
		mTravelTimes[i].push_back(dist);
		vec.push_back(dist);
	}
	vec.push_back(0);

	mTravelTimes.push_back(vec);
}

void OP::PrintTravelTimes(std::string label) {
	size_t pointsSize = mTravelTimes.size();
	std::cout << label << std::endl;
	std::cout << "\t";
	for (int i = 0; i < pointsSize; ++i) {
		std::cout << i << "\t";
	}
	std::cout << std::endl;
	for (size_t i = 0; i < pointsSize; ++i) {
		std::cout << i << ":\t";
		for (size_t j = 0; j < pointsSize; ++j) {
			std::cout << mTravelTimes[i][j] << "\t";
		}
		std::cout << std::endl;
	}
}
