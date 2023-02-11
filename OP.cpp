#include "OP.h"

OP::OP() {
	
}

OP::~OP() {
}

OP::OP(std::vector<TA> attractions, std::vector<Point> points, 
	TA startDepot, TA endDepot, int walksNum, 
	double startTime, double endTime) : mAttractions(attractions), mPoints(points), m_walks_num(walksNum) {
	std::tuple<std::vector<std::vector<double>>, double> tuple = calcTravelTimes(points); //TODO: delete pointer
	mTravelTimes = std::get<0>(tuple);
	mStartDepot = startDepot;
	mEndDepot = endDepot;
	mTimeWindow = TimeWindow{startTime, endTime};
}

OP::OP(std::vector<TA> attractions, std::vector<Point> points, 
	TA startDepot, TA endDepot, int walksNum, 
	double startTime, double endTime, 
	std::vector<std::vector<double>> travelTimes) : mAttractions(attractions), mPoints(points), m_walks_num(walksNum), mTravelTimes(travelTimes) {
	mStartDepot = startDepot;
	mEndDepot = endDepot;
	mTimeWindow = TimeWindow{startTime, endTime};
}

std::tuple<std::vector<std::vector<double>>, double> OP::calcTravelTimes(std::vector<Point>& points){
	size_t pointsSize = points.size();
	std::vector<std::vector<double>> ttMatrix;
	double totalTravelTime = 0;
	double meanTravelTime;
	double val;
	std::cout << std::endl;
	for (int i = 0; i < pointsSize; ++i) {
		std::vector<double> vec;
		for (int j = 0; j < pointsSize; j++) {
			val = points.at(i).euclidean_distance(points.at(j));
			// val = GetEuclideanDistance(points.at(i).pos.lat, points.at(i).pos.lon, points.at(j).pos.lat, points.at(j).pos.lon);
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
