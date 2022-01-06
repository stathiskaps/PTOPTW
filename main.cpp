#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <chrono>
#include <cmath>
#include <variant>
#include <cassert>
#include "Controller.h"
#include "Definitions.h"
#include "Divide.h"

int S::ta_id;

std::vector<std::string> split(const std::string& line) {
	std::string buf;                 // Have a buffer string
	std::stringstream ss(line);       // Insert the string into a stream

	std::vector<std::string> tokens; // Create vector to hold our words

	while (ss >> buf)
		tokens.push_back(buf);

	return tokens;
}

double GetEuclideanDistance2(int x1, int y1, int x2, int y2)	//returns euclidean distance
{
	double eu_dist;
	double nearest;
	eu_dist = pow(x2 - x1, 2) + pow(y2 - y1, 2);
	eu_dist = sqrt(eu_dist);
	nearest = (double)roundf(eu_dist * 100) / 100;
	return nearest;
}

std::tuple<std::vector<std::vector<double>>, double> calcTravelTimesMatrix2(std::vector<Point>& points)
{
	size_t pointsSize = points.size();
	std::vector<std::vector<double>> ttMatrix;
	double totalTravelTime = 0;
	double meanTravelTime;
	double val;
	for (int i = 0; i < pointsSize; ++i) {
		std::vector<double> vec;
		for (int j = 0; j < pointsSize; j++) {
			val = GetEuclideanDistance2(points.at(i).pos.lat, points.at(i).pos.lon, points.at(j).pos.lat, points.at(j).pos.lon);
			vec.push_back(val);
			totalTravelTime += val;
		}
		ttMatrix.push_back(vec);
	}

	meanTravelTime = totalTravelTime / ((double)pointsSize * pointsSize);

	return std::make_tuple(ttMatrix, meanTravelTime);
}





double calcMeanVisitTime(std::vector<TA*> touristAttractions) {
	double totalVisitDuration = 0;
	for (auto& ta : touristAttractions) {
		totalVisitDuration += ta->visitDuration;
	}
	return totalVisitDuration / (touristAttractions.size() - 1); //don't count the depot
}

void init(std::string filename, int numRoutes) {

	std::vector<TA*> touristAttractions; //TODO:delete pointers
	std::vector<Point> points;
	std::string path = "./instances/Cordeau/";
	S::ta_id = 0;

	std::ifstream infile(path.append(filename + ".txt"));

	std::string line;

	std::getline(infile, line);
	std::getline(infile, line);
	std::vector<std::string> poi_data;

	int pointId = 0;

	while (std::getline(infile, line))
	{
		if (line.front() == '#') {
			continue;
		}
		poi_data = split(line);
		if (line.front() == 'R') {
			Point p = Point(pointId++, std::stod(poi_data[2]), std::stod(poi_data[3]));
			Point edge1 = Point(pointId++, std::stod(poi_data[4]), std::stod(poi_data[5]));
			Point edge2 = Point(pointId++, std::stod(poi_data[6]), std::stod(poi_data[7]));

			points.push_back(p);
			points.push_back(edge1);
			points.push_back(edge2);

			touristAttractions.push_back(new Route(
				id::generate(),
				p,
				edge1,
				edge2,
				std::stoi(poi_data[8]),
				std::stoi(poi_data[9]),
				std::stoi(poi_data[12 + std::stoi(poi_data[11])]),
				std::stoi(poi_data[13 + std::stoi(poi_data[11])])
			));
		}
		else {
			if (poi_data[0].empty()) {
				Point p = Point(pointId++, std::stod(poi_data[2]), std::stod(poi_data[3]));
				points.push_back(p);
				touristAttractions.push_back(new Sight(
					id::generate(),
					p,
					std::stoi(poi_data[4]),
					std::stoi(poi_data[5]),
					std::stoi(poi_data[8 + std::stoi(poi_data[7])]),
					std::stoi(poi_data[9 + std::stoi(poi_data[7])])
				));
			} else {
				Point p = Point(pointId++, std::stod(poi_data[1]), std::stod(poi_data[2]));
				points.push_back(p);
				touristAttractions.push_back(new Sight(
					id::generate(),
					p,
					std::stoi(poi_data[3]),
					std::stoi(poi_data[4]),
					std::stoi(poi_data[7 + std::stoi(poi_data[6])]),
					std::stoi(poi_data[8 + std::stoi(poi_data[6])])
				));
			}
		}
		// Strip of the comments.
		
	}

	//double meanShift = std::get<1>(tuple) + calcMeanVisitTime(touristAttractions);

	// Divide divide(touristAttractions);
	// divide.exec();

	// raise(SIGINT);


#if 1
	std::tuple<std::vector<std::vector<double>>, double> tuple = calcTravelTimesMatrix2(points); //TODO: delete pointer
	std::vector<std::vector<double>> ttMatrix = std::get<0>(tuple);

	TA* depot = touristAttractions.at(0);
	touristAttractions.erase(touristAttractions.begin());

	OPTW optw(touristAttractions, ttMatrix, depot, OPEN_DAY_TIME, CLOSE_DAY_TIME);
	Solution sol = optw.solve();
	bool valid = optw.validate();
	std::string msg = valid ? "yes" : "no";
	std::cout << "valid solution? " << msg << std::endl;
	sol.print();
#endif
	
#if 0
	Controller controller(OPEN_DAY_TIME, CLOSE_DAY_TIME, touristAttractions, points, numRoutes);
	controller.exec();
#endif

	for (auto p : touristAttractions) {
		delete p;
	}

	touristAttractions.clear();

}

int main(int argc, char** argv)
{
	int k, numRoutes;
	std::string filename;

	std::cout << "Enter the file's name(without txt): ";
	std::cin >> filename;

	std::cout << "Enter the number of the routes: ";
	std::cin >> numRoutes;

	init(filename, numRoutes);

	std::cin >> k;

	return 0;
}



