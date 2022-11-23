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
#include "ILS.h"
#include "OP.h"

int S::ta_id;

std::vector<std::string> split(const std::string& line) {
	std::string buf;                 // Have a buffer string
	std::stringstream ss(line);       // Insert the string into a stream

	std::vector<std::string> tokens; // Create vector to hold our words

	while (ss >> buf)
		tokens.push_back(buf);

	return tokens;
}







double calcMeanVisitTime(std::vector<TA*> touristAttractions) {
	double totalVisitDuration = 0;
	for (auto& ta : touristAttractions) {
		totalVisitDuration += ta->visitDuration;
	}
	return totalVisitDuration / (touristAttractions.size() - 1); //don't count the depot
}

void init(std::string filename, int numRoutes, int numIntervals) {

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

	// raise(SIGINT);


#if 1
	std::cout << std::endl;
	TA* depot = touristAttractions[0];
	TA* start_depot = depot->clone();
	start_depot->id = START_DEPOT_ID;
	TA* end_depot = depot->clone();
	end_depot->id = END_DEPOT_ID;
	touristAttractions.erase(touristAttractions.begin());

	OP op = OP(touristAttractions, points, start_depot, end_depot, numRoutes);

	ILS_TOPTW ilstoptw = ILS_TOPTW(numRoutes, numIntervals);
	ilstoptw.SolveNew(op);

	delete start_depot;
	delete end_depot;

	/*OPTW optw(touristAttractions, ttMatrix, depot, OPEN_DAY_TIME, CLOSE_DAY_TIME);
	Solution sol = optw.solve();
	bool valid = optw.validate();
	std::string msg = valid ? "yes" : "no";
	std::cout << "valid solution? " << msg << std::endl;
	sol.print();*/
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
	int k, numRoutes, numIntervals;
	std::string filename;

	std::cout << "Enter the file's name(without txt): ";
	std::cin >> filename;

	std::cout << "Enter the number of the routes: ";
	std::cin >> numRoutes;

	std::cout << "Enter the number of intervals that the problem will be divided: ";
	std::cin >> numIntervals;

	init(filename, numRoutes, numIntervals);

	return 0;
}



