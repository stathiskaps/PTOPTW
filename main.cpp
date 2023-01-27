#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <chrono>
#include <unistd.h>
#include <getopt.h>
#include <cmath>
#include <variant>
#include <cassert>
#include "Definitions.h"
#include "ILS.h"
#include "OP.h"

enum InstanceType { fixed, custom };

struct RouteDuration{
	int id, from, to;
	double duration;
};

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

void init(std::string folder, std::string filename, int numRoutes, int numIntervals, InstanceType instance_type) {

	std::vector<TA*> touristAttractions; //TODO:delete pointers
	std::vector<Point> points;

	std::string filepath = "./instances/"+folder+"/"+filename+".txt";
	std::ifstream infile(filepath);
	std::string line;

	std::getline(infile, line);
	std::getline(infile, line);
	std::vector<std::string> poi_data;

	int pointId = 0;
	std::vector<RouteDuration> durations;

	while (std::getline(infile, line))
	{
		poi_data = split(line);

		switch (line.front()) {
		case 'R': {
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
			break;
		}
		case 'D': {
			int id = std::stoi(poi_data[1]), from = std::stoi(poi_data[2]), to = std::stoi(poi_data[3]);
			double duration = std::stod(poi_data[4]);
			durations.push_back(RouteDuration{id, from , to, duration});
			break;
		}
		case '#': {
			continue;
			break;
		}
		default: {
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
			break;
		}
		}

	}

#if 1
	std::cout << std::endl;
	TA* depot = touristAttractions[0];
	TA* start_depot = depot->clone();
	start_depot->id = START_DEPOT_ID;
	TA* end_depot = depot->clone();
	end_depot->id = END_DEPOT_ID;
	touristAttractions.erase(touristAttractions.begin());

	OP op;
	if(instance_type == custom){
		const size_t length = std::sqrt(durations.size());
		std::vector<std::vector<double>> travel_times (length);
		for(size_t i = 0; i < length; ++i){
			travel_times[i] = std::vector<double> (length);
		}
		for(auto& d : durations){
			travel_times[d.from][d.to] = d.duration;
		}
		op = OP(touristAttractions, points, start_depot, end_depot, numRoutes, end_depot->timeWindow.openTime, end_depot->timeWindow.closeTime, travel_times);
	} else {
		op = OP(touristAttractions, points, start_depot, end_depot, numRoutes, end_depot->timeWindow.openTime, end_depot->timeWindow.closeTime);
	}
	

	ILS_TOPTW ilstoptw = ILS_TOPTW(numIntervals);
	ilstoptw.Solve(op);

#endif

	for (auto p : touristAttractions) {
		delete p;
	}

	touristAttractions.clear();

}


int main(int argc, char** argv) {

	InstanceType instance_type = fixed; 

	std::string folder;
	std::string instance;
	int num_of_walks;
	int num_of_intervals;

	int c;
	int option_index = 0;
	bool folder_option_provided = false, instance_option_provided = false, 
		walks_option_provided = false, subs_option_provided = false;
	static struct option long_options[] = {
		{"help", no_argument, 0, 'h'},
		{"custom", no_argument, 0, 'c'},
		{"folder", required_argument, 0, 'f'},
		{"instance", required_argument, 0, 'i'},
		{"walks", required_argument, 0, 'w'},
		{"subs", required_argument, 0, 's'},
		{0, 0, 0, 0}
	};

	while ((c = getopt_long(argc, argv, "hcf:i:w:s:", long_options, &option_index)) != -1) {
		switch (c) {
		case 'h': {
			std::cout << "Please run the program with the following options: "<< std::endl <<
			"./AMTOPTW -f <folder> -i <instance> -w <number of walks> -i <number of sub-intervals> " << std::endl <<
			"or" << std::endl <<
			"./AMTOPTW -folder=\"<folder>\" --instance=\"<instance>\" --walks=<num_of_walks> --intervals=<num_of_intervals>" << std::endl << 
			"e.g. ./AMTOPTW -f Cordeau -i pr11 -w 4 -s 4" << std::endl << 
			"Also, you can add the option -c without any arguments if you want to use custom travel times from a custom made topology" << std::endl;

			return 0;
		}
		case 'c':{
			std::cout << "Will use custom times" << std::endl;
			instance_type = custom;
			break;
		}
		case 'f': {
			if(optarg) {
				folder = optarg;
				folder_option_provided = true;
			} else {
				std::cerr << "Error: argument for option -f is required" << std::endl;
				exit(EXIT_FAILURE);
			}
			break;
		}
		case 'i': {
			if(optarg) {
				instance = optarg;
				instance_option_provided = true;
			} else{
				std::cerr << "Error: argument for option -i is required" << std::endl;
				exit(EXIT_FAILURE);
			}
			break;
		}
		case 'w': {
			if(optarg) {
				num_of_walks = std::stoi(optarg);
				walks_option_provided = true;
			} else{
				std::cerr << "Error: argument for option -w is required" << std::endl;
				exit(EXIT_FAILURE);
			}
			break;
		}
		case 's': {
			if(optarg) {
				num_of_intervals = std::stoi(optarg);
				subs_option_provided = true;
			} else{
				std::cerr << "Error: argument for option -s is required" << std::endl;
				exit(EXIT_FAILURE);
			}
			break;
		}
		default:
			break;
		}
	}

	if (!folder_option_provided) {
		std::cerr << "Error: Option -f is required" << std::endl;
		exit(EXIT_FAILURE);
	}

	if (!instance_option_provided) {
		std::cerr << "Error: Option -i is required" << std::endl;
		exit(EXIT_FAILURE);
	}

	if (!walks_option_provided) {
		std::cerr << "Error: Option -w is required" << std::endl;
		exit(EXIT_FAILURE);
	}

	if (!subs_option_provided) {
		std::cerr << "Error: Option -s is required" << std::endl;
		exit(EXIT_FAILURE);
	}

	// glutInit(&argc, argv);
	// glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE );
	// glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	// glutCreateWindow( "TSP Opengl" ) ;

	// glutDisplayFunc(displayOther ) ;
	// glutIdleFunc(displayOther);
  	// glutReshapeFunc( onResize );
	// glutKeyboardFunc( onKeyDown ) ;
	
	// glutMainLoop();

	glutInit(&argc, argv);

	init(folder, instance, num_of_walks, num_of_intervals, instance_type);

	return 0;
}



