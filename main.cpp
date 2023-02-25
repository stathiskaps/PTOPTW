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
#include <filesystem>
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

std::pair<int, double> init(std::string filepath, std::string filename, int numRoutes, int numIntervals, InstanceType instance_type, ILS::Configuration conf) {

	std::vector<TA> touristAttractions;
	std::vector<Point> points;
	
	std::ifstream infile(filepath);
	std::string line;

	std::getline(infile, line);
	std::getline(infile, line);
	std::vector<std::string> poi_data;

	int pointId = 0;
	std::vector<RouteDuration> durations;

	if(instance_type == custom){

		std::ifstream i("../topology/topology.json");
		json j;
		i >> j;
		
		json preferences = j["preferences"];
		json nodes = j["nodes"];
		json routes = j["routes"];

		for(auto& node : nodes){
			Point p = Point(pointId++, node["lat"], node["lon"]);
			points.push_back(p);
			json timeWindow = node["time_window"];
			touristAttractions.push_back(Sight(
				id::generate(),
				p,
				node["visit_time"],
				preferences[node["category"]],
				timeWindow["start_time"],
				timeWindow["end_time"]
			));
		}

		for(auto& route : routes){
			int id = route["id"], from = route["from"], to = route["to"];
			double duration = route["duration"];
			durations.push_back(RouteDuration{id, from , to, duration});
		}
		
	} else {

		while (std::getline(infile, line))
		{
			poi_data = split(line);
			if(poi_data.size() == 0) break;

			switch (line.front()) {
			case 'R': {
				Point p = Point(pointId++, std::stod(poi_data[2]), std::stod(poi_data[3]));
				Point edge1 = Point(pointId++, std::stod(poi_data[4]), std::stod(poi_data[5]));
				Point edge2 = Point(pointId++, std::stod(poi_data[6]), std::stod(poi_data[7]));

				points.push_back(p);
				points.push_back(edge1);
				points.push_back(edge2);

				touristAttractions.push_back(Route(
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
					touristAttractions.push_back(Sight(
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
					touristAttractions.push_back(Sight(
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
	}

	TA start_depot(touristAttractions[0]);
	TA end_depot(touristAttractions[0]);
	start_depot.id = START_DEPOT_ID;
	end_depot.id = END_DEPOT_ID;
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
		op = OP(touristAttractions, points, start_depot, end_depot, numRoutes, end_depot.timeWindow.openTime, end_depot.timeWindow.closeTime, travel_times);
	} else {
		op = OP(touristAttractions, points, start_depot, end_depot, numRoutes, end_depot.timeWindow.openTime, end_depot.timeWindow.closeTime);
	}
	

	ILS ils = ILS(numIntervals, filename, conf);

	return ils.Solve(op);

}


int main(int argc, char** argv) {

	InstanceType instance_type = fixed; 

	std::string folder;
	std::string instance;
	double execution_time_limit = 0;
	int num_of_walks, num_of_intervals, total_score = 0;
	ILS::Configuration conf;
	bool run_all = false, run_all_cases = false;

	int c;
	int option_index = 0;
	bool folder_option_provided = false, instance_option_provided = false, 
		walks_option_provided = false, subs_option_provided = false;
	static struct option long_options[] = {
		{"all", no_argument, 0, 'a'},
		{"all-cases", no_argument, 0, 'r'},
		{"help", no_argument, 0, 'h'},
		{"write", no_argument, 0, 'w'},
		{"json", no_argument, 0, 'j'},
		{"python", no_argument, 0, 'p'},
		{"custom", no_argument, 0, 'c'},
		{"graphics", no_argument, 0, 'g'},
		{"folder", required_argument, 0, 'f'},
		{"instance", required_argument, 0, 'i'},
		{"walks", required_argument, 0, 'm'},
		{"solutions", required_argument, 0, 's'},
		{"time", required_argument, 0, 't'},
		{0, 0, 0, 0}
	};

	while ((c = getopt_long(argc, argv, "ahrwcjgf:i:m:s:t:", long_options, &option_index)) != -1) {
		switch (c) {
		case 'h': {
			std::cout << "Please run the program with the following options: "<< std::endl <<
			"./AMTOPTW -f <folder> -i <instance> -m <number of walks> -s <number of sub-intervals> " << std::endl <<
			"or" << std::endl <<
			"./AMTOPTW -folder=\"<folder>\" --instance=\"<instance>\" --walks=<num_of_walks> --intervals=<num_of_intervals>" << std::endl << 
			"e.g. ./AMTOPTW -f Cordeau -i pr11 -m 4 -s 4" << std::endl << 
			"Also, you can add the option -c without any arguments if you want to use custom travel times from a custom made topology" << std::endl;

			return 0;
		}
		case 'a':{
			run_all = true;
			std::cout << "Will run all instances" << std::endl;
			break;
		}
		case 'w': {
			conf.write_results = true;
			break;
		}
		case 'j': {
			conf.write_solution = true;
			break;
		}
		case 'r': {
			run_all_cases = true;
			break;
		}
		case 'g':{
			conf.graphics = true;
			break;
		}
		case 't':{
			if(optarg) {
				conf.execution_time_limit = std::atof(optarg);
				conf.time_limited_execution = true;
				std::cout << "Will run ILS for " << conf.execution_time_limit << " seconds" << std::endl;
			} else {
				std::cerr << "Error: argument for option -t is required" << std::endl;
				exit(EXIT_FAILURE);
			}
			break;
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
		case 'm': {
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

	if (!folder_option_provided && !run_all) {
		std::cerr << "Error: Option -f is required" << std::endl;
		exit(EXIT_FAILURE);
	}

	if (!instance_option_provided && !run_all) {
		std::cerr << "Error: Option -i is required" << std::endl;
		exit(EXIT_FAILURE);
	}

	if (!walks_option_provided && !run_all) {
		std::cerr << "Error: Option -m is required" << std::endl;
		exit(EXIT_FAILURE);
	}

	if (!subs_option_provided && !run_all && !run_all_cases) {
		std::cerr << "Error: Option -s is required" << std::endl;
		exit(EXIT_FAILURE);
	}

	if(subs_option_provided && run_all_cases){
		std::cerr << "Error: Options -r and -s can't be used together" << std::endl;
		exit(EXIT_FAILURE);
	}

	if(conf.graphics){
		glutInit(&argc, argv);
	}

	if(run_all && run_all_cases){
		std::cerr << "Error: Options -r and -a cannot be used together, please pick one" << std::endl;
		exit(EXIT_FAILURE);
	}

	if(!run_all && !run_all_cases){
		std::string filepath = "./instances/"+folder+"/"+instance+".txt";
		init(filepath, instance, num_of_walks, num_of_intervals, instance_type, conf);
	} else {
		if(run_all_cases){
			std::string filepath = "./instances/"+folder+"/"+instance+".txt";
			std::ofstream outfile("./output/"+instance+"_"+std::to_string(num_of_walks)+".csv", std::ofstream::out);
			outfile << "s,score,time" << std::endl;

			for(size_t k = 1; k < 5; ++k){
				auto [score, time] = init(filepath, instance, num_of_walks, k, instance_type, conf);
				outfile << k << "," << score << "," << time << std::endl;
			}
			outfile.close();
		} else if(run_all){
			std::string directories[2] = {"./instances/Solomon/", "./instances/Cordeau/"};
			for(size_t i = 0; i < 2; ++i){
				const std::filesystem::path directory_path(directories[i]);
				std::vector<std::filesystem::path> entries;
				for (const auto &entry : std::filesystem::directory_iterator(directory_path)) {
					if (entry.is_regular_file()) {
						entries.push_back(entry);
					}
				}
				std::sort(entries.begin(), entries.end(), [](const auto &a, const auto &b) {
					return a.filename().string() < b.filename().string();
				});

				for (const auto &entry : entries){
					auto filename = entry.filename();
					for(size_t j = 1; j < 5; ++j){
						for(size_t k = 1; k < 5; ++k){
							auto [score, time] = init(entry.string(), filename.replace_extension().string(), j, k, instance_type, conf);
							total_score += score;
						}
					}
				}
			}
			std::cout << "Total score: " << total_score << std::endl;
		}
	}
	

	return 0;
}



