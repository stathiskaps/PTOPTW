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
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

enum InstanceType { fixed, custom };


std::vector<Point> points;

/* Initialize OpenGL Graphics */
void initGL() {
   // Set "clearing" or background color
   glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black and opaque
}

/* Handler for window re-size event. Called back when the window first appears and
   whenever the window is re-sized with its new width and height */
void reshape(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
   // Compute aspect ratio of the new window
   if (height == 0) height = 1;                // To prevent divide by 0
   GLfloat aspect = (GLfloat)width / (GLfloat)height;
 
   // Set the viewport to cover the new window
   glViewport(0, 0, width, height);
 
   // Set the aspect ratio of the clipping area to match the viewport
   glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
   glLoadIdentity();             // Reset the projection matrix
   if (width >= height) {
     // aspect >= 1, set the height from -1 to 1, with larger width
      gluOrtho2D(-1.0 * aspect, 1.0 * aspect, -1.0, 1.0);
   } else {
      // aspect < 1, set the width to -1 to 1, with larger height
     gluOrtho2D(-1.0, 1.0, -1.0 / aspect, 1.0 / aspect);
   }
}

void createGraph(){
	struct Vertex { int foo; };
    struct Edge { std::string blah; };

    using namespace boost;
    using graph_t  = adjacency_list<listS, vecS, directedS, Vertex, Edge >;
    using vertex_t = graph_traits<graph_t>::vertex_descriptor;
    using edge_t   = graph_traits<graph_t>::edge_descriptor;

    //Instantiate a graph
    graph_t g;

    // Create two vertices in that graph
    vertex_t u = boost::add_vertex(Vertex{123}, g);
    vertex_t v = boost::add_vertex(Vertex{456}, g);

    // Create an edge conecting those two vertices
    boost::add_edge(u, v, Edge{"Hello"}, g);

    boost::write_graphviz(std::cout, g, [&] (auto& out, auto v) {
       out << "[label=\"" << g[v].foo << "\"]";
      },
      [&] (auto& out, auto e) {
       out << "[label=\"" << g[e].blah << "\"]";
    });
    std::cout << std::flush;
}

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

	std::string filepath = "./instances/"+folder+"/"+filename+".txt";

	std::ifstream infile(filepath);

	std::string line;

	std::getline(infile, line);
	std::getline(infile, line);
	std::vector<std::string> poi_data;

	int pointId = 0;
	std::vector<double> travel_times;

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
			std::string from = poi_data[2], to = poi_data[3];
			double duration = std::stod(poi_data[4]);
			travel_times.push_back(duration);
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
				double open_time, close_time, visit_dur, profit;
				if(instance_type == custom){
					visit_dur = std::stoi(poi_data[3]);
					profit = std::stoi(poi_data[4]);
					open_time = std::stoi(poi_data[5]);
					close_time = std::stoi(poi_data[6]);
				} else {
					visit_dur = std::stoi(poi_data[4]);
					profit = std::stoi(poi_data[5]);
					open_time = std::stoi(poi_data[8 + std::stoi(poi_data[7])]);
					close_time = std::stoi(poi_data[9 + std::stoi(poi_data[7])]);
				}

				touristAttractions.push_back(new Sight(
					id::generate(),
					p,
					visit_dur,
					profit,
					open_time,
					close_time

				));
			} else {
				Point p = Point(pointId++, std::stod(poi_data[1]), std::stod(poi_data[2]));
				points.push_back(p);
				double open_time, close_time;
				if(instance_type == custom){
					open_time = std::stoi(poi_data[5]);
					close_time = std::stoi(poi_data[6]);
				} else {
					open_time = std::stoi(poi_data[7 + std::stoi(poi_data[6])]);
					close_time = std::stoi(poi_data[8 + std::stoi(poi_data[6])]);
				}

				touristAttractions.push_back(new Sight(
					id::generate(),
					p,
					std::stoi(poi_data[3]),
					std::stoi(poi_data[4]),
					open_time,
					close_time
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

	OP op = OP(touristAttractions, points, start_depot, end_depot, numRoutes, end_depot->timeWindow.openTime, end_depot->timeWindow.closeTime);

	ILS_TOPTW ilstoptw = ILS_TOPTW(numIntervals);
	ilstoptw.SolveNew(op);

	/*OPTW optw(touristAttractions, ttMatrix, depot, OPEN_DAY_TIME, CLOSE_DAY_TIME);
	Solution sol = optw.solve();
	bool valid = optw.validate();
	std::string msg = valid ? "yes" : "no";
	std::cout << "valid solution? " << msg << std::endl;
	sol.print();*/
#endif

	for (auto p : touristAttractions) {
		delete p;
	}

	touristAttractions.clear();

}





int main(int argc, char** argv)
{

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
	// Init(n);
	
	// glutMainLoop();

	// glutInit(&argc, argv);

	init(folder, instance, num_of_walks, num_of_intervals, instance_type);

	return 0;
}



