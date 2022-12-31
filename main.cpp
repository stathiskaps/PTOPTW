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
#include "Definitions.h"
#include "ILS.h"
#include "OP.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>


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

void init(std::string author, std::string filename, int numRoutes, int numIntervals) {

	std::vector<TA*> touristAttractions; //TODO:delete pointers

	std::string filepath = "./instances/"+author+"/"+filename+".txt";

	std::ifstream infile(filepath);

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
	if(argc < 5) {
		std::cout << "Please run the program with the following arguments: "<< std::endl <<
		 "./AMTOPTW <author> <filename> <number of walks> <number of sub-intervals> " << std::endl <<
		 "e.g. ./AMTOPTW Cordeau pr11 4 4" << std::endl;
		return 0;
	}
	// int n;
	// printf("Please insert the number of cities [0,10000]: ");
	// scanf("%d", &n);
	// printf("\n");
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

	std::string author = argv[1];
	std::string filename = argv[2];
	int num_of_walks = std::stoi(argv[3]);
	int num_of_intervals = std::stoi(argv[4]);

	init(author, filename, num_of_walks, num_of_intervals);

	return 0;
}



