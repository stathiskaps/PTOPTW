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

// Define a simple struct to represent a point in 2D space
struct Node
{
    double x, y;
};

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

struct Bounds{
	double minLat, minLon, maxLat, maxLon;
};


void drawPoints(std::vector<Point> points){

	Bounds bounds;
	for(auto p: points){
		if(p.pos.lat < bounds.minLat){
		bounds.minLat = p.pos.lat;
		}
		if(p.pos.lat > bounds.maxLat){
			bounds.maxLat = p.pos.lat;
		}
		if(p.pos.lon < bounds.minLon){
			bounds.minLon = p.pos.lon;
		}
		if(p.pos.lon > bounds.maxLon){
			bounds.maxLon = p.pos.lon;
		}
	}

	// Create a vector of nodes to represent the graph

	const double radius = 2;

    std::vector<Node> nodes = {
        {1, 2},
        {3, 4},
        {5, 6},
        {7, 8}
    };

	std::cout << "Points.size = " << points.size() << std::endl;

    // Create a vector of pairs of indices to represent the routes between nodes
    std::vector<std::pair<int, int>> routes = {
        {0, 1},
        {1, 2},
        {2, 3}
    };

    // Open an output stream to write the SVG file
    std::ofstream out("graph.svg");



    // Write the SVG file header
    out << "<svg viewBox=\"" << bounds.minLat - radius << " " << bounds.minLon - radius << " " << bounds.maxLat - bounds.minLat + radius*2 << " " << bounds.maxLon - bounds.minLon + radius*2 << "\" xmlns=\"http://www.w3.org/2000/svg\">" << std::endl;


	// Write the routes as lines in the SVG file
    for (const auto& [n1, n2] : routes)
    {
        const Point& point1 = points[n1];
        const Point& point2 = points[n2];
        out << "<line x1=\"" << point1.pos.lat << "\" y1=\"" << point1.pos.lon << "\" x2=\"" << point2.pos.lat << "\" y2=\"" << point2.pos.lon << "\" style=\"stroke:rgb(66,66,66);stroke-width:0.5\" />" << std::endl;
    }

    // Write the nodes as circles in the SVG file
    for (const Point& p : points)
    {
        out << "<circle cx=\"" << p.pos.lat << "\" cy=\"" << p.pos.lon  << "\" r=\"" << radius << "\" />" << std::endl;
		out << "<text x=\""<< p.pos.lat << "\" y=\""<< p.pos.lon << "\" text-anchor=\"middle\" font-size=\"2px\" fill=\"white\" alignment-baseline=\"middle\">" << p.id  << "</text>" << std::endl;
    }



    // Write the SVG file footer
    out << "</svg>" << std::endl;

    // Close the output stream
    out.close();
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

void init(std::string filename, int numRoutes, int numIntervals) {

	std::vector<TA*> touristAttractions; //TODO:delete pointers
	std::vector<Point> points;
	std::string path = "./instances/Cordeau/";

	Bounds bounds;

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
	drawPoints(points);

#if 1
	std::cout << std::endl;
	TA* depot = touristAttractions[0];
	TA* start_depot = depot->clone();
	start_depot->id = START_DEPOT_ID;
	TA* end_depot = depot->clone();
	end_depot->id = END_DEPOT_ID;
	touristAttractions.erase(touristAttractions.begin());

	OP op = OP(touristAttractions, points, start_depot, end_depot, numRoutes, end_depot->timeWindow.openTime, end_depot->timeWindow.closeTime);

	ILS_TOPTW ilstoptw = ILS_TOPTW(numRoutes, numIntervals);
	ilstoptw.SolveNew(op);

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



