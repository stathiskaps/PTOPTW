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

#include <GL/glut.h>

#define WINDOW_WIDTH  400
#define WINDOW_HEIGHT 400

struct point{
	int x, y;
};

// Define a list of nodes to render
std::vector<std::pair<double, double>> nodes = {{0.5, 0.5}, {0.75, 0.75},
                                                 {1.0, 1.0},  {1.25, 1.25}};

/* Initialize OpenGL Graphics */
void initGL() {
   // Set "clearing" or background color
   glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black and opaque
}

void display() {
   glClear(GL_COLOR_BUFFER_BIT);   // Clear the color buffer with current clearing color
 
   // Define shapes enclosed within a pair of glBegin and glEnd
   glBegin(GL_QUADS);              // Each set of 4 vertices form a quad
      glColor3f(1.0f, 0.0f, 0.0f); // Red
      glVertex2f(-0.8f, 0.1f);     // Define vertices in counter-clockwise (CCW) order
      glVertex2f(-0.2f, 0.1f);     //  so that the normal (front-face) is facing you
      glVertex2f(-0.2f, 0.7f);
      glVertex2f(-0.8f, 0.7f);
 
      glColor3f(0.0f, 1.0f, 0.0f); // Green
      glVertex2f(-1.0f, -0.6f);
      glVertex2f(-0.1f, -0.6f);
      glVertex2f(-0.1f,  0.0f);
      glVertex2f(-0.7f,  0.0f);
 
      glColor3f(0.2f, 0.2f, 0.2f); // Dark Gray
      glVertex2f(-0.9f, -0.7f);
      glColor3f(1.0f, 1.0f, 1.0f); // White
      glVertex2f(-0.5f, -0.7f);
      glColor3f(0.2f, 0.2f, 0.2f); // Dark Gray
      glVertex2f(-0.5f, -0.3f);
      glColor3f(1.0f, 1.0f, 1.0f); // White
      glVertex2f(-0.9f, -0.3f);
   glEnd();
 
   glBegin(GL_TRIANGLES);          // Each set of 3 vertices form a triangle
      glColor3f(0.0f, 0.0f, 1.0f); // Blue
      glVertex2f(0.1f, -0.6f);
      glVertex2f(0.7f, -0.6f);
      glVertex2f(0.4f, -0.1f);
 
      glColor3f(1.0f, 0.0f, 0.0f); // Red
      glVertex2f(0.3f, -0.4f);
      glColor3f(0.0f, 1.0f, 0.0f); // Green
      glVertex2f(0.9f, -0.4f);
      glColor3f(0.0f, 0.0f, 1.0f); // Blue
      glVertex2f(0.6f, -0.9f);
   glEnd();
 
   glBegin(GL_POLYGON);            // These vertices form a closed polygon
      glColor3f(1.0f, 1.0f, 0.0f); // Yellow
      glVertex2f(0.4f, 0.2f);
      glVertex2f(0.6f, 0.2f);
      glVertex2f(0.7f, 0.4f);
      glVertex2f(0.6f, 0.6f);
      glVertex2f(0.4f, 0.6f);
      glVertex2f(0.3f, 0.4f);
   glEnd();
 
   glFlush();  // Render now
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

point* arrayCities = NULL;	//to store the cities positions
int nCities;								//number of cities
int* position = NULL;				//order to visit the cities
double** distances;					//a pre-computed array of distances between cities (to speedup)
double initialDistance;			//a value to print at the end of the running

//
// To display onto window using OpenGL commands
//
void displayOther()
{
	//
	// clear window to black
	//
	glClearColor( 0, 0 , 0 , 0 );
	glClear( GL_COLOR_BUFFER_BIT );
	glLoadIdentity();
	//orange color
	glColor3ub(200, 100, 15);
	
	//draw path
	glBegin(GL_LINE_LOOP);
	for (int k = 0; k < nCities; k++)
		glVertex2i(arrayCities[position[k]].x, arrayCities[position[k]].y);
	glEnd();

	//draw points, aka cities
	glPointSize(3);
	//light-gray color
	glColor3ub(200, 200, 200);
	for (int k = 0; k < nCities; k++)
	{
		glBegin(GL_POINTS);
			glVertex2i(arrayCities[k].x, arrayCities[k].y);
		glEnd();
	}
	glutSwapBuffers();
}

//return the total distance from initial city to the last one + last one to the initial one
double computeDistance()
{
	double distancePath = 0.0;
	for (int k = 0; k < nCities - 1; k++)
		distancePath += distances[position[k]][position[k + 1]];
	distancePath += distances[position[nCities - 1]][position[0]];
	return distancePath;
}

//swap 2 different random position in the path
// if distance is less, then swap, otherwise keep it
void Swap2Elements(int* e1, int *e2)
{
	do
	{
		*e1 = rand() % nCities;
		*e2 = rand() % nCities;
	} while (*e1 == *e2);
	int temp = position[*e1];
	position[*e1] = position[*e2];
	position[*e2] = temp;
}

// performs the iterations (are 25.000)
// a limit is added to avoid finish all iterations 9
void getShortestPath()
{
	int k = 0;
	int limit = 10000;
	while (k < 25000 && limit > 0)
	{
		int e1, e2;
		double d1 = computeDistance();
		Swap2Elements(&e1, &e2);
		double d2 = computeDistance();
		if (d2 > d1)
		{
			//revert the swaps
			int temp = position[e1];
			position[e1] = position[e2];
			position[e2] = temp;
			//std::cout << "No swap" << std::endl;
			limit--;
		}
		else
			printf("Distance iteration-%d : %lf\n", k,computeDistance());
		k++;
	}
}

//
// key function for ASCII charachters like ESC, a,b,c..,A,B,..Z
//
void onKeyDown(unsigned char key, int x, int y )
{
   // exit when ESC is pressed.
	if (key == 27)
		exit(0);
	 else if (key == 'a')
	 {
		 getShortestPath();
		 printf("Initial distance: %lf; Shortest path: %lf\n", initialDistance, computeDistance());
	 }
	 else if (key == 's')
	 {
			int e1, e2;
			double d1 = computeDistance();
			Swap2Elements(&e1, &e2);
			double d2 = computeDistance();
			if (d2 > d1)
			{
				//revert the swaps
				int temp = position[e1];
				position[e1] = position[e2];
				position[e2] = temp;
			}
			printf("Distance: %lf \n", computeDistance());

	 }
    // to refresh the window it calls display() function
   glutPostRedisplay() ;
}

//
// This function is called when the window size changes.
// w : is the new width of the window in pixels.
// h : is the new height of the window in pixels.
//
void onResize( int w, int h )
{
	glViewport( 0,0,w,h) ;
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	glOrtho( -w/2, w/2, -h/2, h/2, -1, 1);
	glMatrixMode( GL_MODELVIEW);
	glLoadIdentity();
}

// compute the euclidean distance between two points
double distanceFinder(point a, point b)
{
	return sqrt(pow(a.x - b.x, 2.0) + pow(a.y - b.y, 2.0));
}

//initialize the opengl setup and variables
void Init(int size) {
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
	glEnable(GL_POINT_SMOOTH);
	
	nCities = size;			// set the number of cities given by the user
	srand(time(NULL));	//to start each execution as different

	//dynamically create the arrays
	//arrayCities = new point[size];
	arrayCities = (point*)malloc(size * sizeof(point));
	position = (int*)malloc(size * sizeof(int));
	distances = (double**)malloc(size * sizeof(double*));

	// place the cities on the screen
	for (int k = 0; k < size; k++)
	{
		distances[k] = (double*)malloc(size * sizeof(double));
		int ri = rand() % (WINDOW_WIDTH + 1);
		ri -= WINDOW_WIDTH / 2;
		arrayCities[k].x = ri;
		ri = rand() % (WINDOW_HEIGHT + 1);
		ri -= WINDOW_HEIGHT / 2;
		arrayCities[k].y = ri;
		position[k] = k;
	}

	//computes the distance from all cities to all cities
	for (int y = 0; y < size; y++)
		for (int x = 0; x < size; x++)
			distances[x][y] = distanceFinder(arrayCities[x], arrayCities[y]);	//from city y to each x

	initialDistance = computeDistance();
	printf("Initial distance: %lf\n", initialDistance);
}

// Define a simple struct to represent a point in 2D space
struct Node
{
    double x, y;
};

void displayTriangle() {
  // Clear the screen to black
  glClearColor(0.0, 0.0, 0.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT);

  // Set the drawing color to white
  glColor3f(1.0, 1.0, 1.0);

  // Set the triangle vertices
  glBegin(GL_TRIANGLES);
  glVertex2f(-0.5, -0.5);
  glVertex2f(0.5, -0.5);
  glVertex2f(0.0, 0.5);
  glEnd();

  // Flush the OpenGL buffers to the screen
  glFlush();
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
    for (const Point& p : points){
        out << "<circle cx=\"" << p.pos.lat << "\" cy=\"" << p.pos.lon  << "\" r=\"" << radius << "\" />" << std::endl;
    }

	for(const Point& p : points){
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
	int n;
	printf("Please insert the number of cities [0,10000]: ");
	scanf("%d", &n);
	printf("\n");
	glutInit(&argc, argv );
	glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE );
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	glutCreateWindow( "TSP Opengl" ) ;

	glutDisplayFunc(displayOther ) ;
	glutIdleFunc(displayOther);
  	glutReshapeFunc( onResize );
	//
	// keyboard registration
	//
	glutKeyboardFunc( onKeyDown ) ;
	Init(n);
	
	glutMainLoop();

	// glutInit(&argc, argv);          // Initialize GLUT
	// glutInitWindowSize(640, 480);   // Set the window's initial width & height - non-square
	// glutInitWindowPosition(50, 50); // Position the window's initial top-left corner
	// glutCreateWindow("Viewport Transform");  // Create window with the given title
	// glutDisplayFunc(display);       // Register callback handler for window re-paint event
	// glutReshapeFunc(reshape);       // Register callback handler for window re-size event
	// initGL();                       // Our own OpenGL initialization

	// glEnable(GL_BLEND);
	// glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	// glEnable(GL_LINE_SMOOTH);
	// glEnable(GL_POLYGON_SMOOTH);
	// glEnable(GL_POINT_SMOOTH);

	// // Enter the GLUT main loop
	// glutMainLoop();

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



