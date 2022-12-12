#include <GL/glut.h>
#include <thread>

#define WINDOW_WIDTH  800
#define WINDOW_HEIGHT 800

void resize(int w, int h) {	
	int windowWidth = glutGet(GLUT_WINDOW_WIDTH);
	int windowHeight = glutGet(GLUT_WINDOW_HEIGHT);

	int minSize = 500;
    if (windowWidth < minSize) windowWidth = minSize;
    if (windowHeight > minSize) windowHeight = minSize;
	
	// Calculate the aspect ratio of the window
	float aspectRatio = (float) windowWidth / (float) windowHeight;

	glViewport(0, 0, windowWidth, windowHeight);


	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();


	glOrtho(-100 * aspectRatio, 100 * aspectRatio, -100, 100, -1, 1);
	glMatrixMode( GL_MODELVIEW);
	glLoadIdentity();
}

// void display() {
// 	// Clear the screen to black
// 	glClearColor(0.0, 0.0, 0.0 , 1.0);
// 	glClear(GL_COLOR_BUFFER_BIT);
// 	glLoadIdentity();

// 	const float pointSize = 20;

// 	// //orange color
// 	// glColor3ub(200, 200, 200);
// 	// glLineWidth(2.0);

// 	// glBegin(GL_LINES);
// 	// for(int i = 0; i != points.size() - 2; ++i){
// 	// 	glVertex2d(points[i].pos.lat, points[i].pos.lon);
// 	// 	glVertex2d(points[i+1].pos.lat, points[i+1].pos.lon);
// 	// }
// 	// glEnd();
	

// 	// Use glPointSize() to set the size of the points
// 	glPointSize(pointSize);
// 	// Set the drawing color to orange
// 	glColor3ub(200, 103, 51);

// 	// Loop through the list of points and render each one
// 	for (const auto& point : points) {
// 		// Extract the coordinates and label from the point_label pair
// 		double lat = point.pos.lat;
// 		double lon = point.pos.lon;
// 		std::string id = std::to_string(point.id);
// 		int length = id.size();

// 		// Use glBegin() and glEnd() to render the point
// 		glBegin(GL_POINTS);
// 		glVertex2d(lat, lon);
// 		glEnd();

		
// 		glColor3ub(255, 255, 255); //white
// 		glRasterPos2d(lat-length, lon-1);
// 		// Use a loop to draw each character of the text string
// 		for (const auto& c : id) {
// 			glutBitmapCharacter(GLUT_BITMAP_8_BY_13, c);
// 		}
// 		glColor3ub(200, 103, 51); //orange
// 	}

// 	// Flush the OpenGL buffers to the screen
// 	glFlush();
// }

void myInit(){
    glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
	glEnable(GL_POINT_SMOOTH);
}