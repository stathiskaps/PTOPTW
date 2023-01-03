#include "Graphics.h"

void onResize(int w, int h) {	
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
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void myInit(){
	glShadeModel(GL_SMOOTH);
    glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
	glEnable(GL_POINT_SMOOTH);
}