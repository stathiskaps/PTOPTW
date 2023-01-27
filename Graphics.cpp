#include "Graphics.h"

// Graphics::Graphics(){
// 	Graphics::zoom = 1.0;
// }

float Graphics::zoom = 1.0f;
float Graphics::x_offset = 0.0f;
float Graphics::y_offset = 0.0f;
int Graphics::startX = 0;
int Graphics::startY = 0;
bool Graphics::isDragging = false;

void Graphics::onResize(int w, int h) {	
	int windowWidth = glutGet(GLUT_WINDOW_WIDTH);
	int windowHeight = glutGet(GLUT_WINDOW_HEIGHT);

	// int minSize = 500;
    // if (windowWidth < minSize) windowWidth = minSize;
    // if (windowHeight > minSize) windowHeight = minSize;
	
	// Calculate the aspect ratio of the window
	float aspectRatio = (float) windowWidth / (float) windowHeight;

	glViewport(0, 0, windowWidth, windowHeight);

	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	
	glOrtho(-400 * aspectRatio + x_offset, 400 * aspectRatio + x_offset, -400 + y_offset, 400 + y_offset, -1, 1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void Graphics::mouseButton(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON) {
        if (state == GLUT_DOWN) {
            Graphics::startX = x;
            Graphics::startY = y;
            Graphics::isDragging = true;
        } else if (state == GLUT_UP) {
            Graphics::isDragging = false;
        }
    } else if(button == 3) {
		Graphics::zoom *= 1.1;
	} else if(button == 4){
		Graphics::zoom /= 1.1;
	}
	glutPostRedisplay();
}

void Graphics::mouseMove(int x, int y) {
    if (isDragging) {
        int dx = x - startX;
        int dy = y - startY;
		x_offset += dx;
		y_offset += dy;
        // move the viewbox left/right or up/down based on dx and dy
    }
	glutPostRedisplay();
}

void Graphics::drawFilledCircle(float x, float y, float radius) {
    int num_segments = 36; // The number of segments used to draw the circle. More segments will result in a smoother circle
    float theta = 2 * 3.1415926 / float(num_segments); 
    float c = cosf(theta);
    float s = sinf(theta);
    float t;

    float x_ = radius;
    float y_ = 0;

    glBegin(GL_TRIANGLE_FAN);
    glVertex2f(x, y); // Center of circle
    for(int ii = 0; ii <= num_segments; ii++) {
        glVertex2f(x + x_, y + y_);
        t = x_;
        x_ = c * x_ - s * y_;
        y_ = s * t + c * y_;
    }
    glEnd();
}

void Graphics::arrowKeys(int key, int x, int y) {
  switch(key) {
    case GLUT_KEY_UP:
      y_offset += 0.1;
      break;
    case GLUT_KEY_DOWN:
      y_offset -= 0.1;
      break;
    case GLUT_KEY_LEFT:
      x_offset -= 0.1;
      break;
    case GLUT_KEY_RIGHT:
      x_offset += 0.1;
      break;
  }
  glutPostRedisplay();
}

void Graphics::myInit(){
	glShadeModel(GL_SMOOTH);
    glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
	glEnable(GL_POINT_SMOOTH);

}