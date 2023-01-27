#pragma once
#include <GL/glut.h>
#include <thread>
#include <vector>
#include <iostream>
#include <math.h>

#ifndef GRAPHICS_H
#define GRAPHICS_H

#define WINDOW_WIDTH  800
#define WINDOW_HEIGHT 800

class Graphics{
     struct Color {
        uint8_t r;
        uint8_t g;
        uint8_t b;
    };
public:

    static float zoom; // zoom level
    static float x_offset; // x offset for navigation
    static float y_offset; // y offset for navigation
    static int startX, startY;
    static bool isDragging;


    static void onResize(int w, int h);
    static void mouseButton(int button, int dir, int x, int y);
    static void myInit();
    static void arrowKeys(int, int, int);
    static void mouseMove(int, int);
    static void drawFilledCircle(float cx, float cy, float r);
    static constexpr Color colors[] = {
        {236, 110, 46},
        {46, 158, 236},
        {112, 236, 46},
        {158, 46, 236}
    };
    static constexpr int color_count = sizeof(colors) / sizeof(colors[0]);
};

#endif