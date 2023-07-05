//
// Created by PC on 8/10/2022.
//

#ifndef DEF_1705077_RTX_METHODS_H
#define DEF_1705077_RTX_METHODS_H
#include "1705077_classes.h"
#include <bits/stdc++.h>
#include "bitmap_image.hpp"

extern Camera * cam;
extern Scene * scene;
extern Screen * screen;
extern Floor * sceneFloor;
extern SpotLight * spl;
extern double windowWidth, windowHeight;
extern Vec3 captureCamPos;
extern std::vector<Vec3> capturePixelEndpoints;
extern std::vector<Vec3> intersectionPoints;
extern std::vector<Vec3> viewPort;
extern std::vector<std::pair<Vec3, Object *>> testRayIntersections;

int fileNo = 1;
/// Capture screen via raytrace and write to a bitmap
void capture();

#endif //DEF_1705077_RTX_METHODS_H
