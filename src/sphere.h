#ifndef SPHERE_H
#define SPHERE_H

#include "sphere_data_utility.h"
#include "search.h"
#include "sphere_operation.h"
#include"sphere_RMSUtils.h"

#include "sphere_lp.h"
#include"point.h"

// The complete Sphere algorithm
point_set_t* sphereWSImpLP(point_set_t* point_set, int k);
point_set_t* sphereWSImpLP1(point_set_t* point_set, int k);
void runSphere(vector<Point> dataP, int r, int k, vector<Point> &result, vector<Point> curSky, double &time);

#endif