#ifndef ROBOT_MATH_H_
#define ROBOT_MATH_H_

#include <iostream>
#include <math.h>
#include <Eigen/Dense>

using namespace Eigen;

#define sign(x) ((x>0.0 || x==0.0) ? 1.0 : -1.0)

#define random(a,b) (1.0*rand()/RAND_MAX*(b-a)+a)

#endif