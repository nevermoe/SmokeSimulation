#ifndef SIMPLE_VECTOR_H
#define SIMPLE_VECTOR_H

#include "core.h"

class SimpleVector3d {
public:
	SimpleVector3d(GLdouble x = 0, GLdouble y = 0, GLdouble z = 0) {x_ = x; y_ = y; z_ = z;}
	GLdouble x_;
	GLdouble y_;
	GLdouble z_;
};

typedef SimpleVector3d Velocity;
typedef SimpleVector3d Position;
typedef SimpleVector3d Normal;

#endif
