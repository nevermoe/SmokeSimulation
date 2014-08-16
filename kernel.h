#ifndef KERNEL_H
#define KERNEL_H

#include "core.h"
class Kernel {
public:
	static GLdouble SmoothKernel(GLdouble r2, GLdouble h);
	static GLdouble SharpKernel(GLdouble r2, GLdouble h);
};


#endif
