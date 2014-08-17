#include "core.h"
#include "kernel.h"

GLdouble Kernel::SmoothKernel(GLdouble r2, GLdouble h)
{
	return std::max(1.0-r2/(h*h), 0.0);
}

GLdouble Kernel::SharpKernel(GLdouble r2, GLdouble h)
{
	return std::max(h*h/std::max(r2,1.0e-5)-1.0, 0.0);
}
