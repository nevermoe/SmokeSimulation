#ifndef CORE_H
#define CORE_H

#define PI 3.14159
#define INF 99999999

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctype.h>
#include <png.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
using namespace std;

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>
#include <GL/glew.h>
#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include <GL/gl.h>
#include <GL/glu.h>


typedef Eigen::Triplet<GLdouble> T; 
#include "arcball.h"
#include "event_listener.h"
//#include "quaternion.h"
//#include "my_color.h"
//#include "my_traits.h"
//typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;
#include "camera.h"
//#include "my_math.h"
//#include "parameterization.h"
#include "object.h"
#include "controller.h"

////////////////////////////////////////////////////////////////////////////////
/*
This file just has a bunch of common stuff used by all objects. It mainly just
includes GL and some other standard headers.
*/
#endif
