/*
 *  smoke3D.h
 *  smoke3D
 *
 */

#include "object.h"
#include "types.h"
#include "utility.h"
#include "utility.h"
#include "solver.h"
#include "advect.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "types.h"

class Smoke3D: public Object {
public:
	Smoke3D();
	void SimulateStep();
	void enforce_boundary();
	void advection();
	void project();
	void Show();
private:
	FLOAT *** u[3];
	FLOAT *** b;
	FLOAT *** c;
	int frame;
};
