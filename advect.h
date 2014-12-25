/*
 *  advect.h
 *  smoke3D
 *
 */

#include "types.h"

namespace advect {
	void advect( FLOAT ****u, FLOAT ***c, int n, FLOAT dt );
}