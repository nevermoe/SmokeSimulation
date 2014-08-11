#include "core.h"

class Cell {
public:
	ParticleList particles_;
	GLdouble press_;
	int type_;
};

class Grid {
public:
	Cell*** cells_;
	int nX_, nY_, nZ_;

	Velocity*** velX_;	//for MAC grid, vel is stored at the edge of cells
	Velocity*** velY_;	//so it cannot be stored in Cell struct;
	Velocity*** velZ_;

	GLdouble cellSize_;	//side length of one cell
};
