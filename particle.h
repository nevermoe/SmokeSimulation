#ifndef PARTICLE_H
#define PARTICLE_H

#include "core.h"

class Particle {
public: 
	GLdouble density_;	
	GLdouble mass_;
	Position position_;
	Velocity velocity_;
};


typedef std::vector<Particle* > ParticleList;

#endif
