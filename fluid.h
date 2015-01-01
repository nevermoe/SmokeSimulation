/* Adapted from: Johannes Schmid, 2006 
 * https://graphics.ethz.ch/teaching/former/imagesynthesis_06/miniprojects/p3 */

#ifndef _FLUID_H
#define _FLUID_H

#include "core.h"
#include "object.h"
#include "renderer.h"

#define DT  0.1
#define RES 40			// box resolution
#define N ((RES)-2)			// valid simulation area
#define SIZE ((RES)*(RES)*(RES))
#define _I(x,y,z) (((x)*(RES)*(RES))+((y)*(RES))+(z))

#define FOR_ALL_CELL for (int i=1; i<=(N); i++) {\
	for (int j=1; j<=(N); j++) {\
		for (int k=1; k<=(N); k++) {

#define END_FOR }}}

#define ALMOST_EQUAL(a, b) ((fabs(a-b)<0.00001f)?true:false)

class Fluid: public Object
{
protected:
	float _buffers[10][SIZE];
public:
	float *_density, *_densityTmp;			// density
	float *_velX, *_velXTmp;			// velocity in x direction
	float *_velY, *_velYTmp;			// velocity in y direction
	float *_velZ, *_velZTmp;			// velocity in z direction
	float _dt;

	Renderer * _renderer;				//register a renderer

protected:
	// simulation methods
	// beware: i changed stam's implementation from a destiation, source ordering
	// of the parameters to source, destiation
	void AddSource(float* src, float *dst);
	void AddBuoyancy();
	void EnforceBoundary(int b, float* quantity);
	void Diffuse(int b, float* velTmp, float* vel, float diff);
	void Advect(int b, float* quantityTmp, float* quantity, float* velX, float* velY, float* velZ);
	void Project(void);
	void VorticityConfinement();

	void VelocityStep();
	void DensityStep();

	// utility methods
	void ClearBuffer(float* buf);
	void ClearSources(void);

	void _GenerateSmoke();
public:
	float sd[SIZE], su[SIZE], sv[SIZE], sw[SIZE], sT[SIZE];	// sources for density and velocities
	float diffusion, viscosity, buoyancy, vc_eps;

	Fluid(void);
	~Fluid(void);

	virtual void SimulateStep();
	virtual void Show();
};

#endif

