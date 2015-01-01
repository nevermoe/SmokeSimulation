/* Adapted from: Johannes Schmid, 2006 
 * https://graphics.ethz.ch/teaching/former/imagesynthesis_06/miniprojects/p3 */

#ifndef _RENDERER_H
#define _RENDERER_H

#include "core.h"

class Fluid;	// forward definition

#define ALMOST_EQUAL(a, b) ((fabs(a-b)<0.00001f)?true:false)

class Renderer
{
private:
	int _sx, _sy;			// screen width and height
	int _ax, _ay;			// mouse movement anchor

	float _quat[4];			// view rotation quaternion
	float _lquat[4];		// light orientation quaternion (bit of a hack...)
	float _dist;			// viewer distance to origin

	FILE* _fp;				// data file
	int _N;					// data resolution
	int _nframes;			// number of frames in data file
	int _cur_frame;

		// OpenGL variables
	unsigned char* _texture_data;
	unsigned int _txt[3];				// texture handles
	unsigned int _prog[2];				// program handles
	unsigned int _font_base;			// first display list for font
	double _persp_m[16], _ortho_m[16];	// projection matrices

	float _light_dir[3];
	int _ray_templ[4096][3];

		// draw the outline of the cube
	void draw_cube(void);

		// draw the slices. m must be the current rotation matrix.
		// if frame==true, the outline of the slices will be drawn as well
	void draw_slices(GLdouble m[][4], bool frame);
		// intersect a plane with the cube, helper function for draw_slices()
	std::vector<Eigen::Vector3f> intersect_edges(float a, float b, float c, float d);

	void gen_ray_templ(int edgelen);
	void cast_light(int edgelen, float* dens, unsigned char* intensity);
	inline void light_ray(int x, int y, int z, int n, float decay, float* dens, unsigned char* intensity);


	void init_GL(void);

public:
	Renderer();
	~Renderer();

	bool _draw_slice_outline;
	char* _dispstring;

	void FillTexture(Fluid* fluid);		// generate data from simulation
	void Render(void);						// draw the volume
};


#endif
