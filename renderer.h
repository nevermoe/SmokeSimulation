/* Adapted from: Johannes Schmid, 2006 
 * https://graphics.ethz.ch/teaching/former/imagesynthesis_06/miniprojects/p3 */

#ifndef _RENDERER_H
#define _RENDERER_H

#include "core.h"
#include "fluid.h"

class Fluid;	// forward definition

#ifndef ALMOST_EQUAL
#define ALMOST_EQUAL(a, b) ((fabs(a-b)<0.00001f)?true:false)
#endif
#define SLICE_NUM			64.0f


class Renderer
{
private:

	// texture data
	unsigned char* _textureData;
	// texture handle
	unsigned int _hTexture;				
	//lighting infomations
	Eigen::Vector3f _lightDir;
	Eigen::Vector3f _lightPos;
	int _rayTemplate[4096][3];

	GLfloat _cubeVertices[8][3];
	GLfloat _cubeEdges[12][2][3];

	// draw the outline of the cube
	void DrawCube(void);

	// draw the slices. mvMatrix must be the MODELVIEW_MATRIX
	void DrawSlices(GLdouble mvMatrix[16]);

	// intersect a plane with the cube, helper function for DrawSlices()
	// plane equation is Ax + By + Cz + D = 0
	std::vector<Eigen::Vector3f> IntersectEdges(float A, float B, float C, float D);

	void GenerateRayTemplate(int edgeLen);
	void CastLight(int edgelen, const float* dens, unsigned char* intensity);
	inline void LightRay(int x, int y, int z, int n, float decay, 
			const float* dens, unsigned char* intensity);


	void InitGL();

	// if _isDrawSliceOutline==true, the outline of the slices will be drawn as well
	bool _isDrawSliceOutline;
public:
	Renderer();
	~Renderer();
	void SetLightPostion(Eigen::Vector3f &pos);


	void FillTexture(Fluid* fluid);		// generate texture from smoke density 
	void Render(void);					// draw the volume
};


#endif
