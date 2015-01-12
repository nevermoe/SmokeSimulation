#ifndef _RENDERER_H
#define _RENDERER_H

#include "core.h"

#ifndef ALMOST_EQUAL
#define ALMOST_EQUAL(a, b) ((fabs(a-b)<0.00001f)?true:false)
#endif

#ifndef _I
#define _I(x,y,z) (((x)*(_RES)*(_RES))+((y)*(_RES))+(z))	//FIXME
#endif

#ifndef FOR_ALL_CELL
#define FOR_ALL_CELL for (int i=1; i<=(_N); i++) {\
	for (int j=1; j<=(_N); j++) {\
		for (int k=1; k<=(_N); k++) {
#define END_FOR }}}
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
	float *_volumeData;

	GLfloat _cubeVertices[8][3];
	GLfloat _cubeEdges[12][2][3];

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
	bool _isRendering;

	int _SIZE;		//size of volume data
	int _N;			//
	int _RES;		//
public:
	Renderer(float* volumeData, int RES);
	~Renderer();
	void SetLightPostion(Eigen::Vector3f &pos);
	void SetRendering(bool isRendering);
	void SetSliceOutline(bool isDrawSliceOutline);


	void FillTexture();		// generate texture from smoke density 
	void Render();					// draw the volume
	// draw the outline of the cube
	void DrawCube();
	void DrawLight();
	void DrawVolumeData();

};


#endif
