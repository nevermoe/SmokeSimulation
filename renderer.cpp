#include "renderer.h"



Renderer::Renderer()
{
	_textureData = NULL;
	_isDrawSliceOutline = false;

	// cube vertices
	GLfloat cv[][3] = {
		{1.0f, 1.0f, 1.0f}, {-1.0f, 1.0f, 1.0f}, {-1.0f, -1.0f, 1.0f}, {1.0f, -1.0f, 1.0f},
		{1.0f, 1.0f, -1.0f}, {-1.0f, 1.0f, -1.0f}, {-1.0f, -1.0f, -1.0f}, {1.0f, -1.0f, -1.0f}
	};

	// cube edges have the form edges[n][0][xyz] + t*edges[n][1][xyz]
	GLfloat ce[12][2][3] = {
		{{1.0f, 1.0f, -1.0f}, {0.0f, 0.0f, 1.0f}},
		{{-1.0f, 1.0f, -1.0f}, {0.0f, 0.0f, 1.0f}},
		{{-1.0f, -1.0f, -1.0f}, {0.0f, 0.0f, 1.0f}},
		{{1.0f, -1.0f, -1.0f}, {0.0f, 0.0f, 1.0f}},

		{{1.0f, -1.0f, 1.0f}, {0.0f, 1.0f, 0.0f}},
		{{-1.0f, -1.0f, 1.0f}, {0.0f, 1.0f, 0.0f}},
		{{-1.0f, -1.0f, -1.0f}, {0.0f, 1.0f, 0.0f}},
		{{1.0f, -1.0f, -1.0f}, {0.0f, 1.0f, 0.0f}},

		{{-1.0f, 1.0f, 1.0f}, {1.0f, 0.0f, 0.0f}},
		{{-1.0f, -1.0f, 1.0f}, {1.0f, 0.0f, 0.0f}},
		{{-1.0f, -1.0f, -1.0f}, {1.0f, 0.0f, 0.0f}},
		{{-1.0f, 1.0f, -1.0f}, {1.0f, 0.0f, 0.0f}}
	};

	memcpy(_cubeVertices, cv, sizeof(_cubeVertices) );
	memcpy(_cubeEdges, ce, sizeof(_cubeEdges) );

	InitGL();
}

Renderer::~Renderer()
{
	if (_textureData)
		free(_textureData);
}

void Renderer::SetLightPostion(Eigen::Vector3f &pos)
{
	_lightPos = pos;
	_lightDir = -_lightPos;

	GenerateRayTemplate(RES);
}

void Renderer::InitGL()
{
	glEnable(GL_TEXTURE_3D);
	glDisable(GL_DEPTH_TEST);
	glCullFace(GL_FRONT);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


	glGenTextures(2, &_hTexture);

	glActiveTextureARB(GL_TEXTURE0_ARB);
	glBindTexture(GL_TEXTURE_3D, _hTexture);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S,GL_CLAMP);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T,GL_CLAMP);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R,GL_CLAMP);

}

void Renderer::Render(void)
{
	GLdouble mvMatrix[16];
	glGetDoublev(GL_MODELVIEW_MATRIX, mvMatrix);

	glDisable(GL_DEPTH_TEST);	//Important in this rendering

	DrawCube();
	
	//core rendering parts
	DrawSlices(mvMatrix);

}


void Renderer::DrawCube(void)
{
	glDisable(GL_TEXTURE_3D);
	glDisable(GL_FRAGMENT_PROGRAM_ARB);

	glEnable(GL_CULL_FACE);
	float coeff = 0.5;
	glColor4f(0.9f*coeff, 0.8f*coeff, 0.4f*coeff, 1.0f);
	GLfloat (*cv)[3] = _cubeVertices;

	glBegin(GL_POLYGON);
	glVertex3fv(cv[0]); glVertex3fv(cv[1]); glVertex3fv(cv[2]); glVertex3fv(cv[3]);
	glEnd();
	glBegin(GL_POLYGON);
	glVertex3fv(cv[0]); glVertex3fv(cv[4]); glVertex3fv(cv[5]); glVertex3fv(cv[1]);
	glEnd();
	glBegin(GL_POLYGON);
	glVertex3fv(cv[0]); glVertex3fv(cv[3]); glVertex3fv(cv[7]); glVertex3fv(cv[4]);
	glEnd();

	glBegin(GL_POLYGON);
	glVertex3fv(cv[7]); glVertex3fv(cv[6]); glVertex3fv(cv[5]); glVertex3fv(cv[4]);
	glEnd();
	glBegin(GL_POLYGON);
	glVertex3fv(cv[2]); glVertex3fv(cv[6]); glVertex3fv(cv[7]); glVertex3fv(cv[3]);
	glEnd();
	glBegin(GL_POLYGON);
	glVertex3fv(cv[1]); glVertex3fv(cv[5]); glVertex3fv(cv[6]); glVertex3fv(cv[2]);
	glEnd();

	glDisable(GL_CULL_FACE);

	glBegin(GL_LINES);
	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
	glVertex3fv(cv[0]); glVertex3fv(cv[1]);
	glVertex3fv(cv[1]); glVertex3fv(cv[2]);
	glVertex3fv(cv[2]); glVertex3fv(cv[3]);
	glVertex3fv(cv[3]); glVertex3fv(cv[0]);

	glVertex3fv(cv[4]); glVertex3fv(cv[5]);
	glVertex3fv(cv[5]); glVertex3fv(cv[6]);
	glVertex3fv(cv[6]); glVertex3fv(cv[7]);
	glVertex3fv(cv[7]); glVertex3fv(cv[4]);

	glVertex3fv(cv[0]); glVertex3fv(cv[4]);
	glVertex3fv(cv[1]); glVertex3fv(cv[5]);
	glVertex3fv(cv[2]); glVertex3fv(cv[6]);
	glVertex3fv(cv[3]); glVertex3fv(cv[7]);
	glEnd();

	
	//draw light
	glPointSize(13.0f);
	glBegin(GL_POINTS);
	glColor4f(0.0f, 1.0f, 1.0f, 1.0f);
	glVertex3f(_lightPos[0], _lightPos[1], _lightPos[2]);
	glEnd();


}

class Convexcomp
{
	private:
		const Eigen::Vector3f &p0, &up;
	public:
		Convexcomp(const Eigen::Vector3f& p0, const Eigen::Vector3f& up) : p0(p0), up(up) {}

		bool operator()(const Eigen::Vector3f& a, const Eigen::Vector3f& b) const
		{
			Eigen::Vector3f va = a-p0, vb = b-p0;
			//return dot(up, cross(va, vb)) >= 0;
			return up.dot(va.cross(vb)) >= 0;
		}
};

void Renderer::DrawSlices(GLdouble mvMatrix[16])
{
	int i;
	Eigen::Vector3f viewdir(-mvMatrix[2], -mvMatrix[6], -mvMatrix[10]);	//FIXME

	viewdir.normalize();
	// find cube vertex that is closest to the viewer
	GLfloat (*cv)[3] = _cubeVertices;
	for (i=0; i<8; i++) {
		float x = cv[i][0] + viewdir[0];
		float y = cv[i][1] + viewdir[1];
		float z = cv[i][2] + viewdir[2];
		if ((x>=-1.0f)&&(x<=1.0f)
				&&(y>=-1.0f)&&(y<=1.0f)
				&&(z>=-1.0f)&&(z<=1.0f))
		{
			break;
		}
	}
	assert(i != 8);

	// our slices are defined by the plane equation A*x + B*y +C*z + D = 0
	// (a,b,c), the plane normal, are given by viewdir
	// d is the parameter along the view direction. the first d is given by
	// inserting previously found vertex into the plane equation
	float d0 = -(viewdir[0]*cv[i][0] + viewdir[1]*cv[i][1] + viewdir[2]*cv[i][2]);
	float dStep = 2*d0/SLICE_NUM;
	int n = 0;
	for (float d = -d0; d < d0; d += dStep) {
		// IntersectEdges returns the intersection points of all cube edges with
		// the given plane that lie within the cube
		std::vector<Eigen::Vector3f> pt = IntersectEdges(viewdir[0], viewdir[1], viewdir[2], d);

		if (pt.size() > 2) {
			// sort points to get a convex polygon
			std::sort(pt.begin()+1, pt.end(), Convexcomp(pt[0], viewdir));

			glEnable(GL_TEXTURE_3D);
			glBegin(GL_POLYGON);
			for (i=0; i<pt.size(); i++){
				glColor3f(1.0, 1.0, 1.0);
				glTexCoord3d((pt[i][0]+1.0)/2.0, (pt[i][1]+1)/2.0, (pt[i][2]+1.0)/2.0);//FIXME
				glVertex3f(pt[i][0], pt[i][1], pt[i][2]);
			}
			glEnd();

			if (_isDrawSliceOutline)
			{
				glDisable(GL_TEXTURE_3D);
				glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
				glBegin(GL_POLYGON);
				for (i=0; i<pt.size(); i++) {
					glColor3f(0.0, 0.0, 1.0);
					glVertex3f(pt[i][0], pt[i][1], pt[i][2]);
				}
				glEnd();
				glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			}
		}
		n++;
	}
}

std::vector<Eigen::Vector3f> Renderer::IntersectEdges(float A, float B, float C, float D)
{
	float t;
	Eigen::Vector3f p;
	std::vector<Eigen::Vector3f> res;
	GLfloat (*edges)[2][3] = _cubeEdges;
	

	for (int i=0; i<12; i++) {
		t = -(A*edges[i][0][0] + B*edges[i][0][1] + C*edges[i][0][2] + D)
			/ (A*edges[i][1][0] + B*edges[i][1][1] + C*edges[i][1][2]);
		if ((t>0)&&(t<2)) {
			p[0] = edges[i][0][0] + edges[i][1][0]*t;
			p[1] = edges[i][0][1] + edges[i][1][1]*t;
			p[2] = edges[i][0][2] + edges[i][1][2]*t;
			res.push_back(p);
		}
	}

	return res;
}


void Renderer::FillTexture(Fluid* fluid)
{
	if (_textureData == NULL) {
		_textureData = new unsigned char[SIZE*4];
	}
	memset(_textureData, 0, SIZE*4);

	const float *density = fluid->GetDensity();
	unsigned char* intensity = new unsigned char[SIZE];
	memset(intensity, 0, SIZE);
	CastLight(RES, density, intensity);

#if 1
	//FIXME: It is important to beware that texture coordinate
	//is in reverse order  of the simulation coordinate
	for (int i = 0 ; i < RES ; i++) {
		for (int j = 0 ; j < RES ; j++) {
			for (int k = 0 ; k < RES ; k++) {
				//unsigned char c = 200;
				int texIndex = i*RES*RES + j*RES + k;	/*reverse order*/
				int densIndex = k*RES*RES + j*RES + i;
				unsigned char c = intensity[densIndex];
				_textureData[texIndex*4] = c;
				_textureData[texIndex*4+1] = c;
				_textureData[texIndex*4+2] = c;
				_textureData[texIndex*4+3] = (density[densIndex]>0.1f) ? 
					255 : (unsigned char) (density[densIndex]*2550.0f);
			}
		}
	}
#else
	for (int i=0; i<SIZE; i++) {
		unsigned char c = intensity[i];
		//unsigned char c = 200;
		_textureData[(i<<2)] = c;
		_textureData[(i<<2)+1] = c;
		_textureData[(i<<2)+2] = c;
		_textureData[(i<<2)+3] = (density[i]>0.1f) ? 255 : (unsigned char) (density[i]*2550.0f);
	}
#endif

	delete []intensity;

	glActiveTextureARB(GL_TEXTURE0_ARB);
	glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA, RES, RES, RES, 0, GL_RGBA, GL_UNSIGNED_BYTE, _textureData);
}

void Renderer::GenerateRayTemplate(int edgeLen)
{
	_rayTemplate[0][0] = _rayTemplate[0][2] = _rayTemplate[0][2] = 0;
	float fx = 0.0f, fy = 0.0f, fz = 0.0f;
	int x = 0, y = 0, z = 0;
	float lx = _lightDir[0] + 0.000001f, ly = _lightDir[1] + 0.000001f, lz = _lightDir[2] + 0.000001f;
	int xinc = (lx > 0) ? 1 : -1;
	int yinc = (ly > 0) ? 1 : -1;
	int zinc = (lz > 0) ? 1 : -1;
	float tx, ty, tz;
	int i = 1;
	int len = 0;
	int maxlen = 3*edgeLen*edgeLen;
	while (len <= maxlen)
	{
		// fx + t*lx = (x+1)   ->   t = (x+1-fx)/lx
		tx = (x+xinc-fx)/lx;
		ty = (y+yinc-fy)/ly;
		tz = (z+zinc-fz)/lz;

		if ((tx<=ty)&&(tx<=tz)) {
			_rayTemplate[i][0] = _rayTemplate[i-1][0] + xinc;
			x =+ xinc;
			fx = x;

			if (ALMOST_EQUAL(ty,tx)) {
				_rayTemplate[i][1] = _rayTemplate[i-1][1] + yinc;
				y += yinc;
				fy = y;
			} else {
				_rayTemplate[i][1] = _rayTemplate[i-1][1];
				fy += tx*ly;
			}

			if (ALMOST_EQUAL(tz,tx)) {
				_rayTemplate[i][2] = _rayTemplate[i-1][2] + zinc;
				z += zinc;
				fz = z;
			} else {
				_rayTemplate[i][2] = _rayTemplate[i-1][2];
				fz += tx*lz;
			}
		} else if ((ty<tx)&&(ty<=tz)) {
			_rayTemplate[i][0] = _rayTemplate[i-1][0];
			fx += ty*lx;

			_rayTemplate[i][1] = _rayTemplate[i-1][1] + yinc;
			y += yinc;
			fy = y;

			if (ALMOST_EQUAL(tz,ty)) {
				_rayTemplate[i][2] = _rayTemplate[i-1][2] + zinc;
				z += zinc;
				fz = z;
			} else {
				_rayTemplate[i][2] = _rayTemplate[i-1][2];
				fz += ty*lz;
			}
		} else {
			assert((tz<tx)&&(tz<ty));
			_rayTemplate[i][0] = _rayTemplate[i-1][0];
			fx += tz*lx;
			_rayTemplate[i][1] = _rayTemplate[i-1][1];
			fy += tz*ly;
			_rayTemplate[i][2] = _rayTemplate[i-1][2] + zinc;
			z += zinc;
			fz = z;
		}

		len = _rayTemplate[i][0]*_rayTemplate[i][0]
			+ _rayTemplate[i][1]*_rayTemplate[i][1]
			+ _rayTemplate[i][2]*_rayTemplate[i][2];
		i++;
	}
}

#define DECAY 0.06f
void Renderer::CastLight(int n /*edgelen*/, const float* dens, unsigned char* intensity)
{

	int i,j;
	int sx = (_lightDir[0]>0) ? 0 : n-1;
	int sy = (_lightDir[1]>0) ? 0 : n-1;
	int sz = (_lightDir[2]>0) ? 0 : n-1;

	float decay = 1.0f/(n*DECAY);

	for (i=0; i<n; i++)
		for (j=0; j<n; j++) {
			if (!ALMOST_EQUAL(_lightDir[0], 0) )
				LightRay(sx, i, j, n, decay, dens, intensity);
			if (!ALMOST_EQUAL(_lightDir[1], 0) )
				LightRay(i, sy, j, n, decay, dens, intensity);
			if (!ALMOST_EQUAL(_lightDir[2], 0) )
				LightRay(i, j, sz, n, decay, dens, intensity);
		}
}


#define AMBIENT 100
inline void Renderer::LightRay(int x, int y, int z, int n, float decay, const float* dens, unsigned char* intensity)
{
	int xx = x, yy = y, zz = z, i = 0;
	int offset;

	int l = 200;
	float d;

	do {
		offset = ((xx*n) + yy)*n + zz;//FIXME
		if (intensity[offset] > 0)
			intensity[offset] = (unsigned char) ((intensity[offset] + l)*0.5f);
		else
			intensity[offset] = (unsigned char) l;
		d = dens[offset]*255.0f;
		if (l > AMBIENT) {
			l -= d*decay;
			if (l < AMBIENT)
				l = AMBIENT;
		}

		i++;
		xx = x + _rayTemplate[i][0];
		yy = y + _rayTemplate[i][1];
		zz = z + _rayTemplate[i][2];
	} while ((xx>=0)&&(xx<n)&&(yy>=0)&&(yy<n)&&(zz>=0)&&(zz<n));
}

