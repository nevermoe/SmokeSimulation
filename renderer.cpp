#include "renderer.h"
#include <GL/glut.h>

int const fbo_width = 1280;
int const fbo_height = 960;

GLuint fb, color, depth;

void CHECK_FRAMEBUFFER_STATUS()
{                                                         
	std::cout << __FUNCTION__ << " " << __LINE__ << std::endl;
	GLenum status;
	status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER); 
	switch(status) {
	case GL_FRAMEBUFFER_COMPLETE:
		break;

	case GL_FRAMEBUFFER_UNSUPPORTED:
	/* choose different formats */
		std::cout << "unsupport" << std::endl;
		break;

	default:
		/* programming error; will fail on all hardware */
		fputs("Framebuffer Error\n", stderr);
		exit(-1);
	}
}

float const light_dir[]={1,1,1,0};
float const light_color[]={1,0.95,0.9,1};

void Renderer::prepare()
{
	static float a=0, b=0, c=0;

	//glBindTexture(GL_TEXTURE_2D, 0);
	//_showTex->UnBind();
	glEnable(GL_TEXTURE_2D);
#if 0
	glBindFramebuffer(GL_FRAMEBUFFER, fb);
#else
	_fbo->Bind();
#endif

	glViewport(0,0, fbo_width, fbo_height);

	glClearColor(1,1,1,0);
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45, 1, 1, 10);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);

	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);

	glLightfv(GL_LIGHT0, GL_POSITION, light_dir);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_color);

	glTranslatef(0,0,-5);

	glRotatef(a, 1, 0, 0);
	glRotatef(b, 0, 1, 0);
	glRotatef(c, 0, 0, 1);

	glutSolidTeapot(0.75);

	a=fmod(a+0.1, 360.);
	b=fmod(b+0.5, 360.);
	c=fmod(c+0.25, 360.);
#if 0
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
#else
	_fbo->Disable();
#endif
}

void Renderer::final1()
{
	static float a=0, b=0, c=0;

	const int win_width  = 1280;
	const int win_height = 960;
	const float aspect = (float)win_width/(float)win_height;

	
	glViewport(0,0, win_width, win_height);

	glClearColor(1.,1.,1.,0.);
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45, aspect, 1, 10);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0,0,-5);

	//glRotatef(b, 0, 1, 0);

	//b=fmod(b+0.5, 360.);

	glEnable(GL_TEXTURE_2D);
#if 0
	glBindTexture(GL_TEXTURE_2D, color);
#else
	_showTex->Bind();
#endif

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glDisable(GL_LIGHTING);

#if 0
	float cube[][5]=
	{
		{-1, -1, -1,  0,  0},
		{ 1, -1, -1,  1,  0},
		{ 1,  1, -1,  1,  1},
		{-1,  1, -1,  0,  1},

		{-1, -1,  1, -1,  0},
		{ 1, -1,  1,  0,  0},
		{ 1,  1,  1,  0,  1},
		{-1,  1,  1, -1,  1},
	};
	unsigned int faces[]=
	{
		0, 1, 2, 3,
		1, 5, 6, 2,
		5, 4, 7, 6,
		4, 0, 3, 7,
		3, 2, 6, 7,
		4, 5, 1, 0
	};

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_TEXTURE_COORD_ARRAY);

	glVertexPointer(3, GL_FLOAT, 5*sizeof(float), &cube[0][0]);
	glTexCoordPointer(2, GL_FLOAT, 5*sizeof(float), &cube[0][3]);

	glCullFace(GL_BACK);
	glDrawElements(GL_QUADS, 24, GL_UNSIGNED_INT, faces);

	//glCullFace(GL_FRONT);
	//glDrawElements(GL_QUADS, 24, GL_UNSIGNED_INT, faces);

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);
#else

	glBegin(GL_QUADS);

	glVertex3f(-1.0f, -1.0f, -1.0f); // The bottom left corner
	glTexCoord2f(0.0f, 0.0f);

	glVertex3f(1.0f, -1.0f, -1.0f); // The top left corner
	glTexCoord2f(1.0f, 0.0f);

	glVertex3f(1.0f, 1.0f, -1.0f); // The top right corner
	glTexCoord2f(1.0f, 1.0f);

	glVertex3f(-1.0f, 1.0f, -1.0f); // The bottom right corner
	glTexCoord2f(0.0f, 1.0f);

	glEnd();
#endif

	_showTex->UnBind();

}



//#################################################################################
Renderer::Renderer(float* volumeData, int RES)
{

	_isDrawSliceOutline = false;
	_isRendering = true;

	_RES = RES;
	_N = _RES - 2;
	_SIZE = _RES*_RES*_RES;
	_volumeData = volumeData;

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

	//create a volume texture
	_textureData = NULL;
	_volumeTex3D = new GLTexture(_RES, _RES, _RES, GL_RGBA, GL_RGBA,
			GL_UNSIGNED_BYTE, GL_LINEAR, GL_CLAMP);

#if 0
	//create framebuffer and bind to a texture
	_showTex = new GLTexture(1280, 960, 0, GL_RGBA, GL_RGBA, 
			GL_UNSIGNED_BYTE, GL_LINEAR, GL_CLAMP);
	_showTex->LoadToGPU();
	
	_fbo = new FramebufferObject();
	_fbo->Bind();
	_fbo->AttachTexture(GL_TEXTURE_2D, _showTex->GetTextureID(), GL_COLOR_ATTACHMENT0_EXT);
	_fbo->IsValid();
	_fbo->Disable();
#endif

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

	GenerateRayTemplate(_RES);
}

void Renderer::InitGL()
{
	//glEnable(GL_TEXTURE_3D);
	//glDisable(GL_DEPTH_TEST);
	//glCullFace(GL_FRONT);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	//############################################
	glGenFramebuffers(1, &fb);
	glGenTextures(1, &color);
	glGenRenderbuffers(1, &depth);

	glBindFramebuffer(GL_FRAMEBUFFER, fb);

#if 0
	glBindTexture(GL_TEXTURE_2D, color);
	glTexImage2D(	GL_TEXTURE_2D, 
			0, 
			GL_RGBA, 
			fbo_width, fbo_height,
			0, 
			GL_RGBA, 
			GL_UNSIGNED_BYTE, 
			NULL);
	
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, color, 0);

#else
	_showTex = new GLTexture(1280, 960, 0, GL_RGBA, GL_RGBA, 
			GL_UNSIGNED_BYTE, GL_LINEAR, GL_CLAMP);
	_showTex->LoadToGPU(NULL);
#endif

	//glBindRenderbuffer(GL_RENDERBUFFER, depth);
	//glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, fbo_width, fbo_height);
	//glFramebufferRenderbuffer(GL_DRAW_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depth);

	_fbo = new FramebufferObject();
	_fbo->Bind();
	_fbo->AttachTexture(GL_TEXTURE_2D, _showTex->GetTextureID(), GL_COLOR_ATTACHMENT0_EXT);
	_fbo->IsValid();
	_fbo->Disable();

	CHECK_FRAMEBUFFER_STATUS();

}

void Renderer::drawSlice(float z)
{
    glBegin(GL_QUADS);
    glTexCoord3f(0.0f, 0.0f, z); glVertex2f(-1.0f, -1.0f);
    glTexCoord3f(1.0f, 0.0f, z); glVertex2f(1.0f, -1.0f);
    glTexCoord3f(1.0f, 1.0f, z); glVertex2f(1.0f, 1.0f);
    glTexCoord3f(0.0f, 1.0f, z); glVertex2f(-1.0f, 1.0f);
    glEnd();
}

void Renderer::Render(void)
{

#if 1

	//glEnable(GL_TEXTURE_2D);
	_fbo->Bind();
#if 1
	glClearColor(1,1,1,0);
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	//glPushAttrib(GL_VIEWPORT_BIT);	//save current state

	//glViewport(0,0, fbo_width, fbo_height);
	GLdouble mvMatrix[16];
	glGetDoublev(GL_MODELVIEW_MATRIX, mvMatrix);

	glDisable(GL_DEPTH_TEST);	//Important in this rendering

	DrawCube();
	DrawLight();
	if(_isRendering) {
		//core rendering parts
		DrawSlices(mvMatrix);
	}
	else 
		DrawVolumeData();

	glEnable(GL_DEPTH_TEST);

	//glPopAttrib();
#else
	prepare();
#endif

	_fbo->Disable();


	//###################
#if 1
	glEnable(GL_TEXTURE_2D);
#else
	glEnable(GL_TEXTURE_3D);
#endif
	//glPushMatrix();
	//const int win_width  = 1280;
	//const int win_height = 960;
	//const float aspect = (float)win_width/(float)win_height;
	//glMatrixMode(GL_PROJECTION);
	//glLoadIdentity();
	//gluPerspective(45, aspect, 1, 10);

	//glMatrixMode(GL_MODELVIEW);
	//glLoadIdentity();
	//glTranslatef(0,0,-5);

	_showTex->Bind();
	glBegin(GL_QUADS);


#if 0
    for(int z=0; z<_RES; z++) {
        // attach texture slice to FBO
        _fbo->AttachTexture(GL_TEXTURE_3D, _showTex->GetTextureID(), GL_COLOR_ATTACHMENT0_EXT, 0, z);
        // render
        drawSlice(z + 0.5);
		//std::cout << z << std::endl;
    }

#else

	
	glVertex3f(-1.0f, -1.0f, -1.0f); // The bottom left corner
	glTexCoord3f(-1.0f, -1.0f, -1.0f); // The bottom left corner
	//glTexCoord2f(0.0f, 0.0f);

	glVertex3f(1.0f, -1.0f, -1.0f); // The top left corner
	glTexCoord3f(1.0f, -1.0f, -1.0f); // The top left corner
	//glTexCoord2f(1.0f, 0.0f);

	glVertex3f(1.0f, 1.0f, -1.0f); // The top right corner
	glTexCoord3f(1.0f, 1.0f, -1.0f); // The top right corner
	//glTexCoord2f(1.0f, 1.0f);

	glVertex3f(-1.0f, 1.0f, -1.0f); // The bottom right corner
	glTexCoord3f(-1.0f, 1.0f, -1.0f); // The bottom right corner
	//glTexCoord2f(0.0f, 1.0f);

	glEnd();
#endif


	_showTex->UnBind();
	//glPopMatrix();

//	glDisable(GL_TEXTURE_3D);
#else
	prepare();
	final1();
#endif


}

void Renderer::DrawLight()
{
	
	glPushAttrib(GL_POINT_BIT | GL_CURRENT_BIT);	//save current state
	//draw light
	glPointSize(13.0f);
	glBegin(GL_POINTS);
	glColor4f(0.0f, 1.0f, 1.0f, 1.0f);
	glVertex3f(_lightPos[0], _lightPos[1], _lightPos[2]);
	glEnd();
	glPopAttrib();		//restore painting state
}

void Renderer::DrawVolumeData()
{
	glBegin(GL_POINTS);
	FOR_ALL_CELL {
		if(!ALMOST_EQUAL(_volumeData[_I(i,j,k)], 0) ) {
			glVertex3f(((float)i/_N-0.5)*2, ((float)j/_N-0.5)*2, ((float)k/_N-0.5)*2 );
		}
	} END_FOR
	glEnd();

	glPointSize(13.0f);
	glBegin(GL_POINTS);
	glColor4f(0.0f, 1.0f, 1.0f, 1.0f);
	glVertex3f(_lightPos[0], _lightPos[1], _lightPos[2]);
	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
	glEnd();
	glPointSize(1.0f);

}


void Renderer::SetRendering(bool isRendering)
{
	_isRendering = isRendering;
}

void Renderer::SetSliceOutline(bool isDrawSliceOutline)
{
	_isDrawSliceOutline = isDrawSliceOutline;
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

	glDisable(GL_TEXTURE_3D);
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


void Renderer::FillTexture()
{
	if (_textureData == NULL) {
		_textureData = new unsigned char[_SIZE*4];
	}
	memset(_textureData, 0, _SIZE*4);

	const float *density = _volumeData;
	unsigned char* intensity = new unsigned char[_SIZE];
	memset(intensity, 0, _SIZE);
	CastLight(_RES, density, intensity);

#if 1
	//FIXME: It is important to beware that texture coordinate
	//is in reverse order  of the simulation coordinate
	for (int i = 0 ; i < _RES ; i++) {
		for (int j = 0 ; j < _RES ; j++) {
			for (int k = 0 ; k < _RES ; k++) {
				//unsigned char c = 200;
				int texIndex = i*_RES*_RES + j*_RES + k;	/*reverse order*/
				int densIndex = k*_RES*_RES + j*_RES + i;
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
	for (int i=0; i<_SIZE; i++) {
		unsigned char c = intensity[i];
		//unsigned char c = 200;
		_textureData[(i<<2)] = c;
		_textureData[(i<<2)+1] = c;
		_textureData[(i<<2)+2] = c;
		_textureData[(i<<2)+3] = (density[i]>0.1f) ? 255 : (unsigned char) (density[i]*2550.0f);
	}
#endif

	delete []intensity;

	_volumeTex3D->LoadToGPU(_textureData);
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

