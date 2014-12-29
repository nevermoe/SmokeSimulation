/* Author: Johannes Schmid and Ingemar Rask, 2006, johnny@grob.org */
#include "viewer.h"
#include "fluid.h"

bool canRunFragmentProgram(const char* programString);	// from FragmentUtils.cpp

// cube vertices
GLfloat cv[][3] = {
	{1.0f, 1.0f, 1.0f}, {-1.0f, 1.0f, 1.0f}, {-1.0f, -1.0f, 1.0f}, {1.0f, -1.0f, 1.0f},
	{1.0f, 1.0f, -1.0f}, {-1.0f, 1.0f, -1.0f}, {-1.0f, -1.0f, -1.0f}, {1.0f, -1.0f, -1.0f}
};

// edges have the form edges[n][0][xyz] + t*edges[n][1][xyz]
float edges[12][2][3] = {
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

Viewer::Viewer()
{
	std::cout << __FILE__ << " " << __FUNCTION__ << std::endl;

	_dist = 5.0f;
	_fp = NULL;
	_texture_data = NULL;

	_draw_slice_outline = true;

	init_GL();

	/*
	   _light_dir[0] = -1.0f;
	   _light_dir[1] = 0.50f;
	   _light_dir[2] = 0.0f;
	   */
	_light_dir[0] = 0.0f;
	_light_dir[1] = 1.0f;
	_light_dir[2] = -1.0f;
	gen_ray_templ(N+2);
}

Viewer::~Viewer()
{
	if (_texture_data)
		free(_texture_data);
	if (_fp)
		fclose(_fp);
}

void Viewer::init_GL(void)
{
	glEnable(GL_TEXTURE_3D);
	glDisable(GL_DEPTH_TEST);
	glCullFace(GL_FRONT);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


	glGenTextures(1, &_txt[0]);

	glActiveTextureARB(GL_TEXTURE0_ARB);
	glBindTexture(GL_TEXTURE_3D, _txt[0]);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S,GL_CLAMP);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T,GL_CLAMP);
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R,GL_CLAMP);

}

void Viewer::draw(void)
{
	GLdouble m[4][4], mt[4][4];
	glGetDoublev(GL_MODELVIEW_MATRIX, &m[0][0]);
	

	glDisable(GL_DEPTH_TEST);	//Important in this rendering

	draw_cube();

	draw_slices(m, _draw_slice_outline);

}


void Viewer::draw_cube(void)
{
	glDisable(GL_TEXTURE_3D);
	glDisable(GL_FRAGMENT_PROGRAM_ARB);

	glEnable(GL_CULL_FACE);
	float l = 0.5;
	glColor4f(0.9f*l, 0.8f*l, 0.4f*l, 1.0f);

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

	/*
	   glPointSize(3.0f);
	   glBegin(GL_POINTS);
	   glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
#define PDIST 2.0f
glVertex3f(-_light_dir[0]*PDIST,-_light_dir[1]*PDIST,-_light_dir[2]*PDIST);
#undef PDIST
glEnd();
*/

}

class Convexcomp
{
	private:
		const Vec3 &p0, &up;
	public:
		Convexcomp(const Vec3& p0, const Vec3& up) : p0(p0), up(up) {}

		bool operator()(const Vec3& a, const Vec3& b) const
		{
			Vec3 va = a-p0, vb = b-p0;
			return dot(up, cross(va, vb)) >= 0;
		}
};

void Viewer::draw_slices(GLdouble m[][4], bool frame)
{
	int i;

	//printf("draw_slices\n");
	Vec3 viewdir(m[0][2], m[1][2], m[2][2]);
	viewdir.Normalize();
	// find cube vertex that is closest to the viewer
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

	// our slices are defined by the plane equation a*x + b*y +c*z + d = 0
	// (a,b,c), the plane normal, are given by viewdir
	// d is the parameter along the view direction. the first d is given by
	// inserting previously found vertex into the plane equation
	float d0 = -(viewdir[0]*cv[i][0] + viewdir[1]*cv[i][1] + viewdir[2]*cv[i][2]);
	float dd = 2*d0/64.0f;
	int n = 0;
	for (float d = -d0; d < d0; d += dd) {
		// intersect_edges returns the intersection points of all cube edges with
		// the given plane that lie within the cube
		std::vector<Vec3> pt = intersect_edges(viewdir[0], viewdir[1], viewdir[2], d);

		if (pt.size() > 2) {
			// sort points to get a convex polygon
			std::sort(pt.begin()+1, pt.end(), Convexcomp(pt[0], viewdir));

			glEnable(GL_TEXTURE_3D);
			glBegin(GL_POLYGON);
			for (i=0; i<pt.size(); i++){
				glColor3f(1.0, 1.0, 1.0);
				glTexCoord3d((pt[i][0]+1.0)/2.0, (-pt[i][1]+1)/2.0, (pt[i][2]+1.0)/2.0);
				glVertex3f(pt[i][0], pt[i][1], pt[i][2]);
			}
			glEnd();

#if 1
			frame = false;
			if (frame)
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
#endif
		}
		n++;
	}

	//	printf("Sliced: %d\n", n);
}

std::vector<Vec3> Viewer::intersect_edges(float a, float b, float c, float d)
{
	int i;
	float t;
	Vec3 p;
	std::vector<Vec3> res;

	for (i=0; i<12; i++) {
		t = -(a*edges[i][0][0] + b*edges[i][0][1] + c*edges[i][0][2] + d)
			/ (a*edges[i][1][0] + b*edges[i][1][1] + c*edges[i][1][2]);
		if ((t>0)&&(t<2)) {
			p[0] = edges[i][0][0] + edges[i][1][0]*t;
			p[1] = edges[i][0][1] + edges[i][1][1]*t;
			p[2] = edges[i][0][2] + edges[i][1][2]*t;
			res.push_back(p);
		}
	}

	return res;
}


void Viewer::frame_from_sim(Fluid* fluid)
{
	if (_texture_data == NULL)
		_texture_data = (unsigned char*) malloc((N+2)*(N+2)*(N+2)*4);

	unsigned char* l = (unsigned char*) malloc((N+2)*(N+2)*(N+2));
	memset(l,0,(N+2)*(N+2)*(N+2));
	cast_light(N+2, fluid->d, l);

	for (int i=0; i<(N+2)*(N+2)*(N+2); i++) {
#if 0
		unsigned char c = l[i];
		//printf("%f ",l[i]);
#else
		unsigned char c = 200-fluid->d[i];
		//printf("%d ",c);
#endif
		_texture_data[(i<<2)] = c;
		_texture_data[(i<<2)+1] = c;
		_texture_data[(i<<2)+2] = c;
		_texture_data[(i<<2)+3] = (fluid->d[i]>0.1f) ? 255 : (unsigned char) (fluid->d[i]*2550.0f);
	}

	free(l);

	glActiveTextureARB(GL_TEXTURE0_ARB);
	glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA, N+2, N+2, N+2, 0, GL_RGBA, GL_UNSIGNED_BYTE, _texture_data);
}

#define ALMOST_EQUAL(a, b) ((fabs(a-b)<0.00001f)?true:false)
void Viewer::gen_ray_templ(int edgelen)
{
	_ray_templ[0][0] = _ray_templ[0][2] = _ray_templ[0][2] = 0;
	float fx = 0.0f, fy = 0.0f, fz = 0.0f;
	int x = 0, y = 0, z = 0;
	float lx = _light_dir[0] + 0.000001f, ly = _light_dir[1] + 0.000001f, lz = _light_dir[2] + 0.000001f;
	int xinc = (lx > 0) ? 1 : -1;
	int yinc = (ly > 0) ? 1 : -1;
	int zinc = (lz > 0) ? 1 : -1;
	float tx, ty, tz;
	int i = 1;
	int len = 0;
	int maxlen = 3*edgelen*edgelen;
	while (len <= maxlen)
	{
		// fx + t*lx = (x+1)   ->   t = (x+1-fx)/lx
		tx = (x+xinc-fx)/lx;
		ty = (y+yinc-fy)/ly;
		tz = (z+zinc-fz)/lz;

		if ((tx<=ty)&&(tx<=tz)) {
			_ray_templ[i][0] = _ray_templ[i-1][0] + xinc;
			x =+ xinc;
			fx = x;

			if (ALMOST_EQUAL(ty,tx)) {
				_ray_templ[i][1] = _ray_templ[i-1][1] + yinc;
				y += yinc;
				fy = y;
			} else {
				_ray_templ[i][1] = _ray_templ[i-1][1];
				fy += tx*ly;
			}

			if (ALMOST_EQUAL(tz,tx)) {
				_ray_templ[i][2] = _ray_templ[i-1][2] + zinc;
				z += zinc;
				fz = z;
			} else {
				_ray_templ[i][2] = _ray_templ[i-1][2];
				fz += tx*lz;
			}
		} else if ((ty<tx)&&(ty<=tz)) {
			_ray_templ[i][0] = _ray_templ[i-1][0];
			fx += ty*lx;

			_ray_templ[i][1] = _ray_templ[i-1][1] + yinc;
			y += yinc;
			fy = y;

			if (ALMOST_EQUAL(tz,ty)) {
				_ray_templ[i][2] = _ray_templ[i-1][2] + zinc;
				z += zinc;
				fz = z;
			} else {
				_ray_templ[i][2] = _ray_templ[i-1][2];
				fz += ty*lz;
			}
		} else {
			assert((tz<tx)&&(tz<ty));
			_ray_templ[i][0] = _ray_templ[i-1][0];
			fx += tz*lx;
			_ray_templ[i][1] = _ray_templ[i-1][1];
			fy += tz*ly;
			_ray_templ[i][2] = _ray_templ[i-1][2] + zinc;
			z += zinc;
			fz = z;
		}

		len = _ray_templ[i][0]*_ray_templ[i][0]
			+ _ray_templ[i][1]*_ray_templ[i][1]
			+ _ray_templ[i][2]*_ray_templ[i][2];
		i++;
	}
}

#define DECAY 0.06f
void Viewer::cast_light(int n /*edgelen*/, float* dens, unsigned char* intensity)
{

	int i,j;
	int sx = (_light_dir[0]>0) ? 0 : n-1;
	int sy = (_light_dir[1]>0) ? 0 : n-1;
	int sz = (_light_dir[2]>0) ? 0 : n-1;

	float decay = 1.0f/(n*DECAY);

	for (i=0; i<n; i++)
		for (j=0; j<n; j++) {
			if (!ALMOST_EQUAL(_light_dir[0], 0))
				light_ray(sx,i,j,n,decay,dens,intensity);
			if (!ALMOST_EQUAL(_light_dir[1], 0))
				light_ray(i,sy,j,n,decay,dens,intensity);
			if (!ALMOST_EQUAL(_light_dir[2], 0))
				light_ray(i,j,sz,n,decay,dens,intensity);
		}
}


#define AMBIENT 100
inline void Viewer::light_ray(int x, int y, int z, int n, float decay, float* dens, unsigned char* intensity)
{
	int xx = x, yy = y, zz = z, i = 0;
	int offset;

	int l = 200;
	float d;

	do {
		offset = ((zz*n) + yy)*n + xx;
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
		xx = x + _ray_templ[i][0];
		yy = y + _ray_templ[i][1];
		zz = z + _ray_templ[i][2];
	} while ((xx>=0)&&(xx<n)&&(yy>=0)&&(yy<n)&&(zz>=0)&&(zz<n));
}

