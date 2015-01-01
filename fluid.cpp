/* Adapted from: Johannes Schmid, 2006 
 * https://graphics.ethz.ch/teaching/former/imagesynthesis_06/miniprojects/p3 */

#include "fluid.h"

Fluid::Fluid()
{
	diffusion = 0.00001f;
	viscosity = 0.000000f;
	buoyancy  = 4.0f;
	vc_eps    = 5.0f;
	_dt = DT;

	for (int i=0; i<10; i++)
		ClearBuffer(_buffers[i]);

	int i=0;
	_density=_buffers[i++]; _densityTmp=_buffers[i++];
	_velX=_buffers[i++]; _velXTmp=_buffers[i++];
	_velY=_buffers[i++]; _velYTmp=_buffers[i++];
	_velZ=_buffers[i++]; _velZTmp=_buffers[i++];

	ClearSources();

	_renderer = new Renderer();
}

Fluid::~Fluid()
{
}

#define BOUNDARY
#undef BOUNDARY
void Fluid::EnforceBoundary(int b, float* quantity)
{
#ifdef BOUNDARY
	int i, j;
	for (i=1; i<=N; i++)
	{
		for (j=1; j<=N; j++) {
			quantity[_I(0,i,j)]    = (b==1) ? -quantity[_I(1,i,j)] : quantity[_I(1,i,j)];
			quantity[_I(N+1,i,j)]  = (b==1) ? -quantity[_I(N,i,j)] : quantity[_I(N,i,j)];
			quantity[_I(i,0,j)]    = (b==2) ? -quantity[_I(i,1,j)] : quantity[_I(i,1,j)];
			quantity[_I(i,N+1,j)]  = (b==2) ? -quantity[_I(i,N,j)] : quantity[_I(i,N,j)];
			quantity[_I(i,j,0)]    = (b==3) ? -quantity[_I(i,j,1)] : quantity[_I(i,j,1)];
			quantity[_I(i,j,N+1)]  = (b==3) ? -quantity[_I(i,j,N)] : quantity[_I(i,j,N)];
		}
	}

	quantity[_I(0,0,0)]       = (quantity[_I(1,0,0)]    + quantity[_I(0,1,0)]     + quantity[_I(0,0,1)])    / 3;
	quantity[_I(0,N+1,0)]     = (quantity[_I(1,N+1,0)]  + quantity[_I(0,N,0)]     + quantity[_I(0,N+1,1)])  / 3;
	quantity[_I(N+1,0,0)]     = (quantity[_I(N,0,0)]    + quantity[_I(N+1,1,0)]   + quantity[_I(N+1,0,1)])  / 3;
	quantity[_I(N+1,N+1,0)]   = (quantity[_I(N,N+1,0)]  + quantity[_I(N+1,N,0)]   + quantity[_I(N+1,N+1,1)])/ 3;
	quantity[_I(0,0,N+1)]     = (quantity[_I(1,0,N+1)]  + quantity[_I(0,1,N+1)]   + quantity[_I(0,0,N)])    / 3;
	quantity[_I(0,N+1,N+1)]   = (quantity[_I(1,N+1,N+1)]+ quantity[_I(0,N,N+1)]   + quantity[_I(0,N+1,N)])  / 3;
	quantity[_I(N+1,0,N+1)]   = (quantity[_I(N,0,N+1)]  + quantity[_I(N+1,1,N+1)] + quantity[_I(N+1,0,N)])  / 3;
	quantity[_I(N+1,N+1,N+1)] = (quantity[_I(N,N+1,N+1)]+ quantity[_I(N+1,N,N+1)] + quantity[_I(N+1,N+1,N)])/ 3;
#endif
}

void Fluid::AddSource(float* src, float *dst)
{
	int i;

	for (i=0; i<SIZE; i++)
		dst[i] += src[i]*_dt;
}

void Fluid::AddBuoyancy()
{
	int i;

	for (i=0; i<SIZE; i++)
		_velY[i] += _density[i]*buoyancy*_dt;//FIXME
}

inline void Fluid::Diffuse(int b, float* velTmp, float* vel, float diff)
{
	int l;
	float a=_dt*diff*N*N*N;
	for (l=0; l<20; l++) 
	{
		FOR_ALL_CELL {
			vel[_I(i,j,k)] = (velTmp[_I(i,j,k)] + a*(
						vel[_I(i-1,j,k)]+vel[_I(i+1,j,k)]+
						vel[_I(i,j-1,k)]+vel[_I(i,j+1,k)]+
						vel[_I(i,j,k-1)]+vel[_I(i,j,k+1)]))/(1+6*a);
		} END_FOR
		EnforceBoundary(b,vel);
	}
}

inline void Fluid::Advect(int b, float* quantityTmp, float* quantity, float* velX, float* velY, float* velZ)
{
	int i0, j0, k0, i1, j1, k1;
	float sx0, sx1, sy0, sy1, sz0, sz1, v0, v1;
	float xx, yy, zz, dt0;
	dt0 = _dt*N;
	FOR_ALL_CELL {
		xx = i-dt0*velX[_I(i,j,k)];
		yy = j-dt0*velY[_I(i,j,k)];
		zz = k-dt0*velZ[_I(i,j,k)];
		if (xx<0.5) xx=0.5f; if (xx>N+0.5) xx=N+0.5f; i0=(int)xx; i1=i0+1;
		if (yy<0.5) yy=0.5f; if (yy>N+0.5) yy=N+0.5f; j0=(int)yy; j1=j0+1;
		if (zz<0.5) zz=0.5f; if (zz>N+0.5) zz=N+0.5f; k0=(int)zz; k1=k0+1;
		sx1 = xx-i0; sx0 = 1-sx1;
		sy1 = yy-j0; sy0 = 1-sy1;
		sz1 = zz-k0; sz0 = 1-sz1;
		v0 = sx0 * ( sy0*quantityTmp[_I(i0,j0,k0)] + sy1*quantityTmp[_I(i0,j1,k0)] ) +
			sx1 * ( sy0*quantityTmp[_I(i1,j0,k0)] + sy1*quantityTmp[_I(i1,j1,k0)] );
		v1 = sx0 * ( sy0*quantityTmp[_I(i0,j0,k1)] + sy1*quantityTmp[_I(i0,j1,k1)] ) +
			sx1 * ( sy0*quantityTmp[_I(i1,j0,k1)] + sy1*quantityTmp[_I(i1,j1,k1)] );
		quantity[_I(i,j,k)] = sz0*v0 + sz1*v1;
	} END_FOR
	EnforceBoundary(b,_density);
}

void Fluid::Project(void)
{
	float* p = _velXTmp;	
	float* div = _velYTmp;	// temporary buffers, use old velocity buffers
	int l;
	float h;
	h = 1.0f/N;
	FOR_ALL_CELL {
				div[_I(i,j,k)] = -h*(
					_velX[_I(i+1,j,k)]-_velX[_I(i-1,j,k)]+
					_velY[_I(i,j+1,k)]-_velY[_I(i,j-1,k)]+
					_velZ[_I(i,j,k+1)]-_velZ[_I(i,j,k-1)])/3;
				p[_I(i,j,k)] = 0;
	} END_FOR
	EnforceBoundary(0,div); EnforceBoundary(0,p);
	for (l=0; l<20; l++) 
	{
		FOR_ALL_CELL {
					p[_I(i,j,k)] = (div[_I(i,j,k)]+
						p[_I(i-1,j,k)]+p[_I(i+1,j,k)]+
						p[_I(i,j-1,k)]+p[_I(i,j+1,k)]+
						p[_I(i,j,k-1)]+p[_I(i,j,k+1)])/6;
		} END_FOR
		EnforceBoundary(0,p);
	}
	FOR_ALL_CELL {
				_velX[_I(i,j,k)] -= (p[_I(i+1,j,k)]-p[_I(i-1,j,k)])/3/h;
				_velY[_I(i,j,k)] -= (p[_I(i,j+1,k)]-p[_I(i,j-1,k)])/3/h;
				_velZ[_I(i,j,k)] -= (p[_I(i,j,k+1)]-p[_I(i,j,k-1)])/3/h;
	} END_FOR
	EnforceBoundary(1,_velX); EnforceBoundary(2,_velY);
}

void Fluid::VorticityConfinement()
{
	int ijk;
	//temp buffers
	float *curlX = _velXTmp, *curlY = _velYTmp, *curlZ=_velZTmp, *curl=_densityTmp;
	float dt0 = _dt * vc_eps;

	FOR_ALL_CELL {
				ijk = _I(i,j,k);
				// curlx = dw/dy - dv/dz
				curlX[ijk] = (_velZ[_I(i,j+1,k)] - _velZ[_I(i,j-1,k)]) * 0.5f -
					(_velY[_I(i,j,k+1)] - _velY[_I(i,j,k-1)]) * 0.5f;

				// curly = du/dz - dw/dx
				curlY[ijk] = (_velX[_I(i,j,k+1)] - _velX[_I(i,j,k-1)]) * 0.5f -
					(_velZ[_I(i+1,j,k)] - _velZ[_I(i-1,j,k)]) * 0.5f;

				// curlz = dv/dx - du/dy
				curlZ[ijk] = (_velY[_I(i+1,j,k)] - _velY[_I(i-1,j,k)]) * 0.5f -
					(_velX[_I(i,j+1,k)] - _velX[_I(i,j-1,k)]) * 0.5f;

				// curl = |curl|
				curl[ijk] = sqrtf(curlX[ijk]*curlX[ijk] +
						curlY[ijk]*curlY[ijk] +
						curlZ[ijk]*curlZ[ijk]);
	} END_FOR

	FOR_ALL_CELL {
				ijk = _I(i,j,k);
				float nX = (curl[_I(i+1,j,k)] - curl[_I(i-1,j,k)]) * 0.5f;
				float nY = (curl[_I(i,j+1,k)] - curl[_I(i,j-1,k)]) * 0.5f;
				float nZ = (curl[_I(i,j,k+1)] - curl[_I(i,j,k-1)]) * 0.5f;
				float len1 = 1.0f/(sqrtf(nX*nX+nY*nY+nZ*nZ)+0.0000001f);
				nX *= len1;
				nY *= len1;
				nZ *= len1;
				_velX[ijk] += (nY*curlZ[ijk] - nZ*curlY[ijk]) * dt0;
				_velY[ijk] += (nZ*curlX[ijk] - nX*curlZ[ijk]) * dt0;
				_velZ[ijk] += (nX*curlY[ijk] - nY*curlX[ijk]) * dt0;
	} END_FOR
}

#define DIFFUSE
#define ADVECT

void Fluid::VelocityStep()
{
#if 1
	AddSource(su, _velX);
	AddSource(sv, _velY);
	AddSource(sw, _velZ);
#endif
	AddBuoyancy();
	VorticityConfinement();

#ifdef DIFFUSE
	std::swap(_velXTmp, _velX); std::swap(_velYTmp, _velY); std::swap(_velZTmp, _velZ);
	Diffuse(1, _velXTmp, _velX, viscosity);
	Diffuse(2, _velYTmp, _velY, viscosity);
	Diffuse(3, _velZTmp, _velZ, viscosity);
	Project();
#endif
#ifdef ADVECT
	std::swap(_velXTmp, _velX); std::swap(_velYTmp, _velY); std::swap(_velZTmp, _velZ);
	Advect(1, _velXTmp, _velX, _velXTmp, _velYTmp, _velZTmp);
	Advect(2, _velYTmp, _velY, _velXTmp, _velYTmp, _velZTmp);
	Advect(3, _velZTmp, _velZ, _velXTmp, _velYTmp, _velZTmp);
	Project();
#endif
}

void Fluid::DensityStep()
{
#if 1
	AddSource(sd, _density);
#endif
#ifdef DIFFUSE
	std::swap(_densityTmp, _density);
	Diffuse(0, _densityTmp, _density, diffusion);
#endif
#ifdef ADVECT
	std::swap(_densityTmp, _density);
	Advect(0, _densityTmp, _density, _velX, _velY, _velZ);

#if 1
	//decrease density
	for (int k=1; k<=N; k++) {
		for (int j=1; j<=N; j++) {
			for (int i=1; i<=N; i++) {
				_density[_I(k,j,i)] -= 0.001;
				if(_density[_I(k,j,i)] < 0)
					_density[_I(k,j,i)] = 0;
			}
		}
	}
#endif

#endif
}

void Fluid::_GenerateSmoke()
{
	const int centerY = RES/4;
	const int centerZ = RES/2;
	float dens = (rand()%1000)/1000.0f;
	for (int i = 1 ; i <= N ; i++) {
		for (int j = 1 ; j <= N ; j++) {
			if (hypot(i-centerY, j-centerZ) < RES/10) {
				this->_density[_I(5,i,j)] = dens;
				this->_velX[_I(5,i,j)] = 5.0f;
			}
		}
	}

}

void Fluid::SimulateStep()
{
	_GenerateSmoke();

	VelocityStep();
	DensityStep();

}


void Fluid::Show()
{
#if 1
	_renderer->FillTexture(this);
	_renderer->Render();
#else
	glBegin(GL_POINTS);
	FOR_ALL_CELL {
		if(!ALMOST_EQUAL(_density[_I(i,j,k)], 0) ) {
			glVertex3f(((float)i/N-0.5)*2, ((float)j/N-0.5)*2, ((float)k/N-0.5)*2 );
		}
	} END_FOR
	glEnd();
#endif
}


void Fluid::ClearBuffer(float* buf)
{
	for (int i=0; i<SIZE; i++) {
		buf[i] = 0.0f;
	}
}

void Fluid::ClearSources(void)
{
	for (int i=0; i<SIZE; i++) {
		sd[i] = su[i] = sv[i] = sw[i] = 0.0f;
	}
}

