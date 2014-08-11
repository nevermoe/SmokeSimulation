#ifndef FLIP_H
#define FLIP_H

#include "core.h"

#define _DEBUG_ 1

typedef Eigen::Triplet<GLdouble> DoubleTriplet; 
typedef Eigen::SparseMatrix<GLdouble> DoubleSpaMat;
#define	N               32
#define MAX_STEP        600

#define ALPHA           0.95
#define DT              0.6e-2
//#define DENSITY         0.5
#define RHO				0.5
#define	GRAVITY         9.8

#define	WALL_THICKNESS  (1.0/N)

#define SOLID			0
#define FLUID			1
#define AIR				2

#define PARTICLE_PER_CELL 2	//particle no per cell in 1 dimension

#define LOOP_FOR_CELLS(nX, nY, nZ) for(int i=0;i<(nX);i++) for(int j=0;j<(nY);j++) for(int k=0;k<(nZ);k++) {
#define LOOP_FOR_INNER_CELLS(nX, nY, nZ) for(int i=1;i<(nX)-1;i++) for(int j=1;j<(nY-1);j++) for(int k=1;k<(nZ-1);k++) {
#define LOOP_FOR_PARTICLES(pl) for(int l=0;l<pl.size();l++){Particle* p=particles_[l];

#define LOOP_FOR_XVELS for(int i=1;i<grid_.nX_;i++) \
			for(int j=1;j<grid_.nY_-1;j++) \
			for(int k=1;k<grid_.nZ_-1;k++) {
#define LOOP_FOR_YVELS for(int i=1;i<grid_.nX_-1;i++) \
			for(int j=1;j<grid_.nY_;j++) \
			for(int k=1;k<grid_.nZ_-1;k++) {
#define LOOP_FOR_ZVELS for(int i=1;i<grid_.nX_-1;i++) \
			for(int j=1;j<grid_.nY_-1;j++) \
			for(int k=1;k<grid_.nZ_;k++) {

#define END_LOOP }

#define MAT2LINEAR(i,j,k) (i)*grid_.nX_*grid_.nY_+(j)*grid_.nY_+(k)

class Flip: public Object {
public:
	Flip();
	void Init();
	void ComputeDensity();
	void EnforceBoundary();
	void Project();
	void AddExtForce();
	void AdvectParticles();
	void SaveGridVel();
	void ExtrapolateVelocity();
	void SolvePICFLIP();
	void MatchParticlesToCell();
	void TransferParticleVelToCell();
	void MixPICFLIP();
	void ParticleCollisionDetection();
	void GenerateSurface();

	virtual void Reset();
	virtual void SimulateStep();
	virtual void Show();
private:
	void _GenerateFluidBlock(GLdouble xL, GLdouble xU, GLdouble yL, GLdouble yU, 
			GLdouble zL, GLdouble zU);
	void _MarkFluid();
	void _InsertToCell(Particle* particle);
	void _InitGrid();
	void _InitBoundary();
	ParticleList _GetNeighborParticles(int si, int sj, int sk, int w, int h, int d);
	GLdouble _ReciprocalKernal(GLdouble d);


	int step_;
	ParticleList particles_;
	Grid grid_;
	Grid savedGrid_;
	
	GLdouble simualateAreaX_;
	GLdouble simualateAreaY_;
	GLdouble simualateAreaZ_;

	GridObject *boundary;

};

#endif
