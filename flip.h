#ifndef FLIP_H
#define FLIP_H

#include "core.h"
#include "object.h"
#include "gridobject.h"
#include "simple_vector.h"
#include "kernel.h"
#include "allocator.h"
#include "file_operator.h"

class FileOperator;


typedef Eigen::Triplet<GLdouble> DoubleTriplet; 
typedef Eigen::SparseMatrix<GLdouble> DoubleSpaMat;
#define	N               32
#define MAX_STEP        600

#define VEL_FILE_PATH	"data/vel/"

#define ALPHA           1.0
#define DT              1e-2
#define DENSITY         0.5	//FIXME
#define RHO				1.0
#define	GRAVITY         9.8

#define	WALL_THICKNESS  (1.0/N)	//FIXME
#define usolid			0		//FIXME
#define KERNEL_SHARP	1.4

#define SOLID			0
#define FLUID			1
#define AIR				2

#define PARTICLE_PER_CELL 2	//particle no per cell in 1 dimension

#define LOOP_FOR_CELLS for(int i=0;i<grid_.nX_;i++)\
			for(int j=0;j<grid_.nY_;j++)\
			for(int k=0;k<grid_.nZ_;k++) {

#define LOOP_FOR_INNER_CELLS for(int i=1;i<grid_.nX_-1;i++)\
			for(int j=1;j<grid_.nY_-1;j++)\
			for(int k=1;k<grid_.nZ_-1;k++) {

#define LOOP_FOR_PARTICLES(pl) for(int l=0;l<pl.size();l++)\
					{Particle* p=pl[l];

#define LOOP_FOR_XVELS for(int i=1;i<grid_.nX_;i++)\
			for(int j=1;j<grid_.nY_-1;j++)\
			for(int k=1;k<grid_.nZ_-1;k++) {
#define LOOP_FOR_YVELS for(int i=1;i<grid_.nX_-1;i++)\
			for(int j=1;j<grid_.nY_;j++)\
			for(int k=1;k<grid_.nZ_-1;k++) {
#define LOOP_FOR_ZVELS for(int i=1;i<grid_.nX_-1;i++)\
			for(int j=1;j<grid_.nY_-1;j++)\
			for(int k=1;k<grid_.nZ_;k++) {

#define END_LOOP }

#define MAT2LINEAR(i,j,k) (i)*grid_.nY_*grid_.nZ_+(j)*grid_.nZ_+(k)

class Particle {
public: 
	GLdouble density_;	
	GLdouble mass_;
	Position position_;
	Velocity velocity_;
	int type_;
};


typedef std::vector<Particle* > ParticleList;

class Cell {
public:
	ParticleList particles_;
	GLdouble press_;
	int type_;
};

class Grid {
public:
	Cell*** cells_;
	int nX_, nY_, nZ_;

	GLdouble*** velX_;	//for MAC grid, vel is stored at the edge of cells
	GLdouble*** velY_;	//so it cannot be stored in Cell struct;
	GLdouble*** velZ_;

	GLdouble cellSize_;	//side length of one cell
};

class Flip: public Object {
public:
	Flip();
	virtual ~Flip();
	void Init();
	Grid* GetGrid();
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
	void ParticleCollisionDetection();
	void GenerateSurface();
	void ShowBoundary();

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
	void _InitDensity();
	void _DeleteMemory();
	void _WriteVelocity();

	ParticleList _GetNeighborParticles(int si, int sj, int sk, int w, int h, int d);

	int step_;

	GLdouble maxDensity_; //max density of particles

	ParticleList particles_;
	Grid grid_;
	Grid savedGrid_;
	
	GLdouble simualateAreaX_;
	GLdouble simualateAreaY_;
	GLdouble simualateAreaZ_;

	GridObject *boundary;

	FileOperator *velFileOp_;
};

#endif
