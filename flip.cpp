#include "core.h"
#include "flip.h"

Flip::Flip()
{
	Init();
}

Flip::~Flip()
{
	_DeleteMemory();
}

void Flip::_DeleteMemory()
{
	LOOP_FOR_PARTICLES(particles_) {
		delete p;
	} END_LOOP
	particles_.clear();

	Allocator<Cell> cellFreer;
	Allocator<GLdouble> velFreer;

	cellFreer.Free3D(grid_.cells_);
	velFreer.Free3D(grid_.velX_);
	velFreer.Free3D(grid_.velY_);
	velFreer.Free3D(grid_.velZ_);

	cellFreer.Free3D(savedGrid_.cells_);
	velFreer.Free3D(savedGrid_.velX_);
	velFreer.Free3D(savedGrid_.velY_);
	velFreer.Free3D(savedGrid_.velZ_);

	delete boundary;
}

void Flip::Reset()
{
	Object::Reset();
	/*FIXME: free memory and reinitialize*/
	_DeleteMemory();
	Init();
}

void Flip::ShowBoundary()
{
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	GLdouble vertices[][3] = { 
		{ -simualateAreaX_/2, -simualateAreaY_/2, -simualateAreaZ_/2 },
		{ simualateAreaX_/2, -simualateAreaY_/2, -simualateAreaZ_/2 },
		{ simualateAreaX_/2, simualateAreaY_/2, -simualateAreaZ_/2 },
		{ -simualateAreaX_/2, simualateAreaY_/2, -simualateAreaZ_/2 }, 
		{ -simualateAreaX_/2, -simualateAreaY_/2, simualateAreaZ_/2 },
		{ simualateAreaX_/2, -simualateAreaY_/2, simualateAreaZ_/2 }, 
		{ simualateAreaX_/2, simualateAreaY_/2, simualateAreaZ_/2 },
		{ -simualateAreaX_/2, simualateAreaY_/2, simualateAreaZ_/2 } };

	int faces[][4] = { { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 4, 5, 6, 7 },
		{ 0, 4, 7, 3 }, { 0, 1, 5, 4 }, { 0, 3, 2, 1 } };

	GLdouble colors[][3] = { { 0.0f, 1.0f, 1.0f }, { 1.0f, 0.0f, 1.0f },
		{ 1.0f, 1.0f, 0.0f }, { 0.0f, 0.5f, 0.5f },
		{ 0.5f, 0.0f, 0.5f }, { 0.5f, 0.5f, 0.0f } };

	glBegin(GL_QUADS);
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 4; j++)
			glVertex3dv(vertices[faces[i][j]]);
	}
	glEnd();
}

void Flip::Show()
{
	ShowBoundary();
	
	glBegin(GL_POINTS);
	
	for (int i = 0 ; i < particles_.size() ; i++) {
		GLdouble x = particles_[i]->position_.x_ - simualateAreaX_ / 2;
		GLdouble y = particles_[i]->position_.y_ - simualateAreaY_ / 2;
		GLdouble z = particles_[i]->position_.z_ - simualateAreaZ_ / 2;
		glVertex3d(x, y, z);
	}
	
	glEnd();
}

void Flip::_InsertToCell(Particle* particle)
{
	int i = (int)(particle->position_.x_ / grid_.cellSize_);
	int j = (int)(particle->position_.y_ / grid_.cellSize_);
	int k = (int)(particle->position_.z_ / grid_.cellSize_);

	grid_.cells_[i][j][k].particles_.push_back(particle);
}

void Flip::_GenerateFluidBlock(GLdouble xL, GLdouble xU, GLdouble yL, GLdouble yU, 
			GLdouble zL, GLdouble zU)
{
	GLdouble particleIntervalX = simualateAreaX_ / ((PARTICLE_PER_CELL)*grid_.nX_);
	GLdouble particleIntervalY = simualateAreaY_ / ((PARTICLE_PER_CELL)*grid_.nY_);
	GLdouble particleIntervalZ = simualateAreaZ_ / ((PARTICLE_PER_CELL)*grid_.nZ_);
	GLdouble startX = 1.25 * particleIntervalX;
	GLdouble startY = 1.25 * particleIntervalY;
	GLdouble startZ = 1.25 * particleIntervalZ;

	for (GLdouble i = startX ; i <= simualateAreaX_-startX;
			i += particleIntervalX) {
		for (GLdouble j = startY ; j <= simualateAreaY_-startY; 
				j += particleIntervalY) {
			for (GLdouble k = startZ ; k <= simualateAreaZ_-startZ;
					k += particleIntervalZ) {
				if (i > xL && i < xU && j > yL && j < yU && k > zL && k < zU) {
					Particle* particle = new Particle;
					particle->position_.x_ = i;
					particle->position_.y_ = j;
					particle->position_.z_ = k;
					particles_.push_back(particle);
#if 0
					std::cout << i << " " << j << " " << k << std::endl;
#endif
					_InsertToCell(particle);
				}
			}
		}
	}

	_MarkFluid();
}


void Flip::_MarkFluid()
{
	LOOP_FOR_CELLS {
		if(grid_.cells_[i][j][k].type_ != SOLID) {
			grid_.cells_[i][j][k].type_ = AIR;
		}
		if(!grid_.cells_[i][j][k].particles_.empty() &&
				grid_.cells_[i][j][k].type_ != SOLID) {
			grid_.cells_[i][j][k].type_ = FLUID;
		}
	}END_LOOP
}

void Flip::_InitGrid()
{
	grid_.nX_ = grid_.nY_ = grid_.nZ_ = N;

	Allocator<Cell> cellAlloc;
	Allocator<GLdouble> velAlloc;

	grid_.cells_ = cellAlloc.Alloc3D(grid_.nX_, grid_.nY_, grid_.nZ_);
	grid_.velX_ = velAlloc.Alloc3D(grid_.nX_+1, grid_.nY_, grid_.nZ_);
	grid_.velY_ = velAlloc.Alloc3D(grid_.nX_, grid_.nY_+1, grid_.nZ_);
	grid_.velZ_ = velAlloc.Alloc3D(grid_.nX_, grid_.nY_, grid_.nZ_+1);
	grid_.cellSize_ = simualateAreaX_ / grid_.nX_;

	savedGrid_.cells_ = cellAlloc.Alloc3D(grid_.nX_, grid_.nY_, grid_.nZ_);
	savedGrid_.velX_ = velAlloc.Alloc3D(grid_.nX_+1, grid_.nY_, grid_.nZ_);
	savedGrid_.velY_ = velAlloc.Alloc3D(grid_.nX_, grid_.nY_+1, grid_.nZ_);
	savedGrid_.velZ_ = velAlloc.Alloc3D(grid_.nX_, grid_.nY_, grid_.nZ_+1);
	savedGrid_.cellSize_ = simualateAreaX_ / grid_.nX_;
}

void Flip::_InitBoundary()
{
	//set all cell as air
	LOOP_FOR_CELLS {
		grid_.cells_[i][j][k].type_ = AIR;
	} END_LOOP

	GLdouble rate = 0.01;
	//init wall and wall normals
	boundary = new GridObject;
	Allocator<Normal> normAlloc;
	boundary->normals = normAlloc.Alloc3D(grid_.nX_, grid_.nY_, grid_.nZ_);
	for (int i = 0 ; i < grid_.nX_ ; i++) {
		for (int j = 0 ; j < grid_.nY_ ; j++) {
			grid_.cells_[i][j][0].type_ = SOLID;
			grid_.cells_[i][j][grid_.nZ_-1].type_ = SOLID;

			boundary->normals[i][j][0].x_ += 0;
			boundary->normals[i][j][0].y_ += 0;
			boundary->normals[i][j][0].z_ += rate * grid_.cellSize_;

			boundary->normals[i][j][grid_.nZ_-1].x_ += 0;
			boundary->normals[i][j][grid_.nZ_-1].y_ += 0;
			boundary->normals[i][j][grid_.nZ_-1].z_ += -rate * grid_.cellSize_;
		}
	}

	for (int j = 0 ; j < grid_.nX_ ; j++) {
		for (int k = 0 ; k < grid_.nY_ ; k++) {
			grid_.cells_[0][j][k].type_ = SOLID;
			grid_.cells_[grid_.nX_-1][j][k].type_ = SOLID;

			boundary->normals[0][j][k].x_ += rate * grid_.cellSize_;
			boundary->normals[0][j][k].y_ += 0;
			boundary->normals[0][j][k].z_ += 0;

			boundary->normals[grid_.nX_-1][j][k].x_ += -rate * grid_.cellSize_;
			boundary->normals[grid_.nX_-1][j][k].y_ += 0;
			boundary->normals[grid_.nX_-1][j][k].z_ += 0;
		}
	}

	for (int i = 0 ; i < grid_.nX_ ; i++) {
		for (int k = 0 ; k < grid_.nY_ ; k++) {
			grid_.cells_[i][0][k].type_ = SOLID;
			grid_.cells_[i][grid_.nY_-1][k].type_ = SOLID;

			boundary->normals[i][0][k].x_ += 0;
			boundary->normals[i][0][k].y_ += rate * grid_.cellSize_;
			boundary->normals[i][0][k].z_ += 0;

			boundary->normals[i][grid_.nY_-1][k].x_ += 0;
			boundary->normals[i][grid_.nY_-1][k].y_ += -rate * grid_.cellSize_;
			boundary->normals[i][grid_.nY_-1][k].z_ += 0;
		}
	}
}

void Flip::_InitDensity()
{
	GLfloat h = DENSITY/grid_.nX_;
	LOOP_FOR_INNER_CELLS {
		Particle *p = new Particle;
		p->position_.x_ = (i+0.5)*h;
		p->position_.y_ = (j+0.5)*h;
		p->position_.z_ = (k+0.5)*h;
		p->mass_ = 1.0;
		particles_.push_back(p);
	} END_LOOP

	MatchParticlesToCell();
	maxDensity_ = 1.0;
	ComputeDensity();
	maxDensity_ = 0.0;
	LOOP_FOR_PARTICLES(particles_){
		maxDensity_ = std::max(maxDensity_,p->density_);
		delete p;
	} END_LOOP
	particles_.clear();

	LOOP_FOR_CELLS {
		grid_.cells_[i][j][k].particles_.clear();
	} END_LOOP
}

void Flip::Init()
{
	step_ = 0;
	simualateAreaX_ = simualateAreaY_ = simualateAreaZ_ = 1.0;

	//allocate memory
	_InitGrid();

	//init boundary and mark solid
	_InitBoundary();

	//_InitDensity();

	/*several test cases*/
#if 0
	_GenerateFluidBlock(0.2, 0.4, 0.06, 0.4, 0.2, 0.8);
	_GenerateFluidBlock(WALL_THICKNESS, 1.0-WALL_THICKNESS, WALL_THICKNESS, 
			0.06, WALL_THICKNESS, 1.0-WALL_THICKNESS);
#endif
#if 0
	_GenerateFluidBlock(WALL_THICKNESS, 1.0-WALL_THICKNESS, WALL_THICKNESS, 
			0.6, WALL_THICKNESS, 1.0-WALL_THICKNESS);
#endif
#if 1
	_GenerateFluidBlock(0.6, 0.8, 0.6, 
			0.8, 0.25, 0.75);
	_GenerateFluidBlock(WALL_THICKNESS, 1.0-WALL_THICKNESS, WALL_THICKNESS, 
			0.1, WALL_THICKNESS, 1.0-WALL_THICKNESS);
#endif
#if 0
	_GenerateFluidBlock(0.5, 0.75, 0.5, 
			0.75, 0.25, 0.5);
#endif
#if 0
	_GenerateFluidBlock(0.6, 0.8, 0.6, 
			0.8, 0.25, 0.75);
	_GenerateFluidBlock(0.2, 0.4, 0.2, 
			0.3, 0.25, 0.75);
#endif
}

void Flip::EnforceBoundary()
{
	LOOP_FOR_CELLS {
		if(grid_.cells_[i][j][k].type_ == SOLID) {
			grid_.velX_[i][j][k] = usolid;
			grid_.velX_[i+1][j][k] = usolid;
			grid_.velY_[i][j][k] = usolid;
			grid_.velY_[i][j+1][k] = usolid;
			grid_.velZ_[i][j][k] = usolid;
			grid_.velZ_[i][j][k+1] = usolid;
		}
	} END_LOOP
}

void Flip::Project()
{
	GLdouble dx = grid_.cellSize_;
	GLdouble div = 0;


	//compute right-hand side divergence
	int index;	//linear index
	int rhsLength = MAT2LINEAR(grid_.nX_-1,grid_.nY_-1,grid_.nZ_-1) + 1;
	Eigen::VectorXd rightHand(rhsLength);
	rightHand.setZero();
	LOOP_FOR_CELLS {
		if (grid_.cells_[i][j][k].type_ == FLUID) {
			index = MAT2LINEAR(i,j,k);
			div = - ( (grid_.velX_[i+1][j][k] - grid_.velX_[i][j][k])
				+ (grid_.velY_[i][j+1][k] - grid_.velY_[i][j][k]) 
				+ (grid_.velZ_[i][j][k+1] - grid_.velZ_[i][j][k]) ) / dx;
			rightHand(index) = div;
		}
	} END_LOOP
	
	/*FIXME: modify divergence to account for solid velocities*/
	LOOP_FOR_CELLS {
		index = MAT2LINEAR(i,j,k);
		if (grid_.cells_[i][j][k].type_ == FLUID) {
			if(grid_.cells_[i-1][j][k].type_ == SOLID) {
				rightHand(index) -= (grid_.velX_[i][j][k] - usolid) / dx;
			}
			if(grid_.cells_[i+1][j][k].type_ == SOLID) {
				rightHand(index) += (grid_.velX_[i+1][j][k] - usolid) / dx;
			}

			if(grid_.cells_[i][j-1][k].type_ == SOLID) {
				rightHand(index) -= (grid_.velY_[i][j][k] - usolid) / dx;
			}
			if(grid_.cells_[i][j+1][k].type_ == SOLID) {
				rightHand(index) += (grid_.velY_[i][j+1][k] - usolid) / dx;
			}

			if(grid_.cells_[i][j][k-1].type_ == SOLID) {
				rightHand(index) -= (grid_.velZ_[i][j][k] - usolid) / dx;
			}
			if(grid_.cells_[i][j][k+1].type_ == SOLID) {
				rightHand(index) += (grid_.velZ_[i][j][k+1] - usolid) / dx;
			}
		}
	} END_LOOP


	//set coefficient matrix using sparse matrix
	std::vector<DoubleTriplet> coefficients;
	DoubleSpaMat coeffMat(rhsLength, rhsLength);
	coeffMat.setZero();
	int index0, indexIPlus1, indexJPlus1, indexKPlus1,
		indexIMinus1, indexJMinus1, indexKMinus1;	//linear index
	GLdouble diagIJK, iPlus1, jPlus1, kPlus1, 
			 iMinus1, jMinus1, kMinus1;	//value
	GLdouble scale = DT / (RHO*dx*dx);

	LOOP_FOR_CELLS {

		index0 = MAT2LINEAR(i,j,k);
		indexIPlus1 = MAT2LINEAR(i+1,j,k);
		indexJPlus1 = MAT2LINEAR(i,j+1,k);
		indexKPlus1 = MAT2LINEAR(i,j,k+1);
		indexIMinus1 = MAT2LINEAR(i-1,j,k);
		indexJMinus1 = MAT2LINEAR(i,j-1,k);
		indexKMinus1 = MAT2LINEAR(i,j,k-1);
		diagIJK = iPlus1 = jPlus1 = kPlus1 
			= iMinus1 = jMinus1 = kMinus1 = 0.0;

		if(grid_.cells_[i][j][k].type_ != FLUID)
			continue;

		if(grid_.cells_[i][j][k].type_ == FLUID && 
				grid_.cells_[i+1][j][k].type_ == FLUID) {
			diagIJK += scale;
			iPlus1 = -scale;
		}
		else if (grid_.cells_[i][j][k].type_ == FLUID &&
				grid_.cells_[i+1][j][k].type_ == AIR) {
			diagIJK += scale;
		}


		if(grid_.cells_[i][j][k].type_ == FLUID && 
				grid_.cells_[i][j+1][k].type_ == FLUID) {
			diagIJK += scale;
			jPlus1 = -scale;
		}
		else if (grid_.cells_[i][j][k].type_ == FLUID &&
				grid_.cells_[i][j+1][k].type_ == AIR) {
			diagIJK += scale;
		}


		if(grid_.cells_[i][j][k].type_ == FLUID && 
				grid_.cells_[i][j][k+1].type_ == FLUID) {
			diagIJK += scale;
			kPlus1 = -scale;
		}
		else if (grid_.cells_[i][j][k].type_ == FLUID &&
				grid_.cells_[i][j][k+1].type_ == AIR) {
			diagIJK += scale;
		}


		if(grid_.cells_[i][j][k].type_ == FLUID && 
				grid_.cells_[i-1][j][k].type_ == FLUID) {
			diagIJK += scale;
			iMinus1 = -scale;
		}
		else if (grid_.cells_[i][j][k].type_ == FLUID &&
				grid_.cells_[i-1][j][k].type_ == AIR) {
			diagIJK += scale;
		}

		
		if(grid_.cells_[i][j][k].type_ == FLUID && 
				grid_.cells_[i][j-1][k].type_ == FLUID) {
			diagIJK += scale;
			jMinus1 = -scale;
		}
		else if (grid_.cells_[i][j][k].type_ == FLUID &&
				grid_.cells_[i][j-1][k].type_ == AIR) {
			diagIJK += scale;
		}


		if(grid_.cells_[i][j][k].type_ == FLUID && 
				grid_.cells_[i][j][k-1].type_ == FLUID) {
			diagIJK += scale;
			kMinus1 = -scale;
		}
		else if (grid_.cells_[i][j][k].type_ == FLUID &&
				grid_.cells_[i][j][k-1].type_ == AIR) {
			diagIJK += scale;
		}

		coefficients.push_back(DoubleTriplet(index0, index0, diagIJK) );
		coefficients.push_back(DoubleTriplet(index0, indexIPlus1, iPlus1) );
		coefficients.push_back(DoubleTriplet(index0, indexJPlus1, jPlus1) );
		coefficients.push_back(DoubleTriplet(index0, indexKPlus1, kPlus1) );
		coefficients.push_back(DoubleTriplet(index0, indexIMinus1, iMinus1) );
		coefficients.push_back(DoubleTriplet(index0, indexJMinus1, jMinus1) );
		coefficients.push_back(DoubleTriplet(index0, indexKMinus1, kMinus1) );
	} END_LOOP
	
	//set coefficient matrix
	coeffMat.setFromTriplets(coefficients.begin(), coefficients.end() );
	
	Eigen::VectorXd pressure(rhsLength);
	Eigen::ConjugateGradient<DoubleSpaMat> solver;
	solver.compute(coeffMat);
	if(solver.info() != Eigen::Success) {
		std::cerr << "ConjugateGradient solver error!" << std::endl;
		exit(-1);
	}
	pressure = solver.solve(rightHand);

	//save pressure to grid
	for(int l = 0 ; l < pressure.size() ; l++) {
		//change linear index to 3D index
		int i = l / (grid_.nY_*grid_.nZ_);
		int j = (l - i * grid_.nY_ * grid_.nZ_) / grid_.nZ_;
		int k = (l - i * grid_.nY_ * grid_.nZ_) % grid_.nZ_;

		grid_.cells_[i][j][k].press_ = pressure[l];
	}
	
	
	//update velocities
	//FIXME:set boundary vel to zero
	scale = DT/(RHO*dx);

	//update x vel;
	LOOP_FOR_XVELS {
		if(grid_.cells_[i][j][k].type_ == FLUID) {
			grid_.velX_[i][j][k] -= scale*grid_.cells_[i][j][k].press_;
			grid_.velX_[i+1][j][k] += scale*grid_.cells_[i][j][k].press_;
		}
	} END_LOOP

	//update y vel;
	LOOP_FOR_YVELS {
		if(grid_.cells_[i][j][k].type_ == FLUID) {
			grid_.velY_[i][j][k] -= scale*grid_.cells_[i][j][k].press_;
			grid_.velY_[i][j+1][k] += scale*grid_.cells_[i][j][k].press_;
		}
	} END_LOOP

	//update z vel;
	LOOP_FOR_ZVELS {
		if(grid_.cells_[i][j][k].type_ == FLUID) {
			grid_.velZ_[i][j][k] -= scale*grid_.cells_[i][j][k].press_;
			grid_.velZ_[i][j][k+1] += scale*grid_.cells_[i][j][k].press_;
		}
	} END_LOOP

}

void Flip::AddExtForce()
{
	LOOP_FOR_PARTICLES(particles_) {
		p->velocity_.y_ += -DT*GRAVITY;
	} END_LOOP
}

void Flip::SaveGridVel()
{
	LOOP_FOR_XVELS {
		savedGrid_.velX_[i][j][k] = grid_.velX_[i][j][k];
	} END_LOOP

	LOOP_FOR_YVELS {
		savedGrid_.velY_[i][j][k] = grid_.velY_[i][j][k];
	} END_LOOP

	LOOP_FOR_ZVELS {
		savedGrid_.velZ_[i][j][k] = grid_.velZ_[i][j][k];
	} END_LOOP
}

void Flip::ExtrapolateVelocity()
{
}

void Flip::SolvePICFLIP()
{
	Allocator<GLdouble> velAlloc;
	GLdouble*** deltaVelX = velAlloc.Alloc3D(grid_.nX_+1, grid_.nY_, grid_.nZ_);
	GLdouble*** deltaVelY = velAlloc.Alloc3D(grid_.nX_, grid_.nY_+1, grid_.nZ_);
	GLdouble*** deltaVelZ = velAlloc.Alloc3D(grid_.nX_, grid_.nY_, grid_.nZ_+1);

	LOOP_FOR_XVELS {
		deltaVelX[i][j][k] = grid_.velX_[i][j][k] - savedGrid_.velX_[i][j][k];
	} END_LOOP

	LOOP_FOR_YVELS {
		deltaVelY[i][j][k] = grid_.velY_[i][j][k] - savedGrid_.velY_[i][j][k];
	} END_LOOP

	LOOP_FOR_ZVELS {
		deltaVelZ[i][j][k] = grid_.velZ_[i][j][k] - savedGrid_.velZ_[i][j][k];
	} END_LOOP

	//trilinear interpolation
	GLdouble tx, ty, tz;
	LOOP_FOR_PARTICLES(particles_) {
		int i = (int)(p->position_.x_ / grid_.cellSize_);
		int j = (int)(p->position_.y_ / grid_.cellSize_);
		int k = (int)(p->position_.z_ / grid_.cellSize_);

		tx = (p->position_.x_ - i*grid_.cellSize_) / grid_.cellSize_;
		ty = (p->position_.y_ - j*grid_.cellSize_) / grid_.cellSize_;
		tz = (p->position_.z_ - k*grid_.cellSize_) / grid_.cellSize_;

		Velocity pic, flip;
		pic.x_ = (tx*grid_.velX_[i+1][j][k] + (1-tx)*grid_.velX_[i][j][k]);
		pic.y_ = (ty*grid_.velY_[i][j+1][k] + (1-ty)*grid_.velY_[i][j][k]);
		pic.z_ = (tz*grid_.velZ_[i][j][k+1] + (1-tz)*grid_.velZ_[i][j][k]);
		flip.x_ = p->velocity_.x_ + (tx*deltaVelX[i+1][j][k] + (1-tx)*deltaVelX[i][j][k]);
		flip.y_ = p->velocity_.y_ + (ty*deltaVelY[i][j+1][k] + (1-ty)*deltaVelY[i][j][k]);
		flip.z_ = p->velocity_.z_ + (tz*deltaVelZ[i][j][k+1] + (1-tz)*deltaVelZ[i][j][k]);
		p->velocity_.x_ = ALPHA * flip.x_ + (1-ALPHA) * pic.x_;
		p->velocity_.y_ = ALPHA * flip.y_ + (1-ALPHA) * pic.y_;
		p->velocity_.z_ = ALPHA * flip.z_ + (1-ALPHA) * pic.z_;
	} END_LOOP
}

void Flip::ComputeDensity()
{

	LOOP_FOR_PARTICLES(particles_) {
		
		int i = (int)(p->position_.x_ / grid_.cellSize_);
		int j = (int)(p->position_.y_ / grid_.cellSize_);
		int k = (int)(p->position_.z_ / grid_.cellSize_);

		// Find Neighbors
		ParticleList neighbors = _GetNeighborParticles(i,j,k,1,1,1);
		GLdouble wsum = 0.0;
		for( int m = 0; m < neighbors.size(); m++ ) {
			Particle* np = neighbors[m];
			Eigen::Vector3d subVec(p->position_.x_-np->position_.x_, 
					p->position_.y_-np->position_.y_,
					p->position_.z_-np->position_.z_);
			GLdouble w = np->mass_*Kernel::SmoothKernel(subVec.norm(), 4.0*DENSITY/N);
			wsum += w;
		}
		p->density_ = wsum / maxDensity_;
	} END_LOOP
}

/*FIXME: Loop for particles is better?*/
void Flip::TransferParticleVelToCell()
{
	int i, j, k;
	GLdouble weightedVel;
	GLdouble weight, totalWeight;
	ParticleList pl;
	GLdouble pX, pY, pZ;
	GLdouble cellCenterX, cellCenterY, cellCenterZ;
	
	//interpolate x velocity to grid
	LOOP_FOR_XVELS {
		pl = _GetNeighborParticles(i, j, k, 1, 2, 2);
		weight = totalWeight = 0;
		weightedVel = 0;
		LOOP_FOR_PARTICLES(pl) {
			cellCenterX = (i)*grid_.cellSize_;
			cellCenterY = (j+0.5)*grid_.cellSize_;
			cellCenterZ = (k+0.5)*grid_.cellSize_;
			pX = p->position_.x_;
			pY = p->position_.y_;
			pZ = p->position_.z_;
			Eigen::Vector3d subVec(pX-cellCenterX, 
					pY-cellCenterY,	pZ-cellCenterZ);
			weight = Kernel::SharpKernel(subVec.norm(), KERNEL_SHARP);
			weightedVel += p->velocity_.x_ * weight;
			totalWeight += weight;
		} END_LOOP
		if(weight != 0) {
			grid_.velX_[i][j][k] = weightedVel / totalWeight;
		}
		else
			grid_.velX_[i][j][k] = 0;
	}END_LOOP

	//interpolate y velocity to grid
	LOOP_FOR_YVELS {
		pl = _GetNeighborParticles(i, j, k, 2, 1, 2);
		weight = totalWeight = 0;
		weightedVel = 0;
		LOOP_FOR_PARTICLES(pl) {
			cellCenterX = (i+0.5)*grid_.cellSize_;
			cellCenterY = (j)*grid_.cellSize_;
			cellCenterZ = (k+0.5)*grid_.cellSize_;
			pX = p->position_.x_;
			pY = p->position_.y_;
			pZ = p->position_.z_;
			Eigen::Vector3d subVec(pX-cellCenterX, 
					pY-cellCenterY,	pZ-cellCenterZ);
			weight = Kernel::SharpKernel(subVec.norm(), KERNEL_SHARP);
			weightedVel += p->velocity_.y_ * weight;
			totalWeight += weight;
		} END_LOOP
		if(weight != 0) {
			grid_.velY_[i][j][k] = weightedVel / totalWeight;
		}
		else
			grid_.velY_[i][j][k] = 0;
	}END_LOOP

		//interpolate z velocity to grid
	LOOP_FOR_ZVELS {
		pl = _GetNeighborParticles(i, j, k, 2, 2, 1);
		weight = totalWeight = 0;
		weightedVel = 0;
		LOOP_FOR_PARTICLES(pl) {
			cellCenterX = (i+0.5)*grid_.cellSize_;
			cellCenterY = (j+0.5)*grid_.cellSize_;
			cellCenterZ = (k)*grid_.cellSize_;
			pX = p->position_.x_;
			pY = p->position_.y_;
			pZ = p->position_.z_;
			Eigen::Vector3d subVec(pX-cellCenterX, 
					pY-cellCenterY,	pZ-cellCenterZ);
			weight = Kernel::SharpKernel(subVec.norm(), KERNEL_SHARP);
			weightedVel += p->velocity_.z_ * weight;
			totalWeight += weight;
		} END_LOOP
		if(weight != 0)
			grid_.velZ_[i][j][k] = weightedVel / totalWeight;
		else
			grid_.velZ_[i][j][k] = 0;

	} END_LOOP

}

ParticleList Flip::_GetNeighborParticles(int i, int j, int k, int w, int h, int d)
{
	ParticleList res;
	for( int si=i-w; si<=i+w-1; si++ )
		for( int sj=j-h; sj<=j+h-1; sj++ ) 
			for( int sk=k-d; sk<=k+d-1; sk++ ) {
		if(si < 0 || si > grid_.nX_-1 || 
				sj < 0 || sj > grid_.nY_-1 || sk < 0 || sk > grid_.nZ_-1) 
			continue;
		ParticleList pl = grid_.cells_[si][sj][sk].particles_;
		for( int l=0; l<pl.size(); l++ ) { 
			Particle *p = pl[l];
			res.push_back(p);
		}
	}
	return res;
}


void Flip::MatchParticlesToCell()
{
	LOOP_FOR_CELLS {
		grid_.cells_[i][j][k].particles_.clear();
	} END_LOOP

	LOOP_FOR_PARTICLES(particles_) {
		_InsertToCell(p);
	} END_LOOP

	_MarkFluid();

}

void Flip::AdvectParticles()
{
	LOOP_FOR_PARTICLES(particles_) {
		p->position_.x_ += DT*p->velocity_.x_;
		p->position_.y_ += DT*p->velocity_.y_;
		p->position_.z_ += DT*p->velocity_.z_;
	} END_LOOP
}

void Flip::ParticleCollisionDetection()
{

#if 0
	LOOP_FOR_PARTICLES(particles_) {
		int i = (int)(p->position_.x_ / grid_.cellSize_);
		int j = (int)(p->position_.y_ / grid_.cellSize_);
		int k = (int)(p->position_.z_ / grid_.cellSize_);
		ParticleList neighbors = _GetNeighborParticles(i,j,k,1,1,1);
		// First, find the normalized vector n from the center of 
		// circle1 to the center of circle2
		LOOP_FOR_PARTICLES(neighbors) {
			Particle* np = neighbors[l];
			Eigen::Vector3d n(p->position_.x_-np->position_.x_,
					p->position_.y_-np->position_.y_,
					p->position_.z_-np->position_.z_);
#if 0
			std::cout << n.norm() << " " << p->density_ << std::endl;
#endif
			if (n.norm() >= p->density_/8) 
				continue;
			n.normalize();
			// Find the length of the component of each of the movement
			// vectors along n. 
			// a1 = v1 . n
			// a2 = v2 . n
			Eigen::Vector3d v1(p->velocity_.x_, p->velocity_.y_, p->velocity_.z_);
			Eigen::Vector3d v2(p->velocity_.x_, p->velocity_.y_, p->velocity_.z_);
			GLdouble a1 = v1.dot(n);
			GLdouble a2 = v2.dot(n);

			// Using the optimized version, 
			// optimizedP =  2(a1 - a2)
			//              -----------
			//                m1 + m2
			GLdouble optimizedP = (2.0 * (a1 - a2)) / (p->mass_ + np->mass_);

			// Calculate v1', the new movement vector of circle1
			// v1' = v1 - optimizedP * m2 * n
			v1 = v1 - optimizedP * np->mass_ * n;

			// Calculate v1', the new movement vector of circle1
			// v2' = v2 + optimizedP * m1 * n
			v2 = v2 + optimizedP * p->mass_ * n;

			p->velocity_.x_ = v1.x();
			p->velocity_.y_ = v1.y();
			p->velocity_.z_ = v1.z();
			
			np->velocity_.x_ = v2.x();
			np->velocity_.y_ = v2.y();
			np->velocity_.z_ = v2.z();
		}END_LOOP

	} END_LOOP
#endif

	/*FIXME: using more flexible collision detection*/
	LOOP_FOR_PARTICLES(particles_) {
		int i, j, k;
		p->position_.x_ = std::max(p->position_.x_,0.0);
		p->position_.y_ = std::max(p->position_.y_,0.0);
		p->position_.z_ = std::max(p->position_.z_,0.0);
		i = (int)(p->position_.x_ / grid_.cellSize_);
		j = (int)(p->position_.y_ / grid_.cellSize_);
		k = (int)(p->position_.z_ / grid_.cellSize_);
		i = (i >= grid_.nX_ ? grid_.nX_-1 : i);
		j = (j >= grid_.nY_ ? grid_.nY_-1 : j);
		k = (k >= grid_.nZ_ ? grid_.nZ_-1 : k);

		while(i < 1 || j < 1 || k < 1 || 
				i > grid_.nX_-2 || j > grid_.nY_-2 || k > grid_.nZ_-2) {
			if(grid_.cells_[i][j][k].type_ == SOLID) {
				p->position_.x_ += boundary->normals[i][j][k].x_;
				p->position_.y_ += boundary->normals[i][j][k].y_;
				p->position_.z_ += boundary->normals[i][j][k].z_;
				if (boundary->normals[i][j][k].x_ != 0)
					p->velocity_.x_ = 0;
				if (boundary->normals[i][j][k].y_ != 0)
					p->velocity_.y_ = 0;
				if (boundary->normals[i][j][k].z_ != 0)
					p->velocity_.z_ = 0;
			}
			p->position_.x_ = std::max(p->position_.x_,0.0);
			p->position_.y_ = std::max(p->position_.y_,0.0);
			p->position_.z_ = std::max(p->position_.z_,0.0);
			i = (int)(p->position_.x_ / grid_.cellSize_);
			j = (int)(p->position_.y_ / grid_.cellSize_);
			k = (int)(p->position_.z_ / grid_.cellSize_);
			i = (i >= grid_.nX_ ? grid_.nX_-1 : i);
			j = (j >= grid_.nY_ ? grid_.nY_-1 : j);
			k = (k >= grid_.nZ_ ? grid_.nZ_-1 : k);
		}
	} END_LOOP
}

void Flip::GenerateSurface()
{

}

void Flip::SimulateStep()
{
	//timestep++;
	step_++;

	//add externaal force, e.g. gravity
	AddExtForce();

	//match particles to corresponding cells
	MatchParticlesToCell();
	
	//ComputeDensity();
	
	//transfer particle velocities to cell
	TransferParticleVelToCell();

	SaveGridVel();

	//solve pressure, enforce boundaries, extrapolate velocities
	Project();
	EnforceBoundary();

	//solve pic, flip and transfer velocity to particles
    SolvePICFLIP();
	
	AdvectParticles();

	ParticleCollisionDetection();
	
	//genenrate surface mesh
	GenerateSurface();
}
