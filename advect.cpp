/*
 *  advect.cpp
 *  smoke3D
 *
 */

#include "advect.h"
#include "utility.h"

static FLOAT spline_cubic(const FLOAT a[4], FLOAT x) {
	int i, j;
	FLOAT alpha[4], l[4], mu[4], z[4];
	FLOAT b[4], c[4], d[4];
	for(i = 1; i < 3; i++) {
		alpha[i] = 3.0 * (a[i+1] - a[i]) - 3.0 * (a[i] - a[i-1]);
	}
	l[0] = 1.0;
	mu[0] = 0.0;
	z[0] = 0.0;
	for(i = 1; i < 3; i++) {
		l[i] = 4.0 - mu[i-1];
		mu[i] = 1.0 / l[i];
		z[i] = (alpha[i] - z[i-1]) / l[i];
	}
	l[3] = 1.0;
	z[3] = 0.0;
	c[3] = 0.0;
	for(j = 2; 0 <= j; j--) {
		c[j] = z[j] - mu[j] * c[j+1];
		b[j] = a[j+1] - a[j] - (c[j+1] + 2.0 * c[j]) / 3.0;
		d[j] = (c[j+1] - c[j]) / 3.0;
	}
	
	FLOAT minv = min(a[1],a[2]);
	FLOAT maxv = max(a[2],a[1]);
	return min(maxv,max(minv,(a[1] + b[1] * x + c[1] * x * x + d[1] * x * x * x )));
}

static FLOAT spline( FLOAT ***d, int width, int height, int depth, FLOAT x, FLOAT y, FLOAT z ) {
	FLOAT f[16];
	FLOAT xn[4];
	FLOAT zn[4];
	
	x = max(0.0,min(width,x));
	y = max(0.0,min(height,y));
	z = max(0.0,min(depth,z));
	
	for( int k=0; k<4; k++ ) {
		for( int j=0; j<4; j++ ) {
			for( int i=0; i<4; i++ ) {
				int h = max(0,min(width-1,(int)x - 1 + i));
				int v = max(0,min(height-1,(int)y - 1 + j));
				int g = max(0,min(depth-1,(int)z - 1 + k));
				f[4*j+i] = d[h][v][g];
			}
		}
		
		for( int j=0; j<4; j++ ) {
			xn[j] = spline_cubic( &f[4*j], x - (int)x );
		}
		zn[k] = spline_cubic( xn, y - (int)y );
	}
	
	return spline_cubic( zn, z - (int)z );
}

FLOAT interp( FLOAT ***d, int width, int height, int depth, FLOAT x, FLOAT y, FLOAT z ) {
	return spline( d, width, height, depth, x, y, z );
}

static FLOAT u_ref( FLOAT ****u, int dir, int i, int j, int k, int N ) {
	if( dir == 0 )
		return u[0][max(0,min(N,i))][max(0,min(N-1,j))][max(0,min(N-1,k))];
	else if( dir == 1 )
		return u[1][max(0,min(N-1,i))][max(0,min(N,j))][max(0,min(N-1,k))];
	else
		return u[2][max(0,min(N-1,i))][max(0,min(N-1,j))][max(0,min(N,k))];
}
		   
void semiLagrangian( FLOAT ***d, FLOAT ***d0, int width, int height, int depth, FLOAT ****u, int N, FLOAT dt ) {
	OPENMP_FOR for( int n=0; n<width*height*depth; n++ ) {
		int i = (n%(width*height))%width;
		int j = (n%(width*height))/width;
		int k = n/(width*height);
		
		d[i][j][k] = interp( d0, width, height, depth, i-N*u[0][i][j][k]*dt, j-N*u[1][i][j][k]*dt, k-N*u[2][i][j][k]*dt );
	}
}

void advect::advect( FLOAT ****u, FLOAT ***c, int N, FLOAT dt ) {
	
	// Compute Fluid Velocity At Each Staggered Faces And Concentration Cell Centers
	static FLOAT ***ux[3] = { alloc3D(N+1,N,N), alloc3D(N+1,N,N), alloc3D(N+1,N,N) };
	static FLOAT ***uy[3] = { alloc3D(N,N+1,N), alloc3D(N,N+1,N), alloc3D(N,N+1,N) };
	static FLOAT ***uz[3] = { alloc3D(N,N,N+1), alloc3D(N,N,N+1), alloc3D(N,N,N+1) };
	static FLOAT ***uc[3] = { alloc3D(N,N,N), alloc3D(N,N,N), alloc3D(N,N,N) };
	static FLOAT ***out[4] = { alloc3D(N+1,N,N), alloc3D(N,N+1,N), alloc3D(N,N,N+1), alloc3D(N,N,N) };
	
	FOR_EVERY_X_FLOW {
		ux[0][i][j][k] = u[0][i][j][k];
		ux[1][i][j][k] = (u_ref(u,1,i-1,j,k,N)+u_ref(u,1,i-1,j+1,k,N)+u_ref(u,1,i,j,k,N)+u_ref(u,1,i,j+1,k,N))/4.0;
		ux[2][i][j][k] = (u_ref(u,2,i-1,j,k,N)+u_ref(u,2,i-1,j,k+1,N)+u_ref(u,2,i,j,k,N)+u_ref(u,2,i,j,k+1,N))/4.0;
	} END_FOR
	
	FOR_EVERY_Y_FLOW {
		uy[0][i][j][k] = (u_ref(u,0,i,j-1,k,N)+u_ref(u,0,i+1,j-1,k,N)+u_ref(u,0,i,j,k,N)+u_ref(u,0,i+1,j,k,N))/4.0;
		uy[1][i][j][k] = u[1][i][j][k];
		uy[2][i][j][k] = (u_ref(u,2,i,j-1,k,N)+u_ref(u,2,i,j-1,k+1,N)+u_ref(u,2,i,j,k,N)+u_ref(u,2,i,j,k+1,N))/4.0;
	} END_FOR
	
	FOR_EVERY_Z_FLOW {
		uz[0][i][j][k] = (u_ref(u,0,i,j,k-1,N)+u_ref(u,0,i+1,j,k-1,N)+u_ref(u,0,i,j,k,N)+u_ref(u,0,i+1,j,k,N))/4.0;
		uz[1][i][j][k] = (u_ref(u,1,i,j,k-1,N)+u_ref(u,1,i,j+1,k-1,N)+u_ref(u,1,i,j,k,N)+u_ref(u,1,i,j+1,k,N))/4.0;
		uz[2][i][j][k] = u[2][i][j][k];
	} END_FOR
	
	FOR_EVERY_CELL {
		uc[0][i][j][k] = 0.5*u[0][i][j][k]+0.5*u[0][i+1][j][k];
		uc[1][i][j][k] = 0.5*u[1][i][j][k]+0.5*u[1][i][j+1][k];
		uc[2][i][j][k] = 0.5*u[2][i][j][k]+0.5*u[2][i][j][k+1];
	} END_FOR
	
	// BackTrace X Flow
	semiLagrangian( out[0], u[0], N+1, N, N, ux, N, dt );
	
	// BackTrace Y Flow
	semiLagrangian( out[1], u[1], N, N+1, N, uy, N, dt );
	
	// BackTrace Z Flow
	semiLagrangian( out[2], u[2], N, N, N+1, uz, N, dt );
	
	// BackTrace Concentration
	semiLagrangian( out[3], c, N, N, N, uc, N, dt );
	
	// Copy Back To The Result
	copy3D(u[0],out[0],N);
	copy3D(u[1],out[1],N);
	copy3D(u[2],out[2],N);
	copy3D(c,out[3],N);
}
