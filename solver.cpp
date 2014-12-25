/*
 *  solver.cpp
 *  smoke3D
 *
 */

#include "solver.h"
#include "utility.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define FOR_EVERY_COMP(N)	for( int i=0; i<N; i++ ) for( int j=0; j<N; j++ ) for( int k=0; k<N; k++ ) {

// Clamped Fetch
static FLOAT x_ref( FLOAT ***x, int fi, int fj, int fk, int i, int j, int k, FLOAT ***w, int n ) {
	i = min(max(0,i),n-1);
	j = min(max(0,j),n-1);
	k = min(max(0,k),n-1);
	if( w[i][j][k] > 0.5) {
		return x[fi][fj][fk];
	} else {
		return x[i][j][k];
	}
}

// Ans = Ax
static void compute_Ax( FLOAT ***x, FLOAT ***w, FLOAT ***ans, int n ) {
	FLOAT h2 = 1.0/(n*n);
	FOR_EVERY_COMP(n) {
		ans[i][j][k] = (x_ref(x,i,j,k,i+1,j,k,w,n)+x_ref(x,i,j,k,i-1,j,k,w,n)+
						x_ref(x,i,j,k,i,j+1,k,w,n)+x_ref(x,i,j,k,i,j-1,k,w,n)+
						x_ref(x,i,j,k,i,j,k+1,w,n)+x_ref(x,i,j,k,i,j,k-1,w,n)-6.0*x[i][j][k])/h2;
	} END_FOR
}

// Gauss-Seidel Iteration
static void gaussseidel( FLOAT ***x, FLOAT ***b, FLOAT ***w, int n, int t ) {
	FLOAT h2 = 1.0/(n*n);
	for( int k=0; k<t; k++ ) {
		FOR_EVERY_COMP(n) {
			x[i][j][k] = (x_ref(x,i,j,k,i+1,j,k,w,n)+x_ref(x,i,j,k,i-1,j,k,w,n)+
						  x_ref(x,i,j,k,i,j+1,k,w,n)+x_ref(x,i,j,k,i,j-1,k,w,n)+
						  x_ref(x,i,j,k,i,j,k+1,w,n)+x_ref(x,i,j,k,i,j,k-1,w,n)-h2*b[i][j][k]) / 6.0;
		} END_FOR
	}
}

// ans = x^T * x
static FLOAT product( FLOAT ***x, FLOAT ***y, int n ) {
	FLOAT ans = 0.0;
	FOR_EVERY_COMP(n) {
		ans += x[i][j][k]*y[i][j][k];
	} END_FOR
	return ans;
}

// x = 0
static void clear( FLOAT ***x, int n ) {
	FOR_EVERY_COMP(n) {
		x[i][j][k] = 0.0;
	} END_FOR
}

// x <= y
static void copy( FLOAT ***x, FLOAT ***y, int n ) {
	FOR_EVERY_COMP(n) {
		x[i][j][k] = y[i][j][k];
	} END_FOR
}

// Ans = x + a*y
static void op( FLOAT ***x, FLOAT ***y, FLOAT ***ans, FLOAT a, int n ) {
	static FLOAT ***tmp = alloc3D(n,n,n);
	FOR_EVERY_COMP(n) {
		tmp[i][j][k] = x[i][j][k]+a*y[i][j][k];
	} END_FOR
	copy(ans,tmp,n);
}

static void smooth( FLOAT ***x, FLOAT ***b, FLOAT ***w, int n, int t ) {
	// Smooth Using Gaus-Seidel Method
	gaussseidel( x, b, w, n, t );
}

// r = b - Ax
static void residual( FLOAT ***x, FLOAT ***b, FLOAT ***w, FLOAT ***r, int n ) {
	compute_Ax(x,w,r,n);
	op( b, r, r, -1.0, n );
}

// Shrink the image
static void shrink( FLOAT ***fine, FLOAT ***coarse, int fn ) {
	FOR_EVERY_COMP(fn/2) {
		// TODO: Interpolate Smoothly.
		FLOAT sum = 0.0;
		for( int i2=0; i2<2; i2++ ) for( int j2=0; j2<2; j2++ ) for( int k2=0; k2<2; k2++ ) {
			sum += fine[2*i+i2][2*j+j2][2*k+k2];
		}
		coarse[i][j][k] = sum/8.0;
	} END_FOR
}

// Expand the image
static void expand( FLOAT ***coarse, FLOAT ***fine, int fn ) {
	FOR_EVERY_COMP(fn) {
		// TODO: Interpolate Smoothly
		fine[i][j][k] = coarse[i/2][j/2][k/2];
	} END_FOR
}

// (V-Cycle Only) Multigrid Method
// TODO: Might Be Better Implement Full Multigrid Method
#define MAX_LAYER		8
static void mgv( FLOAT ***x, FLOAT ***b, FLOAT ***w, int n, int recr=0 ) {
	
	// Memory Saving Part
	static FLOAT ***fine_r[MAX_LAYER];
	static FLOAT ***fine_e[MAX_LAYER];
	static FLOAT ***coarse_r[MAX_LAYER];
	static FLOAT ***coarse_e[MAX_LAYER];
	static FLOAT ***coarse_w[MAX_LAYER];
	static FLOAT initialized = false;
	if( ! initialized  ) {
		for( int n=0; n<MAX_LAYER; n++ ) {
			fine_r[n] = NULL;
			fine_e[n] = NULL;
			coarse_r[n] = NULL;
			coarse_e[n] = NULL;
			coarse_w[n] = NULL;
		}
		initialized = true;
	}
	
	if( ! fine_r[recr] ) fine_r[recr] = alloc3D(n,n,n);
	if( ! fine_e[recr] ) fine_e[recr] = alloc3D(n,n,n);
	if( ! coarse_r[recr] ) coarse_r[recr] = alloc3D(n/2,n/2,n/2);
	if( ! coarse_e[recr] ) coarse_e[recr] = alloc3D(n/2,n/2,n/2);
	if( ! coarse_w[recr] ) coarse_w[recr] = alloc3D(n/2,n/2,n/2);
	
	clear(fine_r[recr],n);
	clear(fine_e[recr],n);
	clear(coarse_r[recr],n/2);
	clear(coarse_e[recr],n/2);
	clear(coarse_w[recr],n/2);
	
	///////////// Beginning of V-Cycle
	
	// Pre-smoothing
	smooth( x, b, w, n, 4 );
	
	// Compute Residual
	residual( x, b, w, fine_r[recr], n );
	
	// Restrict
	shrink(fine_r[recr],coarse_r[recr],n);
	shrink(w,coarse_w[recr],n);
	
	if( n <= 2 ) {
		// TODO: Should Be Solved Exactly
		smooth( coarse_e[recr], coarse_r[recr], coarse_w[recr], n/2, 10 );
	} else {
		// Recursively Call Itself
		mgv(coarse_e[recr],coarse_r[recr],coarse_w[recr],n/2,recr+1);
	}
	
	// Interpolate
	expand(coarse_e[recr],fine_e[recr],n);
	
	// Correct ( x = x + e )
	op( x, fine_e[recr], x, 1.0, n );
	
	// Post-smoothing
	smooth( x, b, w, n, 4 );
}

FLOAT solver::solve( FLOAT ***x, FLOAT ***b, FLOAT ***w, int n ) {
	static FLOAT ***r = alloc3D(n,n,n);
	clear(x,n);
	mgv(x,b,w,n);
	smooth( x, b, w, n, 4 );
	residual( x, b, w, r, n );
	return sqrt(product( r, r, n ))/(n*n);
}
