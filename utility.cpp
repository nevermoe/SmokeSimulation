/*
 *  utility.cpp
 *  smoke3D
 *
 */

#include "utility.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

FLOAT *** alloc3D( int w, int h, int d ) {
	FLOAT *** field = new FLOAT **[w+1];
	for( int i=0; i<w; i++ ) {
		field[i] = new FLOAT*[h+1];
		for( int j=0; j<h; j++ ) {
			field[i][j] = new FLOAT[d];
			for( int k=0; k<d; k++ ) {
				field[i][j][k] = 0.0;
			}
		}
		field[i][h] = NULL;
	}
	field[w] = NULL;	
	return field;
}

void free3D( FLOAT ***ptr ) {
	for( int i=0; ptr[i]!=NULL; i++ ) {
		for( int j=0; ptr[i][j]!=NULL; j++ ) delete [] ptr[i][j];
		delete [] ptr[i];
	}
	delete [] ptr;
}

void copy3D( FLOAT ***dst, FLOAT ***src, int N ) {
	FOR_EVERY_CELL {
		dst[i][j][k] = src[i][j][k];
	} END_FOR
}