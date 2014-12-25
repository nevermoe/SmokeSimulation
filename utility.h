/*
 *  utility.h
 *  smoke3D
 *
 */

#include "types.h"

#define max(i,j) (i>j?i:j)
#define min(i,j) (i>j?j:i)

#define FOR_EVERY_X_FLOW	for( int i=0; i<N+1; i++ ) for( int j=0; j<N; j++ ) for( int k=0; k<N; k++ ) {
#define FOR_EVERY_Y_FLOW	for( int i=0; i<N; i++ ) for( int j=0; j<N+1; j++ ) for( int k=0; k<N; k++ ) {
#define FOR_EVERY_Z_FLOW	for( int i=0; i<N; i++ ) for( int j=0; j<N; j++ ) for( int k=0; k<N+1; k++ ) {
#define FOR_EVERY_CELL		for( int i=0; i<N; i++ ) for( int j=0; j<N; j++ ) for( int k=0; k<N; k++ ) {
#define END_FOR }

#ifdef _OPENMP
#include <omp.h>
#define OPENMP_FOR	_Pragma("omp parallel for" )
#define OPENMP_SECTION  _Pragma("omp section" )
#define OPENMP_BEGIN	_Pragma("omp parallel" ) {
#define OPENMP_END		}
#define OPENMP_FOR_P	_Pragma("omp for" )
#else
#define OPENMP_FOR
#define OPENMP_SECTION
#define OPENMP_BEGIN
#define OPENMP_END
#define OPENMP_FOR_P
#endif

FLOAT *** alloc3D( int w, int h, int d );
void free3D( FLOAT ***ptr );
void copy3D( FLOAT ***dst, FLOAT ***src, int n );
