#ifndef ALLOC3D_H
#define ALLOC3D_H

#include "core.h"

template <class T>
class Allocator{
public:
	T*** Alloc3D(int w, int h, int d);
	void Free3D(T*** ptr);
};

template <class T>
T*** Allocator<T>::Alloc3D(int w, int h, int d)
{
	T *** field = new T **[w+1];
	for( int i=0; i<w; i++ ) {
		field[i] = new T*[h+1];
		for( int j=0; j<h; j++ ) {
			field[i][j] = new T[d];
			memset(field[i][j], 0, sizeof(T)*d);
		}
		field[i][h] = NULL;
	}
	field[w] = NULL;	
	return field;
}

template <class T>
void Allocator<T>::Free3D(T*** ptr)
{
	for( int i=0; ptr[i] != NULL; i++ ) {
		for( int j=0; ptr[i][j] != NULL; j++ ) {
			delete []ptr[i][j];
		}
		delete []ptr[i];
	}

	delete ptr;
}
#endif
