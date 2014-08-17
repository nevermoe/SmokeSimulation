#include "core.h"
#include "gridobject.h"

GridObject::~GridObject() 
{
	_DeleteMemory();
}

void GridObject::_DeleteMemory()
{
	Allocator<Normal> normFreer;
	normFreer.Free3D(normals);
}
