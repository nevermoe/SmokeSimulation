#ifndef GRIDOBJECT_H
#define GRIDOBJECT_H

#include "core.h"
#include "object.h"
#include "allocator.h"
#include "simple_vector.h"

class GridObject: public Object {
public:
	Normal*** normals;
	virtual ~GridObject();
protected:
	void _DeleteMemory();
};

#endif
