#include "core.h"
#include "arcball.h"

/**
 * \ingroup GLVisualization
 * Default constructor, it sets the _ballRadius to 600
 **/
Arcball::Arcball()
{
	_ballRadius = 600;
	_isRotating = false;
	_isDragging = false;
	_zoomRate = 0.3;
	_winX = _winY = 0;
	Reset();
}

/**
 * \ingroup GLVisualization
 * Set _winX and _winY of the current windows, it's needed every time you resize the window
 * \param w _winX of the rendering window
 * \param h _winY of the rendering window
 **/
void Arcball::SetWidthHeight(int w, int h)
{  
	_winX = w;
	_winY = h;
	_ballRadius = std::min((int)(w/2), (int)(h/2));
}

/**
 * \ingroup GLVisualization
 * Set the radius of the ball (a typical radius for a 1024x768 window is 600
 * \param newRadius The radius of the spherical dragging area
 **/
void Arcball::SetRadius(float newRadius)
{  
	_ballRadius = newRadius;
}

void Arcball::StartDragging(int x, int y)
{
	_startDragX = x;
	_startDragY = y;
	_isDragging = true;
	_upDir = Eigen::Vector3f(0.0, 1.0, 0.0);

	GLdouble mvMatrix[16];
	glGetDoublev(GL_MODELVIEW_MATRIX, mvMatrix);
	_viewDir[0] = -mvMatrix[2];
	_viewDir[1] = -mvMatrix[6];
	_viewDir[2] = -mvMatrix[10];
#if 0
	_upDir[0] = mvMatrix[1];
	_upDir[1] = mvMatrix[5];
	_upDir[2] = mvMatrix[9];
	_rightDir[0] = mvMatrix[0];
	_rightDir[1] = mvMatrix[4];
	_rightDir[2] = mvMatrix[8];
#else
	_upDir = Eigen::Vector3f(0.0, 1.0, 0.0);
	_rightDir = _viewDir.cross(_upDir);
#endif
}

Eigen::Vector3f Arcball::UpdateDragging(int nx, int ny)
{
	_dragRate = 0.003 * _winX/1024;
	Eigen::Vector3f offset;
#if 1
	offset[0] = _dragRate * _rightDir[0] * (nx - _startDragX);
	offset[1] = _dragRate * _upDir[1] * (_startDragY - ny);
	offset[2] = _dragRate * _rightDir[2]*(nx-_startDragX) + _upDir[2]* (_startDragY-ny);
#else
	_viewDir.normalize();
	offset[0] = _dragRate * _viewDir[0] * (nx - _startDragX);
	offset[1] = _dragRate * _rightDir[0] * (_startDragY - ny);
	offset[2] = _dragRate * _upDir[0]*(nx-_startDragX) + _upDir[0]* (_startDragY-ny);
#endif

	_startDragX = nx;
	_startDragY = ny;
	return offset;
}

void Arcball::StopDragging()
{
	_isDragging = false;
}

void Arcball::StartZooming(int x, int y)
{
	//store original transform matrix
	glGetFloatv(GL_MODELVIEW_MATRIX, _startMatrix);

	_startZoomX = x;
	_startZoomY = y;
	_isZooming = true;
}

void Arcball::UpdateZooming(int x, int y)
{
	if (_isZooming) {
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glTranslatef(0.0f, 0.0f, _zoomRate * (x - _startZoomX) );
	}
	ApplyRotationMatrix();
}

void Arcball::StopZooming()
{
	_isZooming = false;
}

/**
 * \ingroup GLVisualization
 * Start the rotation. Use this method in association with the left click.
 * Here you must give directly the coordinates of the mouse as the glut functions extract. This method supposes that the 0,0 is in the upper-left part of the screen
 * \param _x Horizontal position of the mouse (0,0) = upperleft corner (w,h) = lower right
 * \param _y Vertical position of the mouse (0,0) = upperleft corner (w,h) = lower right
 *
 **/
void Arcball::StartRotation(int _x, int _y)
{
	int x = ( (_x)-(_winX/2) );
	int y = ((_winY/2)-_y);

	//store original transform matrix
	glGetFloatv(GL_MODELVIEW_MATRIX, _startMatrix);

	_startRotationVector = ConvertXY(x,y);
	_startRotationVector.normalize();

	_currentRotationVector=  _startRotationVector;
	_isRotating = true;

}

/**
 * \ingroup GLVisualization
 * Update the rotation. Use this method in association with the drag event.
 * Here you must give directly the coordinates of the mouse as the glut functions extract. This method supposes that the 0,0 is in the upper-left part of the screen
 * \param _x Horizontal position of the mouse (0,0) = upperleft corner (w,h) = lower right
 * \param _y Vertical position of the mouse (0,0) = upperleft corner (w,h) = lower right
 **/
void Arcball::UpdateRotation(int _x, int _y)
{  
	int x = ( (_x)-(_winX/2) );
	int y = ((_winY/2)-_y);

	_currentRotationVector = ConvertXY(x,y);

	_currentRotationVector.normalize();

	//Fixed by MY
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	ApplyRotationMatrix();
}

/**
 * \ingroup GLVisualization
 * Apply the computed rotation matrix
 * This method must be invoked inside the \code glutDisplayFunc() \endcode
 *
 **/
void Arcball::ApplyRotationMatrix()
{  
	//recover original matrix
	glMultMatrixf(_startMatrix);

	if (_isRotating) { 
		// Do some rotation according to start and current rotation vectors
		//cerr << _currentRotationVector.transpose() << " " << _startRotationVector.transpose() << endl;
		if ( ( _currentRotationVector - _startRotationVector).norm() > 1E-6 )	{
			Eigen::Vector3d rotationAxis = _currentRotationVector.cross(_startRotationVector);
			rotationAxis.normalize();

			//FIXED by MY
			Eigen::Matrix3d sm;
			for(int i = 0 ; i < 3 ; i++)
				for(int j = 0 ; j < 3 ; j++)
					sm(i,j) = (double)_startMatrix[4*i+j];
			rotationAxis = sm * rotationAxis;

			double val = _currentRotationVector.dot(_startRotationVector);
			val > (1-1E-10) ? val=1.0 : val=val ;
			double rotationAngle = acos(val) * 180.0f/(float)M_PI;

			// rotate around the current position
			glRotatef(rotationAngle * 2, -rotationAxis.x(),  -rotationAxis.y(),-rotationAxis.z());
		}
	}
}

/**
 * \ingroup GLVisualization
 * Stop the current rotation and prepare for a new click-then-drag event
 *
 **/
void Arcball::StopRotation()
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	ApplyRotationMatrix();
	_isRotating = false;
}


Eigen::Vector3d Arcball::ConvertXY(int x, int y)
{

	int d = x*x+y*y;
	float radiusSquared = _ballRadius * _ballRadius;
	if (d > radiusSquared) {
		return Eigen::Vector3d((float)x,(float)y, 0 );
	}
	else {  
		return Eigen::Vector3d((float)x,(float)y, sqrt(radiusSquared - d));
	}
}

/**
 * \ingroup GLVisualization
 * Reset the current transformation to the identity
 **/
void Arcball::Reset()
{  
	// reset matrix
	memset(_startMatrix, 0, sizeof(_startMatrix));
	_startMatrix[0] = 1;
	_startMatrix[1] =0;
	_startMatrix[2] = 0;
	_startMatrix[3] = 0;
	_startMatrix[4] = 0;
	_startMatrix[5] =1;
	_startMatrix[6] = 0;
	_startMatrix[7] = 0;
	_startMatrix[8] = 0;
	_startMatrix[9] =0;
	_startMatrix[10] = 1;
	_startMatrix[11] = 0;
	_startMatrix[12] = 0;
	_startMatrix[13] =0;
	_startMatrix[14] = 0;
	_startMatrix[15] = 1;

	_startZoomX = _startZoomY = 0;
}

