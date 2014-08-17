#define _USE_MATH_DEFINES

#include "core.h"
#include "arcball.h"

/**
 * \ingroup GLVisualization
 * Default constructor, it sets the ballRadius to 600
 **/
Arcball::Arcball()
{
	this->ballRadius=600;
	isRotating=false;
	zoomRate = 0.003;
	width=height=0;
	Reset();
}

/**
 * \ingroup GLVisualization
 * Set width and height of the current windows, it's needed every time you resize the window
 * \param w Width of the rendering window
 * \param h Height of the rendering window
 **/
void Arcball::SetWidthHeight(int w, int h)
{  
	width=w;
	height=h;
	ballRadius = std::min((int)(w/2), (int)(h/2));
}

/**
 * \ingroup GLVisualization
 * Set the radius of the ball (a typical radius for a 1024x768 window is 600
 * \param newRadius The radius of the spherical dragging area
 **/
void Arcball::SetRadius(float newRadius)
{  
	ballRadius = newRadius;
}

void Arcball::StartZooming(int x, int y)
{
	//store original transform matrix
	glGetFloatv(GL_MODELVIEW_MATRIX, startMatrix);

	startZoomX = x;
	startZoomY = y;
	isZooming = true;
}

void Arcball::UpdateZooming(int x, int y)
{
	if (isZooming) {
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glTranslatef(0.0f, 0.0f, zoomRate * (x - startZoomX) );
	}
	ApplyRotationMatrix();
}

void Arcball::StopZooming()
{
	isZooming = false;
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
	int x = ( (_x)-(width/2) );
	int y = ((height/2)-_y);

	//store original transform matrix
	glGetFloatv(GL_MODELVIEW_MATRIX, startMatrix);

	startRotationVector = ConvertXY(x,y);
	startRotationVector.normalize();

	currentRotationVector=  startRotationVector;
	isRotating = true;

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
	int x = ( (_x)-(width/2) );
	int y = ((height/2)-_y);

	currentRotationVector = ConvertXY(x,y);

	currentRotationVector.normalize();

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
	glMultMatrixf(startMatrix);

	if (isRotating) { 
		// Do some rotation according to start and current rotation vectors
		//cerr << currentRotationVector.transpose() << " " << startRotationVector.transpose() << endl;
		if ( ( currentRotationVector - startRotationVector).norm() > 1E-6 )	{
			Eigen::Vector3d rotationAxis = currentRotationVector.cross(startRotationVector);
			rotationAxis.normalize();

			//FIXED by MY
			Eigen::Matrix3d sm;
			for(int i = 0 ; i < 3 ; i++)
				for(int j = 0 ; j < 3 ; j++)
					sm(i,j) = (double)startMatrix[4*i+j];
			rotationAxis = sm * rotationAxis;

			double val = currentRotationVector.dot(startRotationVector);
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
	isRotating = false;
}


Eigen::Vector3d Arcball::ConvertXY(int x, int y)
{

	int d = x*x+y*y;
	float radiusSquared = ballRadius*ballRadius;
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
	fov = INITIAL_FOV;
	// reset matrix
	memset(startMatrix, 0, sizeof(startMatrix));
	startMatrix[0] = 1;
	startMatrix[1] =0;
	startMatrix[2] = 0;
	startMatrix[3] = 0;
	startMatrix[4] = 0;
	startMatrix[5] =1;
	startMatrix[6] = 0;
	startMatrix[7] = 0;
	startMatrix[8] = 0;
	startMatrix[9] =0;
	startMatrix[10] = 1;
	startMatrix[11] = 0;
	startMatrix[12] = 0;
	startMatrix[13] =0;
	startMatrix[14] = 0;
	startMatrix[15] = 1;

	transX = transY = 0;
	startZoomX = startZoomY = currentTransX = currentTransY = 0;
}


const float Arcball::INITIAL_FOV = 30;
const float Arcball::TRANSLATION_FACTOR = 0.01f;

