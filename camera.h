#ifndef _CAMERA_H
#define _CAMERA_H

#include "core.h"

class Camera {
public:
	Camera(GLFWwindow* windowHandle);

	void Reset();
	void SetLight();
	void RegisterParentWindow(GLFWwindow* windowHandle);

	
	// Access functions
	void SetAspect(GLfloat aspect);

private:
	// Perspective controls
	GLfloat _FOV;		// Field of View Angle
	GLfloat _aspect;	// Aspect Ratio
	GLfloat _nearClip;	// Near clipping plane distance
	GLfloat _farClip;	// Far clipping plane distance
	int _winX;
	int _winY;


	GLFWwindow* _windowHandle;
};

/*
The Camera class provides a simple means to controlling the 3D camera. It could
be extended to support more interactive controls. Ultimately. the camera sets the
GL projection and viewing matrices.
*/

#endif
