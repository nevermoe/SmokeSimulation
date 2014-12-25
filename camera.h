#ifndef CAMERA_H
#define CAMERA_H

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
	GLfloat FOV_;		// Field of View Angle
	GLfloat aspect_;	// Aspect Ratio
	GLfloat nearClip_;	// Near clipping plane distance
	GLfloat farClip_;	// Far clipping plane distance
	int winX_;
	int winY_;


	GLFWwindow* windowHandle_;
};

/*
The Camera class provides a simple means to controlling the 3D camera. It could
be extended to support more interactive controls. Ultimately. the camera sets the
GL projection and viewing matrices.
*/

#endif
