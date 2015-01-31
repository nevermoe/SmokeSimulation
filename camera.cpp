#include "core.h"
#include "camera.h"

Camera::Camera(GLFWwindow* windowHandle) {
#ifdef DEBUG_LEVEL
	std::cout << __FILE__ << " " << __FUNCTION__ << std::endl;
#endif
	_windowHandle = windowHandle;

	//SetLight();
	Reset();
}

void Camera::SetLight()
{
	//set light
	static const GLfloat light_position[] = {2.0f, 1.0f, 10.0f, 1.0f};
	static const GLfloat light_ambient[]  = {0.2f, 0.2f, 0.2f, 1.0f};
	static const GLfloat light_diffuse[]  = {1.0f, 1.0f, 0.0f, 1.0f};
	static const GLfloat light_specular[] = {1.0f, 0.0f, 1.0f, 1.0f};

	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);


	/*
	const static GLfloat mat_diffuse[] = {1.0f, 0.0f, 0.0f, 1.0f};
	const static GLfloat green_color[] = {0.0f, 1.0f, 0.0f, 0.3333f};
	const static GLfloat blue_color[] = {0.0f, 0.0f, 1.0f, 0.5f};

	
	GLfloat mat_shininess = 30.0;
	static const GLfloat mat_specular[] = {0.0f, 1.0f, 0.0f, 1.0f};
	static const GLfloat mat_emission[] = {0.0f, 0.0f, 0.0f, 1.0f};

	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR,  mat_specular);
	glMaterialfv(GL_FRONT, GL_EMISSION,  mat_emission);
	glMaterialf (GL_FRONT, GL_SHININESS, mat_shininess);
	*/


	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
}

void Camera::RegisterParentWindow(GLFWwindow* windowHandle)
{
	_windowHandle = windowHandle;
}

void Camera::Reset()
{

#if 1
	glfwGetWindowSize(_windowHandle, &_winX, &_winY);
	_FOV = 60.0f;
	_aspect = (GLfloat)_winX / _winY;
	_nearClip = 0.1f;
	_farClip = 100.0f;

#if 1
	glEnable(GL_DEPTH_TEST);
#else
	glDisable(GL_DEPTH_TEST);
#endif
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	// Set perspective projection
	gluPerspective(_FOV, _aspect, _nearClip, _farClip);

	glViewport(0, 0, _winX, _winY);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

#endif
}

void Camera::SetAspect(GLfloat aspect)		
{
	_aspect = aspect;
}
