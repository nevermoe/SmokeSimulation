#include "core.h"

Camera::Camera(GLFWwindow* windowHandle) {
	windowHandle_ = windowHandle;
	Reset();
}

void Camera::SetLight()
{
	//set light
	static const GLfloat light_position[] = {2.0f, 1.0f, 10.0f, 1.0f};
	static const GLfloat light_ambient[]  = {0.2f, 0.2f, 0.2f, 1.0f};
	static const GLfloat light_diffuse[]  = {1.0f, 1.0f, 1.0f, 1.0f};
	static const GLfloat light_specular[] = {1.0f, 1.0f, 1.0f, 1.0f};

	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
/*
	const static GLfloat mat_diffuse[] = {1.0f, 0.0f, 0.0f, 1.0f};
	const static GLfloat green_color[] = {0.0f, 1.0f, 0.0f, 0.3333f};
	const static GLfloat blue_color[] = {0.0f, 0.0f, 1.0f, 0.5f};

	
	GLfloat mat_shininess = 30.0;
	static const GLfloat mat_specular[] = {0.0f, 0.0f, 0.0f, 1.0f};
	static const GLfloat mat_emission[] = {0.0f, 0.0f, 0.0f, 1.0f};

	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR,  mat_specular);
	glMaterialfv(GL_FRONT, GL_EMISSION,  mat_emission);
	glMaterialf (GL_FRONT, GL_SHININESS, mat_shininess);
*/
	
}

void Camera::SetParentWindow(GLFWwindow* windowHandle)
{
	windowHandle_ = windowHandle;
}

void Camera::Reset()
{
	glfwGetWindowSize(windowHandle_, &winX_, &winY_);
	FOV_ = 60.0f;
	aspect_ = (GLfloat)winX_ / winY_;
	nearClip_ = 0.1f;
	farClip_ = 100.0f;

	glEnable(GL_DEPTH_TEST);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	// Set perspective projection
	gluPerspective(FOV_, aspect_, nearClip_, farClip_);

	glViewport(0, 0, winX_, winY_);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void Camera::SetAspect(GLfloat aspect)		
{
	aspect_ = aspect;
}
