#ifndef DRAWER_H 
#define DRAWER_H

#include "core.h"

////////////////////////////////////////////////////////////////////////////////

class Object: public EventListener{
public:
	Object();


	void ComputeNormal(GLfloat v[3][3], GLfloat normal[]);

	void Cube();

	virtual void SimulateStep() {};
	virtual void Show();
	virtual bool LoadFile(char* fileName);
	virtual void ComputeVertexNormal() {}
	virtual void Reset();


	virtual void MouseButton(GLFWwindow *window, int button,int action,int mods);
	virtual void MouseMotion(GLFWwindow *window, double nx, double ny);
	virtual void Keyboard(GLFWwindow * window, int key, int scancode, int action, int mods) {}
	virtual void Resize(GLFWwindow *window, int x, int y);

	void SetParentWindow(GLFWwindow* windowHandle);

protected:
	std::fstream file_;

	GLfloat rotX_;
	GLfloat rotY_;
	GLfloat depth_;

	//for event handling
	bool isLeftKeyPressed_, isCtrlPressed_, isRightKeyPressed_, isMiddleKeyPressed_;

	Arcball arcball_;

	GLFWwindow* windowHandle_;

#define SMOOTH_SHADING 0
#define FLAT_SHADING 1
#define CREASE_SHADING 2

#define GAUSSIAN_CURVATURE 1
#define MEAN_CURVATURE 2
};

////////////////////////////////////////////////////////////////////////////////

#endif
