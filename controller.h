#ifndef CONTROLLER_H 
#define CONTROLLER_H

#include "core.h"
#include "object.h"
#include "arcball.h"
#include "camera.h"
#include "timer.h"

class Controller {
public:
	Controller(int argc,char **argv, const char *windowName);
	~Controller();

	void BeginLoop();
	void RegisterObject(Object* object);

	void Reset();
	void Render();
	void InitCamera();

	int GetActiveObject(int mx, int my);

	void Quit();

	// Event handlers
	void Resize(GLFWwindow *window, int x, int y);
	void Keyboard(GLFWwindow * window, int key, int scancode,int action, int mods);
	void MouseButton(GLFWwindow *window, int btn, int action, int mods);
	void MouseMotion(GLFWwindow *window, double x, double y);

private:
	// Window management
	std::string windowName_;
	std::stringstream titleInfo_;
	GLFWwindow *windowHandle_;
	int winX_, winY_;

	// Input
	bool isLeftKeyPressed_, isCtrlPressed_, isRightKeyPressed_, isMiddleKeyPressed_;
	bool isLightOn_;
	
	double prevX_, prevY_;

	Object** objects_;
	int objectNo_;
	int maxObjectNo_;
	int activeObj_;

	Camera* camera_;

	void _ComputeFPS();
	Timer timer_;
	double elapsedTime_;
	int numOfFrame_;
	double fps_;
};

#endif
