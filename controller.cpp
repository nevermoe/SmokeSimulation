#include "core.h"
#include "controller.h"

extern Controller* g_controller;

//set callback functions
static void error_callback(int error, const char* description)
{
    fputs(description, stderr);
}

static void resize(GLFWwindow *window, int x,int y)
{
	::g_controller->Resize(window, x, y);
}

static void keyboard(GLFWwindow * window, int key, int scancode,int action, int mods)		
{
	::g_controller->Keyboard(window, key, scancode, action, mods);
}

static void mousebutton(GLFWwindow *window,int button,int action,int mods)	
{
	::g_controller->MouseButton(window, button, action, mods);
}

static void mousemotion(GLFWwindow *window, double x, double y)
{
	::g_controller->MouseMotion(window, x, y);
}

static void mousescroll(GLFWwindow *window, double x, double y)
{
	::g_controller->MouseScroll(window, x, y);
}

void Controller::RegisterObject(Object* object)
{
	object->RegisterParentWindow(_windowHandle);
	object->Reset();
	_objects[_objectNo++] = object;
	_activeObj++;
}


Controller::Controller(int argc,char **argv, const char *windowName) 
{
	_windowName = windowName;

	//for timer
	_elapsedTime = 0;
	_numOfFrame = 0;
	_fps = 0;

	_winX = 1280;
	_winY = 960;

	_activeObj = -1;
	_objectNo = 0;
	_maxObjectNo = 100;
	_objects = new Object*[_maxObjectNo];

	_isCtrlPressed = _isLeftKeyPressed = _isMiddleKeyPressed = _isRightKeyPressed = false;
	_isLightOn = true;
	_prevCursorX = _prevCursorY = 0;

	// Initialize components
	if (!glfwInit()) {
		std::cerr << "glfwInit() failed!" << std::endl;
		exit(-1);
	}

	// Create the window
	_windowHandle = glfwCreateWindow(_winX, _winY, windowName, NULL, NULL);
	if (!_windowHandle) {
		std::cerr << "Create Window failed"  << std::endl;
		exit(-1);
	}
	glfwMakeContextCurrent(_windowHandle);
	glfwSetWindowPos(_windowHandle, 0, 0);

	// Background color
	glClearColor( 0., 0., 0., 1. );

	// Callbacks
	glfwSetErrorCallback(error_callback);
	glfwSetMouseButtonCallback(_windowHandle, mousebutton);
	glfwSetScrollCallback(_windowHandle, mousescroll);
	glfwSetCursorPosCallback(_windowHandle, mousemotion);
	glfwSetKeyCallback(_windowHandle, keyboard);
	glfwSetWindowSizeCallback(_windowHandle, resize);

	if(glewInit() != GLEW_OK) {
		std::cerr << "glewInit() failed!" << std::endl;
		exit(-1);
	}

	//InitCamera();
}

void Controller::InitCamera()
{
	_camera = new Camera(_windowHandle);
}

void Controller::BeginLoop()
{
	while (!glfwWindowShouldClose(_windowHandle))
    {
        /* Render here */
		::g_controller->Render();

        /* Swap front and back buffers */
        glfwSwapBuffers(_windowHandle);

        /* Poll for and process events */
        glfwPollEvents();
    }
}

//////////////////////////////////////////////////////////////////////////////

Controller::~Controller() {
	glFinish();
	glfwDestroyWindow(_windowHandle);
	glfwTerminate();
}

void Controller::Reset() {
#ifdef DEBUG_LEVEL
	std::cout << __FILE__ << " " << __FUNCTION__ << std::endl;
#endif
	_camera->Reset();

	for(int i = 0 ; i < _objectNo ; i++)
		_objects[i]->Reset();
}

void Controller::_ComputeFPS()
{
	_numOfFrame++;
	_elapsedTime += _timer.StopTimer();
	_timer.StartTimer();
	if(_elapsedTime > 1) {
		_fps = _numOfFrame / _elapsedTime;
		//std::cout << _fps  << " "<< _numOfFrame << std::endl;
		_elapsedTime = 0;
		_numOfFrame = 0;
	}
	
	//clear stringstream buffer
	_titleInfo.str("");
	_titleInfo << _windowName;
	_titleInfo << "     FPS: ";
	_titleInfo << std::setprecision(4) << _fps;
	_titleInfo.width(2);
}

void Controller::Render() {

	//compute fps and display 
	_ComputeFPS();
	glfwSetWindowTitle(_windowHandle, _titleInfo.str().c_str() );


	//Begin drawing scene
	_camera->Reset();
	for(int i = 0 ; i < _objectNo ; i++) {
		_objects[i]->SimulateStep();
		_objects[i]->Show();
	}


	//Finish drawing scene
	//glFinish();
}

int Controller::GetActiveObject(int mx, int my)
{
	//FIXME: do some judgement;
	return _activeObj;
}

void Controller::Quit() {
	glFinish();
	glfwDestroyWindow(_windowHandle);
	exit(0);
}


void Controller::Resize(GLFWwindow *window, int x, int y) {
	_winX = x;
	_winY = y;
	
	for(int i = 0 ; i < _objectNo ; i++)
		_objects[i]->Resize(window, x, y);
}

void Controller::Keyboard(GLFWwindow * window, int key, int scancode, int action, int mods)
{
	if (action == GLFW_PRESS) {
		switch(key) {
			case GLFW_KEY_ESCAPE:		// Escape
				Quit();
				break;
			case GLFW_KEY_R:			//reset
				Reset();
				break;
			case GLFW_KEY_LEFT_CONTROL:
				_isCtrlPressed = true;
				break;
		}
	}
	else if(action == GLFW_RELEASE)
		switch (key) {
			case GLFW_KEY_LEFT_CONTROL:
				_isCtrlPressed = false;
				break;
		}

	//route message to active object
	_objects[_activeObj]->Keyboard(window, key, scancode, action, mods);
}


void Controller::MouseButton(GLFWwindow *window, int button,int action,int mods) 
{
	//get active object and then transfer the message to the object
	if(action == GLFW_PRESS) {
		glfwGetCursorPos(window, &_prevCursorX, &_prevCursorY);
	}
	_activeObj = GetActiveObject(_prevCursorX, _prevCursorY);
	_objects[_activeObj]->MouseButton(window, button, action, mods);
}


void Controller::MouseMotion(GLFWwindow *window, double nx, double ny) 
{
	_objects[_activeObj]->MouseMotion(window, nx, ny);
}

void Controller::MouseScroll(GLFWwindow *window, double nx, double ny)
{
	_objects[_activeObj]->MouseScroll(window, nx, ny);
}
