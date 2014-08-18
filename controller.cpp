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

void Controller::RegisterObject(Object* object)
{
	object->SetParentWindow(windowHandle_);
	object->Reset();
	objects_[objectNo_++] = object;
	activeObj_++;
}


Controller::Controller(int argc,char **argv, const char *windowName) {
	windowName_ = windowName;

	//for timer
	elapsedTime_ = 0;
	numOfFrame_ = 0;
	fps_ = 0;

	winX_ = 1280;
	winY_ = 960;

	activeObj_ = -1;
	objectNo_ = 0;
	maxObjectNo_ = 100;
	objects_ = new Object*[maxObjectNo_];

	isCtrlPressed_ = isLeftKeyPressed_ = isMiddleKeyPressed_ = isRightKeyPressed_ = false;
	isLightOn_ = true;
	prevX_ = prevY_ = 0;

	// Initialize components
	if (!glfwInit()) {
		std::cout << "glfwInit() failed!" << std::endl;
		exit(0);
	}

	// Create the window
	windowHandle_ = glfwCreateWindow(winX_, winY_, windowName, NULL, NULL);
	if (!windowHandle_) {
		std::cout << "Create Window failed"  << std::endl;
		exit(0);
	}
	glfwMakeContextCurrent(windowHandle_);
	glfwSetWindowPos(windowHandle_, 0, 0);

	// Background color
	glClearColor( 0., 0., 0., 1. );

	// Callbacks
	glfwSetErrorCallback(error_callback);
	glfwSetMouseButtonCallback(windowHandle_, mousebutton);
	glfwSetCursorPosCallback(windowHandle_, mousemotion);
	glfwSetKeyCallback(windowHandle_, keyboard);
	glfwSetWindowSizeCallback(windowHandle_, resize);

	InitCamera();
}

void Controller::InitCamera()
{
	camera_ = new Camera(windowHandle_);
}

void Controller::BeginLoop()
{
	while (!glfwWindowShouldClose(windowHandle_))
    {
        /* Render here */
		::g_controller->Render();

        /* Swap front and back buffers */
        glfwSwapBuffers(windowHandle_);

        /* Poll for and process events */
        glfwPollEvents();
    }
}

//////////////////////////////////////////////////////////////////////////////

Controller::~Controller() {
	glFinish();
	glfwDestroyWindow(windowHandle_);
	glfwTerminate();
}

void Controller::Reset() {
	camera_->Reset();

	for(int i = 0 ; i < objectNo_ ; i++)
		objects_[i]->Reset();
}

void Controller::_ComputeFPS()
{
	numOfFrame_++;
	elapsedTime_ += timer_.StopTimer();
	timer_.StartTimer();
	if(elapsedTime_ > 1) {
		fps_ = numOfFrame_ / elapsedTime_;
		//std::cout << fps_  << " "<< numOfFrame_ << std::endl;
		elapsedTime_ = 0;
		numOfFrame_ = 0;
	}
	
	//clear stringstream buffer
	titleInfo_.str("");
	titleInfo_ << windowName_;
	titleInfo_ << "     FPS: ";
	titleInfo_ << std::setprecision(4) << fps_;
	titleInfo_.width(2);
}

void Controller::Render() {

	//compute fps and display 
	_ComputeFPS();
	glfwSetWindowTitle(windowHandle_, titleInfo_.str().c_str() );

	//Begin drawing scene
	camera_->Reset();
	for(int i = 0 ; i < objectNo_ ; i++) {
		objects_[i]->SimulateStep();
		objects_[i]->Show();
	}


	//Finish drawing scene
	//glFinish();
}

int Controller::GetActiveObject(int mx, int my)
{
	//do some judgement;
	return activeObj_;
}

void Controller::Quit() {
	glFinish();
	glfwDestroyWindow(windowHandle_);
	exit(0);
}


void Controller::Resize(GLFWwindow *window, int x, int y) {
	winX_ = x;
	winY_ = y;
	
	for(int i = 0 ; i < objectNo_ ; i++)
		objects_[i]->Resize(window, x, y);
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
				isCtrlPressed_ = true;
				break;
		}
	}
	else if(action == GLFW_RELEASE)
		switch (key) {
			case GLFW_KEY_LEFT_CONTROL:
				isCtrlPressed_ = false;
				break;
		}
}


void Controller::MouseButton(GLFWwindow *window, int button,int action,int mods) 
{

	//get active object and then transfer the message to the object
	if(action == GLFW_PRESS) {
		glfwGetCursorPos(window, &prevX_, &prevY_);
	}

	activeObj_ = GetActiveObject(prevX_, prevY_);
	objects_[activeObj_]->MouseButton(window, button, action, mods);
}


void Controller::MouseMotion(GLFWwindow *window, double nx, double ny) 
{
	objects_[activeObj_]->MouseMotion(window, nx, ny);
}

