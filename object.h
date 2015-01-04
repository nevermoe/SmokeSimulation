#ifndef DRAWER_H 
#define DRAWER_H

#include "core.h"
#include "arcball.h"

////////////////////////////////////////////////////////////////////////////////

class Object: public EventListener{
public:
	Object();
	~Object();

	void RegisterShader(const char* vertexProgName, GLenum shaderType);
	int LoadShaderFile(const char* fileName, GLuint shader);
	void RebindShader(GLuint hVertexShader, GLuint hfragShader);
	void CreateShaderProgram();
	void EnableShader();

	void rcSetUinforms();


	void ComputeNormal(GLfloat v[3][3], GLfloat normal[]);

	void Cube();

	virtual void SimulateStep() {};
	virtual void Show();
	virtual bool LoadFile(char* fileName);
	virtual void ComputeVertexNormal() {}
	virtual void Reset();


	virtual void MouseButton(GLFWwindow *window, int button,int action,int mods);
	virtual void MouseMotion(GLFWwindow *window, double nx, double ny);
	virtual void MouseScroll(GLFWwindow *window, double nx, double ny);
	virtual void Keyboard(GLFWwindow * window, int key, int scancode, int action, int mods) {}
	virtual void Resize(GLFWwindow *window, int x, int y);

	void RegisterParentWindow(GLFWwindow* windowHandle);

	//for rendering
	GLuint initTFF1DTex(const char* filename);
	GLuint initFace2DTex(GLuint bfTexWidth, GLuint bfTexHeight);
	GLuint initVol3DTex(const char* filename, GLuint w, GLuint h, GLuint d);
	void checkFramebufferStatus();
	void initFrameBuffer(GLuint texObj, GLuint texWidth, GLuint texHeight);
	void render(GLenum cullFace);

protected:
	std::fstream _file;

	GLfloat _rotX;
	GLfloat _rotY;
	GLfloat _depth;

	//for event handling
	bool _isLeftKeyPressed, _isCtrlPressed, _isRightKeyPressed, _isMiddleKeyPressed;

	Arcball _arcball;

	GLFWwindow* _windowHandle;
	int _winX, _winY;

	GLuint _shaderProg;	 //shader
	GLuint _hShaders[10]; //support upto 10 shaders;
	int _shaderNum;			//current number of shaders

	//for rendering
	float _stepSize;
	GLuint _vao;
	GLuint _frameBuffer;
	// transfer function
	GLuint _tffTexObj;
	GLuint _bfTexObj;
	GLuint _volTexObj;

};

////////////////////////////////////////////////////////////////////////////////

#endif
