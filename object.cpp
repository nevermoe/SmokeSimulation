#include "core.h"
#include "object.h"

Object::Object() 
{
#ifdef DEBUG_LEVEL
	std::cout << __FILE__ << " " << __FUNCTION__ << std::endl;
#endif
	//Reset();
	_shaderProg = 0;
	memset(_hShaders, 0, sizeof(_hShaders) );
	_shaderNum = 0;

	//for rendering
	_stepSize = 0.001f;
}

Object::~Object()
{
}


int Object::LoadShaderFile(const char* fileName, GLuint shader)
{
	std::fstream file;
	GLint compiled;
	int size = 0;
	char* buffer = NULL;

	file.open(fileName, std::ios::in | std::ios::binary);
	if(!file) {
		std::cerr << "Open file " << fileName << " failed!" << std::endl;
		return -1;
	}

	//get file size and allocate memory
	file.seekg(0, std::ios::end);
	size = file.tellg();
	file.seekg(0, std::ios::beg);
	buffer = new char[size];

	file.read(buffer, size);
	//load to shader 
	glShaderSource(shader, 1, (const GLchar **)&buffer, &size);

	delete []buffer;
	file.close();

	return 0;
}

void Object::EnableShader()
{
	CreateShaderProgram();
}

void Object::RebindShader(GLuint hVertexShader, GLuint hfragShader)
{

	GLint testVal;
	GLsizei maxCount = 2;
	GLsizei count;
	GLuint shaders[maxCount];

	//get exists shaders
    glGetAttachedShaders(_shaderProg, maxCount, &count, shaders);
	//and detach then
	for (int i = 0; i < count; i++) {
		glDetachShader(_shaderProg, shaders[i]);
	}

    // Bind index 0 to the shader input variable "VerPos"
    glBindAttribLocation(_shaderProg, 0, "VerPos");
    // Bind index 1 to the shader input variable "VerClr"
    glBindAttribLocation(_shaderProg, 1, "VerClr");

	glAttachShader(_shaderProg, hVertexShader);
	glAttachShader(_shaderProg, hfragShader);
	
	// Attempt to link
	glLinkProgram(_shaderProg);

	// Make sure link worked too
	glGetProgramiv(_shaderProg, GL_LINK_STATUS, &testVal);
	if(testVal == GL_FALSE) {
		char infoLog[1024];
		glGetProgramInfoLog(_shaderProg, 1024, NULL, infoLog);
		std::cerr << "The program " << _shaderProg
			<< " failed to link with the following error:" << std::endl
			<< infoLog << std::endl;
		glDeleteProgram(_shaderProg);
		exit(-1);
	}

    glUseProgram(_shaderProg);
}

void Object::CreateShaderProgram()
{
	GLint testVal;
	GLuint hReturn;

	// Create the final program object, and attach the shaders
	hReturn = glCreateProgram();

	// All done, return our ready to use shader program
	_shaderProg = hReturn;
}


void Object::RegisterShader(const char* progName, GLenum shaderType)
{
	// Temporary Shader objects
	GLuint hShader;
	GLuint hReturn = 0;
	GLint testVal;
	
	// Create shader objects
	hShader = glCreateShader(shaderType);
	
	// Load, If fail clean up and return null
	if(LoadShaderFile(progName, hShader) < 0) {
		glDeleteShader(hShader);
		std::cerr << "The shader at " << progName 
			<< " could not be found." << std::endl;
		exit(-1);
	}

	// Compile 
	glCompileShader(hShader);
	// Check for errors in shader
	glGetShaderiv(hShader, GL_COMPILE_STATUS, &testVal);
	if(testVal == GL_FALSE) {
		char infoLog[1024];
		glGetShaderInfoLog(hShader, 1024, NULL, infoLog);
		std::cerr << "The shader at " << progName 
			<< " failed to compile with the following error:" << std::endl
			<< infoLog << std::endl;
		glDeleteShader(hShader);
		exit(-1);
	}
	// Check for errors in shader
	glGetShaderiv(hShader, GL_COMPILE_STATUS, &testVal);
	if(testVal == GL_FALSE) {
		char infoLog[1024];
		glGetShaderInfoLog(hShader, 1024, NULL, infoLog);
		std::cerr << "The shader at " << progName 
			<< " failed to compile with the following error:" << std::endl
			<< infoLog << std::endl;
		glDeleteShader(hShader);
		exit(-1);
	}

	_hShaders[_shaderNum++] = hShader;


#if 0
	// Now, we need to bind the attribute names to their specific locations
	// List of attributes
	va_list attributeList;
	va_start(attributeList, fragmentProgName);
	// Iterate over this argument list
	char *szNextArg;
	int iArgCount = va_arg(attributeList, int);
	// Number of attributes
	for(int i = 0; i < iArgCount; i++)
	{
		int index = va_arg(attributeList, int);
		szNextArg = va_arg(attributeList, char*);
		glBindAttribLocation(hReturn, index, szNextArg);
	}
	va_end(attributeList);
#endif
}

void Object::RegisterParentWindow(GLFWwindow* windowHandle)
{
	_windowHandle = windowHandle;
}

void Object::Resize(GLFWwindow* windowHandle, int x, int y)
{
	_arcball.SetWidthHeight(x, y);
}

void Object::MouseButton(GLFWwindow *window, int button,int action,int mods) 
{
	double mouseX, mouseY;

	if(button == GLFW_MOUSE_BUTTON_LEFT) {
		::glfwGetCursorPos(window, &mouseX, &mouseY);
		if(action == GLFW_PRESS) {
			_isLeftKeyPressed = true;
			_arcball.StartRotation(mouseX, mouseY);
		}
		else if(action == GLFW_RELEASE) {
			_isLeftKeyPressed = false;
			_arcball.StopRotation();
		}
	}
	else if(button == GLFW_MOUSE_BUTTON_MIDDLE) {
		_isMiddleKeyPressed = (action == GLFW_PRESS);
	}
	else if(button == GLFW_MOUSE_BUTTON_RIGHT) {
#if 0
		if (action == GLFW_PRESS) {
			_isRightKeyPressed = true;
			::glfwGetCursorPos(window, &mouseX, &mouseY);
			_arcball.StartZooming(mouseX, mouseY);
		}
		else if(action ==GLFW_RELEASE) {
			_isRightKeyPressed = false;
			_arcball.StopZooming();
		}
#endif
	}
}

void Object::MouseMotion(GLFWwindow *window, double nx, double ny) 
{
	if(_isLeftKeyPressed && _isCtrlPressed) {
	}
	else if(_isLeftKeyPressed) {
		_arcball.UpdateRotation(nx, ny);
	}
#if 0
	else if(_isRightKeyPressed) {
		_arcball.UpdateZooming(nx, ny);
	}
#endif
}

void Object::MouseScroll(GLFWwindow *window, double nx, double ny)
{
	_arcball.StartZooming(0, 0);
	_arcball.UpdateZooming(-ny, nx);
	_arcball.StopZooming();
}

void Object::ComputeNormal(GLfloat v[3][3], GLfloat normal[]) 
{ 
	float v1[3],v2[3];
	static const int x = 0;
	static const int y = 1;
	static const int z = 2;

	// Calculate two vectors from the three points
	v1[x] = v[0][x] - v[1][x];
	v1[y] = v[0][y] - v[1][y];
	v1[z] = v[0][z] - v[1][z];

	v2[x] = v[1][x] - v[2][x];
	v2[y] = v[1][y] - v[2][y];
	v2[z] = v[1][z] - v[2][z];

	// Take the cross product of the two vectors to get
	// the normal vector which will be stored in out
	normal[x] = v1[y]*v2[z] - v1[z]*v2[y];
	normal[y] = v1[z]*v2[x] - v1[x]*v2[z];
	normal[z] = v1[x]*v2[y] - v1[y]*v2[x];

	// Normalize the vector (shorten length to one)
}

void Object::Reset() {
#ifdef DEBUG_LEVEL
	std::cout << __FILE__ << " " << __FUNCTION__ << std::endl;
#endif
	int width, height;
	glfwGetWindowSize(_windowHandle, &width, &height);
	_winX = width;
	_winY = height;
	_arcball.SetWidthHeight(width, height);

#if 1
	_depth = 3.8;
	_rotX = 10;
	_rotY = -20;


	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glTranslatef(0, 0, -_depth);
	glRotatef(_rotX, 1, 0, 0);
	glRotatef(_rotY, 0, 1, 0);
#endif
}

bool Object::LoadFile(char* fileName)
{
	_file.open(fileName);
	if (!_file) {
		std::cout << "Open file: " << fileName << " failed!" << std::endl;
		return false;
	}
	std::cout << "Open file: " << fileName << " succeed!" << std::endl;

	return true;
}

//this is a test function, just draw a cube
void Object::Cube() { 
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	GLfloat vertices[][3] = { { -1.0f/2, -1.0f/2, -1.0f/2 },
		{ 1.0f/2, -1.0f/2, -1.0f/2}, { 1.0f/2, 1.0f/2, -1.0f/2 },
		{ -1.0f/2, 1.0f/2, -1.0f/2 }, { -1.0f/2, -1.0f/2, 1.0f/2 },
		{ 1.0f/2, -1.0f/2, 1.0f/2 }, { 1.0f/2, 1.0f/2, 1.0f/2 },
		{ -1.0f/2, 1.0f/2, 1.0f/2 } };              // Vertices
	int faces[][4] = { { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 4, 5, 6, 7 },
		{ 0, 4, 7, 3 }, { 0, 1, 5, 4 }, { 0, 3, 2, 1 } };

	GLfloat colors[][3] = { { 0.0f, 1.0f, 1.0f }, { 1.0f, 0.0f, 1.0f },
		{ 1.0f, 1.0f, 0.0f }, { 0.0f, 0.5f, 0.5f },
		{ 0.5f, 0.0f, 0.5f }, { 0.5f, 0.5f, 0.0f } };

	glBegin(GL_QUADS);  
	for (int i = 0; i < 6; i++) {
		glColor3fv(colors[i]);
		
		Eigen::Vector3f nm;
		Eigen::Vector3f v0(vertices[faces[i][0]][0], vertices[faces[i][0]][1], vertices[faces[i][0]][2]);
		Eigen::Vector3f v1(vertices[faces[i][1]][0], vertices[faces[i][1]][1], vertices[faces[i][1]][2]);
		Eigen::Vector3f v2(vertices[faces[i][2]][0], vertices[faces[i][2]][1], vertices[faces[i][2]][2]);
		Eigen::Vector3f l1 = v1-v0;
		Eigen::Vector3f l2 = v2-v0;
		nm = l1.cross(l2);
		nm.normalize();

		glNormal3f(nm(0),nm(1),nm(2));
		for (int j = 0; j < 4; j++)
			glVertex3fv(vertices[faces[i][j]]);
	}
	glEnd();
}

void Object::Show() 
{
	Cube();
}
