#include "core.h"

Object::Object() {
	//Reset();
}

void Object::SetParentWindow(GLFWwindow* windowHandle)
{
	windowHandle_ = windowHandle;
}

void Object::Resize(GLFWwindow* windowHandle, int x, int y)
{
	arcball_.setWidthHeight(x, y);
}

void Object::MouseButton(GLFWwindow *window, int button,int action,int mods) 
{
	double mouseX, mouseY;

	if(button == GLFW_MOUSE_BUTTON_LEFT) {
		::glfwGetCursorPos(window, &mouseX, &mouseY);
		if(action == GLFW_PRESS) {
			isLeftKeyPressed_ = true;
			arcball_.startRotation(mouseX, mouseY);
		}
		else if(action == GLFW_RELEASE) {
			isLeftKeyPressed_ = false;
			arcball_.stopRotation();
		}
	}
	else if(button == GLFW_MOUSE_BUTTON_MIDDLE) {
		isMiddleKeyPressed_ = (action == GLFW_PRESS);
	}
	else if(button == GLFW_MOUSE_BUTTON_RIGHT) {
		if (action == GLFW_PRESS) {
			isRightKeyPressed_ = true;
			::glfwGetCursorPos(window, &mouseX, &mouseY);
			arcball_.startZooming(mouseX, mouseY);
		}
		else if(action ==GLFW_RELEASE) {
			isRightKeyPressed_ = false;
			arcball_.stopZooming();
		}
	}
}

void Object::MouseMotion(GLFWwindow *window, double nx, double ny) 
{
	if(isLeftKeyPressed_ && isCtrlPressed_) {
	}
	else if(isLeftKeyPressed_) {
		arcball_.updateRotation(nx, ny);
	}
	else if(isRightKeyPressed_) {
		arcball_.updateZooming(nx, ny);
	}
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
	int width, height;
	glfwGetWindowSize(windowHandle_, &width, &height);
	arcball_.setWidthHeight(width, height);

	depth_ = 5.0;
	rotX_ = 25;
	rotY_ = -30;


	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glTranslatef(0, 0, -depth_);
	glRotatef(rotX_, 1, 0, 0);
	glRotatef(rotY_, 0, 1, 0);
}

bool Object::LoadFile(char* fileName)
{
	file_.open(fileName);
	if (!file_) {
		std::cout << "Open file: " << fileName << " failed!" << std::endl;
		return false;
	}
	std::cout << "Open file: " << fileName << " succeed!" << std::endl;

	return true;
}

void Object::Cube() {                               // 立方体の描画
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		GLfloat vertices[][3] = { { -1.0f, -1.0f, -1.0f },
				{ 1.0f, -1.0f, -1.0f }, { 1.0f, 1.0f, -1.0f },
				{ -1.0f, 1.0f, -1.0f }, { -1.0f, -1.0f, 1.0f },
				{ 1.0f, -1.0f, 1.0f }, { 1.0f, 1.0f, 1.0f },
				{ -1.0f, 1.0f, 1.0f } };              // 頂点座標値 
		int faces[][4] = { { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 4, 5, 6, 7 },
				{ 0, 4, 7, 3 }, { 0, 1, 5, 4 }, { 0, 3, 2, 1 } };
		// 各面の頂点番号列
		GLfloat colors[][3] = { { 0.0f, 1.0f, 1.0f }, { 1.0f, 0.0f, 1.0f },
				{ 1.0f, 1.0f, 0.0f }, { 0.0f, 0.5f, 0.5f },
				{ 0.5f, 0.0f, 0.5f }, { 0.5f, 0.5f, 0.0f } };
		// 各面の描画色
		glBegin(GL_QUADS);                  // 四角形描画開始
		for (int i = 0; i < 6; i++) {
		//	glBegin(GL_QUAD_STRIP);                  // 四角形描画開始
			glColor3fv(colors[i]);
			for (int j = 0; j < 4; j++)
				glVertex3fv(vertices[faces[i][j]]);
		//	glEnd();
		}                                         // 四角形頂点列の指定
		glEnd();                               // 四角形描画終了
}

void Object::Show() 
{
	Cube();
}
