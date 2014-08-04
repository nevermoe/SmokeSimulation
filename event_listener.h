class EventListener {
public:
	virtual void MouseButton(GLFWwindow *window, int button,int action,int mods) = 0;
	virtual void MouseMotion(GLFWwindow *window, double nx, double ny) = 0;
	virtual void Keyboard(GLFWwindow * window, int key, int scancode, int action, int mods)	= 0;
	virtual void Resize(GLFWwindow *window, int x, int y) = 0;
};
