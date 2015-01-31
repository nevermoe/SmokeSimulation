#include "core.h"
#include "controller.h"
#include "fluid.h"

#define WINDOW_NAME		 "Smoke3D"

//#define VERTEX_PROG_NAME1   "backface.vert"
//#define FRAGMENT_PROG_NAME1 "backface.frag"
//#define VERTEX_PROG_NAME2   "raycasting.vert"
//#define FRAGMENT_PROG_NAME2 "raycasting.frag"

Controller* g_controller = NULL;

int main(int argc, char **argv) {
#if 1
	glutInit(&argc, argv);
#endif
	::g_controller = new Controller(argc, argv, WINDOW_NAME);
	Object* object = new Fluid;

#if 0
	//enable shader
	object->RegisterShader(VERTEX_PROG_NAME1, GL_VERTEX_SHADER);
	object->RegisterShader(FRAGMENT_PROG_NAME1, GL_FRAGMENT_SHADER);
	object->RegisterShader(VERTEX_PROG_NAME2, GL_VERTEX_SHADER);
	object->RegisterShader(FRAGMENT_PROG_NAME2, GL_FRAGMENT_SHADER);
	object->EnableShader();
#endif

	::g_controller->InitCamera();
	::g_controller->RegisterObject(object);
	::g_controller->BeginLoop();
	return 0;
}
