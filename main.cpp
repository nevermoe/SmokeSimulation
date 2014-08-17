#include "core.h"
#include "controller.h"
#include "flip.h"

#define WINDOW_NAME		 "FLIP"
Controller* g_controller = NULL;

int main(int argc, char **argv) {
	//Object* object = new Object;
	Object* object = new Flip;
	::g_controller = new Controller(argc, argv, WINDOW_NAME);
	::g_controller->RegisterObject(object);
	::g_controller->BeginLoop();
	return 0;
}
