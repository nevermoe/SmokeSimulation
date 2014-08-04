#include "core.h"

Controller* g_controller = NULL;

int main(int argc, char **argv) {
	Object* object = new Object;
	::g_controller = new Controller(argc, argv);
	::g_controller->RegisterObject(object);
	::g_controller->BeginLoop();
	return 0;
}
