#include "core.h"
#include "object.h"

class Renderer: public Object {
public:
	

	//for rendering
	GLuint initTFF1DTex(const char* filename);
	GLuint initFace2DTex(GLuint bfTexWidth, GLuint bfTexHeight);
	GLuint initVol3DTex(const char* filename, GLuint w, GLuint h, GLuint d);
	void initFrameBuffer(GLuint texObj, GLuint texWidth, GLuint texHeight);
	void initVBO();
	void render(GLenum cullFace);
	void checkFramebufferStatus();
	void rcSetUinforms();

	virtual void Reset();
	virtual void Show();
	void drawBox(GLenum glFaces);

private:
	//for rendering
	float g_stepSize;
	GLuint g_vao;
	GLuint g_frameBuffer;
	// transfer function
	GLuint g_tffTexObj;
	GLuint g_bfTexObj;
	GLuint g_volTexObj;

};
