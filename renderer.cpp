#include "renderer.h"



#define GL_ERROR() checkForOpenGLError(__FILE__, __LINE__)
int checkForOpenGLError(const char* file, int line)
{
	// return 1 if an OpenGL error occured, 0 otherwise.
	GLenum glErr;
	int retCode = 0;

	glErr = glGetError();
	while(glErr != GL_NO_ERROR) {
		std::cerr << "glError in file " << file
			<< "@line " << line << gluErrorString(glErr) << std::endl;
		retCode = 1;
		exit(EXIT_FAILURE);
	}
	return retCode;
}

// init the vertex buffer object
void Renderer::initVBO()
{
    GLfloat vertices[24] = {
	0.0, 0.0, 0.0,
	0.0, 0.0, 1.0,
	0.0, 1.0, 0.0,
	0.0, 1.0, 1.0,
	1.0, 0.0, 0.0,
	1.0, 0.0, 1.0,
	1.0, 1.0, 0.0,
	1.0, 1.0, 1.0
    };
// draw the six faces of the boundbox by drawwing triangles
// draw it contra-clockwise
// front: 1 5 7 3
// back: 0 2 6 4
// left??0 1 3 2
// right:7 5 4 6    
// up: 2 3 7 6
// down: 1 0 4 5
    GLuint indices[36] = {
	1,5,7,
	7,3,1,
	0,2,6,
    6,4,0,
	0,1,3,
	3,2,0,
	7,5,4,
	4,6,7,
	2,3,7,
	7,6,2,
	1,0,4,
	4,5,1
    };
    GLuint gbo[2];
    
    glGenBuffers(2, gbo);
    GLuint vertexdat = gbo[0];
    GLuint veridxdat = gbo[1];
    glBindBuffer(GL_ARRAY_BUFFER, vertexdat);
    glBufferData(GL_ARRAY_BUFFER, 24*sizeof(GLfloat), vertices, GL_STATIC_DRAW);
    // used in glDrawElement()
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, veridxdat);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 36*sizeof(GLuint), indices, GL_STATIC_DRAW);

    GLuint vao;
    glGenVertexArrays(1, &vao);
    // vao like a closure binding 3 buffer object: verlocdat vercoldat and veridxdat
    glBindVertexArray(vao);
    glEnableVertexAttribArray(0); // for vertexloc
    glEnableVertexAttribArray(1); // for vertexcol

    // the vertex location is the same as the vertex color
    glBindBuffer(GL_ARRAY_BUFFER, vertexdat);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (GLfloat *)NULL);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (GLfloat *)NULL);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, veridxdat);
    // glBindVertexArray(0);
    g_vao = vao;
}

void Renderer::Reset()
{
	Object::Reset();
	initVBO();
	g_tffTexObj = initTFF1DTex("ttf.dat");
	g_bfTexObj = initFace2DTex(winX_, winY_);
	g_volTexObj = initVol3DTex("head256.raw", 256, 256, 225);
	GL_ERROR();
	initFrameBuffer(g_bfTexObj, winX_, winY_);
	GL_ERROR();
}

// init the 1 dimentional texture for transfer function
GLuint Renderer::initTFF1DTex(const char* filename)
{
    // read in the user defined data of transfer function
	std::ifstream inFile(filename, std::ifstream::in);
	if (!inFile) {
		std::cerr << "Error openning file: " << filename << std::endl;
		exit(EXIT_FAILURE);
    }
    
    const int MAX_CNT = 10000;
    GLubyte *tff = (GLubyte *) calloc(MAX_CNT, sizeof(GLubyte));
    inFile.read(reinterpret_cast<char *>(tff), MAX_CNT);
	if (inFile.eof()) {
		size_t bytecnt = inFile.gcount();
		*(tff + bytecnt) = '\0';
		std::cout << "bytecnt " << bytecnt << std::endl;
	}
	else if(inFile.fail()) {
		std::cout << filename << "read failed " << std::endl;
	}
	else {
		std::cout << filename << "is too large" << std::endl;
	}    
	GLuint tff1DTex;
	glGenTextures(1, &tff1DTex);
    glBindTexture(GL_TEXTURE_1D, tff1DTex);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA8, 256, 0, GL_RGBA, GL_UNSIGNED_BYTE, tff);
    free(tff);    
    return tff1DTex;
}

// init the 2D texture for render backface 'bf' stands for backface
GLuint Renderer::initFace2DTex(GLuint bfTexWidth, GLuint bfTexHeight)
{
    GLuint backFace2DTex;
    glGenTextures(1, &backFace2DTex);
    glBindTexture(GL_TEXTURE_2D, backFace2DTex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, bfTexWidth, bfTexHeight, 0, GL_RGBA, GL_FLOAT, NULL);
    return backFace2DTex;
}
// init 3D texture to store the volume data used fo ray casting
GLuint Renderer::initVol3DTex(const char* filename, GLuint w, GLuint h, GLuint d)
{
    
    FILE *fp;
    size_t size = w * h * d;
    GLubyte *data = new GLubyte[size];			  // 8bit
    if (!(fp = fopen(filename, "rb"))) {
		std::cerr << "Error: opening .raw file failed" << std::endl;
        exit(EXIT_FAILURE);
    }
    else
    {
		std::cerr << "OK: open .raw file succeed" << std::endl;
    }
    if ( fread(data, sizeof(char), size, fp)!= size) 
    {
		std::cerr << "Error: read .raw file failed" << std::endl;
        exit(1);
    }
    else
    {
		std::cerr << "OK: read .raw file succeed" << std::endl;
    }
    fclose(fp);

    glGenTextures(1, &g_volTexObj);
    // bind 3D texture target
    glBindTexture(GL_TEXTURE_3D, g_volTexObj);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);	
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_REPEAT);
    // pixel transfer happens here from client to OpenGL server
    glPixelStorei(GL_UNPACK_ALIGNMENT,1);
    glTexImage3D(GL_TEXTURE_3D, 0, GL_INTENSITY, w, h, d, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE,data);

    delete []data;
	std::cerr << "volume texture created" << std::endl;
    return g_volTexObj;
}

void Renderer::checkFramebufferStatus()
{
    GLenum complete = glCheckFramebufferStatus(GL_FRAMEBUFFER);
    if (complete != GL_FRAMEBUFFER_COMPLETE) {
		std::cerr << "framebuffer is not complete" << std::endl;
		exit(EXIT_FAILURE);
    }
}

void Renderer::drawBox(GLenum glFaces)
{
    glEnable(GL_CULL_FACE);
    glCullFace(glFaces);
    glBindVertexArray(g_vao);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, (GLuint *)NULL);
    glDisable(GL_CULL_FACE);
}

void Renderer::render(GLenum cullFace)
{
}

// init the framebuffer, the only framebuffer used in this program
void Renderer::initFrameBuffer(GLuint texObj, GLuint texWidth, GLuint texHeight)
{
    // create a depth buffer for our framebuffer
    GLuint depthBuffer;
    glGenRenderbuffers(1, &depthBuffer);
    glBindRenderbuffer(GL_RENDERBUFFER, depthBuffer);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, texWidth, texHeight);

    // attach the texture and the depth buffer to the framebuffer
    glGenFramebuffers(1, &g_frameBuffer);
    glBindFramebuffer(GL_FRAMEBUFFER, g_frameBuffer);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texObj, 0);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthBuffer);
    checkFramebufferStatus();
    glEnable(GL_DEPTH_TEST);    
}

void Renderer::rcSetUinforms()
{
	// setting uniforms such as
	// ScreenSize 
	// StepSize
	// TransferFunc
	// ExitPoints i.e. the backface, the backface hold the ExitPoints of ray casting
	// VolumeTex the texture that hold the volume data i.e. head256.raw
	GLint screenSizeLoc = glGetUniformLocation(shaderProg_, "ScreenSize");
	if (screenSizeLoc >= 0) {
		glUniform2f(screenSizeLoc, (float)winX_, (float)winY_);
	}
	else {
		std::cerr << "ScreenSize is not bind to the uniform" << std::endl;
	}

	GLint stepSizeLoc = glGetUniformLocation(shaderProg_, "StepSize");
	GL_ERROR();
	if (stepSizeLoc >= 0) {
		glUniform1f(stepSizeLoc, g_stepSize);
	}
	else {
		std::cerr << "StepSize is not bind to the uniform" << std::endl;
	}
	GL_ERROR();
	GLint transferFuncLoc = glGetUniformLocation(shaderProg_, "TransferFunc");
	if (transferFuncLoc >= 0) {
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_1D, g_tffTexObj);
		glUniform1i(transferFuncLoc, 0);
	}
	else {
		std::cerr << "TransferFunc"
			<< "is not bind to the uniform"
			<< std::endl;
	}
	GL_ERROR();    
	GLint backFaceLoc = glGetUniformLocation(shaderProg_, "ExitPoints");
	if (backFaceLoc >= 0) {
		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, g_bfTexObj);
		glUniform1i(backFaceLoc, 1);
	}
	else {
		std::cerr << "ExitPoints"
			<< "is not bind to the uniform"
			<< std::endl;
	}
	GL_ERROR();    
	GLint volumeLoc = glGetUniformLocation(shaderProg_, "VolumeTex");
	if (volumeLoc >= 0) {
		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_3D, g_volTexObj);
		glUniform1i(volumeLoc, 2);
	}
	else {
		std::cerr << "VolumeTex"
			<< "is not bind to the uniform"
			<< std::endl;
	}

}

#undef ENABLE_SHADER

void Renderer::Show() 
{
	//Cube();
	//first pass
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, g_frameBuffer);
	//RebindShader(_hShaders[0], _hShaders[1]);
    GL_ERROR();
    drawBox(GL_FRONT);
    //glUseProgram(0);

    GL_ERROR(); 

	//second pass
	//RebindShader(_hShaders[2], _hShaders[3]);
	//rcSetUinforms();
    drawBox(GL_BACK);

    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    //glUseProgram(0);
    GL_ERROR(); 
}
