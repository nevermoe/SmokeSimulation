#ifndef _GL_TEXTURE_H
#define  _GL_TEXTURE_H

#include "core.h"

class GLTexture
{
public:
	//GL texture element byte size can be different from the input element byte size
	//GLTexture();
	GLTexture(int width, int height=0, int depth=0, GLenum elementFormat = GL_RGBA, 
			GLint  internalFormat=GL_RGBA8, GLenum dataType = GL_UNSIGNED_BYTE, GLenum filterType=GL_LINEAR, GLenum borderType=GL_CLAMP_TO_EDGE); 
	~GLTexture();

	void InitTexture(int width, int height=0, int depth=0, GLenum elementFormat = GL_RGBA, 
			GLint   internalFormat=GL_RGBA8, GLenum filterType=GL_LINEAR, GLenum borderType=GL_CLAMP_TO_EDGE); 

	//accessors 
	//set xx
	GLuint GetTextureID();
	GLuint GetHeight();
	GLuint GetWidth();
	GLuint GetDepth();

	GLenum GetTextureType();
	//setxx
	void SetFilterType(GLenum);

	//bind unbind
	void Bind();
	void UnBind();

	//read to memory from File
	bool ReadTextureFromFile(const char* filename, GLenum elementType, int channelNum=1); //called when GL has been initialized

	//
	void LoadToGPU(void* data);	   //glTexImage1D/2D/3D and then free the CPU memory
	void LoadToGPU();

	//upload to GPU a small chunk at a time using subImage void the large memory usage 
	void LoadToGPUWithGLBufferSmallChunk();
	//void LoadToGPUkeepCPU(); //keep CPU memory unreleased

	//subTexture To GPU
	void SubloadToGPU(int offsetX, int offsetY, int offsetZ, int sizeX, int sizeY, int sizeZ, void* data, GLenum elementType);

	/////////texture buffer


protected:
	//failed when return NULL pointer
	void* loadRawFile(int, int*);
	void* loadPPMFile();

protected:
	std::string _fileName;
	std::string _fileType;
	void* _data;
	int _dataSize;


	GLenum _filterType;
	GLenum _borderType;
	GLenum _textureType;
	GLenum _elementFormat;
	unsigned int _internalFormat;	  //1 2 4 8
	GLenum _elementType; //element type in the CPU side
	unsigned int _elementByteSize;

	GLuint _tex; //texture object ID
	GLuint _dim[3]; // width, height, depth;

	bool _isBufferAllocated;
};

#endif
