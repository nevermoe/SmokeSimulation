#include "GLSL_program.h"

GLSLProgram::GLSLProgram()
{
}

GLSLProgram::~GLSLProgram()
{
	if(buffer)
		delete buffer;
}
void GLSLProgram::CheckInput(const char * filename)
{
	FILE *fp = fopen(filename, "r");
	if(!fp) {
		std::cout<<"Can't open - "<<filename;
		int ret = system("pause");
		exit(0);
	}

	else {
		fclose(fp);
		return;
	}
}

std::string* GLSLProgram::LoadShaderToString(const char* filename)
{
	CheckInput(filename);
	int size = 100000;
	char *cbuffer = new GLchar[size];
	memset(cbuffer, 0, sizeof(GLchar)*size);

	//load from file
	std::ifstream shader;
	shader.open(filename, std::ios::in);
	shader.seekg(0,std::ios::end);
	int length = int(shader.tellg());
	shader.seekg(0, std::ios::beg);		 
	shader.read(cbuffer,length);
	shader.close();
	std::string* str = new std::string(cbuffer);

	return str;
}

bool GLSLProgram::MakeProgramFromString(std::string* vertString, std::string* fragString)
{
	GLuint vertexglsl, fragmentglsl;
	//vertex
	vertexglsl = glCreateShader(GL_VERTEX_SHADER);
	GLchar const* vertC =	 vertString->c_str();
	GLint const vertSize = (GLint const)vertString->size();
	glShaderSource(vertexglsl,1, &vertC, &vertSize);
	glCompileShader(vertexglsl);
	_printlog(vertexglsl);
	std::cout << "\nvertexshader compiled!\n" << std::endl;

	//fragment
	fragmentglsl = glCreateShader(GL_FRAGMENT_SHADER);	
	GLchar const* fragC =	 fragString->c_str();
	GLint const fragSize = (GLint const )fragString->size();
	glShaderSource(fragmentglsl,1, &fragC,  &fragSize);
	glCompileShader(fragmentglsl);
	_printlog(fragmentglsl);
	std::cout << " fragmentshader compiled!\n" << std::endl;

	///clear up
	delete vertString;
	delete fragString;			 

	///////////////////////////////////////
	///verify 
	GLint vertcompiled, fragmentcompiled, linked;

	glGetShaderiv(vertexglsl,GL_COMPILE_STATUS,&vertcompiled);
	glGetShaderiv(fragmentglsl,GL_COMPILE_STATUS,&fragmentcompiled);

	if(!vertcompiled)
		return false;

	//create program
	prog = glCreateProgram();

	//Setup GS
	//glProgramParameteriEXT(prog, GL_GEOMETRY_VERTICES_OUT_EXT, 1000);
	//glProgramParameteriEXT(prog, GL_GEOMETRY_INPUT_TYPE_EXT, GL_POINTS);
	//glProgramParameteriEXT(prog, GL_GEOMETRY_OUTPUT_TYPE_EXT, GL_TRIANGLE_STRIP);

	//Attach Shader
	glAttachShader(prog, vertexglsl);
	glAttachShader(prog, fragmentglsl);

	//Link program
	glLinkProgram(prog);

	//get program info
	glGetProgramiv(prog, GL_LINK_STATUS, &linked);
	_printlog(prog);

	if(!linked)
	{
		std::cout << "Link failed! \n"<< std::endl;
		return false;
	}
	return true;
}


bool GLSLProgram::MakeProgram(const char* vertexsh, const char* fragmentsh, const char* geometrysh )
{
	GLuint vertexglsl, fragmentglsl, geometryglsl;
	vertexglsl = glCreateShader(GL_VERTEX_SHADER);
	fragmentglsl = glCreateShader(GL_FRAGMENT_SHADER);
	geometryglsl = glCreateShader(GL_GEOMETRY_SHADER_EXT);

	//TODO make the length of the char adapted to the length of the shader source file
	buffer = new GLchar[100000];
	for(int i = 0; i< 100000; i++)
		buffer[i] = ' ';
	_loadshader(vertexsh, buffer);
	glShaderSource(vertexglsl,1,(GLchar const**)&buffer,&len);
	glCompileShader(vertexglsl);
	_printlog(vertexglsl);
	std::cout << "\nvertexshader compiled!\n" << std::endl;

	if (fragmentsh)
	{
		for(int i = 0; i< 100000; i++)
			buffer[i] = ' ';
		_loadshader(fragmentsh, buffer);//ftransf.glsl
		glShaderSource(fragmentglsl,1,(GLchar const**)&buffer,&len);
		glCompileShader(fragmentglsl);
		_printlog(fragmentglsl);
		std::cout << "fragmentshader compiled!\n" << std::endl;
	}

	if (geometrysh)
	{
		for(int i = 0; i< 100000; i++)
			buffer[i] = ' ';
		_loadshader(geometrysh, buffer);
		glShaderSource(geometryglsl,1,(GLchar const**)&buffer,&len);
		glCompileShader(geometryglsl);
		_printlog(geometryglsl);
		std::cout << "geometry shader compiled!\n" << std::endl;
	}

	std::cout <<"-------------"<<std::endl;

	delete[] buffer;
	GLint vertcompiled, fragmentcompiled, linked;
	glGetShaderiv(vertexglsl,GL_COMPILE_STATUS,&vertcompiled);
	glGetShaderiv(fragmentglsl,GL_COMPILE_STATUS,&fragmentcompiled);

	if(!vertcompiled)
		return false;

	//create program
	prog = glCreateProgram();

	//Setup GS
	//glProgramParameteriEXT(prog, GL_GEOMETRY_VERTICES_OUT_EXT, 1000);
	//glProgramParameteriEXT(prog, GL_GEOMETRY_INPUT_TYPE_EXT, GL_POINTS);
	//glProgramParameteriEXT(prog, GL_GEOMETRY_OUTPUT_TYPE_EXT, GL_TRIANGLE_STRIP);

	//Attach Shader
	glAttachShader(prog, vertexglsl);
	if(fragmentsh)
		glAttachShader(prog, fragmentglsl);
	if(geometrysh)
		glAttachShader(prog, geometryglsl);


	//Link program
	glLinkProgram(prog);

	//get program info
	glGetProgramiv(prog, GL_LINK_STATUS, &linked);
	_printlog(prog);

	if(!linked)
	{
		std::cout << "Link failed! \n"<< std::endl;
		return false;
	}
	return true;
}

void GLSLProgram::BindTexture(unsigned int unit, GLenum texturetype, GLuint textureid, const char* nameinshader)
{
	GLuint loc = glGetUniformLocation(prog, nameinshader);
	//printf("loc: %d, texID:%d, name:%s\n", loc, textureid, nameinshader);
	glUniform1i(loc,unit);	

	glActiveTexture(GL_TEXTURE0 + unit);
	glBindTexture(texturetype, textureid);
}

void GLSLProgram::SetUniform1f(const char* nameinshader, float value)
{
	GLuint loc = glGetUniformLocation(prog, nameinshader);
	glUniform1f(loc, value);
}

void GLSLProgram::SetUniform3f(const char* nameinshader, float* value)
{
	GLuint loc = glGetUniformLocation(prog, nameinshader);
	glUniform3f(loc, value[0], value[1], value[2]);
}

void GLSLProgram::SetUniformMatf(const char* nameinshader, float* value)
{
	GLuint loc = glGetUniformLocation(prog, nameinshader);
	glUniformMatrix4fv(loc, 1, 0, value);
}

void GLSLProgram::SetUniformMat3f(const char* nameinshader, float* value)
{
	GLuint loc = glGetUniformLocation(prog, nameinshader);
	glUniformMatrix3fv(loc, 1, 0, value);
}

void GLSLProgram::SetUniform1uint(const char* nameinshader, unsigned int value)
{
	GLuint loc = glGetUniformLocation(prog, nameinshader);
	glUniform1uiEXT(loc, value);
}

void GLSLProgram::SetUniform3uint(const char* nameinshader, unsigned int x,unsigned int y, unsigned int z)
{
	GLuint loc = glGetUniformLocation(prog, nameinshader);
	glUniform3uiEXT(loc, x, y, z);
}


void GLSLProgram::BeginProgram()
{
	glUseProgram(prog);
}

void GLSLProgram::EndProgram()
{
	glUseProgram(0);
}

//==========private member===========

void GLSLProgram::_loadshader(const char* filename, char* buffer)
{
	CheckInput(filename);


	std::ifstream shader;
	shader.open(filename, std::ios::in);
	shader.seekg(0,std::ios::end);
	len = int(shader.tellg());
	shader.seekg(0, std::ios::beg);

	shader.read(buffer, len);
	shader.close();

}

void GLSLProgram::_printlog(GLuint obj)
{
	int infologLength = 0;
	int maxLength;

	if(glIsShader(obj))
		glGetShaderiv(obj,GL_INFO_LOG_LENGTH,&maxLength);
	else
		glGetProgramiv(obj,GL_INFO_LOG_LENGTH,&maxLength);

	char *infoLog = new char[maxLength];

	if (glIsShader(obj))
		glGetShaderInfoLog(obj, maxLength, &infologLength, infoLog);
	else
		glGetProgramInfoLog(obj, maxLength, &infologLength, infoLog);

	if (infologLength > 0)
		printf("%s\n",infoLog);
}
