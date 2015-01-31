#ifndef _GLSL_PROGRAM_H
#define _GLSL_PROGRAM_H

#include "core.h"

#define MAXPROGRAM 10

class GLSLProgram
{
public:
  GLSLProgram();
  ~GLSLProgram();
  bool MakeProgram(const char* vertexsh, const char* geometrysh=0, const char* fragmentsh=0);
  bool MakeProgramFromString(std::string* vertString, std::string* fragString);
  void BeginProgram(); //glUseProgram(prog)
  void EndProgram();
  void BindTexture(unsigned int unit,GLenum, GLuint, const char* nameinshader);

  void SetUniform1f(const char*,float);
  void SetUniform3f(const char* nameinshader, float* value);
  void SetUniformMatf(const char*, float*);
  void SetUniformMat3f(const char*, float*);
  void SetUniform1uint(const char*, unsigned int);
  void SetUniform3uint(const char* nameinshader, unsigned int x,unsigned int y, unsigned int z);
  static std::string* LoadShaderToString(const char* filename);


private:
  GLuint prog;
  void _loadshader(const char* shadername, char* buffer);
  void _printlog(GLuint obj);
  static void CheckInput(const char * filename);

  //handle
  char* buffer;
  GLint len;


};


#endif
