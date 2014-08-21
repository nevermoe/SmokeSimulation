#ifndef FILE_OPERATOR_H
#define FILE_OPERATOR_H

#include "core.h"
#include "flip.h"
class Flip;

class FileOperator {
public:
	FileOperator(const char *dirPath);
	void SetDirPath(const char *dirPath);
	void MakeDir();
	void WriteFlipVel(std::string fileName, Flip *object);
private:
	std::string dirPath_;
	std::string fileName_;
	std::string fullPath_;
	std::fstream file;
};

#endif
