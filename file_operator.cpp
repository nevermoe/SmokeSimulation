#include "core.h"
#include "file_operator.h"

FileOperator::FileOperator(const char *dirPath)
{
	SetDirPath(dirPath);
}

void FileOperator::SetDirPath(const char *dirPath)
{
	dirPath_ = dirPath;
	MakeDir();
}

void FileOperator::MakeDir()
{
	std::string command = "mkdir -p ";
	command += dirPath_;
	int ret = system(command.c_str());
	if(ret == -1) {
		std::cerr << "mkdir failed" << std::endl;
		return;
	}
}

void FileOperator::WriteFlipVel(std::string fileName, Flip *object)
{
	fullPath_ = dirPath_ + fileName;
	file.open(fullPath_.c_str(), std::ios::out);
	if(!file) {
		std::cerr << "Open file " << fullPath_ << " failed!" << std::endl;
		return;
	}

	Grid *grid = object->GetGrid();

	//start writing
	for(int i = 1 ; i < grid->nX_ ; i++)
		for(int j = 1 ; j < grid->nY_-1 ; j++)
			for(int k = 1 ; k < grid->nZ_-1 ; k++)
				file << grid->velX_[i][j][k] << " ";
	file << std::endl;

	for(int i=1;i<grid->nX_-1;i++)
		for(int j=1;j<grid->nY_;j++)
			for(int k=1;k<grid->nZ_-1;k++) 
				file << grid->velY_[i][j][k] << " ";
	file << std::endl;

	for(int i=1;i<grid->nX_-1;i++)
		for(int j=1;j<grid->nY_-1;j++)
			for(int k=1;k<grid->nZ_;k++) 
				file << grid->velZ_[i][j][k] << " ";
	file << std::endl;

	file.close();
}
