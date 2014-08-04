TARGET = main
CXXFLAGS = -I./include #-O3
CXXEXTRA_FLAGS = -Wl,-rpath,/usr/local/lib/ #not used
CC = g++
LIBS_PATH = -L/usr/local/lib/OpenMesh
LIBS_PATH += -L/usr/local/lib
LIBS = -lGLEW -lGL -lGLU -lglfw3 -lX11  -lXrandr -lpthread -lXi
LIBS += -lpng -lz
#LIBS += -lOpenMeshCored -lOpenMeshToolsd
LIBS += -lOpenMeshCore -lOpenMeshTools

OBJS = camera.o 
OBJS += object.o 
OBJS += arcball.o
OBJS += controller.o
OBJS += main.o

all: $(OBJS)
	$(CC) $(OBJS) $(LIBS_PATH) $(LIBS) -o $(TARGET) #必须把库的链接放在源文件之后！
	@echo ""
	@echo "\033[31mCompiling Succeed!"
	@echo "\033[32m"
	@echo "\033[0m"

clean:
	rm -f $(OBJS)
	rm -f $(TARGET)
