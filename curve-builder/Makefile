all: bspline-interpolation bspline

FLAGS=-g
ifeq ($(shell uname),Linux)
	LIBS=-lglut
else
	LIBS=-lGL -lGLU -lglut
endif

bspline-interpolation: bspline-interpolation.o interpolate.o lu.o
	g++ $(FLAGS) -o $@ $^ $(LIBS)

bspline: bspline.cc
	g++ $(FLAGS) -o $@ $^ $(LIBS)
