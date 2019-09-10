all: bspline-interpolation bspline

FLAGS=-g
LIBS=-lGL -lGLU -lglut

bspline-interpolation: bspline-interpolation.o interpolate.o lu.o
	g++ $(FLAGS) -o $@ $^ $(LIBS)

bspline: bspline.cc
	g++ $(FLAGS) -o $@ $^ $(LIBS)
