// Simple B-Spline example for 汪艺
// 
// Peter Salvi, 2007
// 
// Time-stamp: <2009.06.08., 21:43:53 (salvi)>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

#include <GL/glut.h>

#include "interpolate.hh"

// ***********
// * Options *
// ***********

int degree = 3;
int const resolution = 100;	// number of points in the curve
int const width = 800;
int const height = 600;
int const tolerance = 5;	// largest allowed hor/ver deviation in pixels

// IDs for the menu elements / keyboard shortcuts
enum MenuCommand { MENU_PRINT, MENU_RESET, MENU_QUIT };

// ********************
// * Global variables *
// ********************

PointVector points;	        // points to be interpolated
DoubleVector knots;		// knot vector
PointVector cpts;		// control points
int dragging = -1;		// -1 means no dragging, otherwise the point#
struct SavedPoints {
  bool save;
  DoubleVector params;
  PointVector points;
} saved_points;

// ***********
// * Display *
// ***********

Point getObjectCoordinates(int x, int y)
{
  double model[16], proj[16];
  int view[4];
  double rx, ry, rz;
  glGetDoublev(GL_MODELVIEW_MATRIX, model);
  glGetDoublev(GL_PROJECTION_MATRIX, proj);
  glGetIntegerv(GL_VIEWPORT, view);
  gluUnProject(x, view[3] - y, 0, model, proj, view, &rx, &ry, &rz);
  return Point(rx, ry);
}

Point getWindowCoordinates(Point p)
{
  double model[16], proj[16];
  int view[4];
  double rx, ry, rz;
  glGetDoublev(GL_MODELVIEW_MATRIX, model);
  glGetDoublev(GL_PROJECTION_MATRIX, proj);
  glGetIntegerv(GL_VIEWPORT, view);
  gluProject(p.x, p.y, p.z, model, proj, view, &rx, &ry, &rz);
  return Point(rx, view[3] - ry);
}

Point interpolate(Point a, double t, Point b)
{
  a.x += (b.x - a.x) * t;
  a.y += (b.y - a.y) * t;
  a.z += (b.z - a.z) * t;
  return a;
}

int findSpan(double t)
{
  int r;
  int const n = cpts.size();
  if(t == knots[n])
    return n - 1;
  for(r = 0; r < knots.size(); ++r)
    if(knots[r] <= t && t < knots[r + 1])
      break;
  return r;
}

Point deBoor(double t)
{
  int r = findSpan(t);
  PointVector d = cpts;

  for(int j = 1; j <= degree; ++j) {
    for(int i = r; i >= r - degree + j; --i) {
      double const alpha =
	(t - knots[i]) / (knots[i + degree + 1 - j] - knots[i]);
      d[i] = interpolate(d[i - 1], alpha, d[i]);
    }
  }

  return d[r];
}

void drawCurve()
{
  // knot vector
  Point min(-0.1, -3.0, 0.0), max(5.0, -3.0, 0.0);
  glColor3d(0.0, 0.0, 1.0);
  glBegin(GL_LINES);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(max.x, max.y, max.z);
  glEnd();
  glPointSize(6.0);
  glBegin(GL_POINTS);
  for(int i = 0, ie = knots.size(); i != ie; ++i) {
    Point const tmp = interpolate(min, knots[i], max);
    glColor3d(0.0, 1.0, 1.0);
    glVertex3d(tmp.x, tmp.y, tmp.z);
  }
  glEnd();

  // interpolated point boxes
  glColor3d(0.0, 1.0, 1.0);
  glPointSize(10.0);
  glBegin(GL_POINTS);
  for(PointVector::const_iterator i = points.begin(); i != points.end(); ++i)
    glVertex3d(i->x, i->y, i->z);
  glEnd();

  if(cpts.empty())
    return;

  // control polygon
  glColor3d(1.0, 1.0, 0.0);
  glBegin(GL_LINE_STRIP);
  for(PointVector::const_iterator i = cpts.begin(); i != cpts.end(); ++i)
    glVertex3d(i->x, i->y, i->z);
  glEnd();

  // control point boxes
  glColor3d(1.0, 0.0, 0.0);
  glPointSize(10.0);
  glBegin(GL_POINTS);
  for(PointVector::const_iterator i = cpts.begin(); i != cpts.end(); ++i)
    glVertex3d(i->x, i->y, i->z);
  glEnd();

  // the curve itself
  glColor3d(0.0, 1.0, 0.0);
  glBegin(GL_LINE_STRIP);
  for(int i = 0; i < resolution; ++i) {
    double const t = (double)i / ((double)resolution - 1.0);
    if(knots[degree] <= t && t <= knots[knots.size() - degree - 1]) {
      Point p = deBoor(t);
      if(saved_points.save) {
	saved_points.params.push_back(t);
	saved_points.points.push_back(p);
      }
      glVertex3d(p.x, p.y, p.z);
    }
  }
  glEnd();
}

void display()
{
  glMatrixMode(GL_MODELVIEW);
  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  glScaled(0.341, 0.455, 1.0);
  glTranslated(-2.45, 1.0, 0.0);
  drawCurve();
  glutSwapBuffers();
}

// ******************
// * User interface *
// ******************

void uniform_knots(int n)
{
  knots.clear();
  for(int i = 0; i < n; ++i)
    knots.push_back((double)i / (double)(n - 1.0));
}

void reset()
{
  saved_points.save = false;
  points.clear();
  cpts.clear();
  uniform_knots(2);
}

void printData()
{
  std::cout << "(make-bspline-curve " << degree << " '( ";
  for(DoubleVector::const_iterator i = knots.begin(); i != knots.end(); ++i)
    std::cout << *i << " ";
  std::cout << ") '( ";
  for(PointVector::const_iterator i = cpts.begin(); i != cpts.end(); ++i)
    std::cout << "(" << i->x << " " << i->y << ") ";
  std::cout << "))" << std::endl;
}

void executeCommand(int command)
{
  switch(static_cast<MenuCommand>(command)) {
  case MENU_PRINT: printData(); break;
  case MENU_RESET: reset(); break;
  case MENU_QUIT: exit(0);
  }
  glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y)
{
  switch(tolower(key)) {
  case 'p' : executeCommand(MENU_PRINT); break;
  case 'r' : executeCommand(MENU_RESET); break;
  case 'q' : executeCommand(MENU_QUIT); break;
  }
}

int nearestPoint(PointVector v, int x, int y)
{
  // Returns the index of the point in v near to (x, y),
  // and -1 if no point is near enough
  int result = -1;

  for(size_t i = 0, ie = v.size(); i < ie; ++i) {
    Point p = getWindowCoordinates(v[i]);
    if(std::fabs(p.x - x) < tolerance && fabs(p.y - y) < tolerance)
      result = i;
  }

  return result;
}

void mouseButton(int button, int state, int x, int y)
{
  if(button == GLUT_LEFT_BUTTON && state == GLUT_UP)
    dragging = -1;
  else if(button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
    if((dragging = nearestPoint(points, x, y)) == -1) {
      // insert new point
      Point p = getObjectCoordinates(x, y);
      points.push_back(p);
      if(points.size() > degree) {
	BSpline bsp = BSplineInterpolate(points, degree);
	cpts = bsp.points;
	knots = bsp.knots;
      }
      glutPostRedisplay();
    }
  } else if(button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN) {
    saved_points.save = true;
    display();
    saved_points.save = false;
    int nearest = nearestPoint(saved_points.points, x, y);
    if(nearest >= 0)
      std::cout << "Selected: " << saved_points.params[nearest] << std::endl;
  }
}

void mouseMotion(int x, int y)
{
  if(dragging < 0)
    return;
  points[dragging] = getObjectCoordinates(x, y);
  if(points.size() > degree) {
    BSpline bsp = BSplineInterpolate(points, degree);
    cpts = bsp.points;
    knots = bsp.knots;
  }
  glutPostRedisplay();
}

int buildPopupMenu()
{
  int menu = glutCreateMenu(executeCommand);
  glutAddMenuEntry("Print curve data\t(p)", MENU_PRINT);
  glutAddMenuEntry("Reset\t(r)", MENU_RESET);
  glutAddMenuEntry("Quit\t(q)", MENU_QUIT);
  return menu;
}

// ******************
// * Initialization *
// ******************

int main(int argc, char *argv[])
{
  if(argc > 1)
    degree = atoi(argv[1]);

  reset();

  // glut window initialization
  glutInit(&argc, argv);
  glutInitWindowSize(width, height);
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
  glutCreateWindow("BSpline Interpolation");

  // register callbacks
  glutDisplayFunc(display);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouseButton);
  glutMotionFunc(mouseMotion);

  // create popup menu
  buildPopupMenu();
  glutAttachMenu(GLUT_RIGHT_BUTTON);

  // turn the flow of control over to glut
  glutMainLoop ();

  return 0;
}
