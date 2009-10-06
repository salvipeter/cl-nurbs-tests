// Simple B-Spline example for 汪艺
// 
// Peter Salvi, 2007
// 
// Time-stamp: <2009.10.06., 11:26:12 (psalvi)>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

#include <GL/glut.h>

// ***********
// * Options *
// ***********

int degree = 3;
int const resolution = 100;	// number of points in the curve
int const width = 800;
int const height = 600;

// *********
// * Types *
// *********

struct Point {
  double x;
  double y;
  double z;
  Point(double a, double b, double c) { x = a; y = b; z = c; }
};

typedef std::vector<Point> PointVector;
typedef std::vector<double> ValueVector;

// IDs for the menu elements / keyboard shortcuts
enum MenuCommand { MENU_EXTEND, MENU_PRINT, MENU_RESET, MENU_QUIT };

// ********************
// * Global variables *
// ********************

PointVector cpts;	        // control points
ValueVector knots;		// knot vector
int dragging = -1;		// -1 means no dragging, otherwise the cpt#
int active_knot = 0;

// ***********
// * Display *
// ***********

Point interpolate(Point a, double t, Point b)
{
  a.x += (b.x - a.x) * t;
  a.y += (b.y - a.y) * t;
  a.z += (b.z - a.z) * t;
  return a;
}

inline Point scalePoint(Point const &p, double const scaling)
{
  return Point(p.x * scaling, p.y * scaling, p.z * scaling);
}

inline Point addPoints(Point const &p1, Point const &p2)
{
  return Point(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
}

double bsp_basis(int i, int n, double t)
{
  if(n == 0) {
    if(knots[i] <= t && t < knots[i+1])
      return 1.0;
    else
      return 0.0;
  }
  double const alpha = (t - knots[i]) / (knots[i + n] - knots[i]);
  double const beta = (knots[i + n + 1] - t) / (knots[i + n + 1] - knots[i + 1]);
  return alpha * bsp_basis(i, n - 1, t) + beta * bsp_basis(i + 1, n - 1, t);
}

Point direct_bspline(double t)
{
  Point result(0.0, 0.0, 0.0);
  for(int i = 0; i < cpts.size(); ++i)
    result = addPoints(result, scalePoint(cpts[i], bsp_basis(i, degree, t)));
  return result;
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
    if(i == active_knot)
      glColor3d(1.0, 1.0, 1.0);
    else
      glColor3d(0.0, 1.0, 1.0);
    glVertex3d(tmp.x, tmp.y, tmp.z);
  }
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
//       Point p = direct_bspline(t);
      Point p = deBoor(t);
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
  if(n < 2 * degree) {
    for(int i = 0; i < n; ++i)
      knots.push_back((double)i / (double)(n - 1.0));
  } else {
    knots.insert(knots.end(), degree, 0.0);
    int k = n - 2 * degree;
    for(int i = 0; i < k; ++i)
      knots.push_back((double)i / (double)(k - 1.0));
    knots.insert(knots.end(), degree, 1.0);
  }
}

void reset()
{
  cpts.clear();
  uniform_knots(degree + 1);
  active_knot = 0;
}

void printData()
{
  std::cout << "(make-bspline-curve " << degree << " '( ";
  for(ValueVector::const_iterator i = knots.begin(); i != knots.end(); ++i)
    std::cout << *i << " ";
  std::cout << ") '( ";
  for(PointVector::const_iterator i = cpts.begin(); i != cpts.end(); ++i)
    std::cout << "(" << i->x << " " << i->y << ") ";
  std::cout << "))" << std::endl;
}

void extendSides()
{
  int const n = knots.size();
  if(n < 2 * (degree + 1))
    return;
  for(int i = 0; i <= degree; ++i) {
    knots[i] = 0.0;
    knots[n - i - 1] = 1.0;
  }
  active_knot = std::min(std::max(active_knot, degree), n - degree - 1);
}

void executeCommand(int command)
{
  switch(static_cast<MenuCommand>(command)) {
  case MENU_EXTEND: extendSides(); break;
  case MENU_PRINT: printData(); break;
  case MENU_RESET: reset(); break;
  case MENU_QUIT: exit(0);
  }
  glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y)
{
  switch(tolower(key)) {
  case 'x' : executeCommand(MENU_EXTEND); break;
  case 'p' : executeCommand(MENU_PRINT); break;    
  case 'l' :
    if(active_knot < knots.size() - 1) {
      ++active_knot;
      glutPostRedisplay();
    } break;
  case 'h' :
    if(active_knot > 0) {
      --active_knot;
      glutPostRedisplay();
    } break;
  case 'j' :
    if(knots[active_knot] >= 0.01) {
      knots[active_knot] -= 0.01;
      std::sort(knots.begin(), knots.end());
      glutPostRedisplay();
    } break;
  case 'k' :
    if(knots[active_knot] <= 0.99) {
      knots[active_knot] += 0.01;
      std::sort(knots.begin(), knots.end());
      glutPostRedisplay();
    } break;
  case 'r' : executeCommand(MENU_RESET); break;
  case 'q' : executeCommand(MENU_QUIT); break;
  }
}

Point getObjectCoordinates(int x, int y)
{
  double model[16], proj[16];
  int view[4];
  double rx, ry, rz;
  glGetDoublev(GL_MODELVIEW_MATRIX, model);
  glGetDoublev(GL_PROJECTION_MATRIX, proj);
  glGetIntegerv(GL_VIEWPORT, view);
  gluUnProject(x, view[3] - y, 0, model, proj, view, &rx, &ry, &rz);
  return Point(rx, ry, rz);
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
  return Point(rx, view[3] - ry, rz);
}

void mouseButton(int button, int state, int x, int y)
{
  if(button == GLUT_LEFT_BUTTON && state == GLUT_UP)
    dragging = -1;
  else if(button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
    for(size_t i = 0, ie = cpts.size(); i < ie; ++i) {
      Point p = getWindowCoordinates(cpts[i]);
      if(std::fabs(p.x - x) < 5 && fabs(p.y - y) < 5)
	dragging = i;
    }
    if(dragging == -1) {	// insert new point
      Point p = getObjectCoordinates(x, y);
      cpts.push_back(p);
      uniform_knots(knots.size() + 1);
      glutPostRedisplay();
    }
  }
}

void mouseMotion(int x, int y)
{
  if(dragging < 0)
    return;
  cpts[dragging] = getObjectCoordinates(x, y);
  glutPostRedisplay();
}

int buildPopupMenu()
{
  int menu = glutCreateMenu(executeCommand);
  glutAddMenuEntry("Extend\t(x)", MENU_EXTEND);
  glutAddMenuEntry("Print\t(p)", MENU_PRINT);
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
  glutCreateWindow("glut-Casteljau");

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
