#include "fx.h"
#include "fx3d.h"

class GLTestWindow : public FXMainWindow {
  FXDECLARE(GLTestWindow)

private:

  FXGLCanvas      *glcanvas;                  // GL Canvas to draw into
  FXGLVisual      *glvisual;                  // OpenGL visual

protected:
  GLTestWindow(){}

public:

  enum{
    ID_CANVAS=FXMainWindow::ID_LAST,
    ID_OPENGL
    };

  // Message handlers
  long onExpose(FXObject*,FXSelector,void*);
  long onConfigure(FXObject*,FXSelector,void*);
  long onCmdOpenGL(FXObject*,FXSelector,void*);

public:

  // GLTestWindow constructor
  GLTestWindow(FXApp* a);

  // Initialize
  void create();

  // Draw scene
  void drawScene();

  // GLTestWindow destructor
  virtual ~GLTestWindow();
  };



// Message Map
FXDEFMAP(GLTestWindow) GLTestWindowMap[]={

  //________Message_Type_________ID_____________________Message_Handler_______
  FXMAPFUNC(SEL_PAINT,     GLTestWindow::ID_CANVAS,   GLTestWindow::onExpose),
  FXMAPFUNC(SEL_CONFIGURE, GLTestWindow::ID_CANVAS,   GLTestWindow::onConfigure),
  FXMAPFUNC(SEL_COMMAND,   GLTestWindow::ID_OPENGL,   GLTestWindow::onCmdOpenGL),
  };

// Implementation
FXIMPLEMENT(GLTestWindow,FXMainWindow,GLTestWindowMap,ARRAYNUMBER(GLTestWindowMap))

// Construct a GLTestApp
GLTestWindow::GLTestWindow(FXApp* a):FXMainWindow(a,"OpenGL Test Application",NULL,NULL,DECOR_ALL,0,0,800,600){

  // A Visual to drag OpenGL
  glvisual=new FXGLVisual(getApp(),VISUAL_DOUBLEBUFFER|VISUAL_STEREO);

  // Drawing glcanvas
  glcanvas=new FXGLCanvas(this,glvisual,this,ID_CANVAS,LAYOUT_FILL_X|LAYOUT_FILL_Y|LAYOUT_TOP|LAYOUT_LEFT);
  }


// Destructor
GLTestWindow::~GLTestWindow(){
  delete glvisual;
  }



// Create and initialize
void GLTestWindow::create(){
  FXMainWindow::create();
  show(PLACEMENT_SCREEN);
  }



// Widget has been resized
long GLTestWindow::onConfigure(FXObject*,FXSelector,void*){
  if(glcanvas->makeCurrent()){
    glViewport(0,0,glcanvas->getWidth(),glcanvas->getHeight());
    glcanvas->makeNonCurrent();
    }
  return 1;
  }



// Widget needs repainting
long GLTestWindow::onExpose(FXObject*,FXSelector,void*){
  drawScene();
  return 1;
  }




// Draws a simple box using the given corners
void drawBox(GLfloat xmin, GLfloat ymin, GLfloat zmin, GLfloat xmax, GLfloat ymax, GLfloat zmax) {
  glBegin(GL_TRIANGLE_STRIP);
    glNormal3f(0.,0.,-1.);
    glVertex3f(xmin, ymin, zmin);
    glVertex3f(xmin, ymax, zmin);
    glVertex3f(xmax, ymin, zmin);
    glVertex3f(xmax, ymax, zmin);
  glEnd();

  glBegin(GL_TRIANGLE_STRIP);
    glNormal3f(1.,0.,0.);
    glVertex3f(xmax, ymin, zmin);
    glVertex3f(xmax, ymax, zmin);
    glVertex3f(xmax, ymin, zmax);
    glVertex3f(xmax, ymax, zmax);
  glEnd();

  glBegin(GL_TRIANGLE_STRIP);
    glNormal3f(0.,0.,1.);
    glVertex3f(xmax, ymin, zmax);
    glVertex3f(xmax, ymax, zmax);
    glVertex3f(xmin, ymin, zmax);
    glVertex3f(xmin, ymax, zmax);
  glEnd();

  glBegin(GL_TRIANGLE_STRIP);
    glNormal3f(-1.,0.,0.);
    glVertex3f(xmin, ymin, zmax);
    glVertex3f(xmin, ymax, zmax);
    glVertex3f(xmin, ymin, zmin);
    glVertex3f(xmin, ymax, zmin);
  glEnd();

  glBegin(GL_TRIANGLE_STRIP);
    glNormal3f(0.,1.,0.);
    glVertex3f(xmin, ymax, zmin);
    glVertex3f(xmin, ymax, zmax);
    glVertex3f(xmax, ymax, zmin);
    glVertex3f(xmax, ymax, zmax);
  glEnd();

  glBegin(GL_TRIANGLE_STRIP);
    glNormal3f(0.,-1.,0.);
    glVertex3f(xmax, ymin, zmax);
    glVertex3f(xmax, ymin, zmin);
    glVertex3f(xmin, ymin, zmax);
    glVertex3f(xmin, ymin, zmin);
  glEnd();
  }


// Draw the GL scene
void GLTestWindow::drawScene(){
  const GLfloat lightPosition[]={15.,10.,5.,1.};
  const GLfloat lightAmbient[]={.1f,.1f,.1f,1.};
  const GLfloat lightDiffuse[]={.9f,.9f,.9f,1.};
  const GLfloat redMaterial[]={1.,0.,0.,1.};
  const GLfloat blueMaterial[]={0.,0.,1.,1.};

  GLdouble width = glcanvas->getWidth();
  GLdouble height = glcanvas->getHeight();
  GLdouble aspect = height>0 ? width/height : 1.0;

  // Make context current
  glcanvas->makeCurrent();

  glViewport(0,0,glcanvas->getWidth(),glcanvas->getHeight());

  glClearColor(1.0,1.0,1.0,1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  glEnable(GL_DEPTH_TEST);

  glDisable(GL_DITHER);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(30.,aspect,1.,100.);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(5.,10.,15.,0.,0.,0.,0.,1.,0.);

  glShadeModel(GL_SMOOTH);
  glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
  glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmbient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHTING);

  glMaterialfv(GL_FRONT, GL_AMBIENT, blueMaterial);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, blueMaterial);

  const float angle = 0.0;

  glPushMatrix();
  glRotated(angle, 0., 1., 0.);
  drawBox(-1, -1, -1, 1, 1, 1);

  glMaterialfv(GL_FRONT, GL_AMBIENT, redMaterial);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, redMaterial);

  glPushMatrix();
  glTranslated(0.,1.75,0.);
  glRotated(angle, 0., 1., 0.);
  drawBox(-.5,-.5,-.5,.5,.5,.5);
  glPopMatrix();

  glPushMatrix();
  glTranslated(0.,-1.75,0.);
  glRotated(angle, 0., 1., 0.);
  drawBox(-.5,-.5,-.5,.5,.5,.5);
  glPopMatrix();

  glPushMatrix();
  glRotated(90., 1., 0., 0.);
  glTranslated(0.,1.75,0.);
  glRotated(angle, 0., 1., 0.);
  drawBox(-.5,-.5,-.5,.5,.5,.5);
  glPopMatrix();

  glPushMatrix();
  glRotated(90., -1., 0., 0.);
  glTranslated(0.,1.75,0.);
  glRotated(angle, 0., 1., 0.);
  drawBox(-.5,-.5,-.5,.5,.5,.5);
  glPopMatrix();

  glPushMatrix();
  glRotated(90., 0., 0., 1.);
  glTranslated(0.,1.75,0.);
  glRotated(angle, 0., 1., 0.);
  drawBox(-.5,-.5,-.5,.5,.5,.5);
  glPopMatrix();

  glPushMatrix();
  glRotated(90., 0., 0., -1.);
  glTranslated(0.,1.75,0.);
  glRotated(angle, 0., 1., 0.);
  drawBox(-.5,-.5,-.5,.5,.5,.5);
  glPopMatrix();

  glPopMatrix();

  // Swap if it is double-buffered
  if(glvisual->isDoubleBuffer()){
    glcanvas->swapBuffers();
    }

  // Make context non-current
  glcanvas->makeNonCurrent();
  }


// Pop a dialog showing OpenGL properties
long GLTestWindow::onCmdOpenGL(FXObject*,FXSelector,void*){
  return 1;
  }

// Here we begin
int main(int argc,char *argv[]){

  // Make application
  FXApp application("GLTest","FoxTest");

  // Open the display
  application.init(argc,argv);

  // Make window
  new GLTestWindow(&application);

  // Create the application's windows
  application.create();

  // Run the application
  return application.run();
  }
