#ifndef OBJECT_HPP
#define OBJECT_HPP

#include <QGLViewer/qglviewer.h>


class Object {
  
protected:
  
public:
Object();

  virtual void animate() = 0;
  virtual void draw() = 0;
};

#endif
