#ifndef SPHERE_HPP
#define SPHERE_HPP

#include "Object.hpp"
#include <list>

class Sphere : public Object {
private :
  float m_size;
  float cr, cg, cb;

  const static int parallels_count = 4;
  const static int meridians_count = 5;
  const static int size_array = 18*meridians_count + 18*meridians_count* (parallels_count-1);

  static GLfloat* vertices;
  static GLfloat* normals;
  
public :
  
  Sphere(float s = 1);
  ~Sphere();

  void setColor(float r, float g, float b);
  void setSize(float s);
  
  void animate();
  void draw();

  static void create_array();
  static void create_vertex(float theta, float phi, int &index);

};
  

#endif
