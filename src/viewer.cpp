/****************************************************************************

 Copyright (C) 2002-2014 Gilles Debunne. All rights reserved.

 This file is part of the QGLViewer library version 2.7.2.

 http://www.libqglviewer.com - contact@libqglviewer.com

 This file may be used under the terms of the GNU General Public License 
 versions 2.0 or 3.0 as published by the Free Software Foundation and
 appearing in the LICENSE file included in the packaging of this file.
 In addition, as a special exception, Gilles Debunne gives you certain 
 additional rights, described in the file GPL_EXCEPTION in this package.

 libQGLViewer uses dual licensing. Commercial/proprietary software must
 purchase a libQGLViewer Commercial License.

 This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

*****************************************************************************/
#include "Grid.hpp"
#include "viewer.hpp"
#include "ui_parameters.hpp"
#include "error.hpp"

#include <QCursor>
#include <QKeyEvent>
#include <QMap>
#include <QMenu>
#include <QMouseEvent>

using namespace std;

Viewer::~Viewer() {
  if (plot_ && stream_plot.is_open()) {
      stream_plot<<"set term pop; set out;";
      stream_plot.close();
    }
  
  _surface.clear();
}

void help_parse() {
  std::cout<<"\n     *** WAVE: Help ***\n"<<std::endl;
  std::cout<<"Synopsis: \n     .\\main <options>\n\nOptions:"<<std::endl;
  std::cout<<"     -l, -load <file>: load configuration file"<<std::endl;
  std::cout<<"     -e, -export <name>: export amplitude grid for each frequencies in the file <name>.obj"<<std::endl;
  std::cout<<"     -i, -import <name>: import amplitude grid for each frequencies in the file <name>.obj"<<std::endl;
  std::cout<<"     -stop <t>: stop animation and exit at time t"<<std::endl;
  std::cout<<"     -em <name>: export heightfields and render files in a set of files <name><frame number>.ong and <name><frame number>.xml"<<std::endl;
  std::cout<<"     -es, -export_step <n>: export every n frames"<<std::endl;
  std::cout<<"     -h, -help: print help\n"<<std::endl;
  exit(0);
}

void Viewer::treatArguments(int argc, char **argv) {
  running_ = false;
  plot_ = false;
  std::cout<<"treat arg "<<argc<<std::endl;
  for (int i = 1;  i < argc; ++i) {
    std::string s(argv[i]);
    if (s == "-l" || s == "-load") {
      if (argc < i + 2) {
  	std::cerr<<"\nERROR: wrong number of arguments\n"<<std::endl;
  	help_parse();
      }
      std::cout<<"Loading configuration file:"<<" "<<argv[i+1]<<std::endl;
      _surface.setImportConf(argv[i+1]);
      ++i;
    } else if (s == "-i" || s == "-import") {
      if (argc < i + 2) {
  	std::cerr<<"\nERROR: wrong number of arguments\n"<<std::endl;
  	help_parse();
      }
      std::cout<<"Importing"<<" "<<argv[i+1]<<std::endl;
      _surface.setImport(argv[i+1]);
      ++i;
    } else if (s == "-e" || s == "-export") {
      if (argc < i + 2) {
  	std::cerr<<"\nERROR: wrong number of arguments\n"<<std::endl;
  	help_parse();
      }
      std::cout<<"Exporting"<<" "<<argv[i+1]<<std::endl;
      _surface.setExport(argv[i+1]);
      ++i;
    } else if (s == "-p" || s == "-plot") {
      if (argc < i + 2) {
  	std::cerr<<"\nERROR: wrong number of arguments\n"<<std::endl;
  	help_parse();
      }
      std::cout<<"plotting in "<<" "<<argv[i+1]<<std::endl;
      export_file_data = argv[i+1];
      plot_ = true;
      ++i;
    } else if (s == "-em") {
      if (argc < i + 2) {
  	std::cerr<<"\nERROR: wrong number of arguments\n"<<std::endl;
  	help_parse();
      }
      std::cout<<"Exporting (mitsuba) "<<" "<<argv[i+1]<<std::endl;
      export_m = true;
      export_file_m = argv[i+1];
      std::stringstream ss_plot;
      // ss_plot <<export_file_m<<"_plot.txt";
      // str_plot = std::string(ss_plot.str());
      ++i;
    } else if (s == "-d" || s == "-data") {
      if (argc < i + 2) {
  	std::cerr<<"\nERROR: wrong number of arguments\n"<<std::endl;
  	help_parse();
      }
      std::cout<<"Saving amplitude data in"<<" "<<argv[i+1]<<std::endl;
      _surface.setData(argv[i+1]);
      ++i;
    } else if (s == "-stop") {
      if (argc < i + 2) {
  	std::cerr<<"\nERROR: wrong number of arguments\n"<<std::endl;
  	help_parse();
      }
      std::cout<<"Stop at t = "<<argv[i+1]<<std::endl;
      // _surface.setStopTime(atoi(argv[i+1]));
      stop_time = atoi(argv[i+1]);
      ++i;
    } else if (s == "-export_step" || s == "-es") {
      if (argc < i + 2) {
  	std::cerr<<"\nERROR: wrong number of arguments\n"<<std::endl;
  	help_parse();
      }
      std::cout<<"Export every "<<argv[i+1]<<" steps"<<std::endl;
      _surface.setExportStep(atoi(argv[i+1]));
      export_step_m = atoi(argv[i+1]);
      ++i;
    } else if (s == "-r" || s == "-run") {
      running_ = true;
    } else if (s == "-h" || s == "-help") {
      std::cout<<"help"<<std::endl;
      help_parse();
    } else {
      std::cerr<<"\nERROR: Unknown option\n"<<std::endl;
      help_parse();
    }
  }
}

void Viewer::animate() {
  try {
    if (time_ > stop_time) {
      std::exit(0);
    }
  _surface.update();

  if (plot_) {
  uint n = time_/export_step_m;
  std::string s0 = "";
  if (n < 10) {
    s0 = "000";
  } else if (n < 100) {
    s0 = "00";
  } else if (n < 1000) {
    s0 = "0";
  }
  
  std::stringstream ss_dat;
  ss_dat <<export_file_data<<s0<<n<<".dat";
  std::string str_dat(ss_dat.str());
  _surface.drawHeighField(str_dat);
  stream_plot<<"set output \""<<export_file_data<<"_2d_"<<s0<<n<<".png\"\n";
  stream_plot<<"splot '"<<str_dat<<"' with pm3d\n";
  }
  ++time_;
   } catch (std::exception& e) {
    std::cerr << "Exception catched : " << e.what() << std::endl;
    _surface.clear();
    throw;
  }
}
  
void Viewer::draw() {
  // Grid grid(100, 100, 0.1);
  // grid.draw();
  // glPushMatrix();
  // sphere.draw();
  // glTranslatef(2, 0, 0);
  //   sphere.draw();
  // glPopMatrix();


float pos[4] = {1.0, 1.0, 1.0, 0.0};
  // Directionnal light
  glLightfv(GL_LIGHT0, GL_POSITION, pos);

   //  // // Spot light
   // qglviewer::Vec pos1 = light1->position();
   // pos[0] = float(pos1.x);
   // pos[1] = float(pos1.y);
   // pos[2] = float(pos1.z);
   // glLightfv(GL_LIGHT1, GL_POSITION, pos);
   // glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION,
   //           light1->inverseTransformOf(qglviewer::Vec(0, 0, 1)));

   // Point light
   qglviewer::Vec pos2 = light2->position();
   pos[0] = float(pos2.x);
   pos[1] = float(pos2.y);
   pos[2] = float(pos2.z);
   glLightfv(GL_LIGHT2, GL_POSITION, pos);
 
  //  drawLight(GL_LIGHT0);
  //  if (light1->grabsMouse())
  //    drawLight(GL_LIGHT1, 1.2f);
  //  else
   //      drawLight(GL_LIGHT1);

  // if (light2->grabsMouse())
  //   drawLight(GL_LIGHT2, 1.2f);
  // else
     drawLight(GL_LIGHT2);
  
   //  sphere.draw();
  _surface.draw();
 // // Draw the intersection line
 //  glBegin(GL_LINES);
 //  glVertex3fv(orig);
 //  glVertex3fv(orig + 100.0 * dir);
 //  glEnd();


   if (found) {
     glPushMatrix();
     glTranslatef(selectedPoint.x, selectedPoint.y, selectedPoint.z);
     sphere.draw();
     glPopMatrix();
   }
  
}


void Viewer::init() {
   glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // Light0 is the default ambient light
  glEnable(GL_LIGHT0);

  // Light1 is a spot light
   // glEnable(GL_LIGHT1);
   // const GLfloat light_ambient[4] = {0.2f, 0.2f, 0.2f, 1.0};//{0.8f, 0.2f, 0.2f, 1.0};
   // const GLfloat light_diffuse[4] = {0.0, 0.0f, 0.0f, 0.0};
   // const GLfloat light_specular[4] = {0.5, 0.5, 0.5, 1.0};

   // glLightf(GL_LIGHT1, GL_SPOT_EXPONENT, 3.0);
   // glLightf(GL_LIGHT1, GL_SPOT_CUTOFF, 20.0);
   // glLightf(GL_LIGHT1, GL_CONSTANT_ATTENUATION, 0.5);
   // glLightf(GL_LIGHT1, GL_LINEAR_ATTENUATION, 1.0);
   // glLightf(GL_LIGHT1, GL_QUADRATIC_ATTENUATION, 1.5);
   // glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient);
   // glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular);
   // glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse);

  // Light2 is a classical directionnal light
  glEnable(GL_LIGHT2);
  const GLfloat light_ambient2[4] = {0.2f, 0.2f, 0.2f, 1.0};
  const GLfloat light_diffuse2[4] = {0.2f, 0.2f, 0.2, 1.0};
  const GLfloat light_specular2[4] = {0.2, 0.2, 0.2, 1.0};

  glLightfv(GL_LIGHT2, GL_AMBIENT, light_ambient2);
  glLightfv(GL_LIGHT2, GL_SPECULAR, light_specular2);
  glLightfv(GL_LIGHT2, GL_DIFFUSE, light_diffuse2);

  //   light1 = new qglviewer::ManipulatedFrame();
   light2 = new qglviewer::ManipulatedFrame();
  // setMouseTracking(true);

  // light1->setPosition(1, 1, 0.5);
  // Align z axis with -position direction : look at scene center
   //   light1->setOrientation(qglviewer::Quaternion(qglviewer::Vec(1, 1, 0), -light1->position()));

    light2->setPosition(0.0, 0.0, 1);


  _surface.reset();
  
  // Restore previous viewer state.
  restoreStateFromFile();

   // Add custom key description (see keyPressEvent).
  setKeyDescription(Qt::Key_W, "Toggles wire frame display");
  setKeyDescription(Qt::Key_S, "Toogles sources and boundary points display");
  setKeyDescription(Qt::Key_I, "Toogles in waves display");

  setMouseTracking(true);

  setSceneRadius(30);

  srand(time(NULL));
  sphere.create_array();
  sphere.setSize(0.1);
  sphere.setColor(0.9f, 0.2f, 0.1f);
  // Opens help window
  help();
  
  glDisable(GL_CULL_FACE);

  if (plot_) {
  
 std::stringstream ss_plot;
 ss_plot <<export_file_data<<"_plot.txt";
 str_plot = std::string(ss_plot.str());
 if (stream_plot.is_open()) {
   stream_plot.close();
 }
 INFO("Writing "<<str_plot);
 stream_plot.open(str_plot);
 stream_plot<<"set view map\n";
 stream_plot<<"unset key\n";
 stream_plot<<"unset tics\n";
 stream_plot<<"unset border\n";
 stream_plot<<"unset colorbox\n";
 stream_plot<<"set terminal png size 600, 600\n";
  }
 if (running_) {
   startAnimation();
 }
}


void Viewer::postSelection(const QPoint &point) {
  // Compute orig and dir, used to draw a representation of the intersecting
  // line
  camera()->convertClickToLine(point, orig, dir);
  FLOAT k;
  if (dir.z == 0) {
    found = false;
    return;
  }
  found = true;
  k  = -orig.z/dir.z;
  selectedPoint.x = orig.x + k*dir.x;
  selectedPoint.y = orig.y + k*dir.y;
  selectedPoint.z = 0;
  // Find the selectedPoint coordinates, using camera()->pointUnderPixel().
  //  selectedPoint = camera()->pointUnderPixel(point, found);
  //  std::cout<<selectedPoint<<std::endl;
  // selectedPoint -= 0.01f * dir; // Small offset to make point clearly visible.
  // Note that "found" is different from (selectedObjectId()>=0) because of the
  // size of the select region.

  _surface.addDrop(selectedPoint.x,selectedPoint.y, ui_parameters::size_drop_);
}

void Viewer::keyPressEvent(QKeyEvent *e) {
  // Get event modifiers key
  const Qt::KeyboardModifiers modifiers = e->modifiers();

  // A simple switch on e->key() is not sufficient if we want to take state key
  // into account. With a switch, it would have been impossible to separate 'F'
  // from 'CTRL+F'. That's why we use imbricated if...else and a "handled"
  // boolean.
  bool handled = false;
  if ((e->key() == Qt::Key_W) && (modifiers == Qt::NoButton)) {
    wireframe_ = !wireframe_;
    if (wireframe_) {
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    } else {
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }
    handled = true;
    update();
  } else if ((e->key() == Qt::Key_Backspace) && (modifiers == Qt::NoButton)) {
    _surface.reset();
    handled = true;
    update();
  } else if ((e->key() == Qt::Key_S) && (modifiers == Qt::NoButton)) {
    _surface.draw_sources = !_surface.draw_sources;
    handled = true;
    update();
  } else if ((e->key() == Qt::Key_O) && (modifiers == Qt::NoButton)) {
    _surface.draw_obs = !_surface.draw_obs;
    handled = true;
    update();
  } else if ((e->key() == Qt::Key_I) && (modifiers == Qt::NoButton)) {
    ui_parameters::show_in_field = !ui_parameters::show_in_field;
    handled = true;
    update();
  } else if ((e->key() == Qt::Key_P) && (modifiers == Qt::NoButton)) {
    ui_parameters::show_scattered_field = !ui_parameters::show_scattered_field;
    handled = true;
    update();
  } else if ((e->key() == Qt::Key_E) && (modifiers == Qt::NoButton)) {
    _surface.exportSurfaceObj("surface.obj");
    _surface.exportObstacleGrid("obstacle.obj");
    handled = true;
    update();
  }
  
  // ... and so on with other else/if blocks.

  if (!handled)
    QGLViewer::keyPressEvent(e);
}

QString Viewer::helpString() const {
  QString text("<h2>S i m p l e V i e w e r</h2>");
  text += "Use the mouse to move the camera around the object. ";
  text += "You can respectively revolve around, zoom and translate with the "
          "three mouse buttons. ";
  text += "Left and middle buttons pressed together rotate around the camera "
          "view direction axis<br><br>";
  text += "Pressing <b>Alt</b> and one of the function keys "
          "(<b>F1</b>..<b>F12</b>) defines a camera keyFrame. ";
  text += "Simply press the function key again to restore it. Several "
          "keyFrames define a ";
  text += "camera path. Paths are saved when you quit the application and "
          "restored at next start.<br><br>";
  text +=
      "Press <b>F</b> to display the frame rate, <b>A</b> for the world axis, ";
  text += "<b>Alt+Return</b> for full screen mode and <b>Control+S</b> to save "
          "a snapshot. ";
  text += "See the <b>Keyboard</b> tab in this window for a complete shortcut "
          "list.<br><br>";
  text += "Double clicks automates single click actions: A left button double "
          "click aligns the closer axis with the camera (if close enough). ";
  text += "A middle button double click fits the zoom of the camera and the "
          "right button re-centers the scene.<br><br>";
  text += "A left button double click while holding right button pressed "
          "defines the camera <i>Revolve Around Point</i>. ";
  text += "See the <b>Mouse</b> tab and the documentation web pages for "
          "details.<br><br>";
  text += "Press <b>Escape</b> to exit the viewer.";
  return text;
}
