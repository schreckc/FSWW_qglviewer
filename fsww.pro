CONFIG       += ordered
TEMPLATE      = subdirs
SUBDIRS       = QGLViewer designerPlugin

QT_VERSION=$$[QT_VERSION]

equals (QT_MAJOR_VERSION, 5) {
	cache()
}

ROOT_DIRECTORY = .
INCLUDE_DIR = $${ROOT_DIRECTORY}
LIB_DIR = $${ROOT_DIRECTORY}/QGLViewer
OBJECTS_DIR = obj/
EIGEN=/usr/include/eigen-3.4.0
OPTIM =/usr/local/include/optim

# Include path, taking QGLViewer path prefix into account
INCLUDEPATH *= $${INCLUDE_DIR}
               DEPENDPATH  *= $${INCLUDE_DIR}

               INCLUDEPATH += $${EIGEN}
               INCLUDEPATH += $${OPTIM}

CONFIG -= debug debug_and_release
CONFIG *= release

TEMPLATE = app
TARGET   = fsww

SRC_DIR = src/
HEADERS  = $${SRC_DIR}*.hpp
HEADERS  += $${SRC_DIR}*.h
SOURCES  = $${SRC_DIR}*.cpp

QMAKE_CXXFLAGS += -g -fopenmp -O3 -D__MODE_DEBUG=3 -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable -DOPTIM_ENABLE_EIGEN_WRAPPERS
LIBS += -L./QGLViewer -lQGLViewer-qt5 -lpthread -lGL -lgomp -fopenmp -lSDL2 -lSDL2_image -L/usr/lib -loptim -lblas
include( ./examples.pri )
