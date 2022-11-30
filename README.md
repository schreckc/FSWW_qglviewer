Fundemental Solutions For Water Wave Animation
------------------------------------------------------------------------------

This code implements a version the method described in

    "Fundemental Solutions For Water Wave Animation"
    by Camille Schreck, Christian Hafner and Chris Wojtan
    published in ACM Transactions on Graphics (SIGGRAPH 2019)


This is an early version, I still haven't tested all the functionality properly. It should be easier to compile and use the viewer than the previous version using the Magnum Engine. The code is in c++ and use Qglviewer to visualize the results.

This code is not particularly well-tested or optimized, but it can be used to
reproduce (at least some) of the examples in the paper.

INSTALL:
qmake fsww.pro
make

<optional> The sources contain a cuda implementation of part of this program. To use it you need to install cuda and make sure the macro USE_CUDA is defined before compiling.

RUN:
./fsww <options>
Options:
     -l, -load <file>: load configuration file (examples of configuation files in folder conf)
     -e, -export <name>: export amplitude grid for each frequencies in the file <name>.obj
     -i, -import <name>: import amplitude grid for each frequencies in the file <name>.obj
     -stop <t>: stop animation and exit at time t
     -em <name>: export heightfields and render files in a set of files <name><frame number>.ong and <name><frame number>.xml
     -es, -export_step <n>: export every n frames
     -h, -help: print help				

Some test configuration files can be found in the folder conf.
An example to try the program would be:
./fsww -l conf/circle_static.conf 
(planar wave reflecting on a circular obstacle)


OPTION:
You can use the macro defined in src/setting.hpp to switch between interactive or static waves (INTERACTIVE_), using or not the Cuda implementation (USE_CUDA), and using a projective grid or a fix one (PROJECTED_GRID).
The number of threads used by the program is defined in the same file (src/setting.hpp) by the macro
NTHREADS_.

