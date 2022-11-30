/* 
 * File: Grid.cpp
 *
 * Copyright (C) 2019  Camille Schreck
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "Grid.hpp"
#include <iostream>
 #include "error.hpp"
#include "settings.hpp"

using namespace settings;

inline int Grid::index(int i, int j) const {
  return n_cols*i + j;
}
  
inline int Grid::row(int ind) const {
  return ind/n_cols;
}

inline int Grid::col(int ind) const {
  return ind - (row(ind)*n_cols);
}


Grid::Grid() {
  n_rows = 0;
  n_cols = 0;
  cell_size = 0;
   setColor(1, 1, 1);
}

Grid::Grid(int rows, int cols, FLOAT cs) {
   setColor(1, 1, 1);
  n_rows = rows;
  n_cols = cols;
  n_nodes = n_rows*n_cols;
  cell_size = cs;

  nodes = std::vector<FLOAT>(n_nodes);
  reset(0);
}

Grid::~Grid() {}

void Grid::animate() {}
void Grid::draw() {
  glBegin(GL_TRIANGLES);
  glColor3f(cr, cg, cb);
  for (int i = 0; i < n_rows - 1; i++) {
    for (int j = 0; j < n_cols - 1; j++) {
      FLOAT nx = nodes[index(i+1,j)] - nodes[index(i,j)];
      FLOAT ny = nodes[index(i,j)] - nodes[index(i,j+1)];
      FLOAT nz = cell_size;
      FLOAT n = sqrt(nx*nx + ny*ny + nz*nz);
      nx /= n;
      ny /= n;
      nz /= n;
      
      glNormal3f(nx, ny, nz);
      glVertex3f(i*cell_size, j*cell_size, nodes[index(i,j)]);
      glNormal3f(nx, ny, nz);
      glVertex3f((i+1)*cell_size, j*cell_size, nodes[index(i+1,j)]);
      glNormal3f(nx, ny, nz);
      glVertex3f((i+1)*cell_size, (j+1)*cell_size, nodes[index(i+1,j+1)]);

      glNormal3f(nx, ny, nz);
      glVertex3f((i+1)*cell_size, (j+1)*cell_size, nodes[index(i+1,j+1)]);
      glNormal3f(nx, ny, nz);
      glVertex3f(i*cell_size, (j+1)*cell_size, nodes[index(i,j+1)]);
      glNormal3f(nx, ny, nz);
      glVertex3f(i*cell_size, j*cell_size, nodes[index(i,j)]);
    }
  }
  glEnd();
  
  //     glBegin(GL_LINES);
  // glColor3f(cr, cg, cb);
  //  for (int i = 0; i < n_rows - 1; i++) {
  // for (int j = 0; j < n_cols - 1; j++) {
  //     glVertex3f(i*cell_size, j*cell_size, nodes[index(i,j)]);
  //     glVertex3f((i+1)*cell_size, j*cell_size, nodes[index(i+1,j)]);
  //     glVertex3f((i+1)*cell_size, j*cell_size, nodes[index(i+1,j)]);
  //     glVertex3f((i+1)*cell_size, (j+1)*cell_size, nodes[index(i+1,j+1)]);
  //     glVertex3f((i+1)*cell_size, (j+1)*cell_size, nodes[index(i+1,j+1)]);
  //     glVertex3f(i*cell_size, j*cell_size, nodes[index(i,j)]);

  //     glVertex3f((i+1)*cell_size, (j+1)*cell_size, nodes[index(i+1,j+1)]);
  //     glVertex3f(i*cell_size, (j+1)*cell_size, nodes[index(i,j+1)]);
  //     glVertex3f(i*cell_size, (j+1)*cell_size, nodes[index(i,j+1)]);
  //     glVertex3f(i*cell_size, j*cell_size, nodes[index(i,j)]);
  //     glVertex3f(i*cell_size, j*cell_size, nodes[index(i,j)]);
  //     glVertex3f((i+1)*cell_size, (j+1)*cell_size, nodes[index(i+1,j+1)]);
  // }
  // glEnd();

  // glBegin(GL_LINES);
  //  glColor3f(0, 0, 1);
  // for (int i = 0; i < n_rows - 1; i++) {
  //   for (int j = 0; j < n_cols - 1; j++) {
  //        FLOAT nx = nodes[index(i+1,j)] - nodes[index(i,j)];
  //     FLOAT ny = nodes[index(i,j)] - nodes[index(i,j+1)];
  //     FLOAT nz = cell_size;
  //     FLOAT n = 10*sqrt(nx*nx + ny*ny + nz*nz);
  //     nx /= n;
  //     ny /= n;
  //     nz /= n;
      
  //     glVertex3f(i*cell_size, j*cell_size, nodes[index(i,j)]);
  //     glVertex3f(i*cell_size+nx, j*cell_size+ny, nodes[index(i,j)]+nz);
  //   }
  // }
  //  glEnd();
}

bool Grid::isEmpty() const {
  return n_nodes == 0;
}

void Grid::setColor(float r, float g, float b) {
   cr = r;
   cg = g;
   cb = b;
}


int Grid::getNbRows() const {
  return n_rows;
}

int Grid::getNbCols() const {
  return n_cols;
}

FLOAT Grid::getCellSize() const {
  return cell_size;
}

void Grid::setCellSize(FLOAT cs) {
  cell_size = cs;
}

int Grid::getHeight() const {
  return n_cols*cell_size;
}
int Grid::getWidth() const {
  return n_rows*cell_size;
}

VEC2 Grid::toWorld(int i, int j) const {
  VEC2 out(i*cell_size, j*cell_size);
  return out;
}

FLOAT Grid::operator()(int i, int j) const {
  if (i < 0 || i >= n_rows || j < 0 || j >= n_cols) {
    return 0;
  }
  int ind = index(i, j);
  return nodes[ind];
}
FLOAT & Grid::operator()(int i, int j) {
  if (i < 0 || i >= n_rows || j < 0 || j >= n_cols) {
    nodes[0] = 0;
    return nodes[0];
  }
  int ind = index(i, j);
  return nodes[ind];
}

FLOAT Grid::interpolatedValue(FLOAT x, FLOAT y) const {
  ERROR(false, "TODO: Grid::interpolatedValue", "");
}

void Grid::reset(FLOAT val) {
  #pragma omp for
  for (int i = 0; i < n_nodes; ++i) {
    nodes[i] = val;
  }
}
  
void Grid::setValues(const SDL_Surface * texture) {
  int w = texture->w;
  int h = texture->h;
  FLOAT pix_per_cell_x = (FLOAT)h/(FLOAT)n_rows;
  FLOAT pix_per_cell_y = (FLOAT)w/(FLOAT)n_cols;
  FLOAT byte_per_cell_x = texture->format->BytesPerPixel * pix_per_cell_x;
  FLOAT byte_per_cell_y = texture->format->BytesPerPixel * pix_per_cell_y;
  unsigned char* pixels = (unsigned char*) texture->pixels;
  
  #pragma omp for
  for (int i = 0; i < n_rows; ++i) {
    int x = i*byte_per_cell_x/texture->format->BytesPerPixel;//
    x *= texture->format->BytesPerPixel;
    for (int j = 0; j < n_cols; ++j) {
      int y = j*byte_per_cell_y/texture->format->BytesPerPixel;
      y *= texture->format->BytesPerPixel;
      
      nodes[i*n_cols+j] = (int)pixels[(int)(x*w+y)]/256.0;
       if (nodes[i*n_cols+j] < 0.1) {
       	nodes[i*n_cols+j] = 0;
       } else if (nodes[i*n_cols+j] > 0.9) {
       	nodes[i*n_cols+j] = 1;
       } 
    }
  }
}

std::ostream& Grid::exportObj(std::ostream& os) const {
  INFO("export Grid");
  os << "# Grid\n";
  for (int i = 0; i < n_rows; i++) {
    for (int j = 0; j < n_cols; j++) {
      VEC2 posv = VEC2(2*i*cell_size/scale_ - 1, 2*(FLOAT)j*cell_size/scale_ - 1);
      FLOAT h = nodes[index(i, j)];
      if (h == 0 && cell_size == cell_size_obs) { 
       	h = -1;
      } 
      os <<"v "<<posv(0)<<" "<<posv(1)<<" "<<h<<"\n";
    }
  }
  for (int i = 0; i < n_rows-1; i++) {
    for (int j = 0; j < n_cols-1; j++) {
      int idx = j + i * n_cols + 1;
      int J   = 1;
      int I   = n_cols;
      os << "f "<<idx<<" "<<idx + J<<" "<<idx + I<<"\n";
      os << "f "<<idx + I<<" "<<idx + J<<" "<<idx + I + J<<"\n";
    }
  }
  return os;
}


Grid& Grid::operator=(const Grid& g) {
  n_rows = g.n_rows;
  n_cols = g.n_cols;
  n_nodes = n_rows*n_cols;
  cell_size = g.cell_size;

  nodes = std::vector<FLOAT>(n_nodes);
  for (int i = 0; i < n_nodes; ++i) {
    nodes[i] = g.nodes[i];
  }
  return *this;
}

std::ostream& operator<<(std::ostream& os, const Grid& g) {
  os << g.n_rows<<" "<<g.n_cols<<" ";
  for (int i = 0; i < g.n_nodes; ++i) {
    os << g.nodes[i]<<" ";
  }
  return os;
}
std::istream& operator >> (std::istream& is, Grid& g) {
  is >> g.n_rows >> g.n_cols;
  g.n_nodes = g.n_rows*g.n_cols;
  g.nodes = std::vector<FLOAT>(g.n_nodes);
  for (int i = 0; i < g.n_nodes; ++i) {
    is >> g.nodes[i];
  }
  return is;
}
