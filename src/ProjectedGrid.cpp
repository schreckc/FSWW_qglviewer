/* 
 * File: ProjectedGrid.cpp
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

#include "ProjectedGrid.hpp"
#include "error.hpp"
#include "settings.hpp"

using namespace settings;

inline int ProjectedGrid::index(int i, int j) const {
  return n_cols*i + j;
}
  
inline int ProjectedGrid::row(int ind) const {
  return ind/n_cols;
}

inline int ProjectedGrid::col(int ind) const {
  return ind - (row(ind)*n_cols);
}

ProjectedGrid::ProjectedGrid() {
  n_rows = 0;
  n_cols = 0;
  n_nodes = 0;
  min_size = -1;
}



ProjectedGrid::ProjectedGrid(int nr, int nc) {
  n_rows = nr;
  n_cols = nc;
  n_nodes = n_rows*n_cols;
  viewer_pos = std::vector<VEC3>(n_nodes);
  displacement = std::vector<VEC2>(n_nodes);
  sizes = std::vector<FLOAT>(n_nodes);
#pragma omp for
  for (int i = 0; i < n_nodes; ++i) {
    viewer_pos[i] = VEC3(0, 0, 0);
  }
}

ProjectedGrid::~ProjectedGrid() {}

VEC3 ProjectedGrid::operator()(int i, int j) const {
  int ind = index(i, j);
  assert(ind >= 0 && ind < n_nodes);
  return viewer_pos[ind];
}

// VEC3 &ProjectedGrid::operator()(int i, int j) {
//  int ind = index(i, j);
//  assert(ind >= 0 && ind < n_nodes);
//  return viewer_pos[ind];
//}

VEC3 ProjectedGrid::operator()(int i) const {
  assert(i >= 0 && i < n_nodes);
  return viewer_pos[i];
}

// VEC3 &ProjectedGrid::operator()(int i) {
//   assert(i >= 0 && i < n_nodes);
//   return viewer_pos[i];
// }

VEC2 ProjectedGrid::getPosWorld(int i, int j) const {
  if (i < 0 || i >= n_rows || j < 0 || j > n_cols) {
    return VEC2(0, 0);
  }
  int ind = index(i, j);
  assert(ind >= 0 && ind < n_nodes);
  VEC2 vp(viewer_pos[ind](0), viewer_pos[ind](1));
  VEC2 wp = viewer2world(vp);
  return wp;
}

VEC2 ProjectedGrid::getPosWorld(int i) const {
  if (i >= 0 && i < n_nodes) {
    return VEC2(0, 0);
  }
  VEC2 vp(viewer_pos[i](0), viewer_pos[i](1));
  VEC2 wp = viewer2world(vp);
  return wp;
}

void ProjectedGrid::setPosOnThePlane(int i, int j, FLOAT x, FLOAT y) {
  int ind = index(i, j);
  assert(ind >= 0 && ind < n_nodes);

  viewer_pos[ind](0) = x;
  viewer_pos[ind](1) = y;
}

void ProjectedGrid::setHeight(int i, int j, FLOAT z) {
  int ind = index(i, j);
  assert(ind >= 0 && ind < n_nodes);
  viewer_pos[ind](2) = z;
}

void ProjectedGrid::setDisplacement(int i, int j, VEC2 d) {
  int ind = index(i, j);
  assert(ind >= 0 && ind < n_nodes);
  displacement[ind] = d*scale_;
}

void ProjectedGrid::setDisplacement(int i, int j, FLOAT x, FLOAT y) {
  int ind = index(i, j);
  assert(ind >= 0 && ind < n_nodes);
  displacement[ind] = VEC2(x, y)*scale_;
}


void ProjectedGrid::setSizes() {
  min_size = -1;
  for (int i = 0; i < n_rows; ++i) {
    for (int j = 0; j < n_cols; ++j) {
      int ind = index(i, j);
      VEC2 p = getPosWorld(i, j);
      int n = 0;
      FLOAT s = 0;
      if (i != 0) {
	VEC2 neigh = getPosWorld(i-1, j);
	s += (neigh - p).norm();
	++n;
      }
      if (j != 0) {
	VEC2 neigh = getPosWorld(i, j-1);
	s += (neigh - p).norm();
	++n;
      }
      if (i != n_rows - 1) {
	VEC2 neigh = getPosWorld(i+1, j);
	s += (neigh - p).norm();
	++n;
      }
      if (j != n_cols - 1) {
	VEC2 neigh = getPosWorld(i, j+1);
	s += (neigh - p).norm();
	++n;
      }
      s /= (FLOAT)n;
      sizes[ind] = s;
      if (min_size < 0 || s < min_size) {
	min_size = s;
      }
    }
  }
  
}

FLOAT ProjectedGrid::minSize() {
  return min_size;
}

FLOAT ProjectedGrid::size(int i) {
  return sizes[i];
}


std::ostream& ProjectedGrid::exportObj(std::ostream& os) const {
  os << "# Projected Grid\n";
  for (int i = 0; i < n_nodes; ++i) {
    VEC3 pos = viewer_pos[i];
    // VEC2 wc = getPosWorld(i);
    // Vector2i oc = world2gridObs(wc);
    // if (displacement[i].norm() < 0.01) {
    //   pos(0) += displacement[i](0);
    //   pos(1) += displacement[i](1);
    // }

    /** todo: create a global parameter to make sure waves inside the 
	obstacle are not higher than the pbstacle */
    if (pos(2) >= 0.019) {
      pos(2) = 0.019;
    }
    os <<"v "<<pos(0)<<" "<<pos(1)<<" "<<pos(2)<<"\n";
  }
  for (int i = 0; i < n_rows-1; i++) {
    for (int j = 0; j < n_cols-1; j++) {
      int idx = j + i * n_cols + 1;
      int J   = 1;
      int I   = n_cols;
      os << "f "<<idx<<" "<<idx + J<<" "<<idx + I<<"\n";
      os << "f "<<idx + J<<" "<<idx + I + J<<" "<<idx + I<<"\n";
    }
  }
  return os;
}

std::ostream& operator<<(std::ostream& os, const ProjectedGrid& g) {
  os << g.n_rows<<" "<<g.n_cols<<" ";
  for (int i = 0; i < g.n_nodes; ++i) {
    os <<"("<<g.viewer_pos[i](0)<<","<<g.viewer_pos[i](1)<<","<<g.viewer_pos[i](2)<<") ";
  }

  return os;
}

// friend std::istream& operator>>(std::istream& is, ProjectedGrid& g) {
//     is >> g.n_rows >> g.n_cols;
//   g.n_nodes = g.n_rows*g.n_cols;
//   g.nodes = std::vector<FLOAT>(g.n_nodes);
//    for (int i = 0; i < g.n_nodes; ++i) {
//      is >> g.nodes[i];
//   }
//    return is;
// }
