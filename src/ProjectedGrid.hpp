/* 
 * File: ProjectedGrid.hpp
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

#ifndef PROJECTED_GRID_HPP
#define PROJECTED_GRID_HPP

#include "settings.hpp"
#include <vector>

class ProjectedGrid {
private:
  int n_rows, n_cols, n_nodes;
  std::vector<VEC3> viewer_pos;
  std::vector<VEC2> displacement;
  std::vector<FLOAT> sizes;
  FLOAT min_size;

  inline int index(int i, int j) const;
  inline int row(int ind) const ;
  inline int col(int ind) const;

public:
  ProjectedGrid();
  ProjectedGrid(int nr, int nc);
  ~ProjectedGrid();

  VEC3 operator()(int i, int j) const;
  //VEC3 &operator()(int i, int j);
  VEC3 operator()(int i) const;
  //VEC3 &operator()(int i);
  VEC2 getPosWorld(int i, int j) const;
  VEC2 getPosWorld(int i) const;
  
  void setPosOnThePlane(int i, int j, FLOAT x, FLOAT y);
  void setHeight(int i, int j, FLOAT z);
  void setDisplacement(int i, int j, VEC2 d);
  void setDisplacement(int i, int j, FLOAT x, FLOAT y);
  
  void setSizes();
  FLOAT minSize();
  FLOAT size(int i);

  std::ostream& exportObj(std::ostream& os) const;
  
  friend std::ostream& operator<<(std::ostream& os, const ProjectedGrid& g);
};

#endif
