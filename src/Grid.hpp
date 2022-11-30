/* 
 * File: Grid.hpp
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

#ifndef GRID_HPP
#define GRID_HPP

#include <vector>
#include "Object.hpp"
#include "definitions.hpp"

#include <SDL2/SDL_image.h>

class Grid: public Object {

private:
  int n_rows, n_cols, n_nodes;
  FLOAT cell_size;
  std::vector<FLOAT> nodes;

  inline int index(int i, int j) const;
  inline int row(int ind) const ;
  inline int col(int ind) const;

  float cr, cg, cb;
  
public:
  Grid();
  Grid(int n_rows, int n_cols, FLOAT cs);
  ~Grid();

  void animate();
  void draw();
  
  bool isEmpty() const;

  void setColor(float r, float g, float b);
  
  int getNbRows() const;
  int getNbCols() const;
  FLOAT getCellSize() const;
  void setCellSize(FLOAT cs);
  int getHeight() const;
  int getWidth() const;
  VEC2 toWorld(int i, int j) const;

  FLOAT operator()(int i, int j) const;
  FLOAT &operator()(int i, int j);
  FLOAT value(FLOAT x, FLOAT y) const;
  FLOAT interpolatedValue(FLOAT x, FLOAT y) const;

  void reset(FLOAT val);
  void setValues(const SDL_Surface *texture);

  void getPixels(unsigned char *pixels) const;
  std::ostream& exportObj(std::ostream& os) const;
  
  Grid& operator=(const Grid& g);
  Grid& operator+=(const Grid& g);
  friend std::ostream& operator<<(std::ostream& os, const Grid& g);
  friend std::istream& operator>>(std::istream& is, Grid& F);
};

#endif
