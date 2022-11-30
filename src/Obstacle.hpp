/* 
 * File: Obstacle.hpp
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

#ifndef OBSTACLE_HPP
#define OBSTACLE_HPP

#include <SDL2/SDL_image.h>
#include "InputPoint.hpp"
#include "EquivalentSource.hpp"
#include "MovingEquivalentSource.hpp"
#include "Grid.hpp"
#include "Wave.hpp"
#include <Eigen/SVD>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>
#include <fstream>

class WaterSurface;
class ComposedObstacle;

class Obstacle {

protected:
  std::vector<FLOAT> wave_lenghts;
  std::vector<std::vector<EquivalentSource*> >sources_l;
  std::vector<std::vector<InputPoint*> >boundaries_l;
  std::vector<std::vector<VEC2C> >normals_l;
  std::vector<Eigen::BDCSVD<MatrixXcf> >svd_l;

  VEC2 pos;
  int nb_wl;

  std::vector<int> nb_old_sources;
  bool movable_;

  virtual void setBoundaries(int w) = 0;
  virtual void setEquivalentSources(int i) = 0;
  
public:

  //interactive
  void setAmpliSources(const std::list<Wave*> &waves, int index, int t);
  void setTransfeMatrix(int index, int t);

  // frequency
  void setAmpliSources(const std::list<Wave*> &waves, int index);
  void setTransfeMatrix(int index);

  Obstacle();
  Obstacle(VEC2 p);
  virtual ~Obstacle();

  void reset(const std::vector<std::list<Wave*> > &waves_l,const std::vector<FLOAT> &wl);
  
  void setPos(VEC2 p);
  void setPos(FLOAT x, FLOAT y);
  VEC2 getPos();
  
  void getBoundariesPos(std::list<VEC2> &boundaries) const;
  void getSourcesPos(std::list<VEC2> &sources) const;

  void update(const WaterSurface *surface);

  FLOAT height(FLOAT x, FLOAT y, FLOAT time) const;
  COMPLEX heightc(FLOAT x, FLOAT y, FLOAT time) const;
  COMPLEX heightc_wl(FLOAT x, FLOAT y, int wl) const;
  COMPLEX heightc_wl(FLOAT x, FLOAT y, int w, FLOAT time) const;

  void getSources(std::list<Wave*> &waves, int w);
  void getSourcesM(std::list<MovingEquivalentSource*> &waves, int w);
 
  std::vector<InputPoint*> &getBoundaries(int w);
  std::vector<VEC2C> &getNormales(int w);
  std::vector<EquivalentSource*> &getSources(int w);
  std::vector<EquivalentSource*> &getSourcesApprox(int w);
  int getBoundariesSize(int w) const;
  int getSourcesSize(int w) const;
  int getSourcesApproxSize(int w) const;

  
  virtual FLOAT getGrid(int i, int j) const = 0;
  virtual std::ofstream& exportMitsuba(std::ofstream &file) const;

  void move(VEC2 dir, FLOAT speed, int time);

  void setMovable(bool mv);
  bool isMovable() const;
  

};

#endif
