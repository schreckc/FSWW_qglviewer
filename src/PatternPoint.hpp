/* 
 * File: PatternPoint.hpp
 *
 * Copyright (C) 2021  Camille Schreck
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

#ifndef PATTERNPOINT_HPP
#define PATTERNPOINT_HPP

#include <list>
#include <vector>
#include "definitions.hpp"

class PatternPoint {
private:
  VEC2 pos;
  COMPLEX ampli;
  
public:
  PatternPoint();
  PatternPoint(VEC2 p);
  PatternPoint(VEC2 p, COMPLEX a);
  ~PatternPoint();
  
  void setPos(VEC2 p);
  void setPos(FLOAT x, FLOAT y);
  void setAmpli(COMPLEX a);
  
  VEC2 getPos() const;
  COMPLEX getAmpli() const;
}; 

#endif
