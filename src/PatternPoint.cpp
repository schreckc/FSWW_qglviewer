/* 
 * File: PatternPoint.cpp
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

#include "PatternPoint.hpp"
#include <fstream>
#include <iostream>

PatternPoint::PatternPoint() {
  pos = VEC2(0, 0);
  ampli = COMPLEX(0, 0);
}

PatternPoint::PatternPoint(VEC2 p) {
  pos = p;
  ampli = COMPLEX(0, 0);
}

PatternPoint::PatternPoint(VEC2 p, COMPLEX a) {
  pos = p;
  ampli = a;
}

PatternPoint::~PatternPoint() {
}

void PatternPoint::setPos(VEC2 p) {
  pos = p;
}

void PatternPoint::setPos(FLOAT x, FLOAT y) {
  pos = VEC2(x, y);
}

void PatternPoint::setAmpli(COMPLEX a) {
  ampli = a;
}

VEC2 PatternPoint::getPos() const {
  return pos;
}

COMPLEX PatternPoint::getAmpli() const {
  return ampli;
}
 
