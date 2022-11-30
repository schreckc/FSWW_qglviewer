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

#ifndef WAVEPATTERN_HPP
#define WAVEPATTERN_HPP

#include "EquivalentSource.hpp"
#include "PatternPoint.hpp"
#include "definitions.hpp"
#include "Grid.hpp"
#include <vector>
#include <Eigen/SVD>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>

class WavePattern {

private:

  std::string pattern_file; 
  VEC2 pos;
  int n_rows, n_cols; // size of the pattern to reproduce
  FLOAT cell_size;
  FLOAT radius; // radius of the circle of sources
  FLOAT wave_length;
  FLOAT wave_number;
  FLOAT ampli;
  FLOAT init_val;
  uint method; // 0 newton, 1 bfgs
  
  std::vector<EquivalentSource*> sources;
  std::vector<PatternPoint> pattern_pts;

  Grid grid;
  Eigen::BDCSVD<MatrixXcf> *svd;
  //  MatrixXcf transfer_mat;
  
public:
  WavePattern();
  WavePattern(FLOAT wl);
  ~WavePattern();

  void setPos(VEC2 p);
  void setPos(FLOAT x, FLOAT y);
  VEC2 getPos() const;

  void setSize(int nr, int nc, FLOAT cs);
  int getNbCols() const;
  int getNbRows() const;
  FLOAT getCellSize() const;
  
  void setAmpli(FLOAT a);
  void setMethod(uint m);
  void setInitVal(FLOAT v);
  
  void createPatternPoints(std::string file);
  void createSources(FLOAT r);
  
  void setTransferMatrix();
  void solve(std::list<Wave*> in);

  const std::vector<PatternPoint> & getPatternPoints() const;
  const std::vector<EquivalentSource*> & getSources() const;

  void getPatternPointsPos(std::list<VEC2> & positions) const;
  void getSourcesPos(std::list<VEC2> & positions) const;

  void energy(std::list<Wave*> in);

  // friend double pattern_energy_hess(const VECX& amplis, VECX* grad_out, MATX* hess_out, void* opt_data);
  // friend double pattern_energy(const VECX& amplis, VECX* grad_out, void* opt_data);
};

#endif
