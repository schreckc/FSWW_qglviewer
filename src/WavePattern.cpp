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

#include "WavePattern.hpp"
#include "error.hpp"
#include "settings.hpp"
#include <SDL2/SDL_image.h>
#include <optim.hpp>
#include <fstream>

using namespace optim;

using namespace settings;

struct wp_data {
  std::vector<EquivalentSource*> sources;
  std::vector<PatternPoint>  pattern_pts;
  std::list<Wave*> in;
  //  ColVec_t pp_ampli;
  MATX transfer_mat;
  FLOAT k;
  // MATX transfer_mat_dx;
  // MATX transfer_mat_dy;
};


double pattern_energy_hess(const ColVec_t& amplis, ColVec_t* grad_out, Mat_t* hess_out, void* opt_data) {
  // amplitudes vals_inp (Re(a_1), Im(a_1), Re(a_2), Im(a_2),...)
  //  size = 2*sources.size()
  //  INFO("pattern_energy");
   wp_data* data = reinterpret_cast<wp_data*>(opt_data);
    
   std::vector<EquivalentSource*> sources = data->sources;
   std::vector<PatternPoint>  pattern_pts = data->pattern_pts;
   std::list<Wave*> in = data->in;
   ColVec_t c(pattern_pts.size());// = data-> pp_ampli;
   for (uint i = 0; i < pattern_pts.size(); ++i) {
     c(i) = abs(pattern_pts[i].getAmpli());
   }
  MATX Phi = data->transfer_mat; //M(2*sources.size(), pattern_pts.size()

  //  INFO("PHI "<<Phi);
  // grad_out = new ColVec_t(Phi.rows());
  // hess_out = new Mat_t(Phi.rows(), Phi.rows());

  if (grad_out && hess_out) {
  for (int i = 0; i < Phi.cols(); ++i) {
    (*grad_out)(i) = 0;
    for (int j = 0; j < Phi.cols(); ++j) {
      (*hess_out)(i, j) = 0;
    }
  }
  }

  double E = 0;
  double c0 = 1, c1 = 1;
  for (int j = 0; j < c.size(); ++j) {
    double win_re = 0, win_im = 0;
    VEC2 ppos = pattern_pts[j].getPos();
    for (auto it: in) {
      win_re += real(it->heightc(ppos,0));
      win_im += imag(it->heightc(ppos,0));
    }
    double aphi_re = 0;
    double aphi_im = 0;
    for (int i = 0; i < Phi.cols()/2; ++i) {
      double ai = amplis(2*i);
      double bi = amplis(2*i+1);
      double mui = Phi(j, 2*i); 
      double nui = Phi(j, 2*i+1);
      aphi_re += ai*mui - bi*nui;
      aphi_im += ai*nui + bi*mui;
      //       INFO("ai:"<<ai<<" bi:"<<bi<<" mui:"<<mui<<" nui"<<nui);
      //  INFO("aphi_re:"<<aphi_re<<" aphi_im:"<<aphi_im<<" ampli:"<<sqrt(aphi_re*aphi_re + aphi_im*aphi_im));
    }
    aphi_re += win_re;
    aphi_im += win_im;
    double epsj = aphi_re*aphi_re + aphi_im*aphi_im - c(j)*c(j);
    if (c(j) == 0) {
    //   INFO("epsj: "<<epsj<<"  cj"<<c(j));
      E += c0*epsj*epsj;
    } else {
      E += c1*epsj*epsj;
    }
    //    E += epsj*epsj;

    if (grad_out) {
    for (int h = 0; h < Phi.cols()/2; ++h) {
      double muh = Phi(j, 2*h);
      double nuh = Phi(j, 2*h+1);
      double d_epsj_da = muh*2*aphi_re + nuh*2*aphi_im;
      // INFO("depsjda: "<<d_epsj_da);
      double d_epsj_db = -nuh*2*aphi_re + muh*2*aphi_im;
      //      INFO("depsjdb: "<<d_epsj_db);
      if (c(j) == 0) {
      	(*grad_out)(2*h) += c0*2*d_epsj_da*epsj;
      	(*grad_out)(2*h+1) += c0*2*d_epsj_db*epsj;
      } else {
      	(*grad_out)(2*h) += c1*2*d_epsj_da*epsj;
      	(*grad_out)(2*h+1) += c1*2*d_epsj_db*epsj;
      }

      // (*grad_out)(2*h) += 2*d_epsj_da*epsj;
      // (*grad_out)(2*h+1) += 2*d_epsj_db*epsj;
    }
    }
    //    INFO("grad \n"<<(*grad_out)(0));

    if (hess_out) {
    for (int h = 0; h < Phi.cols()/2; ++h) {
      	double muh = Phi(j, 2*h);
	double nuh = Phi(j, 2*h+1);
	double d_epsj_dah = muh*2*aphi_re + nuh*2*aphi_im;
	double d_epsj_dbh = -nuh*2*aphi_re + muh*2*aphi_im;
		
      for (int k = 0; k < Phi.cols()/2; ++k) {
	double muk = Phi(j, 2*k);
	double nuk = Phi(j, 2*k+1);
	double d2_epsj_dakdah = 2*muh*muk + 2*nuh*nuk;
	double d2_epsj_dbkdbh = 2*nuh*nuk + 2*muh*muk;
	double d2_epsj_dbkdah = -2*muh*nuk + 2*nuh*muk;
	double d2_epsj_dakdbh = -2*nuh*muk + 2*muh*nuk;
	double d_epsj_dak = muk*2*aphi_re + nuk*2*aphi_im;
	double d_epsj_dbk = -nuk*2*aphi_re + muk*2*aphi_im;

	if (c(j) == 0) {
	  (*hess_out)(2*k, 2*h) += c0*(2*d2_epsj_dakdah*epsj + 2*d_epsj_dah*d_epsj_dak);
	  (*hess_out)(2*k+1, 2*h+1) += c0*(2*d2_epsj_dbkdbh*epsj + 2*d_epsj_dbh*d_epsj_dbk);
	  (*hess_out)(2*k+1, 2*h) += c0*(2*d2_epsj_dbkdah*epsj + 2*d_epsj_dah*d_epsj_dbk);
	  (*hess_out)(2*k, 2*h+1) += c0*(2*d2_epsj_dakdbh*epsj + 2*d_epsj_dbh*d_epsj_dak);
	} else {
	  (*hess_out)(2*k, 2*h) += c1*(2*d2_epsj_dakdah*epsj + 2*d_epsj_dah*d_epsj_dak);
	  (*hess_out)(2*k+1, 2*h+1) += c1*(2*d2_epsj_dbkdbh*epsj + 2*d_epsj_dbh*d_epsj_dbk);
	  (*hess_out)(2*k+1, 2*h) += c1*(2*d2_epsj_dbkdah*epsj + 2*d_epsj_dah*d_epsj_dbk);
	  (*hess_out)(2*k, 2*h+1) += c1*(2*d2_epsj_dakdbh*epsj + 2*d_epsj_dbh*d_epsj_dak);
	}

	// (*hess_out)(2*k, 2*h) += 2*d2_epsj_dakdah*epsj + 2*d_epsj_dah*d_epsj_dak;
	// (*hess_out)(2*k+1, 2*h+1) += 2*d2_epsj_dbkdbh*epsj + 2*d_epsj_dbh*d_epsj_dbk;
	// (*hess_out)(2*k+1, 2*h) += 2*d2_epsj_dbkdah*epsj + 2*d_epsj_dah*d_epsj_dbk;
	// (*hess_out)(2*k, 2*h+1) += 2*d2_epsj_dakdbh*epsj + 2*d_epsj_dbh*d_epsj_dak;


      }

      
    }
    }
    //    INFO("hess\n"<<(*hess_out)(0, 0));
  }
    INFO("energy "<<E);
  // INFO("grad \n"<<*grad_out);
  // INFO("hess\n"<<*hess_out);
  return E;
}


 double pattern_energy(const ColVec_t& amplis, ColVec_t* grad_out, void* opt_data) {
  // amplitudes vals_inp (Re(a_1), Im(a_1), Re(a_2), Im(a_2),...)
  wp_data* data = reinterpret_cast<wp_data*>(opt_data);

     std::vector<EquivalentSource*> sources = data->sources;
  std::vector<PatternPoint>  pattern_pts = data->pattern_pts;
  std::list<Wave*> in = data->in;
  ColVec_t c(pattern_pts.size());// = data-> pp_ampli;
  for (uint i = 0; i < pattern_pts.size(); ++i) {
    c(i) = abs(pattern_pts[i].getAmpli());
  }
  MATX Phi = data->transfer_mat; //M(2*sources.size(), pattern_pts.size()
  
  // ColVec_t c = data-> pp_ampli;
  // MATX Phi = data->transfer_mat; //M(2*sources.size(), pattern_pts.size()

  // grad_out = new ColVec_t(Phi.rows());
  // hess_out = new Mat_t(Phi.rows(), Phi.rows());

  if (grad_out) {
  for (int i = 0; i < Phi.cols(); ++i) {
    (*grad_out)(i) = 0;
  }
  }

   double E = 0;
  for (int j = 0; j < c.size(); ++j) {
       double win_re = 0, win_im = 0;
    VEC2 ppos = pattern_pts[j].getPos();
    for (auto it: in) {
      win_re += real(it->heightc(ppos,0));
      win_im += imag(it->heightc(ppos,0));
    }
    double aphi_re = 0;
    double aphi_im = 0;
    for (int i = 0; i < Phi.cols()/2; ++i) {
      double ai = amplis(2*i);
      double bi = amplis(2*i+1);
      double mui = Phi(j, 2*i); 
      double nui = Phi(j, 2*i+1);
      aphi_re += ai*mui - bi*nui;
      aphi_im += ai*nui + bi*mui;
      // INFO("ai:"<<ai<<" bi:"<<bi<<" mui:"<<mui<<" nui"<<nui);
      //  INFO("aphi_re:"<<aphi_re<<" aphi_im:"<<aphi_im<<" ampli:"<<sqrt(aphi_re*aphi_re + aphi_im*aphi_im));
    }
         aphi_re += win_re;
     aphi_im += win_im;
    double epsj = aphi_re*aphi_re + aphi_im*aphi_im - c(j)*c(j);
    //   INFO("epsj: "<<epsj<<"  cj"<<c(j));
    E += epsj*epsj;

    if (grad_out) {
    for (int h = 0; h < Phi.cols()/2; ++h) {
      double muh = Phi(j, 2*h);
      double nuh = Phi(j, 2*h+1);
      double d_epsj_da = muh*2*aphi_re + nuh*2*aphi_im;
      // INFO("depsjda: "<<d_epsj_da);
      double d_epsj_db = -nuh*2*aphi_re + muh*2*aphi_im;
      //      INFO("depsjdb: "<<d_epsj_db);
      (*grad_out)(2*h) += 2*d_epsj_da*epsj;
      (*grad_out)(2*h+1) += 2*d_epsj_db*epsj;
    }
    }
    //    INFO("grad \n"<<(*grad_out)(0));
  }
     INFO("energy "<<E);
//   // INFO("grad \n"<<*grad_out);
//   // INFO("hess\n"<<*hess_out);
   return E;
 }

void gradPhi(COMPLEX & phi, VEC2C & grad, FLOAT xs, FLOAT ys, FLOAT xp, FLOAT yp, COMPLEX a, FLOAT k) {
  FLOAT rx = xp - xs;
  FLOAT ry = yp - ys;
  FLOAT r  = sqrt(pow(rx, 2.0) + pow(ry, 2.0));
  COMPLEX out_x(0, 0);
  COMPLEX out_y(0, 0);
  phi = 0;
  if (r != 0) {
    // FLOAT der_damp = 0;//damping(r, wave_number);//exp(-damping*pow(wave_number, 2)*r);
    FLOAT cos_phi = rx/r;
    FLOAT sin_phi = ry/r;
    out_x = cos_phi* (COMPLEX(0, -1)/(FLOAT)4.0*k*derHankel(k*r));
	//	-sin_phi/r*der_damp*(-i_/(FLOAT)4.0*Hankel(wave_number*r));
    out_y = sin_phi*(COMPLEX(0, -1)/(FLOAT)4.0*k*derHankel(k*r));
      //+cos_phi/r*der_damp*(-i_/(FLOAT)4.0*Hankel(wave_number*r));
     phi =  COMPLEX(0, -1)/(FLOAT)4.0*Hankel(k*r);
  }
  grad(0) = out_x;
  grad(1) = out_y;
 
  INFO("heightc phi "<<phi<<" "<<k<<" "<<r<<" "<<xp<<" "<<yp<<" "<<xs<<" "<<ys);
  return ;
}

 double pattern_energy_pos_hess(const ColVec_t& pos, ColVec_t* grad_out, Mat_t* hess_out, void* opt_data) {
   // amplitudes vals_inp (Re(a_1), Im(a_1), Re(a_2), Im(a_2),...)
   //  size = 2*sources.size()
   //  INFO("pattern_energy");
    wp_data* data = reinterpret_cast<wp_data*>(opt_data);
    
    std::vector<EquivalentSource*> sources = data->sources;
   std::vector<PatternPoint>  pattern_pts = data->pattern_pts;
   std::list<Wave*> in = data->in;
   ColVec_t c(pattern_pts.size());// = data-> pp_ampli;
   for (uint i = 0; i < pattern_pts.size(); ++i) {
     c(i) = abs(pattern_pts[i].getAmpli());
   }
   // MATX Phi = data->transfer_mat; //M(2*sources.size(), pattern_pts.size()
   FLOAT k = data->k;
   int ns = sources.size();
   int np = pattern_pts.size();
   
   // grad_out = new ColVec_t(Phi.rows());
   // hess_out = new Mat_t(Phi.rows(), Phi.rows());

   if (grad_out && hess_out) {
   for (int i = 0; i < 2*ns; ++i) {
     (*grad_out)(i) = 0;
     for (int j = 0; j < 2*ns; ++j) {
       (*hess_out)(i, j) = 0;
     }
   }
   }

   MATX Phi(np, 2*ns);
   MATX Phi_dx(np, 2*ns);
   MATX Phi_dy(np, 2*ns);

   for (int j = 0; j < np; ++j) {
     VEC2 ppos = pattern_pts[j].getPos();
     for (int i = 0; i < ns; ++i) {
       double sx = pos(2*i);
       double sy = pos(2*i+1);
       COMPLEX phi;
       VEC2C grad;
       gradPhi(phi, grad, sx, sy, ppos(0), ppos(1), sources[i]->getAmpli(), k);
       Phi(j, 2*i) = real(phi);
       Phi(j, 2*i+1) = imag(phi);
       Phi_dx(j, 2*i) = real(grad(0));
       Phi_dx(j, 2*i+1) = imag(grad(0));
       Phi_dy(j, 2*i) = real(grad(1));
       Phi_dy(j, 2*i+1) = imag(grad(1));
     }
   }

   //   INFO("PHI "<<Phi);

   double E = 0;
   double c0 = 1, c1 = 1;
   for (int j = 0; j < np; ++j) {
     double win_re = 0, win_im = 0;
     VEC2 ppos = pattern_pts[j].getPos();
     for (auto it: in) {
       win_re += real(it->heightc(ppos,0));
       win_im += imag(it->heightc(ppos,0));
     }
     double aphi_re = 0;
     double aphi_im = 0;
     for (int i = 0; i < ns; ++i) {
       double ai = real(sources[i]->getAmpli());
       double bi = imag(sources[i]->getAmpli());
       double mui = Phi(j, 2*i); 
       double nui = Phi(j, 2*i+1);
       aphi_re += ai*mui - bi*nui;
       aphi_im += ai*nui + bi*mui;
          INFO("ai:"<<ai<<" bi:"<<bi<<" mui:"<<mui<<" nui"<<nui);
         INFO("aphi_re:"<<aphi_re<<" aphi_im:"<<aphi_im<<" ampli:"<<sqrt(aphi_re*aphi_re + aphi_im*aphi_im));
     }
     aphi_re += win_re;
     aphi_im += win_im;
     double epsj = aphi_re*aphi_re + aphi_im*aphi_im - c(j)*c(j);
     if (c(j) == 0) {
     //   INFO("epsj: "<<epsj<<"  cj"<<c(j));
       E += c0*epsj*epsj;
     } else {
       E += c1*epsj*epsj;
     }
     //    E += epsj*epsj;

     if (grad_out) {
     for (int h = 0; h < ns; ++h) {
       double muh = Phi(j, 2*h);
       double nuh = Phi(j, 2*h+1);
       double ah = real(sources[h]->getAmpli());
       double bh = imag(sources[h]->getAmpli());
       double d_aphi_re_dx = ah*Phi_dx(j, 2*h) - bh*Phi_dx(j, 2*h+1);
       double d_aphi_re_dy = ah*Phi_dy(j, 2*h) - bh*Phi_dy(j, 2*h+1);
       double d_aphi_im_dx = bh*Phi_dx(j, 2*h) + ah*Phi_dx(j, 2*h+1);
       double d_aphi_im_dy = bh*Phi_dy(j, 2*h) + ah*Phi_dy(j, 2*h+1);
       double d_epsj_dx = 2*d_aphi_re_dx*aphi_re + 2*d_aphi_im_dx*aphi_im;
       double d_epsj_dy = 2*d_aphi_re_dy*aphi_re + 2*d_aphi_im_dy*aphi_im; 
//       // INFO("depsjda: "<<d_epsj_da);
//       double d_epsj_db = -nuh*2*aphi_re + muh*2*aphi_im;
//       //      INFO("depsjdb: "<<d_epsj_db);
       if (c(j) == 0) {
       	(*grad_out)(2*h) += c0*2*d_epsj_dx*epsj;
       	(*grad_out)(2*h+1) += c0*2*d_epsj_dy*epsj;
       } else {
       	(*grad_out)(2*h) += c1*2*d_epsj_dx*epsj;
       	(*grad_out)(2*h+1) += c1*2*d_epsj_dy*epsj;
       }

       // (*grad_out)(2*h) += 2*d_epsj_da*epsj;
       // (*grad_out)(2*h+1) += 2*d_epsj_db*epsj;
     }
     }
     //     INFO("grad \n"<<(*grad_out)(0)<<" "<<(*grad_out)(1));

     if (hess_out) {
     for (int h = 0; h < ns; ++h) {
       	 double muh = Phi(j, 2*h);
	 double nuh = Phi(j, 2*h+1);
	 double ah = real(sources[h]->getAmpli());
	 double bh = imag(sources[h]->getAmpli());
	 double d_aphi_re_dxh = ah*Phi_dx(j, 2*h) - bh*Phi_dx(j, 2*h+1);
	 double d_aphi_re_dyh = ah*Phi_dy(j, 2*h) - bh*Phi_dy(j, 2*h+1);
	 double d_aphi_im_dxh = bh*Phi_dx(j, 2*h) + ah*Phi_dx(j, 2*h+1);
	 double d_aphi_im_dyh = bh*Phi_dy(j, 2*h) + ah*Phi_dy(j, 2*h+1);
	 double d_epsj_dxh = 2*d_aphi_re_dxh*aphi_re + 2*d_aphi_im_dxh*aphi_im;
	 double d_epsj_dyh = 2*d_aphi_re_dyh*aphi_re + 2*d_aphi_im_dyh*aphi_im;


       for (int k = 0; k < ns; ++k) {
	 double muk = Phi(j, 2*k);
	 double nuk = Phi(j, 2*k+1);
	 double ak = real(sources[k]->getAmpli());
	 double bk = imag(sources[k]->getAmpli());
	 double d_aphi_re_dxk = ak*Phi_dx(j, 2*k) - bk*Phi_dx(j, 2*k+1);
	 double d_aphi_re_dyk = ak*Phi_dy(j, 2*k) - bk*Phi_dy(j, 2*k+1);
	 double d_aphi_im_dxk = bk*Phi_dx(j, 2*k) + ak*Phi_dx(j, 2*k+1);
	 double d_aphi_im_dyk = bk*Phi_dy(j, 2*k) + ak*Phi_dy(j, 2*k+1);
	 double d_epsj_dxk = 2*d_aphi_re_dxk*aphi_re + 2*d_aphi_im_dxk*aphi_im;
	 double d_epsj_dyk = 2*d_aphi_re_dyk*aphi_re + 2*d_aphi_im_dyk*aphi_im; 
	 
	 double d2_epsj_dxkdxh = 2*d_aphi_re_dxh*d_aphi_re_dxk + 2*d_aphi_im_dxh*d_aphi_im_dxk;
	 double d2_epsj_dykdyh = 2*d_aphi_re_dyh*d_aphi_re_dyk + 2*d_aphi_im_dyh*d_aphi_im_dyk;
	 double d2_epsj_dxkdyh = 2*d_aphi_re_dyh*d_aphi_re_dxk + 2*d_aphi_im_dyh*d_aphi_im_dxk;
	 double d2_epsj_dykdxh = 2*d_aphi_re_dxh*d_aphi_re_dyk + 2*d_aphi_im_dxh*d_aphi_im_dyk;

 	if (c(j) == 0) {
 	  (*hess_out)(2*k, 2*h) += c0*(2*d2_epsj_dxkdxh*epsj + 2*d_epsj_dxh*d_epsj_dxk);
 	  (*hess_out)(2*k+1, 2*h+1) += c0*(2*d2_epsj_dykdyh*epsj + 2*d_epsj_dyh*d_epsj_dyk);
 	  (*hess_out)(2*k+1, 2*h) += c0*(2*d2_epsj_dykdxh*epsj + 2*d_epsj_dxh*d_epsj_dyk);
 	  (*hess_out)(2*k, 2*h+1) += c0*(2*d2_epsj_dxkdyh*epsj + 2*d_epsj_dyh*d_epsj_dxk);
	} else {
	  (*hess_out)(2*k, 2*h) += c1*(2*d2_epsj_dxkdxh*epsj + 2*d_epsj_dxh*d_epsj_dxk);
 	  (*hess_out)(2*k+1, 2*h+1) += c1*(2*d2_epsj_dykdyh*epsj + 2*d_epsj_dyh*d_epsj_dyk);
 	  (*hess_out)(2*k+1, 2*h) += c1*(2*d2_epsj_dykdxh*epsj + 2*d_epsj_dxh*d_epsj_dyk);
 	  (*hess_out)(2*k, 2*h+1) += c1*(2*d2_epsj_dxkdyh*epsj + 2*d_epsj_dyh*d_epsj_dxk);
 	}

 	// (*hess_out)(2*k, 2*h) += 2*d2_epsj_dakdah*epsj + 2*d_epsj_dah*d_epsj_dak;
 	// (*hess_out)(2*k+1, 2*h+1) += 2*d2_epsj_dbkdbh*epsj + 2*d_epsj_dbh*d_epsj_dbk;
 	// (*hess_out)(2*k+1, 2*h) += 2*d2_epsj_dbkdah*epsj + 2*d_epsj_dah*d_epsj_dbk;
 	// (*hess_out)(2*k, 2*h+1) += 2*d2_epsj_dakdbh*epsj + 2*d_epsj_dbh*d_epsj_dak;


       }

      
     }
     }
     //     INFO("hess\n"<<(*hess_out));
   }
     INFO("energy "<<E);
//   // INFO("grad \n"<<*grad_out);
//   // INFO("hess\n"<<*hess_out);
   return E;
 }


 double pattern_energy_pos(const ColVec_t& pos, ColVec_t* grad_out, void* opt_data) {
   // amplitudes vals_inp (Re(a_1), Im(a_1), Re(a_2), Im(a_2),...)
   //  size = 2*sources.size()
   //  INFO("pattern_energy");
    wp_data* data = reinterpret_cast<wp_data*>(opt_data);
    
    std::vector<EquivalentSource*> sources = data->sources;
   std::vector<PatternPoint>  pattern_pts = data->pattern_pts;
   std::list<Wave*> in = data->in;
   ColVec_t c(pattern_pts.size());// = data-> pp_ampli;
   for (uint i = 0; i < pattern_pts.size(); ++i) {
     c(i) = abs(pattern_pts[i].getAmpli());
   }
   // MATX Phi = data->transfer_mat; //M(2*sources.size(), pattern_pts.size()
   FLOAT k = data->k;
   int ns = sources.size();
   int np = pattern_pts.size();
   
   // grad_out = new ColVec_t(Phi.rows());
   // hess_out = new Mat_t(Phi.rows(), Phi.rows());

   if (grad_out) {
   for (int i = 0; i < 2*ns; ++i) {
     (*grad_out)(i) = 0;
   }
   }

   MATX Phi(np, 2*ns);
   MATX Phi_dx(np, 2*ns);
   MATX Phi_dy(np, 2*ns);

   for (int j = 0; j < np; ++j) {
     VEC2 ppos = pattern_pts[j].getPos();
     for (int i = 0; i < ns; ++i) {
       double sx = pos(2*i);
       double sy = pos(2*i+1);
       COMPLEX phi;
       VEC2C grad;
       gradPhi(phi, grad, sx, sy, ppos(0), ppos(1), sources[i]->getAmpli(), k);
       Phi(j, 2*i) = real(phi);
       Phi(j, 2*i+1) = imag(phi);
       Phi_dx(j, 2*i) = real(grad(0));
       Phi_dx(j, 2*i+1) = imag(grad(0));
       Phi_dy(j, 2*i) = real(grad(1));
       Phi_dy(j, 2*i+1) = imag(grad(1));
     }
   }

   //   INFO("PHI "<<Phi);

   double E = 0;
   double c0 = 1, c1 = 1;
   for (int j = 0; j < np; ++j) {
     double win_re = 0, win_im = 0;
     VEC2 ppos = pattern_pts[j].getPos();
     for (auto it: in) {
       win_re += real(it->heightc(ppos,0));
       win_im += imag(it->heightc(ppos,0));
     }
     double aphi_re = 0;
     double aphi_im = 0;
     for (int i = 0; i < ns; ++i) {
       double ai = real(sources[i]->getAmpli());
       double bi = imag(sources[i]->getAmpli());
       double mui = Phi(j, 2*i); 
       double nui = Phi(j, 2*i+1);
       aphi_re += ai*mui - bi*nui;
       aphi_im += ai*nui + bi*mui;
          INFO("ai:"<<ai<<" bi:"<<bi<<" mui:"<<mui<<" nui"<<nui);
         INFO("aphi_re:"<<aphi_re<<" aphi_im:"<<aphi_im<<" ampli:"<<sqrt(aphi_re*aphi_re + aphi_im*aphi_im));
     }
     aphi_re += win_re;
     aphi_im += win_im;
     double epsj = aphi_re*aphi_re + aphi_im*aphi_im - c(j)*c(j);
     if (c(j) == 0) {
     //   INFO("epsj: "<<epsj<<"  cj"<<c(j));
       E += c0*epsj*epsj;
     } else {
       E += c1*epsj*epsj;
     }
     //    E += epsj*epsj;

     if (grad_out) {
     for (int h = 0; h < ns; ++h) {
       double muh = Phi(j, 2*h);
       double nuh = Phi(j, 2*h+1);
       double ah = real(sources[h]->getAmpli());
       double bh = imag(sources[h]->getAmpli());
       double d_aphi_re_dx = ah*Phi_dx(j, 2*h) - bh*Phi_dx(j, 2*h+1);
       double d_aphi_re_dy = ah*Phi_dy(j, 2*h) - bh*Phi_dy(j, 2*h+1);
       double d_aphi_im_dx = bh*Phi_dx(j, 2*h) + ah*Phi_dx(j, 2*h+1);
       double d_aphi_im_dy = bh*Phi_dy(j, 2*h) + ah*Phi_dy(j, 2*h+1);
       double d_epsj_dx = 2*d_aphi_re_dx*aphi_re + 2*d_aphi_im_dx*aphi_im;
       double d_epsj_dy = 2*d_aphi_re_dy*aphi_re + 2*d_aphi_im_dy*aphi_im; 
//       // INFO("depsjda: "<<d_epsj_da);
//       double d_epsj_db = -nuh*2*aphi_re + muh*2*aphi_im;
//       //      INFO("depsjdb: "<<d_epsj_db);
       if (c(j) == 0) {
       	(*grad_out)(2*h) += c0*2*d_epsj_dx*epsj;
       	(*grad_out)(2*h+1) += c0*2*d_epsj_dy*epsj;
       } else {
       	(*grad_out)(2*h) += c1*2*d_epsj_dx*epsj;
       	(*grad_out)(2*h+1) += c1*2*d_epsj_dy*epsj;
       }

       // (*grad_out)(2*h) += 2*d_epsj_da*epsj;
       // (*grad_out)(2*h+1) += 2*d_epsj_db*epsj;
     }
     }
     //     INFO("grad \n"<<(*grad_out)(0)<<" "<<(*grad_out)(1));
   }
     INFO("energy "<<E);
//   // INFO("grad \n"<<*grad_out);
//   // INFO("hess\n"<<*hess_out);
   return E;
 }

WavePattern::WavePattern() {
  pattern_file = "";
  pos = VEC2(0, 0);
  n_rows = 200;
  n_cols = 200;
  cell_size = 0.15;
  wave_length = 0;
  wave_number = 0;
  sources.clear();
  pattern_pts.clear();
  svd = NULL;
  init_val = 1;
  ampli = 0;
  method = 0;
}

WavePattern::WavePattern(FLOAT wl): WavePattern() {
  wave_length = wl;
  wave_number = 2*M_PI/wave_length; 
}


WavePattern::~WavePattern() {
  for (auto it: sources) {
    delete it;
  }
  sources.clear();
  pattern_pts.clear();
  if (svd != NULL) {
    delete svd;
  }
}

void WavePattern::setPos(VEC2 p) {
  pos = p;
}

void WavePattern::setPos(FLOAT x, FLOAT y) {
  pos = VEC2(x, y);
}

VEC2 WavePattern::getPos() const {
  return pos;
}

void WavePattern::setSize(int nr, int nc, FLOAT cs) {
  n_rows = nr;
  n_cols = nc;
  cell_size = cs;
}
int WavePattern::getNbCols() const {
  return n_cols;
}

int WavePattern::getNbRows() const {
  return n_rows;
}
FLOAT WavePattern::getCellSize() const {
  return cell_size;
}

void WavePattern::setAmpli(FLOAT a) {
  ampli = a;
}

void WavePattern::setMethod(uint m) {
  method = m;
}

void WavePattern::setInitVal(FLOAT v) {
  init_val = v;
}

void WavePattern::createPatternPoints(std::string file) {
  pattern_pts.clear();
  
  // SDL_Surface *height_field;
  // height_field = IMG_Load(file.c_str());
  // if(height_field == 0) {
  //   std::cout << "Erreur : " << SDL_GetError() << std::endl;
  //   std::exit(1);
  // }

  // INFO("Grid pattern :"<< n_rows<<" "<< n_cols<<" "<< cell_size);
  // grid = Grid(n_rows, n_cols, cell_size);
  // grid.setCellSize(cell_size);
      
  // grid.setValues(height_field);

  // std::vector<PatternPoint> pts;

  // uint n0 = 0, n1 = 0;
  
  // for (int i = 1; i < grid.getNbRows()-1; ++i) {
  //   for (int j = 1; j < grid.getNbCols()-1; ++j) {
  //     if (grid(i, j) == 0) {
  // 	VEC2 p(((FLOAT)i-(FLOAT)n_rows/2.0)*cell_size + pos(0),
  // 	       ((FLOAT)j-(FLOAT)n_cols/2.0)*cell_size + pos(1));
  // 	PatternPoint pp(p, COMPLEX(0,0));
  // 	pattern_pts.push_back(pp);
  // 	++n0;
  //      } else if (grid(i, j) != 1) {
  //      	VEC2 p(((FLOAT)i-(FLOAT)n_rows/2.0)*cell_size + pos(0),
  //      	       ((FLOAT)j-(FLOAT)n_cols/2.0)*cell_size + pos(1));
  //      	PatternPoint pp(p, COMPLEX(ampli,0));
  //      	pattern_pts.push_back(pp);
  // 	++n1;
  //     }
  //   }
  // }
  // INFO("Nb pattern pts "<< pattern_pts.size()<<"  ampli 0:"<<n0<<"  not ampli 0:"<<n1);
  // SDL_FreeSurface(height_field);

   VEC2 p(15, 15);
   PatternPoint pp(p, COMPLEX(ampli,0));
   pattern_pts.push_back(pp);
   VEC2 p2(16, 15);
   PatternPoint pp2(p2, COMPLEX(ampli,0));
   pattern_pts.push_back(pp2);
}


void WavePattern::createSources(FLOAT r) {
  // radius = r;
  // FLOAT sample_rate_boundary = step_sampling_*wave_length;
  // FLOAT angle_step = sample_rate_boundary/radius;
  // int nb_sources = 2*M_PI/angle_step +1;
  // if (nb_sources < 4) {
  //   nb_sources = 4;
  // }

  // angle_step = 2*M_PI/(FLOAT)nb_sources;
  // VEC2 dir(1, 0);
  // MAT2 rotation;
  // rotation <<
  //   cos(angle_step), -sin(angle_step),
  //   sin(angle_step), cos(angle_step);

  // EquivalentSource* es;
  
  // for (int i = 0; i < nb_sources; ++i) {
  //   es = new EquivalentSource(wave_length, 1); // todo: set right ampli_step
  //   es->setPos(pos + r*dir);
  //   sources.push_back(es);
      
  //   dir = rotation*dir;
  // }
  //   INFO("Nb sources "<< sources.size());

   EquivalentSource* es;
    es = new EquivalentSource(wave_length, 1);
    es->setPos(pos);
    sources.push_back(es);
    // EquivalentSource* es2;
    // es2 = new EquivalentSource(wave_length, 1);
    // es2->setPos(pos+ VEC2(0.1, 0));
    // sources.push_back(es2);
}


void WavePattern::setTransferMatrix() {
  //  INFO("Set Transfer matrix:  "<<index<<"  (wavelength "<<wave_lenghts[index]<<", "<<sources_l[index].size()<<" sources)");
  if (pattern_pts.size() != 0) {
  
    MatrixXcf T = MatrixXcf(pattern_pts.size(), sources.size());
    for (uint i = 0; i < pattern_pts.size(); ++i) {
      for (uint j = 0; j < sources.size(); ++j) {
	T(i, j) = sources[j]->heightc(pattern_pts[i].getPos());
	}
      }    
  //    INFO("SVD..."<<index<<"  (wavelength "<<wave_lenghts[index]<<")");
     svd = new BDCSVD<MatrixXcf>(T,ComputeThinU | ComputeThinV);
     //     transfer_mat = T;
  }
    //    INFO("SVD...done:  "<<index<<"  (wavelength "<<wave_lenghts[index]<<")");
}

void WavePattern::solve(std::list<Wave*> in) {

  //    energy(in);
  
 // VectorXcf c(sources.size());
 // VectorXcf p_in(pattern_pts.size());

 //    //    bool to_solve = false;
 //    for (uint i = 0; i < pattern_pts.size(); ++i) {
 //      p_in(i) = pattern_pts[i].getAmpli();
 //      // VEC2 ppos = pattern_pts[i].getPos();
 //      // for (auto it: in) {
 //      // 	p_in(i) -= it->heightc(ppos,0);
 //      // }
 //      //VEC2 pb = pattern_pts[i]->getPos();
 //    }
 //    // for (uint i = 0; i < boundaries_l[index].size() && !to_solve; ++i) {
 //    //   to_solve = to_solve || (norm(p_in(i)) > 1e-8*wave_lenghts[index]);
 //    // }
 //    //if (to_solve) {
 //    // INFO(p_in);
 //     c = svd->solve(p_in);
  
  ColVec_t pos(2*sources.size());
  ColVec_t amplis(2*sources.size());
 for (int i = 0; i < sources.size(); ++i) {
   pos(2*i) = sources[i]->getPos()(0);
   pos(2*i+1) = sources[i]->getPos()(1);
   amplis(2*i) = init_val;//real(c(i));
   amplis(2*i+1) = init_val;//imag(c(i));
 }

 MATX T(pattern_pts.size(), 2*sources.size());
 for (uint i = 0; i < pattern_pts.size(); ++i) {
   for (uint j = 0; j < sources.size(); ++j) {
     T(i, 2*j) = real(sources[j]->heightc(pattern_pts[i].getPos()));
     T(i, 2*j+1) = imag(sources[j]->heightc(pattern_pts[i].getPos()));
   }
 }

 wp_data opt_data;
 opt_data.pattern_pts = pattern_pts;
 opt_data.sources = sources;
 opt_data.in = in;
 opt_data.k = wave_number;
 opt_data.transfer_mat = T;

    //  ColVec_t grad_out(2*sources.size());
    // Mat_t hess_out(2*sources.size(), 2*sources.size());
    // double E = pattern_energy_hess(amplis, &grad_out, &hess_out, &opt_data);
    // double E2 = pattern_energy_pos_hess(pos, &grad_out, &hess_out, &opt_data);

    // INFO("energy "<<E<<" "<<E2);
    
 
 
 optim::algo_settings_t settings;
 settings.print_level = 4;
 settings.conv_failure_switch = 2;
 settings.iter_max = 2000;

   bool success;
  //    if (method == 0) {
  //      success = optim::newton(amplis,pattern_energy_hess, &opt_data, settings);
  //   } else if (method == 1){
  //     success = optim::bfgs(amplis,pattern_energy, &opt_data, settings);
  //   } else if (method == 2) {
  //     success = optim::nm(amplis,pattern_energy, &opt_data, settings);
  //   } else {
  //     success = optim::de_prmm(amplis,pattern_energy, &opt_data, settings);
  //   }
    
  //   if(success) {
  //     // INFO("pp_amplis \n"<<pp_ampli);
  //        // INFO("trans mat \n"<<T);
  //           INFO("                    SUCCESS amplis \n"<<amplis);
  // 	    //INFO("                    SUCCESS amplis");
  //   } else {
  //     INFO("                    NOPE");
  //   }
  
   for (uint i = 0; i < sources.size(); ++i) {
     sources[i]->setAmplitude(COMPLEX(amplis(2*i), amplis(2*i+1)));
   }
    
    opt_data.pattern_pts = pattern_pts;
    opt_data.sources = sources;
    opt_data.in = in;
    opt_data.k = wave_number;
    opt_data.transfer_mat = T;

    // ColVec_t grad_out(2*sources.size());
    // Mat_t hess_out(2*sources.size(), 2*sources.size());
    // double E = pattern_energy_hess(amplis, &grad_out, &hess_out, &opt_data);
    // double E2 = pattern_energy_pos_hess(pos, &grad_out, &hess_out, &opt_data);

    //  INFO("energy "<<E<<" "<<E2);

    
     if (method == 0) {
        success = optim::newton(pos,pattern_energy_pos_hess, &opt_data, settings);
     } else if (method == 1){
       success = optim::bfgs(pos,pattern_energy_pos, &opt_data, settings);
     } else if (method == 2) {
       success = optim::nm(pos,pattern_energy_pos, &opt_data, settings);
     } else {
       success = optim::de_prmm(pos,pattern_energy_pos, &opt_data, settings);
     }
    
  //  success = optim::newton(pos,pattern_energy_pos_hess, &opt_data, settings);
    if(success) {
   // INFO("pp_amplis \n"<<pp_ampli);
         // INFO("trans mat \n"<<T);
      // INFO("                    SUCCESS amplis \n"<<amplis);
      INFO("                    SUCCESS pos \n"<<pos);
 } else {
   INFO("                    NOPE");
 }

 for (uint i = 0; i < sources.size(); ++i) {
   sources[i]->setPos(pos(2*i), pos(2*i+1));
 }
  
#ifdef OPTIM_AMPLI
  ColVec_t pp_ampli(pattern_pts.size());
  ColVec_t amplis(2*sources.size());
  
   // ColVec_t grad_out(2*sources.size());
   // Mat_t hess_out(2*sources.size(), 2*sources.size());

  for (int i = 0; i < sources.size(); ++i) {
    amplis(2*i) = init_val;//real(c(i));
    amplis(2*i+1) = init_val;//imag(c(i));
  }
  
  MATX T(pattern_pts.size(), 2*sources.size());
  // MATX T_dx(pattern_pts.size(), 2*sources.size());
  // MATX T_dy(pattern_pts.size(), 2*sources.size());
  
  for (uint i = 0; i < pattern_pts.size(); ++i) {
     pp_ampli(i) = real(pattern_pts[i].getAmpli());
      // VEC2 ppos = pattern_pts[i].getPos();
      //  for (auto it: in) {
      // 	pp_ampli(i) -= std::abs(it->heightc(ppos,0));
      //  }
      for (uint j = 0; j < sources.size(); ++j) {
   	T(i, 2*j) = real(sources[j]->heightc(pattern_pts[i].getPos()));
   	T(i, 2*j+1) = imag(sources[j]->heightc(pattern_pts[i].getPos()));
	// VEC2C grad = sources[j]->gradHeightc(pattern_pts[i].getPos());
	// T_dx(i, 2*j) = real(grad(0));
	// T_dx(i, 2*j+1) = imag(grad(0));
	// T_dy(i, 2*j) = real(grad(1));
	// T_dy(i, 2*j+1) = imag(grad(1));
      }
   }
    wp_data opt_data;
    opt_data.pattern_pts = pattern_pts;
    opt_data.sources = sources;
    opt_data.transfer_mat = T;
    // opt_data.transfer_mat_dx = T_dx;
    // opt_data.transfer_mat_dy = T_dy;
    opt_data.in = in;
     // INFO("pp_amplis \n"<<pp_ampli);
    
     // INFO("trans mat \n"<<T);
    
    //     double E = pattern_energy_hess(amplis, &grad_out, &hess_out, &opt_data);
    // //INFO("energy "<<E);

      optim::algo_settings_t settings;
      settings.print_level = 0;
      settings.conv_failure_switch = 2;
      
     bool success;
     if (method == 0) {
       success = optim::newton(amplis,pattern_energy_hess, &opt_data, settings);
    } else if (method == 1){
      success = optim::bfgs(amplis,pattern_energy, &opt_data, settings);
    } else if (method == 2) {
      success = optim::nm(amplis,pattern_energy, &opt_data, settings);
    } else {
      success = optim::de_prmm(amplis,pattern_energy, &opt_data, settings);
    }
    
    if(success) {
      // INFO("pp_amplis \n"<<pp_ampli);
         // INFO("trans mat \n"<<T);
      //      INFO("                    SUCCESS amplis \n"<<amplis);
      INFO("                    SUCCESS amplis");
    } else {
      INFO("                    NOPE");
    }
  
    for (uint i = 0; i < sources.size(); ++i) {
      sources[i]->setAmplitude(COMPLEX(amplis(2*i), amplis(2*i+1)));
    }
#endif

#ifdef OPTIM_RE_IM
    VectorXcf c(sources.size());
    VectorXcf p_in(pattern_pts.size());

    //    bool to_solve = false;
    for (uint i = 0; i < pattern_pts.size(); ++i) {
      p_in(i) = pattern_pts[i].getAmpli();
      VEC2 ppos = pattern_pts[i].getPos();
      for (auto it: in) {
	p_in(i) -= it->heightc(ppos,0);
      }
      //VEC2 pb = pattern_pts[i]->getPos();
    }
    // for (uint i = 0; i < boundaries_l[index].size() && !to_solve; ++i) {
    //   to_solve = to_solve || (norm(p_in(i)) > 1e-8*wave_lenghts[index]);
    // }
    //if (to_solve) {
    // INFO(p_in);
     c = svd->solve(p_in);
     // INFO(c);
     //  MatrixXcf T = MatrixXcf(pattern_pts.size(), sources.size());
     // 	for (uint i = 0; i < pattern_pts.size(); ++i) {
     // 	  for (uint j = 0; j < sources.size(); ++j) {
     // 	T(i, j) = sources[j]->heightc(pattern_pts[i].getPos());
     // 	}
     //   }
     // VectorXcf err(sources.size());
     //   err = T*c - p_in;
     // INFO("err\n"<<T);
    
      for (uint i = 0; i < sources.size(); ++i) {
	if (norm(c(i)) < 1e-8) {
	  sources[i]->setAmplitude(0);
	} else {
	  //	  c(i) = (1.0f - damping_source_)*sources->getAmpli() + damping_source_*c(i);
	  IS_DEF(real(c(i)));
	  IS_DEF(imag(c(i)));
	  sources[i]->setAmplitude(c(i));
	}
      }
      //}
  #endif
}

const std::vector<PatternPoint> & WavePattern::getPatternPoints() const {
  return pattern_pts;
}

const std::vector<EquivalentSource*> & WavePattern::getSources() const {
  return sources;
}

void WavePattern::getPatternPointsPos(std::list<VEC2> & positions) const {
  for (auto it: pattern_pts) {
    positions.push_back(it.getPos());
  }
  
}
void WavePattern::getSourcesPos(std::list<VEC2> & positions) const {
  for (auto it: sources) {
    positions.push_back(it->getPos());
  }
}

void WavePattern::energy(std::list<Wave*> in) {
  std::ofstream file;
  file.open("energy.dat");

  int np = pattern_pts.size();
  INFO("pattern size "<<np);
  VECX c(1);
  for (uint i = 0; i < pattern_pts.size(); ++i) {
    c(i) = abs(pattern_pts[i].getAmpli());
  }

  for (int h = 0; h < 300; ++h) {
    FLOAT sx = h*0.1;
    for (int k = 0; k < 300; ++k) {
      FLOAT sy = k*0.1;
      double E = 0;
      COMPLEX phi;
      double aphi_re, aphi_im;
      for (int j = 0; j < np; ++j) {
	phi = 1;
    	//	double win_re = 0, win_im = 0;
    	VEC2 ppos = pattern_pts[j].getPos();
    	// for (auto it: in) {
    	//   win_re += real(it->heightc(ppos,0));
    	//   win_im += imag(it->heightc(ppos,0));
    	// }
    	double ai = real(sources[0]->getAmpli());
    	double bi = imag(sources[0]->getAmpli());
    	FLOAT rx = ppos(0) - sx;
    	FLOAT ry = ppos(1) - sy;
    	FLOAT r  = sqrt(pow(rx, 2.0) + pow(ry, 2.0));
	
	if (r != 0) {
	  phi =  COMPLEX(0, -1)/(FLOAT)4.0*Hankel(wave_number*r);
	}
    	double mui = real(phi);
    	double nui = imag(phi);
    	aphi_re = ai*mui - bi*nui;
    	aphi_im = ai*nui + bi*mui;
    	double epsj = 100*(aphi_re*aphi_re + aphi_im*aphi_im) - 100*c(j)*c(j);
	 // if (sx == 15 && sy > 12 && sy <18) {
	 //   INFO("sx:"<<sy<<" sy:"<<sx<<" epsj:"<<epsj<<" cj"<<c(j)*c(j));
	 //   INFO("aphi_re:"<<aphi_re<<" aphi_im:"<<aphi_im<<" ampli:"<<aphi_re*aphi_re + aphi_im*aphi_im);
	 // }
    	E += epsj*epsj;
      }
      if (E > 1) {
	E = 0;
      }
      file << h<<" "<<k<<" "<<log(E)<<"\n";
      //  if (sx == 15 && sy > 12 && sy <18) {
      // 	std::cout<< h<<" "<<k<<" "<<E<<"\n";
      // }
    }
       file <<"\n";
  }
  
}
