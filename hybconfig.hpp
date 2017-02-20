  /****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2012 by Emanuel Gull <gull@pks.mpg.de>,
 *                   
 *  based on an earlier version by Philipp Werner and Emanuel Gull
 *
 *
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/
#ifndef HYB_CONFIG_HPP
#define HYB_CONFIG_HPP

#include<alps/ngs.hpp>
#include"hybfun.hpp"
#include"hybsegment.hpp"
#include"hybmatrix.hpp"
//this class knows everything to do with hybridization / bath operators.
//determinants of hybridization matrices, rank one updates, and so on.
//this is where \det \Delta is implemented.
class hybridization_configuration{
public:
  hybridization_configuration(const alps::params &p);

  hybridization_configuration(const hybridization_configuration &rhs) : 
     n_env_(rhs.n_env_), beta_(rhs.beta_), Delta(rhs.Delta), hybmat_(rhs.hybmat_), hyb_strength(rhs.hyb_strength)
  {
    //std::cerr << hybmat_.size() << std::endl;
    //Leo: double loops to accommodate multiple reservoirs 
    for (int i=0; i<hybmat_.size(); i++) //N_Orbital
    {
	for(int j=0; j<hybmat_[i].size(); j++) //N_ENV
  	      hybmat_[i][j].rebuild_hyb_matrix(i, Delta[j]);
    }
  }

  //Leo: new copy-assignment operator using the copy-elision-and-swap idiom
  hybridization_configuration &operator=(hybridization_configuration rhs)        
  {                                     
     std::swap(*this, rhs);
     return *this;         
  }

  //Leo: deprecated copy-assignment operator
//  const hybridization_configuration &operator=(const hybridization_configuration &rhs) 
//  {
//    hybmat_ = rhs.hybmat_;
//    Delta = rhs.Delta;
//    //Leo: double loops to accommodate multiple reservoirs 
//    for (int i=0; i<hybmat_.size(); i++) //N_Orbital
//    {
//	for(int j=0; j<hybmat_[i].size(); j++) //N_ENV
//  	      hybmat_[i][j].rebuild_hyb_matrix(i, Delta[j]);
//    }
//    return *this;
//  }
  
  ~hybridization_configuration() {
//    std::cerr << "Deleting hybconfig\n";
  }
  
  //Leo: although the color of a (anti-)segment is recorded in the segment class,
  //it seems a bit convenient to just pass it as an additional argument
  double hyb_weight_change_insert(const segment &new_segment, int orbital, size_t color);
  void insert_segment(const segment &new_segment, int orbital, size_t color);
  void insert_antisegment(const segment &new_antisegment, int orbital, size_t color);
  double hyb_weight_change_remove(const segment &new_segment, int orbital, size_t color);
  void remove_segment(const segment &new_segment, int orbital, size_t color);
  void remove_antisegment(const segment &new_antisegment, int orbital, size_t color);

  //for worm update
  double hyb_weight_change_worm_creep(const segment &new_segment, double old_worm_head, double worm_tail, int orbital);
  void worm_creep(const segment &new_segment, double old_worm_head, double worm_tail, int orbital);

  //Leo: since the size of reservoir matrices coupled to the same orbital are correlated,
  //     there should be a consistency check after each update to gaurantee we're doing it correctly
  int total_color_matrix_size(int orbital) const;
  //Leo: check overall sign of color matrices
  int overall_color_matrix_sign() const;
  //Leo: initialize the Delta vector 
  void initialize_Delta(const alps::params &p);
  //Leo: TOTALLY EXPERIMENTAL!!!
  void haunt_missing_sign(int orbital);
  double hyb_weight_change_flip(int orbital, size_t color_1, size_t color_2);
  int flip_color(int orbital, size_t color_1, size_t color_2);
    
  void dump();
  void rebuild();
  void rebuild_ordered();
  void rebuild(int orbital);
  void rebuild(std::vector<int> orbital);

  void measure_G(std::vector<std::vector<double> > &G, std::vector<std::vector<double> > &F, const std::vector<std::map<double,double> > &F_prefactor, double sign, double dissipation_weight_ratio) const;
  void measure_Gw(std::vector<std::vector<double> > &Gwr,std::vector<std::vector<double> > &Gwi,std::vector<std::vector<double> > &Fwr,std::vector<std::vector<double> > &Fwi, const std::vector<std::map<double,double> > &F_prefactor, double sign) const;
  void measure_G2w(std::vector<std::vector<std::complex<double> > > &G2w, std::vector<std::vector<std::complex<double> > > &F2w, int N_w2, int N_w_aux, const std::vector<std::map<double,double> > &F_prefactor) const;
  void measure_Gl(std::vector<std::vector<double> > &Gl,std::vector<std::vector<double> > &Fl, const std::vector<std::map<double,double> > &F_prefactor, double sign) const;
  void measure_conductance(std::vector<std::vector<double> > &giwn, double sign) const; //Leo: conductance measurement
  double full_weight() const;

  friend std::ostream &operator<<(std::ostream &os, const hybridization_configuration &hyb_config);

private:
  int n_env_; //Leo: number of reservoirs
  double beta_; //Leo: for conductance measurement
  //the hybridization function Delta
  //Leo: to accomodate multiple reservoirs, this object should be vectorized with dimension N_ENV
  std::vector<hybfun> Delta;
  //Leo: this should be a vector of vector for the same reason
  //The inner vector has dimension N_ENV (new), and the outer vector has dimension N_Orbital (from original code)
  std::vector< std::vector<hybmatrix> > hybmat_;
  std::vector< std::vector<double> > hyb_strength; //Leo: for color flip update
};

std::ostream &operator<<(std::ostream &os, const hybridization_configuration &hyb_config);

#endif

