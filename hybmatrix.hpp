/****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2012 by Emanuel Gull <gull@pks.mpg.de>,
 *                       Hartmut Hafermann <hafermann@cpht.polytechnique.fr>
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
#ifndef HYB_MATRIX
#define HYB_MATRIX
#include<map>
#include<vector>
#include "hybsegment.hpp"
#include "hybfun.hpp"
#include "hybblasmatrix.hpp"

//This is the hybridization matrix class responsible for computing determinant ratios and the like. Derived from a general blas matrix class that can handle elementary blas and lapack operations.

//simple matrix that uses BLAS calls for rank one and matrix vector.
typedef std::map<double,std::size_t> hyb_map_t;
// Leo Fang: a vectorized vesion of the above map for colored operators
//typedef std::vector< std::map<double,std::size_t> > hyb_vector_map_t;

class hybmatrix:public blas_matrix{
public:
  hybmatrix(const alps::params &p){ determinant_=1.; determinant_old_=1.; permutation_sign_=1.; beta_=p["BETA"]; measure_g2w_=p["MEASURE_g2w"]|0; measure_h2w_=p["MEASURE_h2w"]|0; 
  n_env_=p["N_ENV"]|1; //Leo: introduce multiple reservoirs
  time_ordering_sign_=1;      //Leo: 1 means no, -1 means yes
}

  //Leo: copy constructor
  hybmatrix(const hybmatrix &rhs) 
      :blas_matrix(rhs)
      ,cdagger_index_map_(rhs.cdagger_index_map_)
      ,c_index_map_(rhs.c_index_map_)
      ,c_cdagger_map_(rhs.c_cdagger_map_) //Leo: map of operator positions
      ,n_env_(rhs.n_env_)
//      ,Q(rhs.Q)
//      ,R(rhs.R)
//      ,PinvQ(rhs.PinvQ)
      ,S(rhs.S)
      ,S_tilde_inv(rhs.S_tilde_inv)
      ,S_tilde(rhs.S_tilde)
      ,weight_ratio_(rhs.weight_ratio_)
      ,permutation_sign_(rhs.permutation_sign_)
      ,time_ordering_sign_(rhs.time_ordering_sign_)
//      ,disordered_times(rhs.disordered_times)
      ,determinant_(rhs.determinant_)
      ,determinant_old_(rhs.determinant_old_)
      ,beta_(rhs.beta_)
      ,measure_g2w_(rhs.measure_g2w_)
      ,measure_h2w_(rhs.measure_h2w_)
  {}

  //Leo: copy assignment operator
  hybmatrix & operator=(hybmatrix rhs)
  {
     std::swap(*this, rhs);
     return *this;
  }

  //Leo: destructor
  ~hybmatrix() {
//    std::cerr << "Deleting hybmatrix\n";
  }

  double hyb_weight_change_insert(const segment &new_segment, int orbital, const hybfun &Delta);
  double hyb_weight_change_remove(const segment &new_segment, int orbital, const hybfun &Delta);
  void insert_segment(const segment &new_segment, int orbital);
  void remove_segment(const segment &new_segment, int orbital);
  friend std::ostream &operator<<(std::ostream &os, const hybmatrix &hyb_mat);
  void rebuild_hyb_matrix(int orbital, const hybfun &Delta);
  void rebuild_ordered_hyb_matrix(int orbital, const hybfun &Delta);
  double full_weight() const;
  void measure_G(std::vector<double> &G, std::vector<double> &F, const std::map<double,double> &F_prefactor, double sign, int total_size, double dissipation_weight_ratio) const;
  void measure_Gw(std::vector<double> &Gwr, std::vector<double> &Gwi,std::vector<double> &Fwr, std::vector<double> &Fwi, const std::map<double,double> &F_prefactor, double sign) const;
  void measure_Gw_buffer(std::vector<double> &Gwr, std::vector<double> &Gwi,std::vector<double> &Fwr, std::vector<double> &Fwi, const std::map<double,double> &F_prefactor, double sign) const;
  void measure_G2w(std::vector<std::complex<double> > &G2w, std::vector<std::complex<double> >&F2w, int N_w2, int N_w_aux, const std::map<double,double> &F_prefactor) const;
  void measure_Gl(std::vector<double> &Gl, std::vector<double> &Fl, const std::map<double,double> &F_prefactor, double sign) const;
  //void measure_conductance(std::vector<double> &giwn, double sign, int orbital, const hybfun &Delta) const;
  void consistency_check() const;
  //Leo: this seems to be necessary when n_env>1
  //void time_ordering_sign_check(int &time_ordering_sign_, std::set<double> &disordered_times); 
  void time_ordering_sign_check(); 
  int sign() const {return time_ordering_sign_;}    //Leo: access the sign
  inline int head() const {hyb_map_t::const_iterator it=c_cdagger_map_.begin(); return it->second;}  

  inline void access_cdagger_times(std::vector<double> &cdagger_times) const
  {for (hyb_map_t::const_iterator it= cdagger_index_map_.begin(); it != cdagger_index_map_.end(); ++it) cdagger_times.push_back(it->first);}
  inline void access_c_times(std::vector<double> &c_times) const
  {for (hyb_map_t::const_iterator it= c_index_map_.begin(); it != c_index_map_.end(); ++it) c_times.push_back(it->first);}

  //Leo: for swapping two hybmatrix objects
  friend void swap(hybmatrix &A, hybmatrix &B);

private:

  //map of start/end times and their corresponding rows and columns in the matrix.
  hyb_map_t cdagger_index_map_;
  hyb_map_t c_index_map_;

  //Leo: map of start/end times and their corresponding positions in the segment (start or end).
  hyb_map_t c_cdagger_map_;

  //Leo: a vectorized vesion of the above map for colored operators
  //hyb_vector_map_t cdagger_index_vector_map_;
  //hyb_vector_map_t c_index_vector_map_;
 
  // Leo Fang: number of reservoirs
  int n_env_;
   
  //auxiliary column and row vectors for inserting and removing elements.
  std::vector<double> Q;
  std::vector<double> R;
  std::vector<double> PinvQ;
  double S;
  double S_tilde_inv;
  double S_tilde;
  double weight_ratio_;
  double permutation_sign_;
  int time_ordering_sign_;  //Leo: this seems to be necessary when n_env>1; 1 means no, -1 means yes
  std::set<double> disordered_times;
  
  //debug
  double determinant_;
  double determinant_old_;
  
  //physics
  double beta_;
  bool measure_g2w_, measure_h2w_;
};

std::ostream &operator<<(std::ostream   &os, const hybmatrix &hyb_mat); 

#endif
