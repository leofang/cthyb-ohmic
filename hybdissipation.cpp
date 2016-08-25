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

#include "hyblocal.hpp"
#include "hybdissipation.hpp"
#include <algorithm>
#include <iomanip>
#include <boost/math/special_functions/gamma.hpp> 
#include <boost/math/constants/constants.hpp>     


dissipation_configuration::dissipation_configuration(const alps::params &p)
{
  beta_=p["BETA"];
  n_orbitals_=p["N_ORBITALS"];
  n_env_=p["N_ENV"]|1;   //Leo Fang: number of reservoirs

  //Leo: read the dissipaton parameters 
  dissipation_=p["Dissipation"]|false;
  //Leo: checks for adding dissipation. These checks cannot be added to sanity_check 
  //     because dissipation_configuration initializes first...
  if(dissipation_ && !p.defined("r"))
        throw std::invalid_argument("If the ohmic environment is needed, please give a nonzero value for \"r\", otherwise set Dissipation to 0.");
  if(dissipation_ && !p.defined("C0"))
        throw std::invalid_argument("If the ohmic environment is needed, please give a value for the capacitance ratio \"C0\" in [0,1], otherwise set Dissipation to 0.");

  r_ = (dissipation_?(double)(p["r"]):0.);
  C0_= (dissipation_?(double)(p["C0"]):0.);
  //wc_= 100.0/beta_; //Leo: set cutoff to be a hundred times of the temperature
  wc_ = 100.; //Leo: cutoff should be temperature independent; TODO: check it's the largest scale in the system!!!
  kappa_ = 1./(beta_*wc_);
  gamma_ = boost::math::double_constants::euler; //Euler constant
  
  dissipation_coeff_.resize(n_orbitals_, std::vector<double> (2*n_env_, 0.));
  for(int i=0; i<n_orbitals_; i++)
  {
      if(n_env_==1)
      {
         dissipation_coeff_[i][0] = C0_;  //enter from 0-th lead  
         dissipation_coeff_[i][1] = -C0_; //leave to 0-th lead
      }
      else if(n_env_==2)
      {
         dissipation_coeff_[i][0] = C0_;       //enter from 0-th lead   
         dissipation_coeff_[i][1] = -(1.-C0_); //enter from 1-th lead  
         dissipation_coeff_[i][2] = -C0_;      //leave to 0-th lead
         dissipation_coeff_[i][3] = (1.-C0_);  //leave to 1-th lead
      }
      else
         throw std::runtime_error("dissipation_configuration: n_env_ is wrong!");
  }

  extern int global_mpi_rank;
  if(global_mpi_rank==0 && dissipation_)
  {
     std::cout << "dissipation is turned on:\nr = " << r_ << ", C0 = " << C0_;
     if(n_env_>1)
        std::cout << " (so C1 = " << 1.0-C0_ << ")";
     std::cout << std::endl << std::endl; 
  }
}



double dissipation_configuration::dissipation_weight_change(const segment &seg, int orbital, bool insert, 
                                                            const local_configuration &local_config) const
{
   if(!dissipation_) return 1.; //Leo: dissipation is turned off, so no need to compute
  
   double tau;
   double J=0;
 
   //combine special cases: 
   //1. insert a segment or an antisegment into a 0-th order orbital
   //2. remove the only segment (which is also an antisegment) from a 1st order orbital 
   if( (local_config.order(orbital)==0 && insert) || (local_config.order(orbital)==1 && !insert) ) 
   {
       tau = std::abs(seg.t_end_ - seg.t_start_);
       J-=dissipation_coeff_[orbital][seg.c_start_]*dissipation_coeff_[orbital][seg.c_start_+n_env_]*phase_correlator_J(tau);
       return std::exp(J);
   }

   //general case: this is a O(4k+1) operation for insertion into a k-th order orbital 
   //              the design takes the advantage of the fact that phase_correlator_J(0)=0. 
   for(std::set<segment>::const_iterator it=local_config.segments_[orbital].begin();\
       it != local_config.segments_[orbital].end(); ++it)
   {//pair the start and end times of the (anti)segment with other segments
       
       tau = std::abs(seg.t_start_ - it->t_start_);
       J-=dissipation_coeff_[orbital][seg.c_start_]*dissipation_coeff_[orbital][it->c_start_]*phase_correlator_J(tau);
   
       tau = std::abs(seg.t_start_ - it->t_end_);
       J-=dissipation_coeff_[orbital][seg.c_start_]*dissipation_coeff_[orbital][it->c_start_+n_env_]*phase_correlator_J(tau);
  
       tau = std::abs(seg.t_end_ - it->t_start_);
       J-=dissipation_coeff_[orbital][seg.c_start_+n_env_]*dissipation_coeff_[orbital][it->c_start_]*phase_correlator_J(tau); 
  
       tau = std::abs(seg.t_end_ - it->t_end_);
       J-=dissipation_coeff_[orbital][seg.c_start_+n_env_]*dissipation_coeff_[orbital][it->c_start_+n_env_]*phase_correlator_J(tau); 
   }
   //the contribution of the to-be-changed (anti)segment itself
   tau = std::abs(seg.t_start_ - seg.t_end_);
   if(insert)
     J-=dissipation_coeff_[orbital][seg.c_start_]*dissipation_coeff_[orbital][seg.c_start_+n_env_]*phase_correlator_J(tau);
   else //the contribution is subtracted twice, so add it back once
     J+=dissipation_coeff_[orbital][seg.c_start_]*dissipation_coeff_[orbital][seg.c_start_+n_env_]*phase_correlator_J(tau);

//Leo: bad design when the distance happens to be the same as len
//   if(insert) //insert a segment or antisegment 
//   {
//     for(std::set<segment>::const_iterator it=local_config.segments_[orbital].begin();\
//         it != local_config.segments_[orbital].end(); ++it)
//     { //pair the start and end times of the (anti)segment with other segments
//    
//         tau = std::abs(seg.t_start_ - it->t_start_);
//         J-=dissipation_coeff_[orbital][seg.c_start_]*dissipation_coeff_[orbital][it->c_start_]*phase_correlator_J(tau);
//   
//         tau = std::abs(seg.t_start_ - it->t_end_);
//         J-=dissipation_coeff_[orbital][seg.c_start_]*dissipation_coeff_[orbital][it->c_start_+n_env_]*phase_correlator_J(tau);
//  
//         tau = std::abs(seg.t_end_ - it->t_start_);
//         J-=dissipation_coeff_[orbital][seg.c_start_+n_env_]*dissipation_coeff_[orbital][it->c_start_]*phase_correlator_J(tau); 
//  
//         tau = std::abs(seg.t_end_ - it->t_end_);
//         J-=dissipation_coeff_[orbital][seg.c_start_+n_env_]*dissipation_coeff_[orbital][it->c_start_+n_env_]*phase_correlator_J(tau); 
//     }
//     //the contribution of the to-be-moved segment itself
//     tau = std::abs(seg.t_start_ - seg.t_end_);
//     J-=dissipation_coeff_[orbital][seg.c_start_]*dissipation_coeff_[orbital][seg.c_start_+n_env_]*phase_correlator_J(tau);
//   }
//   else //remove a segment or antisegment
//   {
//     double len = std::abs(seg.t_end_-seg.t_start_);
//     for(std::set<segment>::const_iterator it=local_config.segments_[orbital].begin(); 
//         it != local_config.segments_[orbital].end(); ++it)
//     {//Leo: this part is a bit complicated due to the need of avoiding counting the (anti)segment to be removed 
//         int i=0;   
// 
//         tau = std::abs(seg.t_start_ - it->t_start_);
//         if(i<4 && (tau==0||tau==len))  
//            i+=1; //do nothing except adding the counter 
//         else
//            J-=dissipation_coeff_[orbital][seg.c_start_]*dissipation_coeff_[orbital][it->c_start_]*phase_correlator_J(tau);
//   
//         tau = std::abs(seg.t_start_ - it->t_end_);
//         if(i<4 && (tau==0||tau==len))  
//            i+=1; //do nothing except adding the counter 
//         else
//            J-=dissipation_coeff_[orbital][seg.c_start_]*dissipation_coeff_[orbital][it->c_start_+n_env_]*phase_correlator_J(tau);
//  
//         tau = std::abs(seg.t_end_ - it->t_start_);
//         if(i<4 && (tau==0||tau==len))  
//            i+=1; //do nothing except adding the counter 
//         else
//            J-=dissipation_coeff_[orbital][seg.c_start_+n_env_]*dissipation_coeff_[orbital][it->c_start_]*phase_correlator_J(tau); 
//  
//         tau = std::abs(seg.t_end_ - it->t_end_);
//         if(i<4 && (tau==0||tau==len))  
//            i+=1; //do nothing except adding the counter 
//         else
//            J-=dissipation_coeff_[orbital][seg.c_start_+n_env_]*dissipation_coeff_[orbital][it->c_start_+n_env_]*phase_correlator_J(tau); 
//     }
//     //the contribution of the to-be-moved segment itself
//     J-=dissipation_coeff_[orbital][seg.c_start_]*dissipation_coeff_[orbital][seg.c_start_+n_env_]*phase_correlator_J(len);
//   }
   return std::exp(J);
}


//Leo: J(|tau|) is defined such that <T exp(i\phi(tau)) exp(-i\phi(0))> = exp(J(|tau|))
//     Note that the time ordering guarantees that the argument of J is always positive, 
//     so no need to do time wrapping   
double dissipation_configuration::phase_correlator_J(double tau) const
{
//   if(tau<0) tau=std::abs(tau); //TODO: remove this line!
   if(tau<0) 
     throw std::runtime_error("The argument of the phase correlator is negative!");
   else if (tau==0. || tau==beta_)
     return 0.; //skip the function call as J(0)=0
   else
     return -2.0*r_*(std::log(1.0/kappa_*pow(tgamma(1.0+kappa_),2)/(tgamma(1.0+kappa_-tau/beta_)*tgamma(kappa_+tau/beta_))));
   //return -2.0*r_*(std::log(1.0/kappa_*pow(tgamma(1.0+kappa_),2)/(tgamma(1.0+kappa_-tau/beta_)*tgamma(kappa_+tau/beta_)))+gamma_);
}
