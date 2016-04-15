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
#include <algorithm>
#include <iomanip>
// #include <boost/math/special_functions/gamma.hpp> //Leo: for dissipation
// #include <boost/math/constants/constants.hpp>     //Leo: for dissipation

//Leo: disable retarded interaction
local_configuration::local_configuration(const alps::params &p, int crank): crank_(crank), U_(p), mu_(p)//, K_(p)
{
  beta_=p["BETA"];
  n_orbitals_=p["N_ORBITALS"];
//    std::cerr << "Start ...";
  n_env_=p["N_ENV"]|1;   //Leo Fang: number of reservoirs
  segments_.resize(n_orbitals_);
  zero_order_orbital_occupied_.resize(n_orbitals_,false);

  n_segments_.resize(n_orbitals_, std::vector<int> (2*n_env_, 0));

  ////Leo: read the dissipaton parameters (sanity_checks has done necessary checks when this class is initialized!)
  //dissipation_=p["Dissipation"]|0;
  ////Leo: checks for adding dissipation. 
  ////     These checks cannot be added to sanity_check because local_configuration initializes first...
  //if(dissipation_ && !p.defined("r"))
  //      throw std::invalid_argument("If the ohmic environment is needed, please give a nonzero value for \"r\", otherwise set Dissipation to 0.");
  //if(dissipation_ && !p.defined("C0"))
  //      throw std::invalid_argument("If the ohmic environment is needed, please give a value for the capacitance ratio \"C0\" in [0,1], otherwise set Dissipation to 0.");

  //r_ = (dissipation_?(double)(p["r"]):0.);
  //C0_= (dissipation_?(double)(p["C0"]):0.);
  //wc_= 100.0/beta_; //Leo: set cutoff to be a hundred times of the temperature
  //kappa_ = 1./(beta_*wc_);
  //gamma_ = boost::math::double_constants::euler; //Euler constant
  //
  //dissipation_coeff_.resize(n_orbitals_, std::vector<double> (2*n_env_, 0.));
  //for(int i=0; i<n_orbitals_; i++)
  //{
  //    if(n_env_==1)
  //    {
  //       dissipation_coeff_[i][0] = C0_;  //enter from 0-th lead  
  //       dissipation_coeff_[i][1] = -C0_; //leave to 0-th lead
  //    }
  //    else if(n_env_==2)
  //    {
  //       dissipation_coeff_[i][0] = C0_;       //enter from 0-th lead   
  //       dissipation_coeff_[i][1] = -(1.-C0_); //enter from 1-th lead  
  //       dissipation_coeff_[i][2] = -C0_;      //leave to 0-th lead
  //       dissipation_coeff_[i][3] = (1.-C0_);  //leave to 1-th lead
  //    }
  //    else
  //       throw std::runtime_error("local_configuration: n_env_ is wrong!");
  //}

//    std::cerr << " done\n";

  //Leo: disable retarded interaction
//  use_retarded_interaction_=p.defined("RET_INT_K");
//  if(use_retarded_interaction_)
//  {
//    double Kp0=K_.interpolate_deriv(0.0);//K'(0^+)
//    U_.apply_shift(-2.*Kp0); //apply static shift caused by the retarded interaction
//    mu_.apply_shift(+Kp0);
//  }
  
  extern int global_mpi_rank;
  if(global_mpi_rank==0)
  {
    std::cout<<U_<<mu_<<std::endl;
    std::cout<<*this<<std::endl;
  }
}


/* Leo Fang: for test purpose, print out the segment map */
//void local_configuration::print_time_map() const
//{
//       for(std::set<double>::const_iterator it=times_set_.begin(); it!=times_set_.end(); it++)
//       {
//           std::cout << *it << ", ";
//       }
//       std::cout << std::endl << std::endl;
//}


//Leo: for test purpose, print out the segment times and colors
//This function seems redundant, so it is replaced by the overloaded operator<< defined below...
//void local_configuration::print_segments() const
//{
//    for (int i=0; i<n_orbitals_; i++)
//    {
//        std::cout << "Orbital " << i << ":";
//        if(order(i)==0)
//        {
//            std::cout << (zero_order_orbital_occupied(i)?" fully occupied":" empty") << std::endl << std::endl; 
//        }
//        else
//        {
//            std::cout << std::endl;
//            for(std::set<segment>::const_iterator it=segments_[i].begin(); it != segments_[i].end();++it)
//    	    {
//                 std::cout << *it << std::endl;
//    	    }
//    	std::cout << std::endl;
//        }
//    }
//}


//Leo: the overloaded operator is slightly modified to accomodate the colors
std::ostream &operator<<(std::ostream &os, const local_configuration &local_conf)
{
  os<<"local configuration: "<<std::endl;
  for(int i=0; i<local_conf.n_orbitals_; ++i)
  {
    os << "Orbital "<< i << ":" << std::endl; 
    if(local_conf.segments_[i].size()==0)
    {
      os << (local_conf.zero_order_orbital_occupied_[i]?"occupied":"empty") << std::endl;
      if(i!=local_conf.n_orbitals_-1) 
        os << std::endl;
    }
    else
    {
      //os<<i<<" ";
      for(std::set<segment>::const_iterator it=local_conf.segments_[i].begin(); it!=local_conf.segments_[i].end(); ++it)
      {
        //Leo: use the overloaded operator for segments
        //os<<"("<<it->t_start_<<" "<<it->t_end_<<") ";
        os << *it << std::endl;
      }
      if(i!=local_conf.n_orbitals_-1)  os << std::endl;
    }
  }
  //if(local_conf.dissipation_)
  //{
  //   std::cout << "dissipation is turned on:\nr = " << local_conf.r_ << ", C0 = " << local_conf.C0_ << std::endl;
  //   if(local_conf.n_env_>1)
  //      std::cout << "(so C1 = " << 1.0-local_conf.C0_ << ")" << std::endl; 
  //}
  return os;
}


//compute the segment overlap and segment length, and return the weight
//this is an O(k n_orbital) procedure in the expansion order and the number or orbitals.
double local_configuration::local_weight_change(const segment &seg, int orb, bool antisegment) const
{
  //the chemical potential term: just needs a segment length
  double length=seg.t_end_-seg.t_start_;
  double sgn=antisegment?-1.:1.;
  if(length<0) length+=beta_; //wraparound segment
  double weight=std::exp(sgn*mu_[orb]*length);
  //std::cout<<clmagenta<<"chemical potential weight is: "<<weight<<" mu: "<<mu_[orb]<<" length: "<<length<<cblack<<std::endl;
  
  //the interaction term needs the overlap between this orbital and all the other orbitals
  //Leo: why make it static? The lifetime of this overlaps vector should not be infinite!
  static std::vector<double> overlaps(n_orbitals_, 0.); 
  for(int i=0;i<n_orbitals_;++i) { overlaps[i]=0.; }
  for(int i=0;i<n_orbitals_;++i)
  {
    if(i==orb) continue;
    
    if(zero_order_orbital_occupied_[i])
    {
      overlaps[i]=length;
    }
    else
    {
      //find the first segment with time t_start > seg.t_start
      for(std::set<segment>::const_iterator it=segments_[i].begin(); it != segments_[i].end(); ++it)
      {
        overlaps[i]+=segment_overlap(seg, *it);
      }
    }
    //std::cout<<clmagenta<<"weight of orbital: "<<orb<<" wrt orbital: "<<i<<" is: "<<std::exp(U_(orb,i)*overlaps[i])<<" for overlap: "<<overlaps[i]<<cblack<<std::endl;
    weight*=std::exp(-sgn*U_(orb,i)*overlaps[i]);
    /*if(zero_order_orbital_occupied_[i]){
     std::cout<<"weight got an additional factor:"<<std::exp(-sgn*U_(orb,i)*overlaps[i])<<std::endl;
     }*/
  }

  //Leo: disable retarded interaction
  //this is the retarded interaction stuff
  //if(use_retarded_interaction_){
  //  bool is_removal=false;
  //  double retarded_weight=0;
  //  if((seg.t_start_==0) && (seg.t_end_==beta_)){ //not really a segment but a full line
  //    retarded_weight=0.;
  //  }else{
  //    for(int i=0;i<n_orbitals_;++i){
  //      for(std::set<segment>::const_iterator it=segments_[i].begin(); it != segments_[i].end();++it){
  //        retarded_weight+=sgn*K_.interpolate(seg.t_start_-it->t_start_);
  //        retarded_weight-=sgn*K_.interpolate(seg.t_start_-it->t_end_);
  //        retarded_weight-=sgn*K_.interpolate(seg.t_end_-it->t_start_);
  //        retarded_weight+=sgn*K_.interpolate(seg.t_end_-it->t_end_);
  //        //subtract doubly counted contribution if it is equal orbital:
  //        if((it->t_start_==seg.t_start_) || (it->t_start_==seg.t_end_)){ //segment remove
  //          //            std::cout<<"we have an antisegment remove!"<<std::endl;
  //          is_removal=true;
  //        }
  //      }
  //    }
  //    //      std::cout<<"taking care of antisegment remove."<<std::endl;
  //    if(is_removal){
  //      retarded_weight+=K_.interpolate(seg.t_end_-seg.t_start_);
  //    }else{
  //      retarded_weight-=K_.interpolate(seg.t_end_-seg.t_start_);
  //    }
  //  }
  //  //std::cout<<clblue<<"is_removal is: "<<is_removal<<" (true: "<<true<<")"<<std::endl;
  //  weight*=std::exp(retarded_weight);
  //}
  return weight;
}


bool local_configuration::has_overlap(const segment &seg,const int orb) 
{
//  bool overlap = false;
  if(zero_order_orbital_occupied_[orb])
  {
    return true;
  } 
  else 
  {
    //find the first segment with time t_start > seg.t_start
    for(std::set<segment>::const_iterator it=segments_[orb].begin(); it!=segments_[orb].end(); ++it) 
    {
      if (segment_overlap(seg, *it)>0.0) return true;
    }
  }
  return false;
}


//this computes the overlap of two segments. It takes care of the wrapping around zero by splitting the segments up.
double local_configuration::segment_overlap(const segment &seg1, const segment &seg2) const
{
  if(seg1.t_start_>seg1.t_end_) 
    return segment_overlap(segment(seg1.t_start_, beta_), seg2)+segment_overlap(segment(0, seg1.t_end_), seg2);
  if(seg2.t_start_>seg2.t_end_) 
    return segment_overlap(seg1, segment(seg2.t_start_, beta_))+segment_overlap(seg1, segment(0, seg2.t_end_));
  double t1=std::max(seg1.t_start_, seg2.t_start_);
  double t2=std::min(seg1.t_end_, seg2.t_end_);
  return t2-t1<0?0.:t2-t1;
}


//find the distance to the next segment start (the next creation operator)
double local_configuration::find_next_segment_start_distance(double time, int orbital)
{
  if(segments_[orbital].size()==0) return beta_; //no segments present
  std::set<segment>::const_iterator it=segments_[orbital].upper_bound(segment(time, 0.));
  if(it==segments_[orbital].end()) //Leo: "time" falls after the last segment
    return (beta_-time+segments_[orbital].begin()->t_start_); //wrap around
  return it->t_start_-time; //Leo: "time" falls in-between segments
}


//find the distance to the next segment end (the next annihilation operator)
//this is a tiny bit more involved, the end time may be part of the previous segment
double local_configuration::find_next_segment_end_distance(double time, int orbital)
{
  double distance;
  if(segments_[orbital].size()==0) return beta_; //no segments present
  if(segments_[orbital].size()==1)
  {
    distance=segments_[orbital].begin()->t_end_-time;
    return distance<0?distance+beta_:distance; //single segment present
  }
  //first possibility: like start time, the closest end-time is after this segment
  std::set<segment>::const_iterator it=segments_[orbital].upper_bound(segment(time, 0.));
  if(it==segments_[orbital].end()) //Leo: "time" falls after the last segment
	distance=(segments_[orbital].begin()->t_end_-time); //wrap around to end time
  else 
	distance=it->t_end_-time; //Leo: "time" falls in-between segments
  if(distance<0) distance+=beta_;
  //second possibility: the closest end time is part of the previous segment
  if(it==segments_[orbital].begin()) it=segments_[orbital].end(); //wrap around
  it--;
  double distance2=it->t_end_-time; if (distance2<0) distance2+=beta_;
  return std::min(distance, distance2);
}


void local_configuration::insert_segment(const segment &new_segment, int orbital)
{
  segments_[orbital].insert(new_segment);
  if(!times_set_.insert(new_segment.t_start_).second)
  {  
     std::stringstream s; 
     s<<crank_; 
     std::cout<<*this<<std::endl; //Leo: print out the local configuration
     throw std::logic_error("rank "+s.str()+": insert segment start time could not be inserted.");
  }
  if(!times_set_.insert(new_segment.t_end_).second)
  {
     std::stringstream s; 
     s<<crank_; 
     std::cout<<*this<<std::endl; 
     //Leo: with the modified operator<<, this line is redundant
     //std::cout<<"inserted segment: "<<new_segment<<"into orbital: "<<orbital<<std::endl; 
     throw std::logic_error("rank "+s.str()+": insert segment end time could not be inserted.");
  }
}


void local_configuration::insert_antisegment(const segment &new_antisegment, int orbital)
{
  //find segment of which this one is a part
  //full line case
  if(segments_[orbital].size()==0)
  {
    segments_[orbital].insert(new_antisegment);
    zero_order_orbital_occupied_[orbital]=false;
  }
  //general case: need to find a segment, then split it in two.
  else
  {
    std::set<segment>::iterator it=segments_[orbital].upper_bound(new_antisegment);
    if(it==segments_[orbital].begin()) it=segments_[orbital].end(); //wrap around
    it--;
    
    //Leo: this two lines are modified because when spliting a segment we need to paint the colors!
    segment new_later_segment(new_antisegment.t_start_, it->t_end_, new_antisegment.c_start_, it->c_end_);
    segment new_earlier_segment(it->t_start_, new_antisegment.t_end_, it->c_start_, new_antisegment.c_end_);

    segments_[orbital].erase(it);
    segments_[orbital].insert(new_later_segment);
    segments_[orbital].insert(new_earlier_segment);
  }
  if(!times_set_.insert(new_antisegment.t_start_).second)
  {
    std::stringstream s; 
    s<<crank_; 
    std::cout<<*this<<std::endl; //Leo: print out the local configuration
    throw std::logic_error("rank "+s.str()+": insert antisegment start time could not be inserted.");
  }
  if(!times_set_.insert(new_antisegment.t_end_).second)
  {
    std::stringstream s; 
    s<<crank_; 
    std::cout<<*this<<std::endl; //Leo: print out the local configuration
    throw std::logic_error("rank "+s.str()+": insert antisegment start time could not be inserted.");
  }
}


void local_configuration::remove_antisegment(const segment &new_antisegment, int orbital)
{
  //find segment of which this one is a part
  //full line case
  if(segments_[orbital].size()==1)
  {
    segments_[orbital].erase(new_antisegment);
    zero_order_orbital_occupied_[orbital]=true;
  }
  //general case: need to find two segments and merge them
  else
  {
    std::set<segment>::iterator it_later=segments_[orbital].find(new_antisegment);
    std::set<segment>::iterator it_earlier=it_later;
    if(it_earlier==segments_[orbital].begin()) it_earlier=segments_[orbital].end(); //wrap around
    it_earlier--;

    //Leo: this two lines are modified because when merging two segments we need to take care of the colors!
    segment new_segment(it_earlier->t_start_, it_later->t_end_, it_earlier->c_start_, it_later->c_end_);
//    segment new_segment=*it_earlier;
//    new_segment.t_end_=it_later->t_end_;

    segments_[orbital].erase(it_later);
    segments_[orbital].erase(it_earlier);
    segments_[orbital].insert(new_segment);
  }
  if(!times_set_.erase(new_antisegment.t_start_))
  {
    std::cerr<<"in local_configuration::remove_antisegment"<<std::endl;
    std::cerr<<"time to erase was: "<<new_antisegment.t_start_<<std::endl;
    std::cerr<<"new antisegment to remove was: "<<new_antisegment<<std::endl;
    std::cout<<*this<<std::endl;
    throw std::logic_error("did not find start time to remove!");
  }  
  if(!times_set_.erase(new_antisegment.t_end_))
  {
    std::cerr<<"successfully erased new antisegment at: "<<new_antisegment<<std::endl;
    std::cerr<<"successfully erased the segment before that. "<<std::endl;
    std::cerr<<"successfully inserted the segment at: (gone.)"<<std::endl;
    std::cerr<<"the times are: "<<std::endl;
    for (std::set<double>::const_iterator it=times_set_.begin(); it!=times_set_.end();++it)
    { 
	std::cout<<*it<<" ";
    } 
    std::cerr<<std::endl;
    std::cerr<<std::endl;
    std::cerr<<"in local_configuration::remove_antisegment"<<std::endl;
    std::cerr<<"time to erase was: "<<new_antisegment.t_end_<<std::endl;
    std::cerr<<"new antisegment to remove was: "<<new_antisegment<<std::endl;
    std::cout<<*this<<std::endl;
    throw std::logic_error("did not find end time to remove!");
  }
}


segment local_configuration::get_segment(int k, int orbital) const
{
  std::set<segment>::const_iterator it= segments_[orbital].begin();
  if(k>=(int)(segments_[orbital].size())) throw std::logic_error("not enough segments to get this one.");
  advance(it, k);
  return *it;
}


void local_configuration::remove_segment(const segment &new_segment, int orbital)
{
  //std::cout<<clmagenta<<*this<<cblack<<std::endl;
  //std::cout<<clmagenta<<"segment to remove: "<<new_segment<<cblack<<std::endl;
  if(segments_[orbital].erase(new_segment)==0) throw std::logic_error("did not find segment to remove!");
  if(segments_[orbital].size()==0) zero_order_orbital_occupied_[orbital]=false;
  
  if(!times_set_.erase(new_segment.t_start_))
  {
    std::cerr<<"in local_configuration::remove_segment"<<std::endl;
    std::cerr<<"time to erase was: "<<new_segment.t_start_<<std::endl;
    std::cerr<<"new segment to remove was: "<<new_segment<<std::endl;
    std::cout<<*this<<std::endl;
    throw std::logic_error("did not find start time to remove!");
  }  
  if(!times_set_.erase(new_segment.t_end_)) 
  {
    std::cerr<<"in local_configuration::remove_segment"<<std::endl;
    std::cerr<<"time to erase was: "<<new_segment.t_end_<<std::endl;
    std::cerr<<"new segment to remove was: "<<new_segment<<std::endl;
    std::cout<<*this<<std::endl;
    throw std::logic_error("did not find end time to remove!");
  }
}


//compute the segment overlap and segment length, and return the weight
//this is an O(k n_orbital) procedure in the expansion order and the number or orbitals.
double local_configuration::local_energy(const segment &seg, int orb,bool d_mu_only) const
{
    //the chemical potential term: just needs a segment length
    double length=seg.t_end_-seg.t_start_;
    if(length<0) length+=beta_; //wraparound segment
    double energy = mu_[orb]*length;
    if (d_mu_only) return energy;
    //std::cout<<clmagenta<<"chemical potential weight is: "<<weight<<" mu: "<<mu_[orb]<<" length: "<<length<<cblack<<std::endl;
    
    //the interaction term needs the overlap between this orbital and all the other orbitals
    static std::vector<double> overlaps(n_orbitals_, 0.);
    for(int i=0;i<n_orbitals_;++i)  overlaps[i]=0.;
    for(int i=0;i<n_orbitals_;++i)
    {
        if(i==orb) continue;
        
        if(zero_order_orbital_occupied_[i])
        {
            overlaps[i]=length;
        }
        else
        {
            //find the first segment with time t_start > seg.t_start
            for(std::set<segment>::const_iterator it=segments_[i].begin(); it != segments_[i].end();++it)
	    {
                overlaps[i]+=segment_overlap(seg, *it);
            }
        }
        //std::cout<<clmagenta<<"weight of orbital: "<<orb<<" wrt orbital: "<<i<<" is: "<<std::exp(U_(orb,i)*overlaps[i])<<" for overlap: "<<overlaps[i]<<cblack<<std::endl;
        energy -= U_(orb,i)*overlaps[i];
        /*if(zero_order_orbital_occupied_[i]){
         std::cout<<"weight got an additional factor:"<<std::exp(-sgn*U_(orb,i)*overlaps[i])<<std::endl;
         }*/
    }
    //Leo: disable retarded interaction
    //this is the retarded interaction stuff
//    if(use_retarded_interaction_)
//    {
//        bool is_removal=false;
//        double retarded_weight=0;
//        if((seg.t_start_==0) && (seg.t_end_==beta_)){ //not really a segment but a full line
//            retarded_weight=0.;
//        }else{
//            for(int i=0;i<n_orbitals_;++i){
//                for(std::set<segment>::const_iterator it=segments_[i].begin(); it != segments_[i].end();++it){
//                    retarded_weight+=K_.interpolate(seg.t_start_-it->t_start_);
//                    retarded_weight-=K_.interpolate(seg.t_start_-it->t_end_);
//                    retarded_weight-=K_.interpolate(seg.t_end_-it->t_start_);
//                    retarded_weight+=K_.interpolate(seg.t_end_-it->t_end_);
//                }
//            }
//            retarded_weight-=K_.interpolate(seg.t_end_-seg.t_start_);
//        }
//        energy += retarded_weight;
//    }
    return energy;
}


void local_configuration::check_consistency()const 
{
  for(int i=0;i<n_orbitals_;++i)
  {
    if(order(i)<2) continue; //nothing to check if none or only one segment present
    //std::cout<<"testing orbital: "<<i<<" order: "<<order(i)<<std::endl;
    for(std::set<segment>::const_iterator it=segments_[i].begin();it!=segments_[i].end();++it)
    {
      std::set<segment>::const_iterator next_it=it; next_it++;
      if(next_it!=segments_[i].end())
      {
        //std::cout<<"testing it: "<<*it<<" next it: "<<*next_it<<std::endl;
        if(it->t_end_<it->t_start_)
        {
          std::cout<<*this<<std::endl;
          throw std::logic_error("consistency fail: segment does not go the right way!");
        }
        if(it->t_end_>next_it->t_start_)
        {
          std::cout<<*this<<std::endl;
          throw std::logic_error("consistency fail: segment overlaps with next segment");
        }
      }
      else
      {
        next_it=segments_[i].begin();
        //std::cout<<"testing it: "<<*it<<" next it: "<<*next_it<<std::endl;
        if(it->t_start_>it->t_end_)
        { //it wraps around
          if(it->t_end_>next_it->t_start_) 
          {
            std::cout<<*this<<std::endl;
            throw std::logic_error("consistency fail: wraparound segment overlaps with first segment");
          }
        }
        else
        {
          if(it->t_end_<it->t_start_) 
          {
            std::cout<<*this<<std::endl;
            throw std::logic_error("consistency fail: segment does not go the right way");
          }
        }
      }
    }
  }

  //Leo: also check the consistency of number of segments and antisegments
  for (int i=0; i<n_orbitals_; ++i) 
     check_n_segments_consistency(i);
}


//get a vector back that for every creation time has the occupation in all the other orbitals
//we could make this linear in tauprime, right now it is quadratic with lots of searches.
//returns a vector n(\tau') which we need for the improved estimators.
//this actually is n(tau) not n(tau') now - see also comments for get_F_prefactor
//moved determination of density to separate function which is used in get_segment_densities and get_density_vectors
/*
 void local_configuration::get_segment_densities(std::vector<std::vector<std::vector<double> > > &n_tauprime)const{
 n_tauprime.resize(n_orbitals_);
 for(int i=0;i<n_orbitals_;++i){
 n_tauprime[i].resize(segments_[i].size());
 int k=0;
 for(std::set<segment>::const_iterator it=segments_[i].begin();it!=segments_[i].end();++it,++k){
 n_tauprime[i][k].resize(n_orbitals_, 0.);
 //      double tauprime=it->t_start_;
 double tauprime=it->t_end_;
 for(int j=0;j<n_orbitals_;++j){
 if(i==j){
 n_tauprime[i][k][j]=0.;
 continue; //no need to check the same orbital.
 }
 if(segments_[j].size()==0){
 n_tauprime[i][k][j]=zero_order_orbital_occupied_[j]?1.:0.;
 continue;
 }else{
 //                    n_tauprime[i][k][j]=0.;
 //           for(std::set<segment>::const_iterator it2=segments_[j].begin();it2!=segments_[j].end();++it2){
 //           if(it2->t_end_>it2->t_start_){ //regular segment
 //           if(it2->t_start_<tauprime && tauprime <it2->t_end_){
 //           n_tauprime[i][k][j]=1.;
 //           }
 //           }else{
 //           if(tauprime<it2->t_end_ ||tauprime > it2->t_start_){
 //           n_tauprime[i][k][j]=1.;
 //           }
 //           }
 //find first segment after tauprime, call it it_after
 std::set<segment>::const_iterator it_after=segments_[j].upper_bound(segment(tauprime,0.)); //this is the segment that starts after tauprime
 if(it_after==segments_[j].end()) it_after=segments_[j].begin();
 //find the segment which is before it_after. call it it_before
 std::set<segment>::const_iterator it_before=it_after;
 if(it_before==segments_[j].begin()){ it_before=segments_[j].end(); } it_before--; //this is the segment that has the start time before tauprime
 
 //find out the iterator it overlaps with a segment in this orbital. two cases: either it does not wrap and is just in between, or it wraps and is in between.
 bool occuppied=(it_before->t_end_>it_before->t_start_ && (tauprime > it_before->t_start_ && tauprime<it_before->t_end_))
 || (it_before->t_end_<it_before->t_start_ &&(tauprime<it_before->t_end_ || tauprime > it_before->t_start_));
 //std::cout<<i<<" "<<k<<" "<<j<<std::endl;
 n_tauprime[i][k][j]=occuppied?1.:0.;
 //std::cout<<cred<<" segment before: "<<*it_before<<" segment after: "<<*it_after<<" time: "<<tauprime<<" is occuppied: "<<n_tauprime[i][k][j]<<cblack<<std::endl;
 }
 }
 }
 }
 }
 */
/*
 void local_configuration::get_segment_densities(std::vector<std::vector<std::vector<double> > > &n_tauprime)const{//not needed: each element of n_tauprime is requested only once
 n_tauprime.resize(n_orbitals_);
 for(int i=0;i<n_orbitals_;++i){
 n_tauprime[i].resize(segments_[i].size());
 int k=0;
 for(std::set<segment>::const_iterator it=segments_[i].begin();it!=segments_[i].end();++it,++k){
 n_tauprime[i][k].resize(n_orbitals_, 0.);
 //      double tauprime=it->t_start_;
 double tauprime=it->t_end_;
 for(int j=0;j<n_orbitals_;++j){
 if(i==j){
 n_tauprime[i][k][j]=0.;
 continue; //no need to check the same orbital.
 }
 n_tauprime[i][k][j]=density(j,tauprime);
 }
 }
 }
 }
 */


double local_configuration::density(int i, double tauprime) const
{//density on orbital i at time tauprime
  if(segments_[i].size()==0)  return zero_order_orbital_occupied_[i]?1.:0.;
  else
  {
    //find first segment after tauprime, call it it_after
    std::set<segment>::const_iterator it_after=segments_[i].upper_bound(segment(tauprime,0.)); //this is the segment that starts after tauprime
    if(it_after==segments_[i].end())  it_after=segments_[i].begin();
    //find the segment which is before it_after. call it it_before
    std::set<segment>::const_iterator  it_before=it_after;
    if(it_before==segments_[i].begin())
    { 
	it_before=segments_[i].end(); 
    } 
    it_before--; //this is the segment that has the start time before tauprime
    
    //find out the iterator it overlaps with a segment in this orbital. two cases: either it does not wrap and is just in between, or it wraps and is in between.
    bool occuppied=(it_before->t_end_>it_before->t_start_ && (tauprime > it_before->t_start_ && tauprime<it_before->t_end_))
    || (it_before->t_end_<it_before->t_start_ &&(tauprime<it_before->t_end_ || tauprime > it_before->t_start_));
    return occuppied?1.:0.;
    //std::cout<<cred<<" segment before: "<<*it_before<<" segment after: "<<*it_after<<" time: "<<tauprime<<" is occuppied: "<<n_tauprime[i][k][j]<<cblack<<std::endl;
  }
}

//Leo: disable retarded interaction
//double local_configuration::interaction_density_integral(std::set<segment>::const_iterator &it_i) const{
//  //for orbital i, compute \sum_j int_0^beta dt U_ret(tau - t) n_j(t) using the primitive K'(tau) of U(tau)=K''(tau)
//  double integral=0.0; double sgn;
//  for(int j=0; j<n_orbitals_; ++j){
//    for(std::set<segment>::const_iterator it_j=segments_[j].begin(); it_j!=segments_[j].end();++it_j){
//      integral+= K_.interpolate_deriv(it_j->t_end_ - it_i->t_end_) - K_.interpolate_deriv(it_j->t_start_ - it_i->t_end_);
//    }
//  }
//  return integral-2*K_.interpolate_deriv(0.0); //this is the same segment, same time contribution due to the fact that n and c do not commute in this case
//}


//compute the prefactor that will go into the Green's function measurement using Fw
//NOTE: from the equation of motion, we have the correlation function
//H=<T n(tau_1)c(tau_1)c^dagger(tau_2)c(tau_3)c^dagger(tau_4)>
//hence the density has to be attached to the annihilator c(tau_1)
//and we need to evaluate F_prefactor  for the annihilator times
//segments->t_end_. From the equation of motion, we further have the
//correlator F(tau-tau')=-<T c(tau)c^dagger(tau') n(tau')>, where the
//density is coupled to the creator. Strictly we would have to evaluate
//the prefactor for the creator times. However one can show that for a
//function that depends on a single time difference only (or diagonal
//in frequency), one obtains exactly the same result if the prefactor is
//evaluated for the annihilator times. This is *not* the case anymore if F(tau-tau')
//is no longer a diagonal matrix, i.e. for non-diagonal hybridization!
//Note also that for the correlator H it *does* make a difference  whether
//n is coupled the creator or annihilator. To save computations and memory, we
//evaluate F_prefactor for annihilator times only, in favor of the correlator H.
void local_configuration::get_F_prefactor(std::vector<std::map<double,double> > &F_prefactor)const
{
  for(std::size_t i=0; i<n_orbitals_; ++i) F_prefactor[i].clear();
  //this is F_prefactor[orbital i][segment k in orbital i]
  for(std::size_t i=0;i<n_orbitals_;++i)
  {
    std::size_t k=0; //Leo: why do we need k here? It seems useless! 
    for(std::set<segment>::const_iterator it=segments_[i].begin(); it!=segments_[i].end(); ++it,++k)
    {
      F_prefactor[i][it->t_end_] = 0;
      for(std::size_t j=0; j<n_orbitals_; ++j)
      {
        F_prefactor[i][it->t_end_] += 0.5*(U_(i,j)+U_(j,i))*density(j,it->t_end_);
      }
      //Leo: disable retarded interaction
//      if(use_retarded_interaction_)
//      {//also contribute for j==i
//          F_prefactor[i][it->t_end_] += interaction_density_integral(it);
//      }
    }
  }
}


/*
 void local_configuration::measure_density(std::vector<double> &densities, double sign) const{
 for(int i=0;i<n_orbitals_;++i){
 if(segments_[i].size()==0){
 if(zero_order_orbital_occupied_[i]){
 densities[i]+=sign;
 }else{
 densities[i]+=0.;
 }
 }else{
 for(segment_container_t::const_iterator it=segments_[i].begin();it!=segments_[i].end();++it){
 double dist=it->t_end_-it->t_start_;
 if(dist<0) dist+=beta_; //wrap around
 densities[i]+=dist/beta_*sign;
 }
 }
 }
 }
 */


void local_configuration::measure_density(std::vector<double> &densities, double sign) const
{
  for(int i=0;i<n_orbitals_;++i) densities[i]+=segment_density(i)*sign;
}


double local_configuration::segment_density(int i) const
{//density on orbital i (total length of segments)/beta
  double density=0.;
  if(segments_[i].size()==0 && zero_order_orbital_occupied_[i]) 
    density=1.;
  else
  {
    for(segment_container_t::const_iterator it=segments_[i].begin(); it!=segments_[i].end(); ++it)
    {
      double dist=it->t_end_-it->t_start_;
      if(dist<0) dist+=beta_; //wrap around
      density+=dist/beta_;
    }
  }
  return density;
}


double local_configuration::measure_nn(int i, int j) const{
  if(i==j) return segment_density(i); //does not make much sense to call this function for i==j; also correct without this, but slow for i==j
  if(segments_[i].size()==0){
    if(zero_order_orbital_occupied_[i]) return segment_density(j);
    else return 0.0;
  }
  else{//segments in orbital i
    if(segments_[j].size()==0){
      if(zero_order_orbital_occupied_[j]) return segment_density(i);
      else return 0;
    }
    else{
      //segments in i and j -> get overlap
      double overlap=0.;
      for(segment_container_t::const_iterator iti=segments_[i].begin();iti!=segments_[i].end();++iti)//this can be made faster (linear in k)
        for(segment_container_t::const_iterator itj=segments_[j].begin();itj!=segments_[j].end();++itj)
          overlap+=segment_overlap(*iti, *itj);
      return overlap/beta_;
    }
  }
}


/*
 void local_configuration::get_density_vectors(std::vector<std::vector<double> > &n_vector) const{
 for(int i=0;i<n_orbitals_;++i){
 int N_nn=n_vector[i].size()-1;
 for(int n=0;n<=N_nn;++n){
 double tau= n*beta_/static_cast<double>(N_nn);
 n_vector[i][n]=density(i,tau);
 }
 }
 }
 */


void local_configuration::get_density_vectors(std::vector<std::vector<double> > &n_vectors) const{
  //this gives the same result as using local_configuration.density(), but is faster (linear in k instead of k log k for search in the map)
  for(int i=0; i<n_orbitals_; ++i){
    int N_nn=n_vectors[i].size()-1;
    if(segments_[i].size()==0){
      if(zero_order_orbital_occupied_[i]) std::fill(n_vectors[i].begin(), n_vectors[i].end(), 1);
      else std::fill(n_vectors[i].begin(), n_vectors[i].end(), 0);
    }
    else{
      std::fill(n_vectors[i].begin(), n_vectors[i].end(), 1);
      segment_container_t::const_iterator it=segments_[i].end(); it--;
      if(it->t_end_<it->t_start_) n_vectors[i][0]=1;//last segment winds around the circle //n(0)=1
      else n_vectors[i][0]=0;//n(0)=0
      
      int index; // mark segment start and end points
      for(segment_container_t::const_iterator it=segments_[i].begin(); it!=segments_[i].end(); ++it){
        index = (int)(it->t_start_/beta_*N_nn+1);//int conversion always rounds off; assumption: time is never exactly beta
        n_vectors[i][index] *= -1;
        index = (int)(it->t_end_/beta_*N_nn+1);
        n_vectors[i][index] *= -1;
      }
      // fill vector with occupation number
      for(std::size_t n=1; n<n_vectors[i].size(); n++){
        if(n_vectors[i][n]==-1) n_vectors[i][n]=1-n_vectors[i][n-1];//segment starts or ends -> occupation number changes
        else n_vectors[i][n]=n_vectors[i][n-1];//on the same segment -> occupation number identical to that of previous segment
      }
    }
  }
}


void local_configuration::measure_nnw(int i, std::vector<double> &nnw_re, double sign) const{
  int N_W=nnw_re.size();
  //wm=0 has to be treated separately
  nnw_re[0]+= segment_density(i)*beta_;//length of segments
  if(N_W == 1) return;
  for(segment_container_t::const_iterator it=segments_[i].begin();it!=segments_[i].end();++it){//same contribution for winding segments
    double wm=0;
    double dw=2*M_PI/beta_;
    std::complex<double> exp_s=1.;
    std::complex<double> exp_e=1.;
    std::complex<double> dexp_s=std::exp(std::complex<double>(0,dw*it->t_start_) );
    std::complex<double> dexp_e=std::exp(std::complex<double>(0,dw*it->t_end_) );
    for(int m=1;m<N_W;++m){
      wm += dw;
      exp_s*=dexp_s;
      exp_e*=dexp_e;
      nnw_re[m]+=real((exp_e-exp_s)/std::complex<double>(0,wm))*sign;//!!
    }
  }
}


void local_configuration::state_map_segment_insert(state_map &states, const segment &s, int index) const{
  //works also for the case where an operator is inserted exactly at the point of an already present kink
  if(s.t_end_==s.t_start_) return;
  if(s.t_end_<s.t_start_){
    std::cerr<<"fatal logic error inside state_map_segment_insert." << s.t_start_ << " " << s.t_end_ <<std::endl;
    throw std::runtime_error(std::string(__FUNCTION__)+" works only segments where the annihilator comes strictly after the creator");
  }
  state_map::iterator next;
  int current_state=0;

  next=states.upper_bound(s.t_start_);
  if(next==states.end()){ //the new creator is the last operator
    if(!states.empty()){ //there is at least one kink before (can be at exactly the same time) -> get it's state
      next--; current_state=next->second;
    }
    else current_state=0; //no kink before, impurity is in state 0
  }
  else{//not past the end -> there is an operator with a later time
    if(next!=states.begin()){//the new creator is not at the beginning, there is an earlier kink (can be at exactly the same time) -> get it's state
      next--; current_state=next->second;
    }
    else{//it's at the beginning, we'll be inserting the first one
      current_state=0;
    }
  }
  //insert creator
  states[s.t_start_]=current_state+index;//+ for creator //this even works if the kink was exactly at time tau

  for(state_map::iterator it=states.upper_bound(s.t_start_); it!=states.upper_bound(s.t_end_); ++it){//if there is an operator exactly at s.t_end_, it will be updated
    it->second+=index;
  }

  next=states.upper_bound(s.t_end_);

  if(next!=states.begin()){//can never be states.begin(), we have inserted an operator at an earlier time!
    next--; current_state=next->second;
  }
  //insert annihilator
  states[s.t_end_]=current_state-index;//- for annihilator

}


void local_configuration::measure_sector_statistics(std::vector<double> &sector_statistics, double sign) const{//complexity: linear in n_orbitals_
  state_map states; //key is time, value is state
  //state map is organized such that an entry with time t and value s means that the impurity
  //is in state s from time t to the time of the next entry (or beta, for the last one).

  //std::fill(sector_statistics.begin(),sector_statistics.end(),0);
  int full_line_states=0;

  for(int i=0;i<n_orbitals_;++i){
    int index=1<<i;
    if(zero_order_orbital_occupied_[i]){
      full_line_states+=index;//keep track of which time lines (orbitals) are fully occupied
      continue; //no segments->we're done for this orbital
    }
    for(segment_container_t::const_iterator it=segments_[i].begin(); it!=segments_[i].end(); ++it){
      if(it->t_end_<it->t_start_){//winding segment
        state_map_segment_insert(states,segment(0.         ,it->t_end_),index);
        state_map_segment_insert(states,segment(it->t_start_,beta_    ),index);
      }
      else state_map_segment_insert(states,*it,index);
    }
  }
  //if there are no segments, the state did not change:
  if(states.empty()){
    sector_statistics[full_line_states]+=1.*sign;
//    cout << "no segments! " << full_line_states << endl;
  }
  else{
    //otherwise count intervals
    state_map::iterator it=states.end(); it--; int state=it->second; //state of last operator is equal to initial state (interval winds around)
    double tau=0.; int max_state=1<<n_orbitals();
    for(state_map::iterator it=states.begin();it!=states.end();++it){
      sector_statistics[state+full_line_states]+=(it->first-tau)/beta_*sign;
      tau=it->first; state=it->second;
      if(state>=max_state || state<0){ std::cout << "something went wrong! state=" << state << std::endl; exit(1); }
    }
    sector_statistics[state+full_line_states]+=(beta_-tau)/beta_*sign;//don't forget the last interval
  }

}


// return the total weight of the local configuration
//this is an expensive debug operation
double local_configuration::full_weight() const{
  return 1.;
}


//Leo: get the number of segments and antisegments after inserting an antisegment
std::vector<int> local_configuration::get_new_n_segments_insert_antisegment(const segment &new_antisegment, int orbital)
{
  std::vector<int> n_segments_temp = get_n_segments(orbital); //old numbers
  int color = new_antisegment.c_start_;
  
  if(order(orbital)==0) //special case: 0-th order orbital
  {
      n_segments_temp[color]+=1; n_segments_temp[color+n_env_]+=1; 
      return n_segments_temp;
  }
  else  //general case
  {
    std::set<segment>::const_iterator it=segments_[orbital].upper_bound(new_antisegment);
    if(it==segments_[orbital].begin()) it=segments_[orbital].end(); //wrap around
    it--;

//    std::cout << std::endl;
//    std::cout << "new_antisegment = " << new_antisegment << std::endl;
//    std::cout << "it = " << *it << std::endl;
//    std::cout << std::endl;

    n_segments_temp[color+n_env_]+=1; //add one antisegment
    if(it->c_start_ == it->c_end_) 
        n_segments_temp[it->c_start_]-=1; //break one segment
    if(color==it->c_start_) 
       n_segments_temp[color]+=1; //add the left segment 
    if(color==it->c_end_) 
       n_segments_temp[color]+=1; //add the right segment
    return n_segments_temp; 
  }  
}


//Leo: get the number of segments and antisegments after inserting an segment
std::vector<int> local_configuration::get_new_n_segments_insert_segment(const segment &new_segment, int orbital)
{
  std::vector<int> n_segments_temp = get_n_segments(orbital); //old numbers
  int color = new_segment.c_start_;
  
  if(order(orbital)==0) //special case: 0-th order orbital
  {
      n_segments_temp[color]+=1; n_segments_temp[color+n_env_]+=1; 
      return n_segments_temp;
  }
  else  //general case
  {
    std::set<segment>::const_iterator it_later=segments_[orbital].upper_bound(new_segment);
    if(it_later==segments_[orbital].end())  { it_later=segments_[orbital].begin(); }

    std::set<segment>::const_iterator it_earlier=it_later;
    if(it_later==segments_[orbital].begin())  { it_earlier = segments_[orbital].end(); }
    it_earlier--;

    n_segments_temp[color]+=1; //add one segment
    if(it_earlier->c_end_ == it_later->c_start_) 
       n_segments_temp[it_later->c_start_+n_env_]-=1; //break one antisegment
    if(color==it_later->c_start_) 
       n_segments_temp[color+n_env_]+=1; //add the right antisegment 
    if(color==it_earlier->c_end_) 
       n_segments_temp[color+n_env_]+=1; //add the left antisegment
    return n_segments_temp; 
  }  
}



//Leo: get the number of segments and antisegments after removing a segment
std::vector<int> local_configuration::get_new_n_segments_remove_segment(const segment &new_segment, int orbital)
{
  std::vector<int> n_segments_temp = get_n_segments(orbital); //old numbers
  int color = new_segment.c_start_;
  
  if(order(orbital)==1) //special case: 1-th order orbital
  {
      n_segments_temp[color]-=1; n_segments_temp[color+n_env_]-=1; 
      return n_segments_temp;
  }
//  else if(order(orbital)==2) //special case: 2-th order orbital
//  {
//      //in this case it_later = it_earlier, so just do it once. 
//      //Leo: this case seems to very similar to the next one...maybe can be combined?
//      std::set<segment>::const_iterator it_later=segments_[orbital].upper_bound(new_segment);
//      if(it_later==segments_[orbital].end()) { it_later = segments_[orbital].begin(); }
//      n_segments_temp[color]-=1; //substract the segment
//      n_segments_temp[it_later->c_start_+n_env_]+=1; // add one antisegment
//      if(color==it_later->c_start_) 
//         n_segments_temp[color+n_env_]-=1; //subtract the right antisegment  
//      if(color==it_later->c_end_) 
//         n_segments_temp[color+n_env_]-=1; //subtract the end antisegment  
//      return n_segments_temp;
//  }
  else  //general case
  {
      std::set<segment>::const_iterator it_later=segments_[orbital].upper_bound(new_segment);
      if(it_later==segments_[orbital].end())  { it_later = segments_[orbital].begin(); }
  
      std::set<segment>::const_iterator it_earlier=segments_[orbital].find(new_segment);
      if(it_earlier==segments_[orbital].begin())  { it_earlier = segments_[orbital].end(); }
      it_earlier--; 

      if( order(orbital)==2 && (it_later != it_earlier) ) // quick check for n=2
          throw std::logic_error("local_configuration::get_new_n_segments_remove_segment does not capture the correct segments at order n=2. Abort.");
  
      n_segments_temp[color]-=1; //subtract one segment
      if(it_earlier->c_end_ == it_later->c_start_) 
         n_segments_temp[it_later->c_start_+n_env_]+=1; //add one antisegment
      if(color==it_later->c_start_) 
         n_segments_temp[color+n_env_]-=1; //subtract the right antisegment 
      if(color==it_earlier->c_end_) 
         n_segments_temp[color+n_env_]-=1; //subtract the left antisegment 
      return n_segments_temp; 
  }  
}


//Leo: get the number of segments and antisegments after removing an antisegment
std::vector<int> local_configuration::get_new_n_segments_remove_antisegment(const segment &new_antisegment, int orbital)
{
  std::vector<int> n_segments_temp = get_n_segments(orbital); //old numbers
  int color = new_antisegment.c_start_;
  
  if(order(orbital)==1) //special case: 1-th order orbital
  {
      n_segments_temp[color]-=1; n_segments_temp[color+n_env_]-=1; 
      return n_segments_temp;
  }
  else  //general case
  {
    std::set<segment>::const_iterator it_later=segments_[orbital].find(new_antisegment);

    std::set<segment>::const_iterator it_earlier=it_later;
    if(it_later==segments_[orbital].begin()) { it_earlier = segments_[orbital].end(); }
    it_earlier--;

//    if( order(orbital)==2 && (it_later != it_earlier) ) // quick check for n=2
//        throw std::logic_error("local_configuration::get_new_n_segments_remove_segment does not capture the correct segments at order n=2. Abort.");

    n_segments_temp[color+n_env_]-=1; //subtract one antisegment
    if(it_earlier->c_start_ == it_later->c_end_) 
       n_segments_temp[it_later->c_end_]+=1; //add one segment
    if(color==it_later->c_end_) 
       n_segments_temp[color]-=1; //subtract the right segment 
    if(color==it_earlier->c_start_) 
       n_segments_temp[color]-=1; //subtract the left segment

    return n_segments_temp; 
  }  
}



void local_configuration::check_n_segments_consistency(int orbital) const
{  
   if(order(orbital)==0) return; //empty or filled line, no need to check

   std::vector<int> n_segments_count (2*n_env_, 0);
   std::vector<int> n_segments_current = get_n_segments(orbital);
   std::set<segment>::const_iterator it_next;
   for(std::set<segment>::const_iterator it=segments_[orbital].begin(); it != segments_[orbital].end(); ++it)
   {
       if(it->c_start_==it->c_end_)  n_segments_count[it->c_start_]++; // count number of segments
       it_next = it; it_next++;
       if(it_next==segments_[orbital].end())  { it_next=segments_[orbital].begin(); } //hit the last segment, wrap around
       if(it->c_end_==it_next->c_start_) n_segments_count[it->c_end_+n_env_]++; // count number of antisegments
   }

//   std::set<segment>::const_iterator it_begin=segments_[orbital].begin();
//   std::set<segment>::const_iterator it_end=segments_[orbital].end(); it_end--;
//   if(it_end->c_end_==it_begin->c_start_) n_segments_count[it_end->c_end_+n_env_]++; // the wrapping antisegment, if any

   for(int i=0; i<2*n_env_; i++)
   {
       if(n_segments_count[i] != n_segments_current[i])
       {
          std::stringstream temp_stream1, temp_stream2; // stringstream used for the conversion
          for(int j=0; j<2*n_env_; j++)  { temp_stream1 << n_segments_count[j]; }           
          std::string temp1 = "(" + temp_stream1.str() + ")";
          for(int j=0; j<2*n_env_; j++)  { temp_stream2 << n_segments_current[j]; }           
          std::string temp2 = "(" + temp_stream2.str() + ")";
          throw std::logic_error("Error in the number of segments and antisegments! It should be " + temp1 + ", but it was " + temp2 + ". Abort.");
       }    
   }
}



