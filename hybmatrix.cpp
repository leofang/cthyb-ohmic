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

#include "hybmatrix.hpp"
#include "hybblasmatrix.hpp" //Leo: for matrix inversion; TODO: remove this after debugging!
//#include "combinatorial.hpp"   //Leo: for test purpose

//this was changed with update 51!
//nomenclature: the c        is at the segment end  , so the time for c        is new_segment->t_end_
//              the c_dagger is at the segment start, so the time for c_dagger is new_segment->t_start_

//the direct matrix F has the structure F_ij=Delta(\tau_i_end - \tau_j_start) // (check)
//that means that the creation operators enter the rows (the i-s), and the annihilation operators form columns (the j-s) //unchanged (check)

//Careful with the ordering: the ORDERED matrix has creation operators and annihilation operators ordered by start and end time,
//in case of a wraparound segment the wrapping c is the last entry.
//this matrix can always be turned into an ordered matrix by running rebuild_ordered_matrix.
//in general the matrix is not ordered (that would be inefficient), but permutation_sign keeps track of the shift of rows
//and columns.

//Leo: c_dagger & c in the code mean the impurity operator, not bath operator!


//Leo: for swapping two hybmatrix objects
void swap(hybmatrix &A, hybmatrix &B)
{
   using std::swap;

   //first do the blas_matrix swap
   static_cast<blas_matrix&>(A).swap( static_cast<blas_matrix&>(B) );

   //then swap other private members in hybmatrix
   swap(A.cdagger_index_map_  , B.cdagger_index_map_  );
   swap(A.c_index_map_        , B.c_index_map_        );
   swap(A.c_cdagger_map_      , B.c_cdagger_map_      );
   //swap(A.n_env_              , B.n_env_              );
   swap(A.Q                   , B.Q                   );
   swap(A.R                   , B.R                   );
   swap(A.PinvQ               , B.PinvQ               );
   swap(A.S                   , B.S                   );
   swap(A.S_tilde_inv         , B.S_tilde_inv         );
   swap(A.S_tilde             , B.S_tilde             );
   swap(A.weight_ratio_       , B.weight_ratio_       );
   swap(A.permutation_sign_   , B.permutation_sign_   );
   swap(A.time_ordering_sign_ , B.time_ordering_sign_ );
   swap(A.disordered_times    , B.disordered_times    );
   swap(A.inv_Delta_small     , B.inv_Delta_small     );
   swap(A.determinant_        , B.determinant_        );
   swap(A.determinant_old_    , B.determinant_old_    );
   //swap(A.beta_               , B.beta_               );
   //swap(A.measure_g2w_        , B.measure_g2w_        );
   //swap(A.measure_h2w_        , B.measure_h2w_        );
}


//compute the hybridization weight change when an operator pair is inserted
double hybmatrix::hyb_weight_change_insert(const segment &new_segment, int orbital, const hybfun &Delta)
{
  Q.resize(size());
  R.resize(size());
  PinvQ.resize(size());
  //column  Delta_i,last
  for(hyb_map_t::const_iterator it=c_index_map_.begin();it!=c_index_map_.end();++it)
  {
    Q[it->second]=Delta.interpolate(it->first-new_segment.t_start_, orbital); //this is the new column Q
  }
  
  //row Delta_last,i
  for(hyb_map_t::const_iterator it=cdagger_index_map_.begin();it != cdagger_index_map_.end();++it)
  {
    R[it->second]=Delta.interpolate(new_segment.t_end_-it->first, orbital);  //this is the new row R
  }
  S=Delta.interpolate(new_segment.t_end_-new_segment.t_start_, orbital);  //this is the entry S
  //std::cout<<clblue<<"last entry S is: "<<S<<cblack<<std::endl;
  
  S_tilde_inv=S;
  fortran_int_t s=size();
  if(s>0)
  {
    right_multiply(Q, PinvQ); //dgemv
    fortran_int_t inc=1;
    S_tilde_inv-=FORTRAN_ID(ddot)(&s, &(R[0]),&inc,&(PinvQ[0]),&inc);
  }

  //a -1 from the anticommutator from the wraparound segment
  if(new_segment.t_end_<new_segment.t_start_) { weight_ratio_=-S_tilde_inv; }
  else { weight_ratio_=S_tilde_inv; }

  return weight_ratio_;
}


//actually insert an operator pair and change the configuration
void hybmatrix::insert_segment(const segment &new_segment, int orbital)
{
  //std::cout<<clred<<"upon entering insert segment: "<<*this<<cblack<<std::endl;
  //consistency_check();
  //enlarge the M-matrix by one
  int last=size();
  resize(size()+1);
  
  //last element
  operator()(last, last)=1./S_tilde_inv;
  
  fortran_int_t sm1=size()-1;
  if(sm1>0)
  { //this is exactly the content of the loops above, in dger/dgemv blas calls.
   char trans='T', notrans='N';
   double alpha=-1./S_tilde_inv, beta=0.;
   fortran_int_t inc=1;
   fortran_int_t ms=memory_size();
   FORTRAN_ID(dgemv)(&  trans, &sm1, &sm1, &alpha, &(operator()(0,0)), &ms, &(Q[0]), &inc, &beta, &(operator()(0,last)), &ms);
   FORTRAN_ID(dgemv)(&notrans, &sm1, &sm1, &alpha, &(operator()(0,0)), &ms, &(R[0]), &inc, &beta, &(operator()(last,0)), &inc);
   alpha=S_tilde_inv;
   FORTRAN_ID(dger)(&sm1, &sm1, &alpha,&(operator()(last,0)), &inc, &(operator()(0,last)), &ms, &(operator()(0,0)), &ms);
  }
  
  // add the new segment times:
  cdagger_index_map_.insert(std::make_pair(new_segment.t_start_, last));
  c_index_map_      .insert(std::make_pair(new_segment.t_end_  , last));
  c_cdagger_map_    .insert(std::make_pair(new_segment.t_start_, 1)); //Leo: 1 means c_dagger  
  c_cdagger_map_    .insert(std::make_pair(new_segment.t_end_, 0));   //Leo: 0 means c
  
  //keep track of the wraparound sign
  if(new_segment.t_start_>new_segment.t_end_) { permutation_sign_*=-1.; }
  //std::cout<<clred<<*this<<cblack<<std::endl;
  //consistency_check(); 

  //Leo: keep track of the time ordering sign due to operator time ordering
  //TODO: move this to Green function measurements because this sign wouldn't affect other local measurements.
  //time_ordering_sign_check(time_ordering_sign_, disordered_times);
  time_ordering_sign_check();
}


//compute the hybridization weight change when an operator pair is removed
double hybmatrix::hyb_weight_change_remove(const segment &new_segment, int orbital, const hybfun &Delta)
{
  //std::cout<<clgreen<<"proposing to remove the segment: "<<new_segment<<cblack<<std::endl;
  int k1=cdagger_index_map_[new_segment.t_start_];
  int k2=c_index_map_[new_segment.t_end_];
  
  S_tilde = operator()(k1,k2);
  weight_ratio_=1./S_tilde; //Leo: this line is redundant! the weight ratio is simply S_tilde...
  
  // take care of sign changes due to wraparound segments
  if(new_segment.t_start_>new_segment.t_end_) { weight_ratio_ *=-1; }
  
  //std::cout<<"returning removal weight ratio: "<<weight_ratio_<<" S_tilde is: "<<S_tilde<<std::endl;
  //std::cout<<S_tilde<<" weight ratio: "<<weight_ratio_<<std::endl;
  return weight_ratio_;
}


//actually remove an operator pair and change the configuration
void hybmatrix::remove_segment(const segment &new_segment, int orbital)
{
  //std::cout<<clblue<<"upon entering remove segment: "<<*this<<cblack<<std::endl;
  //consistency_check();
  
  //find row and column indices
  std::size_t thisrow=cdagger_index_map_[new_segment.t_start_];
  std::size_t thiscolumn=c_index_map_[new_segment.t_end_];
  //if(thisrow != thiscolumn) throw std::logic_error("row and column got screwed up!");
  std::size_t last=size()-1;
  
  //swap row and column of thisrow and the last row. Det picks up a minus sign for each interchange.
  double row_column_sign=1.;
  if(thisrow != last)
  {
    swap_row(thisrow, last);
    row_column_sign*=-1.;
  }
  if(thiscolumn != last)
  {
    swap_column(thiscolumn, last);
    row_column_sign*=-1.;
  }
  
  //std::cout<<"row column sign is: "<<row_column_sign<<std::endl;
  //std::cout<<"row: "<<thisrow<<" column: "<<thiscolumn<<" size: "<<size()<<std::endl;
  permutation_sign_*=row_column_sign;
  
  //perform rank one update
  /*for (std::size_t i=0; i<last; i++) {
   for (std::size_t j=0; j<last; j++) {
   operator()(i,j) -=operator()(i,last)*operator()(last,j)/operator()(last,last);
   }
   }*/
  fortran_int_t sm1=size()-1;
  if(sm1>0)
  { //this is exactly the content of the loops above, in dger/dgemv blas calls.
    double alpha=-1./operator()(last,last);
    fortran_int_t inc=1;
    fortran_int_t ms=memory_size();
    FORTRAN_ID(dger)(&sm1, &sm1, &alpha,&(operator()(last,0)), &inc, &(operator()(0,last)), &ms, &(operator()(0,0)), &ms);
  }
  
  //adjust index of operator that pointed to last, let it point to thisrow instead
  for(hyb_map_t::iterator it=c_index_map_.begin();it!=c_index_map_.end();++it)
  {
    if(it->second==last) { it->second=thiscolumn; break; }
  }
  for(hyb_map_t::iterator it=cdagger_index_map_.begin();it!=cdagger_index_map_.end();++it)
  {
    if(it->second==last) { it->second=thisrow; break; }
  }
  if(new_segment.t_start_>new_segment.t_end_)
  {
    permutation_sign_*=-1.;
    //std::cout<<"additional permutation sign flip: t_start: "<<new_segment.t_start_<<" t_end: "<<new_segment.t_end_<<std::endl;
  }
  
  //shrink matrix size by one
  resize(size()-1);
  cdagger_index_map_.erase(new_segment.t_start_);
  c_index_map_      .erase(new_segment.t_end_);

  c_cdagger_map_    .erase(new_segment.t_start_);   
  c_cdagger_map_    .erase(new_segment.t_end_);   
  
  //std::cout<<clblue<<*this<<cblack<<std::endl;
  //consistency_check();
  //debug
  /*determinant_/=S_tilde;
   std::cout<<*((blas_matrix*)this)<<std::endl;
   std::cout<<clred<<"incremental determinant: "<<determinant_<<" actual determinant: "<<determinant()<<" prev determinant: "<<determinant_old_<<cblack<<std::endl;
   determinant_old_=determinant_;*/

  //Leo: keep track of the time ordering sign due to operator time ordering
  //TODO: move this to Green function measurements because this sign wouldn't affect other local measurements.
  //time_ordering_sign_check(time_ordering_sign_, disordered_times); 
  time_ordering_sign_check();
}


std::ostream &operator<<(std::ostream &os, const hybmatrix &hyb_mat)
{
  os<< "hyb matrix size: "<<hyb_mat.size()<<", permutation sign: "<<hyb_mat.permutation_sign_ \
    << ", time ordering sign: " << hyb_mat.sign() << std::endl;
  
  if(hyb_mat.size()==0) //don't print the rest if nothing to print
     return os; 
   
  os << "anti-ordered times: ";
  for(std::set<double>::const_iterator it=hyb_mat.disordered_times.begin(); it!=hyb_mat.disordered_times.end(); it++)
      os << *it << ", ";
  os << std::endl;

  os<<"c map: " << std::endl;
  for(hyb_map_t::const_iterator it=hyb_mat.c_index_map_.begin(); it!=hyb_mat.c_index_map_.end(); ++it)
  {
    os<<"("<<it->first<<", "<<it->second<<") ";
  }
  os<<std::endl;

  os<<"cdagger map: " << std::endl;
  for(hyb_map_t::const_iterator it=hyb_mat.cdagger_index_map_.begin(); it!=hyb_mat.cdagger_index_map_.end(); ++it)
  {
    os<<"("<<it->first<<", "<<it->second<<") ";
  }
  os<<std::endl;

  //Leo: for test purpose
  os<<"c & cdagger map: (0 means c, 1 means cdagger; from zero to beta)" << std::endl;
  os << "(";
  for(hyb_map_t::const_iterator it=hyb_mat.c_cdagger_map_.begin(); it!=hyb_mat.c_cdagger_map_.end(); ++it)
  {
    //os<<"("<<it->first<<", "<<it->second<<") "; //Leo: the value (second) is either 0 (c) or 1 (cdagger)
    os << it->second; 
  }
  os << ")" << std::endl;

  //Leo: output the whole blasmatrix object for debug purpose; adapted from ostream of blas_matrix 
  //     Note the format is presented in the Mathematica list structure
  os<<"matrix: " << std::endl << "{";
  for(int i=0; i<hyb_mat.size(); ++i)
  {
    os<<"{";
    for(int j=0; j<hyb_mat.size(); ++j)
    {
      os<< hyb_mat.operator()(i,j) << (j<hyb_mat.size()-1?", ":"");
    }
    if(i<hyb_mat.size()-1)  
       os << "}," << std::endl << " ";
    else
       os << "}";
  }
  os << "}" << std::endl;

  //Leo: output the inverted matrix; make sense only when its dimension is nonzero
  if(hyb_mat.size()>=1)
  {
    blas_matrix temp_matrix(hyb_mat);
    temp_matrix.invert();
    os<<"inverted matrix: " << std::endl << "{";
    for(int i=0; i<temp_matrix.size(); ++i)
    {
      os<<"{";
      for(int j=0; j<temp_matrix.size(); ++j)
      {
        os<< temp_matrix.operator()(i,j) << (j<hyb_mat.size()-1?", ":"");
      }
      if(i<temp_matrix.size()-1)  
         os << "}," << std::endl << " ";
      else
         os << "}";
    }
    os << "}" << std::endl;
  }
  else
    os << "inverted matrix:\n{ }" << std::endl;

  //os << std::endl;
  return os;
}


void hybmatrix::rebuild_hyb_matrix(int orbital, const hybfun &Delta)
{
  //Leo: this line is commented out as it makes no sense --- why we need a fresh blas_matrix copy here?
  //blas_matrix bup(*this);
  
  if(size()<1) return;
  
  //build the matrix inverse:
  for(hyb_map_t::const_iterator it_start=c_index_map_.begin(); it_start != c_index_map_.end(); ++it_start)
  {
    for(hyb_map_t::const_iterator it_end=cdagger_index_map_.begin(); it_end != cdagger_index_map_.end(); ++it_end)
    {
      operator()(it_start->second,it_end->second)=Delta.interpolate(it_start->first-it_end->first, orbital);
    }
  }
  //...then invert it.
  invert();
}


//Leo: this version would not properly work when N_ENV>1, because it relies on the fact
//     that c and c^dagger are alternating, which is not necessarily true in the multi-
//     color cases.
//void hybmatrix::rebuild_ordered_hyb_matrix(int orbital, const hybfun &Delta)
//{
//  if(size()<2) return;
//  //std::cout<<"on entry rebuild orderd: full weight: "<<full_weight()<<" permutation sign: "<<permutation_sign_<<std::endl;
//  //std::cout<<*this<<std::endl;
//  //std::cout<<clblue<<*(blas_matrix*)this<<cblack<<std::endl;
//
//  //order the times properly
//  int k=0;
//  hyb_map_t::iterator it_bup;
//  for(hyb_map_t::iterator it_end=cdagger_index_map_.begin(); it_end != cdagger_index_map_.end();)
//  {
//    it_bup=it_end++;
//    std::pair<double,int> new_entry=*it_bup;
//    new_entry.second=k++;
//    cdagger_index_map_.erase(it_bup);
//    cdagger_index_map_.insert(new_entry);
//  }
//
//  k=0;
//  for(hyb_map_t::iterator it_start=c_index_map_.begin(); it_start != c_index_map_.end();)
//  {
//    it_bup=it_start++;
//    std::pair<double,int> new_entry=*it_bup;
//    new_entry.second=((it_bup==c_index_map_.begin()) && (it_bup->first<cdagger_index_map_.begin()->first))?c_index_map_.size()-1:k++;
//    c_index_map_.erase(it_bup);
//    c_index_map_.insert(new_entry);
//  }
//  
//  //if we have an overlapping segment we need a permutation sign of -1, otherwise it is 1 in the ordered case.
//  if(size()==0)
//  {
//    permutation_sign_=1.;
//  }
//  else
//  {
//    if(c_index_map_.begin()->first<cdagger_index_map_.begin()->first)
//    {
//      permutation_sign_=-1.;
//    }
//    else
//    {
//      permutation_sign_=1.;
//    }
//  }
//  //then rebuild the hybridization matrix
//  rebuild_hyb_matrix(orbital, Delta);
//  //std::cout<<*this<<std::endl;
//  //std::cout<<clred<<*(blas_matrix*)this<<cblack<<std::endl;
//  //std::cout<<"on exit rebuild orderd: full weight: "<<full_weight()<<" permutation sign: "<<permutation_sign_<<std::endl;
//}


//Leo: my version does keep track of the number of times of swapping columns and rows;
//     however, presumably this is a time consuming action...
void hybmatrix::rebuild_ordered_hyb_matrix(int orbital, const hybfun &Delta)
{
  if(size()<2) return; //do nothing to 0*0 and 1*1 matrices

  int counter = 0; //count the number of times of swapping
  int i=0;
  hyb_map_t::iterator it_temp;

  //swap cdagger
  for(hyb_map_t::iterator it=cdagger_index_map_.begin(); it!=cdagger_index_map_.end(); ++it, ++i)
  {
     if(it->second != i)
     {
        it_temp = it; ++it_temp;
        while(it_temp->second != i) { ++it_temp; }

        //swap the values
        it_temp->second = it->second;
        it->second = i;
        counter++;  
     }
  }

  //swap c 
  i=0;
  for(hyb_map_t::iterator it=c_index_map_.begin(); it!=c_index_map_.end(); ++it, ++i)
  {
     if(it->second != i)
     {
        it_temp = it; ++it_temp;
        while(it_temp->second != i) { ++it_temp; }

        //swap the values
        it_temp->second = it->second;
        it->second = i;
        counter++; 
     }
  }
  
  if(counter%2) //get a minus sign if swapped odd number of times
      permutation_sign_ *= -1;

  //then rebuild the hybridization matrix
  rebuild_hyb_matrix(orbital, Delta);

  ////observation: when time is ordered, permutation_sign_ should be equal to time_ordering_sign_ (check!)
  if(permutation_sign_ != time_ordering_sign_)
      throw std::runtime_error("Error in hybmatrix::rebuild_ordered_hyb_matrix: permutation_sign_ != time_ordering_sign_. Abort.");
}


double hybmatrix::full_weight() const
{
  //std::cout<<clcyan<<"det: "<<determinant()<<" ps: "<<permutation_sign_<<cblack<<std::endl;
  return determinant()*permutation_sign_;
}


void hybmatrix::measure_G(std::vector<double> &G, std::vector<double> &F, const std::map<double,double> &F_prefactor, double sign, int total_size, double dissipation_weight_ratio) const
{
//  int debug_counter = 0;

  double N_div_beta=(G.size()-1)/beta_;
  static std::vector<double> cdagger_times(size()); cdagger_times.resize(size());
  static std::vector<double> c_times(size()); c_times.resize(size());
  for (hyb_map_t::const_iterator it= c_index_map_.begin(); it != c_index_map_.end(); ++it) 
  {
    c_times[it->second] = it->first;
  }
  for (hyb_map_t::const_iterator it= cdagger_index_map_.begin(); it != cdagger_index_map_.end(); ++it) 
  {
    cdagger_times[it->second] = it->first;
  }
  //we measure G(tau-tau'):=-<T c(tau) c^dagger(tau')>
  for (int i = 0; i < size(); i++) 
  {
    double f_pref=(F_prefactor.find(c_times[i]))->second;
    for (int j = 0; j < size(); j++) 
    {
      double argument = c_times[i] - cdagger_times[j];
      double bubble_sign = sign;
      if (argument < 0) 
      {
        bubble_sign *=-1.;
        argument += beta_;
      }
      int index = (int) (argument * N_div_beta + 0.5);
      double g = operator() (j, i) * bubble_sign; //Leo: original code

//      //Leo: test...
//      if(total_size > 512)
//          std::runtime_error("The expansion order is too large (>512), abort at measure_G!");
//      if(n_env_ == 2)
//          g*=combinatorial_factor[total_size];

//      double temp = size();
//      double size_ratio = temp*temp/(temp*temp + (total_size-temp)*(total_size-temp));
//      double g = operator() (j, i) * bubble_sign * size_ratio; //Leo: test!!!!!!!

      //TODO: avoid this check for n_env=1 because it's unnecessary in this case
     // if(n_env_>1)
     // {

//      if( disordered_times.find(c_times[i]) != disordered_times.end() )
//      {
//          g *= -1;
//          debug_counter++;
//      }
//      if( disordered_times.find(cdagger_times[j]) != disordered_times.end() )
//      {
//          g *= -1;
//          debug_counter++;
//      }

          //Leo: pick up a minus sign when having two adjacent c or cdagger; see my note. TODO: make it clearer 
//      if( !disordered_times.count(c_times[i]) != !disordered_times.count(cdagger_times[j]) )
//      {
//         g*=-1.;
//       //  std::cout << "A minus sign at c_time " << c_times[i] << " and cdagger_time " << cdagger_times[j] << " !" << std::endl;
//       //  debug_counter++;
//      }
     // }
      
      //g*=permutation_sign_;
      //g*=time_ordering_sign_;
      //int size_dependent_sign = ((size()%2)?-1:1); g*=size_dependent_sign; //Leo: test!
      //if(g<0)  g*=-1;  //Leo: g must be positive so that G=-g is negative. TOTALLY EXPERIMENTAL!
      //double g = operator() (j, i); //Leo: test!
      //double g = operator() (j, i) * bubble_sign * permutation_sign_; //Leo: test!
      //double g = operator() (j, i) * ((i+j)%2?-1:1) * permutation_sign_; //Leo: test!
      //double g = operator() (j, i) * bubble_sign * ((i+j)%2?-1:1) * permutation_sign_; //Leo: test!
      //double g = operator() (j, i) * bubble_sign * time_ordering_sign_; //Leo: due to operator ordering
      //double g = operator() (j, i) * bubble_sign * ((i+j)%2?-1:1); //Leo: test!
      //double g = operator() (j, i) * ((i+j)%2?-1:1); //Leo: test!

      //g*=dissipation_weight_ratio; //Leo: the dissipative environment also contributes to the local Green's function
      //NOTE:  - corresponds to -<T c(tau) c^dag(tau')>
      G[index] -= g; //changed this to-; check consistency with ALPS DMFT loop!
      F[index] -= g*f_pref;
    }
  }

//  if( (disordered_times.size()!=0) && (debug_counter != disordered_times.size()*(size()-disordered_times.size()/2)) )   
//  {
//      std::cout << debug_counter << " minus signs have been added to the estimation of Green's function. ";
//      throw std::runtime_error("BUGGY!");
//  }
}


void hybmatrix::consistency_check() const
{
  for(hyb_map_t::const_iterator it1=c_index_map_.begin(); it1!= c_index_map_.end();++it1)
  {
    for(hyb_map_t::const_iterator it2=c_index_map_.begin(); it2!= c_index_map_.end();++it2)
    {
      if(it1->first != it2->first && it1->second==it2->second)
      {
        std::cout<<clcyan<<"problem; inconsistent c map."<<cblack<<std::endl;
        std::cout<<*this;
        throw std::logic_error("...");
      }
    }
  }
  for(hyb_map_t::const_iterator it1=cdagger_index_map_.begin(); it1!= cdagger_index_map_.end();++it1)
  {
    for(hyb_map_t::const_iterator it2=cdagger_index_map_.begin(); it2!= cdagger_index_map_.end();++it2)
    {
      if(it1->first != it2->first && it1->second==it2->second)
      {
        std::cout<<clcyan<<"problem; inconsistent c map."<<cblack<<std::endl;
        std::cout<<*this;
        throw std::logic_error("...");
      }
    }
  }
}


void hybmatrix::measure_Gl(std::vector<double> &Gl, std::vector<double> &Fl , const std::map<double,double> &F_prefactor, double sign) const{
  static std::vector<double> cdagger_times(size()); cdagger_times.resize(size());
  static std::vector<double> c_times(size()); c_times.resize(size());
  static std::vector<std::complex<double> > cdagger_exp(size()); cdagger_exp.resize(size());
  static std::vector<std::complex<double> > c_exp(size()); c_exp.resize(size());
  int N_l=Gl.size();

  //create map of creator and annihilator times
  for (hyb_map_t::const_iterator it= c_index_map_.begin(); it != c_index_map_.end(); ++it) {
    c_times[it->second] = it->first;
  }
  for (hyb_map_t::const_iterator it= cdagger_index_map_.begin(); it != cdagger_index_map_.end(); ++it) {
    cdagger_times[it->second] = it->first;
  }

  //measures the Legendre coefficients of G(tau-tau'):=-<T c(tau) c^dagger(tau')>
  for (int i = 0; i < size(); i++) {
    double f_pref=(F_prefactor.find(c_times[i]))->second;
    for (int j = 0; j < size(); j++) {
      double M_ji = operator() (j, i) * sign;
      double argument = c_times[i] - cdagger_times[j];
      double bubble_sign = sign;
      if (argument < 0) {
        bubble_sign *=-1.;
        argument += beta_;
      }
      double x=2.0*argument/beta_-1.0;
      double pl_2=1; double pl_1=x; double legendre_p;
      for(int l=0; l<N_l; l++){
        if(l==0) legendre_p=1;
        else if(l==1) legendre_p=x;
        else{ 
          legendre_p=((2*l-1)*x*pl_1-(l-1)*pl_2)/static_cast<double>(l);//l
          pl_2=pl_1; //l-2
          pl_1=legendre_p; //l-1
        }
        double gl = M_ji*legendre_p*bubble_sign;
        double fl = gl*f_pref;
        Gl[l]-=gl/beta_;
        Fl[l]-=fl/beta_;
      }
    }
  }
}


//Leo: When time is ordered, permutation_sign_ should be equal to time_ordering_sign_ (check!)
//     During this check, the operator times at which the operators are not ordered (ex. sequential cdaggers) will
//     be recorded in a set, called disordered_times, for the later usage of measure_G.
//void hybmatrix::time_ordering_sign_check(int &time_ordering_sign_, std::set<double> &disordered_times)
void hybmatrix::time_ordering_sign_check()
{
    //if(n_env_==1) return; //single lead should not have time ordering sign, so do nothing (TODO: check!)

    if(size()==0) 
    {
       time_ordering_sign_ = 1;
       return; 
    }

    int time_ordering_sign_counter = 0;
    int head = c_cdagger_map_.begin()->second;  //determine the head of the sequence should be c or cdagger
    if((head != 0) && (head != 1))
       throw std::runtime_error("Error in hybmatrix::time_ordering_sign_check: meaningless head. Abort.");  
    disordered_times.clear();  //clear the set because we don't know the disordered times at this point yet.

    //construct a sequence {010101...} or {101010...} to be compared with the actual sequence
    std::vector<int> operator_sequence(2*size(), 0); 
    for(int i=(head?0:1); i<operator_sequence.size(); i=i+2)   operator_sequence[i]++;
 
    //comparism
    int i=0;
    for(hyb_map_t::const_iterator it=c_cdagger_map_.begin(); it!=c_cdagger_map_.end(); ++it)
    {
        if(it->second != operator_sequence[i])
        {
            time_ordering_sign_counter++;
            disordered_times.insert(it->first);  //record the time 
        }
        i++; 
    }
    
//    //version1 
//    if(time_ordering_sign_counter!=0 && ((time_ordering_sign_counter/2)%2)) 
//        time_ordering_sign_=-1; 
//    else
//        time_ordering_sign_=1;
   
    //version2 
    //This is the sign coming from (-1)^n, n being the matrix size
    int size_dependent_sign = ((size()%2)?-1:1); 
    if(time_ordering_sign_counter!=0 && ((time_ordering_sign_counter/2)%2)) 
        time_ordering_sign_=-1*(head?1:size_dependent_sign); //size_dependent_sign contributes when the sequence is {010101...}
    else
        time_ordering_sign_=1*(head?1:size_dependent_sign);
 
    //version1 
//    if(time_ordering_sign_counter!=0 && ((time_ordering_sign_counter/2)%2)) 
//        time_ordering_sign_=-1*(head?-1:1); 
//    else
//        time_ordering_sign_=1*(head?-1:1); 
}



//compute the hybridization weight change when the worm is creeping
double hybmatrix::hyb_weight_change_worm_creep(double new_worm_head, double old_worm_head, double worm_tail, int orbital, const hybfun &Delta)
{
  std::size_t row_index = cdagger_index_map_[new_worm_head];
  std::size_t last = size()-1;
  std::size_t column_index = last;

  //find the iterators of the last row
  std::map<double, size_t>::iterator row_it = cdagger_index_map_.begin();
  for( ; row_it != cdagger_index_map_.end(); ++row_it) { if(row_it->second == last) break; }
  
  //create Delta^-1 matrix of size (k-1)*(k-1); this part is adapted from remove_segment
  //first swap the selected row with the last row
  if(row_index != last)
  {
    swap_row(row_index, last);
    permutation_sign_ *= -1.;
    std::swap( cdagger_index_map_[new_worm_head], cdagger_index_map_[row_it->first] );
    assert(cdagger_index_map_[new_worm_head] == last);
    row_index = last;
  }

  //this is the inverse weight of the old configuration
  double S_tilde_old = operator()(row_index, column_index);

  //find the iterators of the last column
  std::map<double, size_t>::iterator column_it = c_index_map_.begin();
  for( ; column_it != c_index_map_.end(); ++column_it) { if(column_it->second == last) break; }

  //then create a copy of the matrix
  inv_Delta_small = new blas_matrix(*this);
  
  fortran_int_t sm1 = last;
  if(sm1>0)
  {
    double alpha = -1./inv_Delta_small->operator()(last, last);
    fortran_int_t inc = 1;
    fortran_int_t ms = inv_Delta_small->memory_size();
    FORTRAN_ID(dger)(&sm1, &sm1, &alpha, &(inv_Delta_small->operator()(last,0)), &inc, &(inv_Delta_small->operator()(0,last)), &ms, &(inv_Delta_small->operator()(0,0)), &ms);
  }

  //now we have the smaller Delta^-1 matrix.
  inv_Delta_small->resize(last);

  //now compute the weight of the new configuration
  Q.resize(last);
  R.resize(last);
  PinvQ.resize(last);

  //column  Delta_i,last
  for(hyb_map_t::const_iterator it = c_index_map_.begin(); it!=c_index_map_.end(); ++it)
  {
     if(it->second != column_index)
        Q[it->second] = Delta.interpolate(it->first - old_worm_head, orbital); //this is the new column Q
  }
  
  //row Delta_last,i
  for(hyb_map_t::const_iterator it = cdagger_index_map_.begin(); it != cdagger_index_map_.end(); ++it)
  {
     if(it->second != row_index)
        R[it->second] = Delta.interpolate(column_it->first - it->first, orbital);  //this is the new row R
  }

  S=Delta.interpolate(column_it->first - old_worm_head, orbital);  //this is the entry S
  
  //S_tilde_inv is the new weight
  S_tilde_inv = S;
  fortran_int_t s = last; //note the size! 
  if(s>0)
  {
    inv_Delta_small->right_multiply(Q, PinvQ); //dgemv
    fortran_int_t inc = 1;
    S_tilde_inv -= FORTRAN_ID(ddot)(&s, &(R[0]), &inc, &(PinvQ[0]),&inc);
  }

  ////pick up a wrapping sign if c^dagger passes through odd number of c
  //std::map<double, size_t>::const_iterator it_low = cdagger_index_map_.lower_bound((old_worm_head>new_worm_head?new_worm_head:old_worm_head));
  //std::map<double, size_t>::const_iterator it_high = cdagger_index_map_.upper_bound((old_worm_head<new_worm_head?new_worm_head:old_worm_head));
  //--it_high;
  //int counter = std::distance(it_low, it_high);
  //weight_ratio_ = S_tilde_old * S_tilde_inv;
  //if( (counter % 2) && !(size() % 2) )
  //   weight_ratio_ *= -1.;

  ////pick up a wrapping sign if c^dagger passes through odd number of c
  //std::map<double, size_t>::const_iterator it_low = c_cdagger_map_.lower_bound((old_worm_head>new_worm_head?new_worm_head:old_worm_head));
  //std::map<double, size_t>::const_iterator it_high = c_cdagger_map_.upper_bound((old_worm_head<new_worm_head?new_worm_head:old_worm_head));
  //int counter = 0;
  //for(std::map<double, size_t>::const_iterator it = it_low; it != it_high; ++it) { if(it->second == 1) counter++; }
  //weight_ratio_ = (counter%2 ? -1. : 1.) * S_tilde_old * S_tilde_inv;

  ////always pick up a minus sign! //TODO: test!
  //weight_ratio_ = S_tilde_old * S_tilde_inv;
  //if( (column_it->first > old_worm_head && column_it->first < new_worm_head) ||
  //    (column_it->first < old_worm_head && column_it->first > new_worm_head) )
  //   weight_ratio_ *= -1.;
  
  //weight_ratio_ = S_tilde_old * S_tilde_inv;
  //if( column_it->first < old_worm_head )
  //   weight_ratio_ *= -1.;
  
  weight_ratio_ = S_tilde_old * S_tilde_inv;
  if( (worm_tail > old_worm_head && worm_tail > new_worm_head) || (worm_tail < old_worm_head && worm_tail < new_worm_head) )
  {
     weight_ratio_ *= -1.;

     if(n_env_ != 1)
     {
        std::map<double, size_t>::const_iterator it_low  = c_cdagger_map_.upper_bound((new_worm_head>worm_tail?worm_tail:new_worm_head));
        std::map<double, size_t>::const_iterator it_high = c_cdagger_map_.lower_bound((new_worm_head<worm_tail?worm_tail:new_worm_head));
        int counter = 0;
        for(std::map<double, size_t>::const_iterator it = it_low; it != it_high; ++it) 
        { 
           if(it->second == 1) 
              counter++; 
           else 
              counter--;
        }
	if( (new_worm_head < worm_tail && new_worm_head < old_worm_head && old_worm_head < worm_tail) ||
	    (worm_tail < new_worm_head && worm_tail < old_worm_head && old_worm_head < new_worm_head) )
	   counter++;

        if(counter != 0)      
           weight_ratio_ *= -1.;
     }
  }

  return weight_ratio_;
}


//actually swap the hybridization line to the old worm head
void hybmatrix::worm_creep(double new_worm_head, double old_worm_head, double worm_tail, int orbital, const hybfun &Delta)
{
  int last = size()-1;
  
  //last element
  operator()(last, last) = 1./S_tilde_inv;
 
  //TODO: see if inv_Delta_small can be reused here
  fortran_int_t sm1 = last;
  if(sm1>0)
  { //this is exactly the content of the loops above, in dger/dgemv blas calls.
   char trans='T', notrans='N';
   double alpha = -1./S_tilde_inv, beta=0.;
   fortran_int_t inc=1;
   fortran_int_t ms = memory_size(); //this line is safe is because inv_Delta_small has the same memory size as "this" during execution
   assert(inv_Delta_small->memory_size() == ms);
   FORTRAN_ID(dgemv)(&  trans, &sm1, &sm1, &alpha, &(inv_Delta_small->operator()(0,0)), &ms, &(Q[0]), &inc, &beta, &(operator()(0,last)), &ms);
   FORTRAN_ID(dgemv)(&notrans, &sm1, &sm1, &alpha, &(inv_Delta_small->operator()(0,0)), &ms, &(R[0]), &inc, &beta, &(operator()(last,0)), &inc);

   alpha=S_tilde_inv;
   FORTRAN_ID(dger)(&sm1, &sm1, &alpha, &(operator()(last,0)), &inc, &(operator()(0,last)), &ms, &(inv_Delta_small->operator()(0,0)), &ms);
   for(fortran_int_t i=0; i<last; ++i)
   { //for each row: copy the entire row (except for the last element).
     memcpy(&(operator()(0,0))+i*memory_size(), &(inv_Delta_small->operator()(0,0))+i*inv_Delta_small->memory_size(), sizeof(double)*(size()-1));
   }
  }
  
  ////TODO: test
  ////pick up a wrapping sign if c^dagger passes through odd number of c
  //std::map<double, size_t>::const_iterator it_low = cdagger_index_map_.lower_bound((old_worm_head>new_worm_head?new_worm_head:old_worm_head));
  //std::map<double, size_t>::const_iterator it_high = cdagger_index_map_.upper_bound((old_worm_head<new_worm_head?new_worm_head:old_worm_head));
  //--it_high;
  //int counter = std::distance(it_low, it_high);
  //if( (counter % 2) && !(size() % 2) )
  //   permutation_sign_ *= -1.;

  if( (worm_tail > old_worm_head && worm_tail > new_worm_head) || (worm_tail < old_worm_head && worm_tail < new_worm_head) )
  {
     permutation_sign_ *= -1.;

     if(n_env_ != 1)
     {
        std::map<double, size_t>::const_iterator it_low  = c_cdagger_map_.upper_bound((new_worm_head>worm_tail?worm_tail:new_worm_head));
        std::map<double, size_t>::const_iterator it_high = c_cdagger_map_.upper_bound((new_worm_head<worm_tail?worm_tail:new_worm_head));
        int counter = 0;
        for(std::map<double, size_t>::const_iterator it = it_low; it != it_high; ++it) 
        { 
           if(it->second == 1) 
              counter++; 
           else 
              counter--;
        }
	if( (new_worm_head < worm_tail && new_worm_head < old_worm_head && old_worm_head < worm_tail) ||
	    (worm_tail < new_worm_head && worm_tail < old_worm_head && old_worm_head < new_worm_head) )
	   counter++;

        if(counter != 0)      
           permutation_sign_ *= -1.;
     }
  }

  ////pick up a wrapping sign if c^dagger passes through odd number of c
  //std::map<double, size_t>::const_iterator it_low = c_cdagger_map_.lower_bound((old_worm_head>new_worm_head?new_worm_head:old_worm_head));
  //std::map<double, size_t>::const_iterator it_high = c_cdagger_map_.upper_bound((old_worm_head<new_worm_head?new_worm_head:old_worm_head));
  //int counter = 0;
  //for(std::map<double, size_t>::const_iterator it = it_low; it != it_high; ++it) { if(it->second == 1) counter++; }
  //if(counter % 2)
  //   permutation_sign_ *= -1.;
  

  //if( column_it->first < old_worm_head )
  //   permutation_sign_ *= -1.;

  //if( (column_it->first > old_worm_head && column_it->first < new_worm_head) ||
  //    (column_it->first < old_worm_head && column_it->first > new_worm_head) )
  //   permutation_sign_ *= -1.;
  
  //if( worm_tail < old_worm_head )
  //   permutation_sign_ *= -1.;

  //erase old times:
  cdagger_index_map_.erase(new_worm_head);
  c_cdagger_map_    .erase(new_worm_head);
  
  // add the new segment times:
  cdagger_index_map_.insert(std::make_pair(old_worm_head, last));
  c_cdagger_map_    .insert(std::make_pair(old_worm_head, 1)); //Leo: 1 means c_dagger
  
  //Leo: keep track of the time ordering sign due to operator time ordering
  //TODO: move this to Green function measurements because this sign wouldn't affect other local measurements.
  //time_ordering_sign_check(time_ordering_sign_, disordered_times);
  time_ordering_sign_check();
}
