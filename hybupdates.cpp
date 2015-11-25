/****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2012 by Emanuel Gull <egull@umich.edu>,
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
#include <iomanip>
#include"hyb.hpp"
#include"hyblocal.hpp"

//this is the heart of the Monte Carlo procedure: we have the following updates:
//1: change the zero order state, swap an empty orbital versus a filled one (8% or 10%)
//2: shift an existing segment start- or end-point (not implemented)
//3: insert or remove a new segment (40% or 50%)
//4: insert or remove an anti-segment (30% or 40%)
//5: perform a segment flip between different orbtials (20%)
//6: perform a global exchange of two orbitals (2%)
//see our review for details of these updates

void hybridization::update()
{
  //one sweep is composed of N_MEAS Monte Carlo updates and one measurement (the latter only if thermalized)
  sweeps++;
  
  double rates[2] = {(spin_flip)?0.5:0.6,(spin_flip)?0.8:1.0};

  for(std::size_t i=0;i<N_meas;++i)
  {
    double update_type=random();
    //Leo: we don't know the update type at this point, so set color_updated = false
    color_updated = false;

    //Leo: Disable global_flip_update because it's too painful to modify it;
    //     besides, this update is removed in the GitHub version (don't know why).
//    if (update_type < 0.02 && global_flip)
//    {
//      global_flip_update();
//    } 
//    else if (update_type<0.1) 
    if (update_type<0.1)
    {
      change_zero_order_state_update();
//    } else if (update_type < 0) {
//      shift_segment_update();
    } 
    else if (update_type < rates[0]) 
    {
      insert_remove_segment_update();
    } 
    else if (update_type < rates[1]) 
    {
      insert_remove_antisegment_update();
    } 
    else 
    {
      insert_remove_spin_flip_update();
    }

    if(is_thermalized())
    {
      measure_order();
      if(MEASURE_time)
      {
        local_config.get_F_prefactor(F_prefactor);//compute segment overlaps in local config
        measure_G(F_prefactor);
      }
    }

    /* Leo Fang: for test purpose, print out the segment map */
    if(sweeps<=debug_number) 
    { 
	std::cout << "|---------------------------------------------------------------------------------|" << std::endl;
	std::cout << "At " << i+(sweeps-1)*N_meas+1 << "-th update:" << std::endl; local_config.print_segments(); 
	std::cout << "|---------------------------------------------------------------------------------|" << std::endl;
        std::cout << std::endl;
    }

  }//N_meas

  if(VERBOSE && sweeps%output_period==0 && crank==0) 
  {
    //  if(VERBOSE && crank==0 && boost::chrono::steady_clock::now() - lasttime > delay) {
    //    lasttime = boost::chrono::steady_clock::now();
    int tot_acc=0,cur_prec = std::cout.precision();
    for (int i=0;i<nacc.size();i++) tot_acc += nacc[i];
    std::cout << std::endl << "|------ Simulation details (master only) after " << sweeps << " sweeps ------|" << std::endl;
    std::cout << "  Total acceptance rate = " << std::setprecision(2) << std::fixed;
    std::cout << (((double)tot_acc)/(sweeps*N_meas))*100 << "%" << std::endl;
    std::cout << "  Individual acceptance rates for update " << std::endl;
    for (int i=0;i<nacc.size();i++) 
    {
      std::cout << "     " << update_type[i] << " = ";
      std::cout << std::setprecision(2) << std::fixed << (((double)nacc[i])/(sweeps*N_meas))*100 << "%";
      std::cout << " (proposal rate = ";
      std::cout << std::setprecision(2) << std::fixed << (((double)nprop[i])/(sweeps*N_meas))*100 << "%)" << std::endl;
    }
    std::cout << "|-----------------------------------------------------------------|" << std::endl;
    std::cout.unsetf(std::ios_base::fixed);
    std::cout.precision(cur_prec);
  }

  //Leo: check the size of colored matrices for each orbital
  for(int i=0; i<n_orbitals; i++)
  { 
    std::stringstream temp_stream; // stringstream used for the conversion
    temp_stream << i;              // add the value of i to the characters in the stream
    if( local_config.order(i)==0 ) continue; //Leo: no segment exists, and so does color, in 0-th order
    if( hyb_config.total_color_matrix_size(i) != local_config.order(i) )
        std::runtime_error("The total size of colored matrices for orbital " + temp_stream.str() + " is incorrect!");
  }

  //if(sweeps%1000==0) {std::cout <<  sweeps << " sweeps has been done...\n";}
}


void hybridization::change_zero_order_state_update()
{
  //choose the orbital in which we do the update
  nprop[0]++;
  int orbital=(int)(random()*n_orbitals);
  
  //changing the zero order state only makes sense if we are at zero order.
  // Leo: Todo: re-organize the code here. A parenthesis can be added anywhere since the outcome is unaffected
  if( (!local_config.order(orbital)) == 0) return;
  
  //propose to change orbital from occuppied to unoccuppied.
  if(local_config.zero_order_orbital_occupied(orbital))
  {
    double local_weight_change=1./local_config.local_weight_change(segment(0,beta), orbital, false);
    if(std::abs(local_weight_change)>random())
    {
      nacc[0]++;
      local_config.set_zero_order_orbital_occupied(orbital, false);
      if(local_weight_change<0)
        sign*=-1.;
    }
  }
  //propose to change from unoccuppied to occuppied.
  else
  {
    double local_weight_change=local_config.local_weight_change(segment(0,beta), orbital, false);
    //std::cout<<cmagenta<<"local weight change is: "<<local_weight_change<<cblack<<std::endl;
    if(std::abs(local_weight_change)>random())
    {
      nacc[0]++;
      local_config.set_zero_order_orbital_occupied(orbital, true);
      if(local_weight_change<0)
        sign*=-1.;
    }
  }
}


//// Perform a complete swap of segments between two orbitals
//// THIS IS TOTALLY EXPERIMENTAL
//// A bare-bone structure for testing
//void hybridization::global_flip_update()
//{
//  int orbital1=0;
//  int orbital2=1;
//
//  // These are the actual orders for each of the orbitals
//  int k1 = local_config.order(orbital1),k2=local_config.order(orbital2);
//  // At present we do nothing if one is empty (can be relaxed, I think)
//  if (k1==0 || k2==0) return;
//  std::cerr << "On entry:" << std::endl;
//  hyb_config.dump();
//  std::vector<segment> seg1(k1),seg2(k2);
//  double total_hyb_weight_change = 1.0,d_e=0.0;
//  for (int k=0;k<k1;k++) {
//    seg1[k] = local_config.get_segment(k,orbital1);
//    d_e -= local_config.local_energy(seg1[k],orbital1);//,true);
//  }
//  for (int k=0;k<k2;k++) {
//    seg2[k] = local_config.get_segment(k,orbital2);
//    d_e -= local_config.local_energy(seg2[k],orbital2);//,true);
//  }
//  
//  for (int k=0;k<k1;k++) {
//    total_hyb_weight_change /= hyb_config.hyb_weight_change_remove(seg1[k],orbital1);
//    hyb_config.remove_segment(seg1[k],orbital1);
//    local_config.remove_segment(seg1[k],orbital1);
//  }
//  for (int k=0;k<k2;k++) {
//    total_hyb_weight_change *= hyb_config.hyb_weight_change_insert(seg2[k],orbital1);
//    hyb_config.insert_segment(seg2[k],orbital1);
//    local_config.remove_segment(seg2[k],orbital2);
//  }
//  for (int k=0;k<k2;k++) {
//    total_hyb_weight_change /= hyb_config.hyb_weight_change_remove(seg2[k],orbital2);
//    hyb_config.remove_segment(seg2[k],orbital2);
//    local_config.insert_segment(seg2[k],orbital1);
//  }
//  for (int k=0;k<k1;k++) {
//    total_hyb_weight_change *= hyb_config.hyb_weight_change_insert(seg1[k],orbital2);
//    hyb_config.insert_segment(seg1[k],orbital2);
//    local_config.insert_segment(seg1[k],orbital2);
//  }
//  for (int k=0;k<k2;k++) d_e += local_config.local_energy(seg2[k],orbital1);
//  for (int k=0;k<k1;k++) d_e += local_config.local_energy(seg1[k],orbital2);
//
//  // This is the total weight change due to the swap. If all orbitals are
//  // equivalent this should be one.
//  double weight_change = exp(d_e)*total_hyb_weight_change;
//  // Since the total expansion order does not change, there should be no
//  // permutation factor appearing here
//  std::cerr << "In between: de = " << d_e << ", total hyb weight change = "<< total_hyb_weight_change << std::endl;
//  // This is the proposed weight. Should be the numbers as before, but for the orbitals exchanged
//  hyb_config.dump();
////  hyb_config.rebuild();
//
//  // Assume it was rejected. We have to restore the old configuration
//  for (int k=0;k<k2;k++) {
//    hyb_config.remove_segment(seg2[k],orbital1);
//  }
//  for (int k=0;k<k1;k++) {
//    hyb_config.remove_segment(seg1[k],orbital2);
//  }
//  for (int k=0;k<k1;k++) {
//    hyb_config.insert_segment(seg1[k],orbital1);
//  }
//  for (int k=0;k<k2;k++) {
//    hyb_config.insert_segment(seg2[k],orbital2);
//  }
////  hyb_config.rebuild();
//  std::cerr << "On exit:" << std::endl;
//  // This should be again the initial configuration
//  hyb_config.dump();
//  exit(-1);
//// Done.
//}
//


// Perform a complete swap of segments between two orbitals
// THIS IS TOTALLY EXPERIMENTAL
// Not (yet) optimized
    //Leo: Disable global_flip_update because it's too painful to modify it;
    //     besides, this update is removed in the GitHub version (don't know why).
//void hybridization::global_flip_update()
//{
//  nprop[6]++;
//  // Pick orbital 1
//  int orbital1=(int)(random()*n_orbitals);
//  // Pick orbital 2 from the rest
//  int orbital2=(int)(random()*(n_orbitals-1));
//  orbital2 = (orbital2<orbital1)?orbital2:1+orbital2;
//  
//  // These are the actual orders for each of the orbitals
//  int k1 = local_config.order(orbital1),k2=local_config.order(orbital2);
//  // At present we do nothing if one is empty (can be relaxed, I think)
//  if (k1==0 || k2==0) {
////    std::cerr << "k1 = " << k1 << " and k2 = " << k2 << std::endl;
//    return;
//  }
//  std::vector<int> orbitals(2);
//  orbitals[0] = orbital1;
//  orbitals[1] = orbital2;
//  //std::cerr << "On entry:" << std::endl;
//  //hyb_config.dump();
////  std::cerr << "Orbitals picked are " << orbital1 << " and " << orbital2 << std::endl;
//  // We need to store the segments for the swap
//  // This is quite clumsy. However, I did not succeed in generating an
//  // intermediate copy of hyb_config. I tried to implement a copy constructor,
//  // but this clashed in a seg-fault when trying to delete it
//  std::vector<segment> seg1(k1),seg2(k2);
//  // I brutally compute the change in hybridization configuration by simply
//  // deleting successivley all segments from orbital 1 and then inserting the
//  // ones from orbital 2; likewise for orbital 2.
//  // The local weight change I compute from the local energy, taking into account
//  // the mu-part only (this is the meaning of the bool in local_energy call;
//  // should be fine for Coulomb only as the segments do not really change, but
//  // may cause trouble when Hund is present or for dynamic Coulomb.
//  double total_hyb_weight_change = 1.0,d_e=0.0;
//  for (int k=0;k<k1;k++) {
//    seg1[k] = local_config.get_segment(k,orbital1);
//    d_e -= local_config.local_energy(seg1[k],orbital1);//,true);
//  }
//  for (int k=0;k<k2;k++) {
//    seg2[k] = local_config.get_segment(k,orbital2);
//    d_e -= local_config.local_energy(seg2[k],orbital2);//,true);
//  }
//  for (int k=0;k<k1;k++) {
//    total_hyb_weight_change /= hyb_config.hyb_weight_change_remove(seg1[k],orbital1);
//    hyb_config.remove_segment(seg1[k],orbital1);
//    local_config.remove_segment(seg1[k],orbital1);
//  }
//  for (int k=0;k<k2;k++) {
//    total_hyb_weight_change *= hyb_config.hyb_weight_change_insert(seg2[k],orbital1);
//    hyb_config.insert_segment(seg2[k],orbital1);
//    local_config.remove_segment(seg2[k],orbital2);
//  }
//  for (int k=0;k<k2;k++) {
//    total_hyb_weight_change /= hyb_config.hyb_weight_change_remove(seg2[k],orbital2);
//    hyb_config.remove_segment(seg2[k],orbital2);
//    local_config.insert_segment(seg2[k],orbital1);
//  }
//  for (int k=0;k<k1;k++) {
//    total_hyb_weight_change *= hyb_config.hyb_weight_change_insert(seg1[k],orbital2);
//    hyb_config.insert_segment(seg1[k],orbital2);
//    local_config.insert_segment(seg1[k],orbital2);
//  }
//  for (int k=0;k<k2;k++) d_e += local_config.local_energy(seg2[k],orbital1);
//  for (int k=0;k<k1;k++) d_e += local_config.local_energy(seg1[k],orbital2);
//
//  // This is the total weight change due to the swap. If all orbitals are
//  // equivalent (and the expansion order is the same) this should be one.
//  double weight_change = exp(d_e)*total_hyb_weight_change;
//  // Since the total expansion order does not change, there should be no
//  // permutation factor appearing here
//  // std::cerr << "In between: de = " << d_e << ", total hyb weight change = "<< total_hyb_weight_change << ", MC weight change = " << weight_change << std::endl;
////  hyb_config.dump();
//  hyb_config.rebuild(orbitals);
//
//  // MC move
//  if(std::abs(weight_change)>random()){
//    nacc[6]++;
//    if(weight_change < 0) sign*=-1.;
//    // Accepted. Since we already changed the configuration, we have nothing to do
//  } else {
//    // Rejected. We have to restore the old configuration
//    for (int k=0;k<k2;k++) {
//      hyb_config.remove_segment(seg2[k],orbital1);
//      local_config.remove_segment(seg2[k],orbital1);
//    }
//    for (int k=0;k<k1;k++) {
//      hyb_config.remove_segment(seg1[k],orbital2);
//      local_config.remove_segment(seg1[k],orbital2);
//    }
//    for (int k=0;k<k1;k++) {
//      hyb_config.insert_segment(seg1[k],orbital1);
//      local_config.insert_segment(seg1[k],orbital1);
//    }
//    for (int k=0;k<k2;k++) {
//      hyb_config.insert_segment(seg2[k],orbital2);
//      local_config.insert_segment(seg2[k],orbital2);
//    }
////    hyb_config.dump();
//    hyb_config.rebuild(orbitals);
//  }
//  local_config.check_consistency();
//  //std::cerr << "On exit:" << std::endl;
//  //Full weight  = " << hyb_config.full_weight() << std::endl;
//  //hyb_config.dump();
////  exit(-1);
//// Done.
//}


void hybridization::shift_segment_update(){
  ///TODO: implement this update!
}


void hybridization::insert_remove_segment_update()
{
  //choose the orbital in which we do the update
  int orbital=(int)(random()*n_orbitals);
   
  //Leo: choose the color in which we do the update
  //Now only two colors (red/1 and blue/0) are considered, but it can be easily changed
  std::size_t color_temp;
  if(n_env == 1) { color_temp = 0; } // Only one color
  else if(n_env==2)
  { 
     if(random()<0.5){ color_temp = 0; } // blue = 0 = R
     else            { color_temp = 1; } // red = 1 = L  
  }
  else 
  {  //Leo: because we set up the sanity_check, this line is redundant here for the moment... 
     throw std::runtime_error("The input N_ENV>2 is currently not supported.");
  }

  if(random()<0.5){ insert_segment_update(orbital, color_temp); }
  else            { remove_segment_update(orbital, color_temp); }
}


void hybridization::insert_remove_antisegment_update()
{
  //choose the orbital in which we do the update
  int orbital=(int)(random()*n_orbitals);

  //Leo: choose the color in which we do the update
  //Now only two colors (red/1 and blue/0) are considered, but it can be easily changed
  std::size_t color_temp;
  if(n_env == 1) { color_temp = 0; } // Only one color
  else if(n_env==2)
  { 
     if(random()<0.5){ color_temp = 0; } // blue = 0 = R
     else            { color_temp = 1; } // red = 1 = L  
  }
  else 
  {  //Leo: because we set up the sanity_check, this line is redundant here for the moment... 
     throw std::runtime_error("The input N_ENV>2 is currently not supported.");
  }

  if(random()<0.5){ insert_antisegment_update(orbital, color_temp); }
  else            { remove_antisegment_update(orbital, color_temp); }
}


void hybridization::insert_remove_spin_flip_update()
{
  //choose the orbital in which we do the update
  int orbital=(int)(random()*n_orbitals);

  //Leo: choose the color in which we do the update
  //Now only two colors (red/1 and blue/0) are considered, but it can be easily changed
  std::size_t color_temp;
  if(n_env == 1) { color_temp = 0; } // Only one color
  else if(n_env==2)
  {  
     if(random()<0.5){ color_temp = 0; } // blue = 0 = R
     else            { color_temp = 1; } // red = 1 = L  
  }
  else 
  {  //Leo: because we set up the sanity_check, this line is redundant here for the moment... 
     throw std::runtime_error("The input N_ENV>2 is currently not supported.");
  }

  spin_flip_update(orbital, color_temp);
}


void hybridization::insert_segment_update(int orbital, std::size_t color_temp)
{
  nprop[1]++;
  //std::cout<<clred<<"starting insertion update."<<cblack<<std::endl;
  if(local_config.order(orbital)==0 && local_config.zero_order_orbital_occupied(orbital)) return; //can't insert segment, orbital is fully occuppied.
  double t_start=random()*beta; //start time of a segment
  if(local_config.exists(t_start)){ /*std::cerr<<"rare event, duplicate: "<<t_start<<std::endl; */ return;} //time already exists.

  double t_next_segment_start=local_config.find_next_segment_start_distance(t_start,orbital);
  double t_next_segment_end=local_config.find_next_segment_end_distance(t_start,orbital);
  //std::cout<<"============================================"<<cblack<<std::endl;
  //std::cout<<clblue<<"orbital: "<<orbital<<" time is: "<<t_start<<" segment start distance: "<<t_next_segment_start<<" end distance: "<<t_next_segment_end<<cblack<<std::endl;
  //std::cout<<*this<<std::endl;
  //std::cout<<"============================================"<<cblack<<std::endl;
  if(t_next_segment_end < t_next_segment_start) return; //we're trying to create a segment on top of another segment. abort.
  
  //draw an end time
  double t_len = random()*t_next_segment_start;
  double t_end=t_start+t_len;
  if(t_end >= beta) t_end-=beta;
  if(local_config.exists(t_end)){ /*std::cerr<<"rare event, duplicate: "<<t_end<<std::endl; */return;} //time already exists.
  if(t_end<=t_start || t_end<=0.0){ /*std::cerr<<"rare event, zero length segment: "<<t_start<<" "<<t_end<<std::endl; */return;} //time already exists.
  
  //Leo: paint the color on the segment
  //compute local weight of the new segment with t_start and t_end
  segment new_segment(t_start, t_end, color_temp, color_temp);
  double local_weight_change=local_config.local_weight_change(new_segment, orbital, false);
  
  //compute hybridization weight change
  double hybridization_weight_change=hyb_config.hyb_weight_change_insert(new_segment, orbital, color_temp);
  
  //compute the proposal probability ratio
  double permutation_factor=beta*t_next_segment_start/(local_config.order(orbital)+1);
  
  //perform metropolis
  double weight_change=local_weight_change*hybridization_weight_change*permutation_factor;
  
  /*std::cout<<" new segment: "<<new_segment<<std::endl;
   std::cout<<clred<<"weight change: "<<weight_change<<" l: "<<local_weight_change<<" h: "<<hybridization_weight_change<<" p: "<<permutation_factor<<cblack<<std::endl;*/
  
  if(std::abs(weight_change)>random())
  {
    nacc[1]++;
    if(weight_change < 0) sign*=-1.;
    local_config.insert_segment(new_segment, orbital);
    hyb_config.insert_segment(new_segment, orbital, color_temp);
    //Leo: record the updated color and set color_updated to true
    color = color_temp;
    color_updated = true;
  }
}


void hybridization::remove_segment_update(int orbital, std::size_t color_temp)
{
  nprop[2]++;
  //std::cout<<clblue<<"starting removal update."<<cblack<<std::endl;
  int k=local_config.order(orbital);
  
  if(k==0) return; //no point, this is an empty orbital //Leo: it could also be occupied, but it doesn't matter here
  
  int segment_nr=(int)(random()*k);
  
  segment segment_to_remove=local_config.get_segment(segment_nr, orbital);

  //Leo: for debug purpose, print the colors of the chosen segment
  std::cout << "remove_segment_update: c_start=" <<segment_to_remove.c_start_<<", "\
            << "c_end=" << segment_to_remove.c_end_ << std::endl;

  //Leo: check if the colors of both ends and the randomly picked color are all the same
  if(color_temp == segment_to_remove.c_start_ && color_temp == segment_to_remove.c_end_) ; //do nothing
  else {return;} //cannot remove because of the different colors
  
  double local_weight_change=1./local_config.local_weight_change(segment_to_remove, orbital, false);
  
  //compute hybridization weight change
  double hybridization_weight_change=1.0/hyb_config.hyb_weight_change_remove(segment_to_remove, orbital, color_temp);
  
  //compute the proposal probability ratio
  double t_next_segment_start=local_config.find_next_segment_start_distance(segment_to_remove.t_start_,orbital);
  double permutation_factor=local_config.order(orbital)/(beta*t_next_segment_start);
  
  //perform metropolis
  double weight_change=local_weight_change*hybridization_weight_change*permutation_factor;
  
  /*if(segment_nr==k-1){
   std::cout<<" segment_to_remove: "<<segment_to_remove<<std::endl;
   std::cout<<"t_next_segment_start: "<<t_next_segment_start<<std::endl;
   std::cout<<clblue<<"weight change: "<<weight_change<<" l: "<<local_weight_change<<" h: "<<hybridization_weight_change<<" p: "<<permutation_factor<<cblack<<std::endl;
   }*/
  if(std::abs(weight_change)>random())
  {
    nacc[2]++;
    if(weight_change < 0) sign*=-1.;
//      double fwo = full_weight();
    local_config.remove_segment(segment_to_remove, orbital);
    hyb_config.remove_segment(segment_to_remove, orbital, color_temp);
    //Leo: record the updated color and set color_updated to true
    color = color_temp;
    color_updated = true;
//      double fwa = full_weight();
//      std::cout << clgreen<<"weight change removal: "<<fwa<<" control: "<<fwo*std::abs(weight_change)<<std::endl;
  }
}


void hybridization::insert_antisegment_update(int orbital, std::size_t color_temp)
{
  nprop[3]++;
  if(local_config.order(orbital)==0 && !local_config.zero_order_orbital_occupied(orbital)) return; //can't insert an antisegment, orbital is empty.
  double t_start=random()*beta; //start time of the anti segment
  if(local_config.exists(t_start)){ /*std::cerr<<"rare event, duplicate: "<<t_start<<std::endl; */return;} //time already exists.
  double t_next_segment_start=local_config.find_next_segment_start_distance(t_start,orbital);
  double t_next_segment_end=local_config.find_next_segment_end_distance(t_start,orbital);
  
  if(t_next_segment_start < t_next_segment_end) return; //we're trying to create an antisegment where there is no segment abort.
  
  //draw an end time
  double t_len = random()*t_next_segment_end;
  double t_end=t_start+t_len; //((t_len<0.1*beta)?t_len:0.1*beta); //random()*t_next_segment_end;
  if(t_end >= beta) t_end-=beta;
  if(local_config.exists(t_end)){ /*std::cerr<<"rare event, duplicate: "<<t_end<<std::endl; */return;} //time already exists.
  if(t_end<=t_start || t_end<=0.0){ /*std::cerr<<"rare event, zero length segment: "<<t_start<<" "<<t_end<<std::endl; */return;  } //time already exists.
  
  //std::cout<<clgreen<<"antisegment insertion update: "<<std::endl<<cblack<<*this<<std::endl;
  //std::cout<<clgreen<<" antisegment start time: (cdagger): "<<t_start<<" end time (c): "<<t_end<<std::endl;
  
  //compute local weight of the removed segment with t_start and t_end
  //Leo: need to check!
  segment new_segment(t_start, t_end, color_temp, color_temp);
  //std::cout<<clred<<"antisegment insert."<<std::endl;
  double local_weight_change=local_config.local_weight_change(new_segment, orbital, true);
  //std::cout<<clred<<"antisegment insert done."<<std::endl;
  
  //compute hybridization weight change //Leo: I don't quite understand...
  segment new_antisegment(t_end, t_start, color_temp, color_temp);
  double hybridization_weight_change=hyb_config.hyb_weight_change_insert(new_antisegment, orbital, color_temp);
  
  //compute the proposal probability ratio
  double permutation_factor=beta*t_next_segment_end/(local_config.order(orbital)+1);
  
  //perform metropolis
  double weight_change=local_weight_change*hybridization_weight_change*permutation_factor;
  
  //std::cout<<clred<<"weight change: "<<weight_change<<" l: "<<local_weight_change<<" h: "<<hybridization_weight_change<<" p: "<<permutation_factor<<cblack<<std::endl;
  
  if(std::abs(weight_change)>random())
  {
    nacc[3]++;
    //std::cout<<cred<<"accepting insert antisegment."<<cblack<<std::endl;
    if(weight_change < 0) sign*=-1.;
    local_config.insert_antisegment(new_antisegment, orbital);
    hyb_config.insert_antisegment(new_antisegment, orbital, color_temp);
    //Leo: record the updated color and set color_updated to true
    color = color_temp;
    color_updated = true;
    //std::cout<<cred<<"done accepting insert antisegment."<<cblack<<std::endl;
  }
}


void hybridization::remove_antisegment_update(int orbital, std::size_t color_temp)
{
  nprop[4]++;
  int k=local_config.order(orbital);
  
  if(k==0) return; //no point, this is an empty orbital
  int segment_nr=(int)(random()*k);
  
  //try to merge segment k and segment k+1
  segment segment_earlier=local_config.get_segment(segment_nr, orbital);
  segment segment_later  =local_config.get_segment(segment_nr==k-1?0:segment_nr+1, orbital);

  //Leo: for debug purpose, print the colors of the chosen segment
  std::cout << "remove_antisegment_update: c_start=" << segment_earlier.c_end_ <<", "\
            << "c_end=" << segment_later.c_start_ << std::endl;

  //Leo: check if the colors of both ends and the randomly picked color are all the same
  if(color_temp == segment_earlier.c_end_ && color_temp == segment_later.c_start_) ; //do nothing
  else {return;} //cannot remove because of the different colors
  
  //std::cout<<clcyan<<"antisegment removal update: "<<cblack<<*this<<std::endl;
  //std::cout<<clcyan<<" antisegment start time: (cdagger) "<<segment_earlier.t_end_<<" end time: (c): "<<segment_later.t_start_<<std::endl;

  //compute local weight of the antisegment. note that time direction here has to be forward
  //Leo: not sure I understand why...
  segment segment_forward(segment_earlier.t_end_, segment_later.t_start_);
  //std::cout<<clgreen<<"antisegment remove."<<std::endl;
  double local_weight_change=1./local_config.local_weight_change(segment_forward, orbital, true);
  //std::cout<<clgreen<<"antisegment remove done."<<std::endl;
  
  //compute hybridization weight change
  //Leo: need to check!
  segment antisegment(segment_later.t_start_, segment_earlier.t_end_, color_temp, color_temp);
  double hybridization_weight_change=hyb_config.hyb_weight_change_remove(antisegment, orbital, color_temp);
  
  //compute the proposal probability ratio
  //Leo: not sure I understand why...
  double t_next_segment_end=local_config.order(orbital)==1?beta:segment_later.t_end_-segment_earlier.t_end_; 
  if(t_next_segment_end<0.) t_next_segment_end+=beta;
  double permutation_factor=local_config.order(orbital)/(beta*t_next_segment_end);
  //perform metropolis
  double weight_change=local_weight_change/hybridization_weight_change*permutation_factor;
  //std::cout<<clblue<<"weight change: "<<weight_change<<" l: "<<local_weight_change<<" h: "<<hybridization_weight_change<<" p: "<<permutation_factor<<cblack<<std::endl;
  
  if(std::abs(weight_change)>random())
  {
    nacc[4]++;
    //std::cout<<cred<<"accepting remove antisegment."<<cblack<<std::endl;
    if(weight_change < 0) sign*=-1.;
    local_config.remove_antisegment(antisegment, orbital);
    hyb_config.remove_antisegment(antisegment, orbital, color_temp);
    //Leo: record the updated color and set color_updated to true
    color = color_temp;
    color_updated = true;
    //std::cout<<cred<<"done accepting remove antisegment."<<cblack<<std::endl;
  }
}


/********************************************************\
 New spin-flip updates
 Idea: take remove_segment_update and insert_segment_update and combine them so one segment is removed from on orbital (spin up or down) and inserted on the corresponding other_orbital (spin down or up), if the other_orbital is not filled.
 \********************************************************/
void hybridization::spin_flip_update(int orbital, std::size_t color_temp)
{  
  nprop[5]++;
  
  int k=local_config.order(orbital);
  
  if(k==0) return; //no point, this is an empty orbital
  //    if (local_config.zero_order_orbital_occupied(orbital)) return;
  int other_orbital=(int)(random()*(n_orbitals-1));
  other_orbital = (other_orbital<orbital)?other_orbital:1+other_orbital;
  int other_orbital_order=local_config.order(other_orbital);
  if (orbital == other_orbital) return;
  if (local_config.zero_order_orbital_occupied(other_orbital)) return;
  int segment_nr=(int)(random()*k);
  segment segment_to_flip=local_config.get_segment(segment_nr, orbital);
  
  //Leo: check if the colors of both ends and the randomly picked color are all the same
  if(color_temp == segment_to_flip.c_start_ && color_temp == segment_to_flip.c_end_) ; //do nothing
  else {return;} //cannot flip because of the different colors
  
  if (local_config.has_overlap(segment_to_flip,other_orbital)) return;
  
  double t_next_segment_start=local_config.find_next_segment_start_distance(segment_to_flip.t_start_,other_orbital);
  double t_next_segment_end=local_config.find_next_segment_end_distance(segment_to_flip.t_start_,other_orbital);
  double seg_length = segment_to_flip.t_end_ - segment_to_flip.t_start_;
  if (seg_length<0.0) seg_length += beta;
  
  //Leo: copied from the cthyb code on GitHub
  // check whether other_orbital is already filled where we want to insert the segment; 
  // confer with insert_antisegment
  if ((t_next_segment_start > t_next_segment_end)||
      (t_next_segment_start<seg_length)||
      (t_next_segment_end-beta>0)) return; 

  //Totally experimential - mistakes here?
  //Leo: this part excludes the possibility of flipping a wrapping segment!
  //     Furthermore, it's impossible for t_end to be larger than beta or smaller than 0
  //     This part is to be fixed, at least to sync with alps-cthyb hosted on GitHub 
  //Leo: fixed by imposing the above part of code and commenting out the third line below
  double t_start = segment_to_flip.t_start_, t_end=segment_to_flip.t_end_;
  if(t_end >= beta) t_end-=beta;
  //if(t_end<=t_start || t_end<=0.0) { /*std::cerr<<"rare (impossible?) event: t_start = t_end."<<std::endl; */return; }
  
  //compute local weight change: As we intend to propose a flip, we can
  //safely ignore the intermediate state and directly compare the energies
  //of the two states involved
  // This energy is associated with the present segment. We give it a negative weight
  // because it is to be removed
  double de = -local_config.local_energy(segment_to_flip,orbital); //Leo: This line is problematic; see GitHub version.
  //double full_weight_orig=full_weight();
  double hybridization_weight_change_1=1./hyb_config.hyb_weight_change_remove(segment_to_flip, orbital, color_temp); // from line 187 - remove_segment_update
  double weight_change_1=hybridization_weight_change_1;
  if (weight_change_1<0) sign*=-1;
  local_config.remove_segment(segment_to_flip, orbital);
  hyb_config.remove_segment(segment_to_flip, orbital, color_temp);
  //double full_weight_removed=full_weight();
  segment new_segment(t_start, t_end, color_temp, color_temp);
  de += local_config.local_energy(new_segment,other_orbital);
  double hybridization_weight_change_2=hyb_config.hyb_weight_change_insert(new_segment, other_orbital, color_temp);   
  double weight_change_2=std::exp(de)*hybridization_weight_change_2;

  //the permutation factor has the probability of proposing this move (the old order in this orbital) divided by the probability of proposing the reverse move (the new order in the new orbital)
  double permutation_factor=k/(double)(other_orbital_order+1.);

  if(std::abs(weight_change_1*weight_change_2*permutation_factor)>random())
  { //Accepted
    nacc[5]++;
    if(weight_change_2 < 0) sign*=-1.;
    local_config.insert_segment(new_segment, other_orbital);
    hyb_config.insert_segment(new_segment, other_orbital, color_temp);
    //Leo: record the updated color and set color_updated to true
    color = color_temp;
    color_updated = true;
  } 
  else 
  { //Not accepted, thus restore old configuration
    double wc = hyb_config.hyb_weight_change_insert(new_segment, orbital, color_temp);
    if (wc<0.0) sign *= -1;
    local_config.insert_segment(new_segment, orbital);
    hyb_config.insert_segment(new_segment, orbital, color_temp);
  }
}
