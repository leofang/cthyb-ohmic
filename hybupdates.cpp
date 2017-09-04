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
#include "hyb.hpp"
#include "hyblocal.hpp"
#include "hybdissipation.hpp"

//#include "combinatorial.hpp" //Leo: test!!!!!!!!!!!!!

//this is the heart of the Monte Carlo procedure: we have the following updates:
//1: change the zero order state, swap an empty orbital versus a filled one (5% or 10%)
//2: shift an existing segment start- or end-point (not implemented)
//3: insert or remove a new segment (40% or 50%)
//4: insert or remove an anti-segment (30% or 40%)
//5: perform a segment flip between different orbtials (20%)
//7. perform a global color flip in a given orbital (5%) (added by Leo)
//see our review for details of these updates
//
//Leo: the global exchange of two orbitals (global_flip_update) is disabled because it's too painful to modify

void hybridization::update()
{
  //one sweep is composed of N_MEAS Monte Carlo updates and one measurement (the latter only if thermalized)
  sweeps++;
  
  double rates[2] = {(spin_flip)?0.55:0.65,(spin_flip)?0.9:1.0};

  for(std::size_t i=0;i<N_meas;++i)
  {
    double update_type=random();

    //if (update_type < 0.02 && global_flip)
    if (update_type < 0.02 && color_flip)
    {
      //global_flip_update();
      color_flip_update();
    } 
    else if (update_type<0.05)
//    if (update_type<0.1)
    {
      change_zero_order_state_update();
//    } else if (update_type < 0) {
//      shift_segment_update();
    }
    else if (update_type < 0.15 && worm_update)
    {
      insert_worm_update();
    }
    else if (update_type < 0.2 && worm_update)
    {
      remove_worm_update();
    }
    else if (update_type < 0.3 && worm_update)
    {
      worm_creep_update();
    }
    else if (update_type < 0.15 && color_swap) //TODO: think a better, cleverer way
    {
      color_swap_update();
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
      //measure_color(); //Leo: for color measurement; TODO: color is measured when accepted, so this is unnecessary!
      if(MEASURE_time && !has_worm)
      {
        local_config.get_F_prefactor(F_prefactor);//compute segment overlaps in local config
        measure_G(F_prefactor);
      }

      if(MEASURE_time_worm && has_worm)
      {
        measure_G_worm();
      }
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
}


void hybridization::change_zero_order_state_update()
{
  nprop[0]++; N_Z++;
  if(has_worm) return;

  //choose the orbital in which we do the update
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
      if(VERY_VERBOSE && sweeps<=debug_number) 
      { 
         debug_output(0, local_weight_change, 1, 1, 1);
      }
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
      if(VERY_VERBOSE && sweeps<=debug_number) 
      { 
         debug_output(0, local_weight_change, 1, 1, 1); 
      }
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

  if(random()<0.5){ insert_segment_update(orbital); }
  else            { remove_segment_update(orbital); }
}


void hybridization::insert_remove_antisegment_update()
{
  //choose the orbital in which we do the update
  int orbital=(int)(random()*n_orbitals);

  if(random()<0.5){ insert_antisegment_update(orbital); }
  else            { remove_antisegment_update(orbital); }
}


void hybridization::insert_remove_spin_flip_update()
{
  //choose the orbital in which we do the update
  int orbital=(int)(random()*n_orbitals);

  spin_flip_update(orbital);
}


void hybridization::color_flip_update()
{
  //choose the orbital in which we do the update
  int orbital=(int)(random()*n_orbitals);

  color_flip_update(orbital);
}


void hybridization::color_swap_update()
{
  //choose the orbital in which we do the update
  int orbital=(int)(random()*n_orbitals);

  color_swap_update(orbital);
}


void hybridization::insert_worm_update()
{
  //choose the orbital in which we do the update
  int orbital=(int)(random()*n_orbitals);

  if(random()<0.5){ insert_worm_segment_update(orbital); }
  else            { insert_worm_antisegment_update(orbital); }
}


void hybridization::remove_worm_update()
{
  //choose the orbital in which we do the update
  int orbital=(int)(random()*n_orbitals);

  if(random()<0.5){ remove_worm_segment_update(orbital); }
  else            { remove_worm_antisegment_update(orbital); }
}


//Leo: experimental global color flipping in a randomly chosen orbital
//Notes:
//1. currently only two colors are supported, so this part should be extended in the future //TODO
//2. local weight will not be affected by this update
void hybridization::color_flip_update(int orbital)
{
  nprop[7]++; N_Z++;
  if(has_worm) return;

  if( local_config.order(orbital)==0 ) return; //no segment for flipping
  
  //Leo: choose the color in which we do the update
  //Now only two colors (red/1 and blue/0) are considered, but it can be easily changed
  std::size_t color_1 = (int)(random()*n_env);
  std::size_t color_2 = (color_1 ? 0 : 1); //TODO: change this line for more colors

  //compute hybridization weight change
  double hybridization_weight_change = hyb_config.hyb_weight_change_flip(orbital, color_1, color_2);

  //compute the dissipation weight change
  double dissipation_weight_change = 1.0;//ohmic_config.color_flip_weight_change(orbital, color_1, color_2, local_config);
  
  //perform metropolis
  double weight_change = hybridization_weight_change * dissipation_weight_change;

  if(std::abs(weight_change)>random())
  {
    nacc[7]++;
    if(weight_change < 0) sign*=-1.;
    local_config.flip_color(orbital, color_1, color_2); 
    int color_diff = hyb_config.flip_color(orbital, color_1, color_2); //return (old) size1-size2
    dissipation_weight_ratio = dissipation_weight_change; //TODO: check if necessary

    //Leo: record the updated color 
    //color = color_temp;
    updated_colors[color_1]++;
    updated_colors[color_2]++;
    ncolor[color_1]++;
    ncolor[color_2]++;
    //updated_colors[color_1] += abs(color_diff);
    //updated_colors[color_2] += abs(color_diff);
    //ncolor[color_1] += abs(color_diff);
    //ncolor[color_2] += abs(color_diff);
    ncolor_diff[color_1] -= color_diff;  //Leo: + for insertion, - for removal 
    ncolor_diff[color_2] += color_diff;  //Leo: + for insertion, - for removal 

    /* Leo Fang: for test purpose, print out a lot of things... */
    if(VERY_VERBOSE && sweeps<=debug_number) 
    { 
       debug_output(7, 1, hybridization_weight_change, dissipation_weight_change, 1); 
    }
  }
}



void hybridization::insert_segment_update(int orbital)
{
  if(!has_worm) //in Z-space
  {
     nprop[1]++; N_Z++;
  }
  else //in G-space
  {
     nprop[13]++; N_G++;
  }

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
  
  //Leo: a length-zero segment is rare event, but I suspect it is still possible
  if(t_len==0) return;

  double t_end=t_start+t_len;
  if(t_end > beta) t_end-=beta;
  if(local_config.exists(t_end)){ /*std::cerr<<"rare event, duplicate: "<<t_end<<std::endl; */return;} //time already exists.

  //Leo: this line seems buggy (and doesn't exist in the GitHub version) because it may
  //     prevent a wrapping segment from being inserted.
  //if(t_end<=t_start || t_end<=0.0){ /*std::cerr<<"rare event, zero length segment: "<<t_start<<" "<<t_end<<std::endl; */return;} //time already exists.
  
  //Leo: choose the color in which we do the update
  //Now only two colors (red/1 and blue/0) are considered, but it can be easily changed
  std::size_t color_temp = (int)(random()*n_env);
  //Leo: this part was older, lenther code 
//  if(n_env == 1) { color_temp = 0; } // Only one color
//  else if(n_env==2)
//  { 
//     if(random()<0.5){ color_temp = 0; } // blue = 0 = R
//     else            { color_temp = 1; } // red = 1 = L  
//  }
//  else 
//  {  //Leo: because we set up the sanity_check, this line is redundant here for the moment... 
//     throw std::runtime_error("The input N_ENV>2 is currently not supported.");
//  }

  //Leo: paint the color on the segment
  segment new_segment(t_start, t_end, color_temp, color_temp);

  //compute local weight of the new segment with t_start and t_end
  double local_weight_change=local_config.local_weight_change(new_segment, orbital, false);
  
  //compute hybridization weight change
  double hybridization_weight_change=hyb_config.hyb_weight_change_insert(new_segment, orbital, color_temp);

  //Leo: compute the dissipation weight change
  double dissipation_weight_change=ohmic_config.dissipation_weight_change(new_segment, orbital, true, local_config);
  
  //compute the proposal probability ratio
  //Leo: the old algorithm must be modified when n_env>1; test!
  double permutation_factor=beta*t_next_segment_start/(local_config.order(orbital)+(has_worm?0:1));
  //double permutation_factor=beta*t_next_segment_start/(n_segments_temp[color_temp]);
 
  //perform metropolis
  double weight_change=local_weight_change*hybridization_weight_change*permutation_factor*dissipation_weight_change;

//  //Leo: test!!!!
//  int k = local_config.order(orbital);
//  //double combinatorial_factor_change = (k+1.)/(4.*k+2.);
//  double combinatorial_factor_change = combinatorial_factor[k]/combinatorial_factor[k+1];
//  double weight_change=local_weight_change*hybridization_weight_change*permutation_factor*dissipation_weight_change*combinatorial_factor_change;
  
  /*std::cout<<" new segment: "<<new_segment<<std::endl;
   std::cout<<clred<<"weight change: "<<weight_change<<" l: "<<local_weight_change<<" h: "<<hybridization_weight_change<<" p: "<<permutation_factor<<cblack<<std::endl;*/
  
  if(std::abs(weight_change)>random())
  {
    if(!has_worm) //in Z-space
    {
       nacc[1]++;
       // //Leo: compute the number of segments and antisegments of the new configuration
       // //only do this in the Z-space
       // std::vector<int> n_segments_temp = local_config.get_new_n_segments_insert_segment(new_segment, orbital);
       // local_config.set_n_segments(orbital, n_segments_temp); //Leo: update the number of segments and antisegments 
    }
    else //in G-space
    {
       nacc[13]++;
    }
    if(weight_change < 0) sign*=-1.;
    local_config.insert_segment(new_segment, orbital);
    hyb_config.insert_segment(new_segment, orbital, color_temp);
    dissipation_weight_ratio = 1.0/dissipation_weight_change; //Leo: keep the weight change (of removal!) for measure_G   

    //Leo: record the updated color and set color_updated to true
    color = color_temp;
    updated_colors[color]++;
    ncolor[color]++;
    ncolor_diff[color]++;  //Leo: + for insertion, - for removal 
    //color_updated = true;

    /* Leo Fang: for test purpose, print out a lot of things... */
    if(VERY_VERBOSE && sweeps<=debug_number) 
    { 
       debug_output( (!has_worm?1:13), local_weight_change, hybridization_weight_change, dissipation_weight_change, permutation_factor); 
    }
  }
}


void hybridization::remove_segment_update(int orbital)
{
  if(!has_worm) //in Z-space
  {
     nprop[2]++; N_Z++;
  }
  else //in G-space
  {
     nprop[14]++; N_G++;
  }

  //std::cout<<clblue<<"starting removal update."<<cblack<<std::endl;
  int k=local_config.order(orbital);
  
  if(k==0) return; //no point, this is an empty orbital //Leo: it could also be occupied, but it doesn't matter here
  
  int segment_nr=(int)(random()*k);
  
  segment segment_to_remove=local_config.get_segment(segment_nr, orbital);

//  //Leo: for debug purpose, print the colors of the chosen segment
//  if(VERY_VERBOSE)
//  {
//      std::cout << "remove_segment_update: c_start=" <<segment_to_remove.c_start_<<", "\
//                << "c_end=" << segment_to_remove.c_end_ << std::endl;
//  }

  std::size_t color_temp = segment_to_remove.c_start_;
  //Leo: check if the colors of both ends and the randomly picked color are all the same
  //Leo: rolling a dice is crucial here!
  if(color_temp != segment_to_remove.c_end_ || color_temp != (int)(random()*n_env)) 
     return; //cannot remove because of the different colors
//  if(color_temp == segment_to_remove.c_end_) ; //do nothing
//  else {return;} //cannot remove because of the different colors
 
  //Leo: get the current number of segments and antisegments
  //std::vector<int> n_segments = local_config.get_n_segments(orbital); 
 
  double local_weight_change=1./local_config.local_weight_change(segment_to_remove, orbital, false);
  
  //compute hybridization weight change
  double hybridization_weight_change=1.0/hyb_config.hyb_weight_change_remove(segment_to_remove, orbital, color_temp);
  
  //Leo: compute the dissipation weight change
  double dissipation_weight_change=1.0/ohmic_config.dissipation_weight_change(segment_to_remove, orbital, false, local_config);
  
  //compute the proposal probability ratio
  double t_next_segment_start=local_config.find_next_segment_start_distance(segment_to_remove.t_start_,orbital);
  //Leo: the old algorithm must be modified when n_env>1
  double permutation_factor = (local_config.order(orbital) + (has_worm?-1:0)) / (beta*t_next_segment_start);
  //double permutation_factor=n_segments[color_temp]/(beta*t_next_segment_start);
  
  //perform metropolis;
  double weight_change=local_weight_change*hybridization_weight_change*permutation_factor*dissipation_weight_change;

//  //Leo: test!!!!!!!!!!!!!!!!
//  double combinatorial_factor_change = combinatorial_factor[k]/combinatorial_factor[k-1];
//  //double combinatorial_factor_change = (4.*k-2.)/k;
//  double weight_change=local_weight_change*hybridization_weight_change*permutation_factor*dissipation_weight_change*combinatorial_factor_change;
  
  /*if(segment_nr==k-1){
   std::cout<<" segment_to_remove: "<<segment_to_remove<<std::endl;
   std::cout<<"t_next_segment_start: "<<t_next_segment_start<<std::endl;
   std::cout<<clblue<<"weight change: "<<weight_change<<" l: "<<local_weight_change<<" h: "<<hybridization_weight_change<<" p: "<<permutation_factor<<cblack<<std::endl;
   }*/
  if(std::abs(weight_change)>random())
  {
    if(!has_worm) //in Z-space
    {
       nacc[2]++;
       // //Leo: compute the number of segments and antisegments of the new configuration
       // std::vector<int> n_segments_temp = local_config.get_new_n_segments_remove_segment(segment_to_remove, orbital);
       // local_config.set_n_segments(orbital, n_segments_temp); //Leo: update the number of segments and antisegments 
    }
    else //in G-space
    {
       nacc[14]++;
    }
    if(weight_change < 0) sign*=-1.;

//      double fwo = full_weight();
    local_config.remove_segment(segment_to_remove, orbital);
    hyb_config.remove_segment(segment_to_remove, orbital, color_temp);
    dissipation_weight_ratio = dissipation_weight_change; //Leo: keep the weight change for measure_G   

    //Leo: record the updated color and set color_updated to true
    color = color_temp;
    updated_colors[color]++;
    ncolor[color]++;
    ncolor_diff[color]--;  //Leo: + for insertion, - for removal 
    //color_updated = true;
//      double fwa = full_weight();
//      std::cout << clgreen<<"weight change removal: "<<fwa<<" control: "<<fwo*std::abs(weight_change)<<std::endl;

    /* Leo Fang: for test purpose, print out the segment map */
    if(VERY_VERBOSE && sweeps<=debug_number) 
    { 
        debug_output( (!has_worm?2:14), local_weight_change, hybridization_weight_change, dissipation_weight_change, permutation_factor); 
    }
  }
}


void hybridization::insert_antisegment_update(int orbital)
{
  if(!has_worm) //in Z-space
  {
     nprop[3]++; N_Z++;
  }
  else //in G-space
  {
     nprop[15]++; N_G++;
  }

  if(local_config.order(orbital)==0 && !local_config.zero_order_orbital_occupied(orbital)) return; //can't insert an antisegment, orbital is empty.
  double t_start=random()*beta; //start time of the anti segment
  if(local_config.exists(t_start)){ /*std::cerr<<"rare event, duplicate: "<<t_start<<std::endl; */return;} //time already exists.
  double t_next_segment_start=local_config.find_next_segment_start_distance(t_start,orbital);
  double t_next_segment_end=local_config.find_next_segment_end_distance(t_start,orbital);
  
  if(t_next_segment_start < t_next_segment_end) return; //we're trying to create an antisegment where there is no segment abort.
  
  //draw an end time
  double t_len = random()*t_next_segment_end;

  //Leo: a length-zero antisegment is rare event, but I suspect it is still possible
  if(t_len==0) return;

  double t_end=t_start+t_len; //((t_len<0.1*beta)?t_len:0.1*beta); //random()*t_next_segment_end;
  if(t_end > beta) t_end-=beta;
  if(local_config.exists(t_end)){ /*std::cerr<<"rare event, duplicate: "<<t_end<<std::endl; */return;} //time already exists.

  //Leo: this line seems buggy (and doesn't exist in the GitHub version) because it may
  //     prevent a wrapping antisegment from being inserted.
//  if(t_end<=t_start || t_end<=0.0){ /*std::cerr<<"rare event, zero length segment: "<<t_start<<" "<<t_end<<std::endl; */return;  } //time already exists.
  
  //std::cout<<clgreen<<"antisegment insertion update: "<<std::endl<<cblack<<*this<<std::endl;
  //std::cout<<clgreen<<" antisegment start time: (cdagger): "<<t_start<<" end time (c): "<<t_end<<std::endl;
  
  //Leo: choose the color in which we do the update
  //Now only two colors (red/1 and blue/0) are considered, but it can be easily changed
  std::size_t color_temp = (int)(random()*n_env);
  //Leo: this part was older, lenther code 
//  if(n_env == 1) { color_temp = 0; } // Only one color
//  else if(n_env==2)
//  { 
//     if(random()<0.5){ color_temp = 0; } // blue = 0 = R
//     else            { color_temp = 1; } // red = 1 = L  
//  }
//  else 
//  {  //Leo: because we set up the sanity_check, this line is redundant here for the moment... 
//     throw std::runtime_error("The input N_ENV>2 is currently not supported.");
//  }

  //compute local weight of the removed segment with t_start and t_end
  //Leo: need to check!
  segment new_segment(t_start, t_end, color_temp, color_temp);
  //std::cout<<clred<<"antisegment insert."<<std::endl;
  double local_weight_change=local_config.local_weight_change(new_segment, orbital, true);
  //std::cout<<clred<<"antisegment insert done."<<std::endl;
  
  //compute hybridization weight change //Leo: I don't quite understand...
  segment new_antisegment(t_end, t_start, color_temp, color_temp);
  double hybridization_weight_change=hyb_config.hyb_weight_change_insert(new_antisegment, orbital, color_temp);

  //Leo: compute the dissipation weight change
  double dissipation_weight_change=ohmic_config.dissipation_weight_change(new_antisegment, orbital, true, local_config);
  
  //compute the proposal probability ratio
  //Leo: the old algorithm must be modified when n_env>1
  double permutation_factor=beta*t_next_segment_end/(local_config.order(orbital)+(has_worm?0:1));
  //double permutation_factor=beta*t_next_segment_end/(n_segments_temp[color_temp+n_env]);

  //perform metropolis
  double weight_change=local_weight_change*hybridization_weight_change*permutation_factor*dissipation_weight_change;

//  //Leo: test!!!!
//  int k = local_config.order(orbital);
//  //double combinatorial_factor_change = (k+1.)/(4.*k+2.);
//  double combinatorial_factor_change = combinatorial_factor[k]/combinatorial_factor[k+1];
//  double weight_change=local_weight_change*hybridization_weight_change*permutation_factor*dissipation_weight_change*combinatorial_factor_change;

  //std::cout<<clred<<"weight change: "<<weight_change<<" l: "<<local_weight_change<<" h: "<<hybridization_weight_change<<" p: "<<permutation_factor<<cblack<<std::endl;
  
  if(std::abs(weight_change)>random())
  {
    if(!has_worm) //in Z-sapce
    {
       nacc[3]++;
       // //Leo: compute the number of segments and antisegments of the new configuration
       // std::vector<int> n_segments_temp = local_config.get_new_n_segments_insert_antisegment(new_antisegment, orbital);
       // local_config.set_n_segments(orbital, n_segments_temp); //Leo: update the number of segments and antisegments 
    }
    else //in G-space
    {
       nacc[15]++;
    }

    //std::cout<<cred<<"accepting insert antisegment."<<cblack<<std::endl;
    if(weight_change < 0) sign*=-1.;
    local_config.insert_antisegment(new_antisegment, orbital);
    hyb_config.insert_antisegment(new_antisegment, orbital, color_temp);
    dissipation_weight_ratio = 1.0/dissipation_weight_change; //Leo: keep the weight change (of removal!) for measure_G   

    //Leo: record the updated color and set color_updated to true
    color = color_temp;
    updated_colors[color]++;
    ncolor[color]++;
    ncolor_diff[color]++;  //Leo: + for insertion, - for removal 
    //color_updated = true;
    //std::cout<<cred<<"done accepting insert antisegment."<<cblack<<std::endl;

    /* Leo Fang: for test purpose, print out the segment map */
    if(VERY_VERBOSE && sweeps<=debug_number) 
    {
        debug_output( (!has_worm?3:15), local_weight_change, hybridization_weight_change, dissipation_weight_change, permutation_factor); 
    }
  }
}


void hybridization::remove_antisegment_update(int orbital)
{
  if(!has_worm) //in Z-space
  {
     nprop[4]++; N_Z++;
  }
  else //in G-space
  {
     nprop[16]++; N_G++;
  }

  int k=local_config.order(orbital);
  
  if(k==0) return; //no point, this is an empty orbital
  int segment_nr=(int)(random()*k);
  
  //try to merge segment k and segment k+1
  segment segment_earlier=local_config.get_segment(segment_nr, orbital);
  segment segment_later  =local_config.get_segment(segment_nr==k-1?0:segment_nr+1, orbital);

//  //Leo: for debug purpose, print the colors of the chosen segment
//  if(VERY_VERBOSE)
//  {
//     std::cout << "remove_antisegment_update: c_start=" << segment_earlier.c_end_ <<", "\
//               << "c_end=" << segment_later.c_start_ << std::endl;
//  }
  
  std::size_t color_temp = segment_earlier.c_end_;
  //Leo: check if the colors of both ends and the randomly picked color are all the same
  //Leo: rolling a dice is crucial here!
  if(color_temp != segment_later.c_start_ || color_temp != (int)(random()*n_env)) 
     return; //cannot remove because of the different colors
//  if(color_temp == segment_later.c_start_) ; //do nothing
//  else {return;} //cannot remove because of the different colors
  
  //std::cout<<clcyan<<"antisegment removal update: "<<cblack<<*this<<std::endl;
  //std::cout<<clcyan<<" antisegment start time: (cdagger) "<<segment_earlier.t_end_<<" end time: (c): "<<segment_later.t_start_<<std::endl;

  //compute local weight of the antisegment. note that time direction here has to be forward
  //Leo: not sure I understand why...
  segment segment_forward(segment_earlier.t_end_, segment_later.t_start_, color_temp, color_temp);
  //std::cout<<clgreen<<"antisegment remove."<<std::endl;
  double local_weight_change=1./local_config.local_weight_change(segment_forward, orbital, true);
  //std::cout<<clgreen<<"antisegment remove done."<<std::endl;
  
  //compute hybridization weight change
  //Leo: need to check!
  segment antisegment(segment_later.t_start_, segment_earlier.t_end_, color_temp, color_temp);
  double hybridization_weight_change=1.0/hyb_config.hyb_weight_change_remove(antisegment, orbital, color_temp);

  //Leo: compute the dissipation weight change
  double dissipation_weight_change=1.0/ohmic_config.dissipation_weight_change(antisegment, orbital, false, local_config);

  //compute the proposal probability ratio
  //Leo: not sure I understand why...
  double t_next_segment_end=local_config.order(orbital)==1?beta:segment_later.t_end_-segment_earlier.t_end_; 
  if(t_next_segment_end<0.) t_next_segment_end+=beta;
  //Leo: the old algorithm must be modified when n_env>1
  double permutation_factor = (local_config.order(orbital) + (has_worm?-1:0)) / (beta*t_next_segment_end);
  //double permutation_factor=n_segments[color_temp+n_env]/(beta*t_next_segment_end);

  //perform metropolis;
  double weight_change=local_weight_change*hybridization_weight_change*permutation_factor*dissipation_weight_change;

//  //Leo: test!!!!!!!!!
//  double combinatorial_factor_change = combinatorial_factor[k]/combinatorial_factor[k-1];
//  //double combinatorial_factor_change = (4.*k-2.)/k;
//  double weight_change=local_weight_change*hybridization_weight_change*permutation_factor*dissipation_weight_change*combinatorial_factor_change;
  //std::cout<<clblue<<"weight change: "<<weight_change<<" l: "<<local_weight_change<<" h: "<<hybridization_weight_change<<" p: "<<permutation_factor<<cblack<<std::endl;
  
  if(std::abs(weight_change)>random())
  {
    if(!has_worm) //in Z-space
    {
       nacc[4]++;
       // //Leo: compute the number of segments and antisegments of the new configuration
       // std::vector<int> n_segments_temp = local_config.get_new_n_segments_remove_antisegment(antisegment, orbital);
       // local_config.set_n_segments(orbital, n_segments_temp); //Leo: update the number of segments and antisegments 
    }
    else //in G-space
    {
       nacc[16]++;
    }
  
    //std::cout<<cred<<"accepting remove antisegment."<<cblack<<std::endl;
    if(weight_change < 0) sign*=-1.;

    local_config.remove_antisegment(antisegment, orbital);
    hyb_config.remove_antisegment(antisegment, orbital, color_temp);
    dissipation_weight_ratio = dissipation_weight_change; //Leo: keep the weight change for measure_G   

    //Leo: record the updated color and set color_updated to true
    color = color_temp;
    updated_colors[color]++;
    ncolor[color]++;
    ncolor_diff[color]--;  //Leo: + for insertion, - for removal 
    //color_updated = true;
    //std::cout<<cred<<"done accepting remove antisegment."<<cblack<<std::endl;

    /* Leo Fang: for test purpose, print out the segment map */
    if(VERY_VERBOSE && sweeps<=debug_number) 
    {
        debug_output( (!has_worm?4:16), local_weight_change, hybridization_weight_change, dissipation_weight_change, permutation_factor); 
    }
  }
}


//remove one segment from a random orbital (spin up or down) and insert it to the other orbital (spin down or up), 
//if the other orbital is not filled.
void hybridization::spin_flip_update(int orbital)
{  
  nprop[5]++;

  //don't do this update if the worm is present
  if(has_worm) return;
  
  int k = local_config.order(orbital);
  if(k==0) return; //no point, this is an empty orbital

  int other_orbital = (int)(random()*(n_orbitals-1));
  other_orbital = (other_orbital<orbital)?other_orbital:1+other_orbital;
  if (local_config.zero_order_orbital_occupied(other_orbital)) return; //the other orbital is completely occupied

  int segment_nr = (int)(random()*k);
  segment segment_to_flip = local_config.get_segment(segment_nr, orbital);
 
  //Leo: check if the color of both ends is the same
  size_t color_temp = segment_to_flip.c_start_; 
  if(color_temp != segment_to_flip.c_end_) return; //cannot flip because of the different colors
  
  //Leo: calculate the distances in the other orbital (NOT the current one!)
  double t_next_segment_start = local_config.find_next_segment_start_distance(segment_to_flip.t_start_, other_orbital);
  double t_next_segment_end   = local_config.find_next_segment_end_distance(segment_to_flip.t_start_, other_orbital);
  double seg_length = segment_to_flip.t_end_ - segment_to_flip.t_start_;
  if (seg_length<0.0) seg_length += beta;
  
  //Leo: check whether other_orbital is already filled where we want to insert the segment 
  if ( (t_next_segment_start > t_next_segment_end) || //overlap with the previous segment
      (t_next_segment_start < seg_length) ) return;   //overlap with the next segment

  //compute local weight change: As we intend to propose a flip, we can safely ignore the 
  //intermediate state and directly compare the energies of the two states involved
  // Leo: the on-site interaction can be ignored because we just guaranteed there is no overlap
  double local_weight_change = std::exp( (local_config.mu(other_orbital)-local_config.mu(orbital)) * seg_length);

  //compute hyb weight change for removal
  double hybridization_weight_change_1 = 1./hyb_config.hyb_weight_change_remove(segment_to_flip, orbital, color_temp);
  if (hybridization_weight_change_1 < 0) sign *= -1;

  //temporarily remove the segment for the chosen orbital
  //Leo: I think it's unnecessary...
  //local_config.remove_segment(segment_to_flip, orbital);
  //hyb_config.remove_segment(segment_to_flip, orbital, color_temp);

  //compute hyb weight change for insertion
  double hybridization_weight_change_2 = hyb_config.hyb_weight_change_insert(segment_to_flip, other_orbital, color_temp);   

  double dissipation_weight_change = 1.0; //TODO: write this part!

  //the permutation factor has the probability of proposing this move (the old order in this orbital) 
  //divided by the probability of proposing the reverse move (the new order in the new orbital)
  //Leo: note that the "next_segment_start_distance" for the two moves are different, and that the 
  //common beta factor is canceled out
  double permutation_factor = k / local_config.find_next_segment_start_distance(segment_to_flip.t_start_, orbital)  //for removal
			      * t_next_segment_start / (local_config.order(other_orbital)+1);                       //for insertion   

  double total_weight_change = local_weight_change * hybridization_weight_change_1 * hybridization_weight_change_2
			       * dissipation_weight_change * permutation_factor;

  if(std::abs(total_weight_change)>random())
  { //Accepted
    nacc[5]++;
    if(hybridization_weight_change_2 < 0) sign *= -1.;

    //remove from the orbital
    local_config.remove_segment(segment_to_flip, orbital);
    hyb_config.remove_segment(segment_to_flip, orbital, color_temp);

    //insert to the other orbital
    local_config.insert_segment(segment_to_flip, other_orbital);
    hyb_config.insert_segment(segment_to_flip, other_orbital, color_temp);

    //Leo: record the updated color and set color_updated to true
    color = color_temp;
    updated_colors[color]++;
    ncolor[color]++;
    //Leo: do nothing to ncolor_diff because the total number of colored segments is not changed.
    //color_updated = true;
  } 
  else 
  { //rejected
    //Leo: nothing to be changed except for the sign
    if(hybridization_weight_change_1 < 0) sign *= -1;
    
    // //Not accepted, thus restore old configuration
    // double wc = hyb_config.hyb_weight_change_insert(segment_to_flip, orbital, color_temp);
    // if (wc<0.0) sign *= -1;
    // local_config.insert_segment(segment_to_flip, orbital);
    // hyb_config.insert_segment(segment_to_flip, orbital, color_temp);
  }

  //Leo: for test purpose, print out a lot of things...
  if(VERY_VERBOSE && sweeps<=debug_number) 
  { 
     debug_output(5, local_weight_change, hybridization_weight_change_1*hybridization_weight_change_2, dissipation_weight_change, permutation_factor); 
  }
}


void hybridization::debug_output(int updatetype, const double &local_weight_change, const double &hybridization_weight_change, const double &dissipation_weight_change, const double &permutation_factor)
{
        static int counter = 1;
        int cur_prec = std::cout.precision();

        //Leo: time ordering the hyb matrices. When time is ordered, permutation_sign_ should be equal to time_ordering_sign_
        //hyb_config.rebuild_ordered();

        std::cout << std::endl;
	std::cout << "|---------------------------------------------------------------------------------|" << std::endl;
//	std::cout << "At " << i+(sweeps-1)*N_meas+1 << "-th update:" << std::endl; local_config.print_segments();
        std::cout << counter << "-th accepted move: " << update_type[updatetype] << std::endl; //local_config.print_segments();
        std::cout << local_config << std::endl;
        std::cout << hyb_config << std::endl;
        if(sign<0)   std::cout << "negative sign!" << std::endl;
        double hyb_full_weight = hyb_config.full_weight();
        if(hyb_full_weight<0)
        {
                std::cout << std::endl << "the weight of hyb_config is negative! (=" << hyb_full_weight \
                          << ")" << std::endl << std::endl;
                //throw std::runtime_error("Abort!"); 
        }
//        if( (hyb_full_weight<0&&hyb_config.overall_color_matrix_sign()>0) || 
//             (hyb_full_weight>0&&hyb_config.overall_color_matrix_sign()<0) )
//             throw std::runtime_error("inconsistent sign, abort!");
        std::cout << std::setprecision(15) << std::fixed;
        std::cout << "local weight change         = " << local_weight_change << std::endl;
        std::cout << "hybridization weight change = " << hybridization_weight_change << std::endl;
        std::cout << "dissipation weight change   = " << dissipation_weight_change << std::endl;
        std::cout << "permutation factor          = " << permutation_factor << std::endl;
        std::cout << "the number configuration    = ";
        for(int i=0; i<n_orbitals; i++)
        {
             std::vector<int> n_segments_temp = local_config.get_n_segments(i);
             std::cout << "(" << local_config.order(i) << "_";
             for (std::vector<int>::const_iterator it=n_segments_temp.begin(); it!=n_segments_temp.end(); it++)
                  std::cout << *it; 
             std::cout << ") ";
        }
        std::cout << std::endl;
	std::cout << "|---------------------------------------------------------------------------------|" << std::endl;

        std::cout.unsetf(std::ios_base::fixed);
        std::cout.precision(cur_prec);

	if(!has_worm)
           check_consistency(counter);

        counter++;
}



void hybridization::check_consistency(const int &counter)
{
  //static int counter = 0;
 
  //Leo: check the size of colored matrices for each orbital
  for(int i=0; i<n_orbitals; i++)
  { 
    std::stringstream temp_stream; // stringstream used for the conversion
    temp_stream << i;              // add the value of i to the characters in the stream
    if( local_config.order(i)==0 ) continue; //Leo: no segment exists, and so does color, in 0-th order
    if( hyb_config.total_color_matrix_size(i) != local_config.order(i) )
        std::runtime_error("The total size of colored matrices for orbital " + temp_stream.str() + " is incorrect!");
  }

  //Leo: check the consistency of the local configuration
  local_config.check_consistency();

  //Leo: time ordering the hyb matrices. When time is ordered, permutation_sign_ should be equal to time_ordering_sign_
  //hyb_config.rebuild_ordered();

  std::cout << "************** consistency checked at " << counter << "-th successful update! **************" << std::endl; 

  //counter++;
}



void hybridization::insert_worm_segment_update(int orbital)
{
  nprop[9]++; N_Z++;

  //cannot have more than one worm at the same time
  if(has_worm) return;

  //can't insert segment, orbital is fully occuppied.
  if(local_config.order(orbital)==0 && local_config.zero_order_orbital_occupied(orbital)) return; 

  double t_start = random()*beta; //start time of a segment
  if(local_config.exists(t_start)) return; //time already exists.

  double t_next_segment_start = local_config.find_next_segment_start_distance(t_start,orbital);
  double t_next_segment_end = local_config.find_next_segment_end_distance(t_start,orbital);

  if(t_next_segment_end < t_next_segment_start) return; //we're trying to create a segment on top of another segment. abort.
  
  //draw an end time
  double t_len = random()*t_next_segment_start;
  
  //Leo: a length-zero segment is rare event, but I suspect it is still possible
  if(t_len==0) return;

  double t_end = t_start + t_len;
  if(t_end > beta) t_end -= beta;
  if(local_config.exists(t_end)) return; //time already exists.

  //Leo: paint the worm color on the segment
  segment new_segment(t_start, t_end, WORM_COLOR, WORM_COLOR);

  //Leo: compute the number of segments and antisegments of the new configuration
  //std::vector<int> n_segments_temp = local_config.get_new_n_segments_insert_segment(new_segment, orbital);

  //compute local weight of the new segment with t_start and t_end
  double local_weight_change=local_config.local_weight_change(new_segment, orbital, false);
  
  //compute hybridization weight change
  //double hybridization_weight_change=hyb_config.hyb_weight_change_insert(new_segment, orbital, color_temp);

  //Leo: compute the dissipation weight change //TODO
  double dissipation_weight_change = 1.0;
  //double dissipation_weight_change=ohmic_config.dissipation_weight_change(new_segment, orbital, true, local_config);
  
  //compute the proposal probability ratio
  double permutation_factor = beta * t_next_segment_start;
 
  //perform metropolis
  double weight_change = eta * local_weight_change * permutation_factor * dissipation_weight_change;

  if(std::abs(weight_change)>random())
  {
    nacc[9]++;
    has_worm = true;
    if(weight_change < 0) sign*=-1.;
    local_config.insert_segment(new_segment, orbital);
    local_config.set_worm(orbital, t_start, t_end);
    //local_config.set_n_segments(orbital, n_segments_temp); //Leo: update the number of segments and antisegments 
    //hyb_config.insert_segment(new_segment, orbital, color_temp);
    dissipation_weight_ratio = 1.0/dissipation_weight_change; //Leo: keep the weight change (of removal!) for measure_G //TODO

    //Leo: record the updated color and set color_updated to true
    //color = color_temp;
    //updated_colors[color]++;
    //ncolor[color]++;
    //ncolor_diff[color]++;  //Leo: + for insertion, - for removal 
    ////color_updated = true;

    //Leo: for test purpose, print out a lot of things...
    if(VERY_VERBOSE && sweeps<=debug_number) 
    { 
       debug_output(9, local_weight_change, 1.0, dissipation_weight_change, permutation_factor); 
    }
  }
}



void hybridization::remove_worm_segment_update(int orbital)
{
  nprop[11]++; N_G++;

  //cannot remove worm if there's none
  if(!has_worm) return;

  //if worm's head and tail are adjacent, get the segment index of the worm, otherwise leave
  int segment_nr = 0;
  if(!local_config.is_worm_segment(segment_nr, orbital)) return;

  segment worm = local_config.get_segment(segment_nr, orbital);

  //Leo: get the current number of segments and antisegments
  //std::vector<int> n_segments = local_config.get_n_segments(orbital); 
 
  double local_weight_change = 1./local_config.local_weight_change(worm, orbital, false);
  
  //Leo: compute the dissipation weight change //TODO
  double dissipation_weight_change = 1.0;
  //double dissipation_weight_change=1.0/ohmic_config.dissipation_weight_change(worm, orbital, false, local_config);
  
  //compute the proposal probability ratio
  double t_next_segment_start = local_config.find_next_segment_start_distance(worm.t_start_,orbital);
  //Leo: the old algorithm must be modified when n_env>1
  double permutation_factor = 1.0/(beta*t_next_segment_start);
  //double permutation_factor=n_segments[color_temp]/(beta*t_next_segment_start);
  
  //perform metropolis;
  double weight_change = local_weight_change * permutation_factor * dissipation_weight_change / eta;

  if(std::abs(weight_change)>random())
  {
    nacc[11]++;
    has_worm = false;
    if(weight_change < 0) sign*=-1.;
    local_config.remove_segment(worm, orbital);
    // //Leo: compute the number of segments and antisegments of the new configuration
    // std::vector<int> n_segments_temp = local_config.rebuild_n_segments(orbital);
    // local_config.set_n_segments(orbital, n_segments_temp); //Leo: update the number of segments and antisegments
//      double fwo = full_weight();
    //hyb_config.remove_segment(worm, orbital, color_temp);
    dissipation_weight_ratio = dissipation_weight_change; //Leo: keep the weight change for measure_G //TODO 

    //Leo: record the updated color and set color_updated to true
    //color = color_temp;
    //updated_colors[color]++;
    //ncolor[color]++;
    //ncolor_diff[color]--;  //Leo: + for insertion, - for removal 
    //color_updated = true;
//      double fwa = full_weight();
//      std::cout << clgreen<<"weight change removal: "<<fwa<<" control: "<<fwo*std::abs(weight_change)<<std::endl;

    /* Leo Fang: for test purpose, print out the segment map */
    if(VERY_VERBOSE && sweeps<=debug_number) 
    { 
        debug_output(11, local_weight_change, 1.0, dissipation_weight_change, permutation_factor); 
    }
  }
}



void hybridization::insert_worm_antisegment_update(int orbital)
{
  nprop[10]++; N_Z++;

  //cannot have more than one worm at the same time
  if(has_worm) return;

  //can't insert an antisegment, orbital is empty.
  if(local_config.order(orbital)==0 && !local_config.zero_order_orbital_occupied(orbital)) return;

  double t_start = random()*beta; //start time of the anti segment
  if(local_config.exists(t_start)) return; //time already exists.
  double t_next_segment_start = local_config.find_next_segment_start_distance(t_start,orbital);
  double t_next_segment_end = local_config.find_next_segment_end_distance(t_start,orbital);
  
  if(t_next_segment_start < t_next_segment_end) return; //we're trying to create an antisegment where there is no segment abort.
  
  //draw an end time
  double t_len = random()*t_next_segment_end;

  //Leo: a length-zero antisegment is rare event, but I suspect it is still possible
  if(t_len==0) return;

  double t_end = t_start+t_len; //((t_len<0.1*beta)?t_len:0.1*beta); //random()*t_next_segment_end;
  if(t_end > beta) t_end -= beta;
  if(local_config.exists(t_end)) return; //time already exists.

  //compute local weight of the removed segment with t_start and t_end
  segment new_segment(t_start, t_end, WORM_COLOR, WORM_COLOR);

  double local_weight_change=local_config.local_weight_change(new_segment, orbital, true);
  
  //compute hybridization weight change //Leo: I don't quite understand...
  segment new_antisegment(t_end, t_start, WORM_COLOR, WORM_COLOR);
  //double hybridization_weight_change=hyb_config.hyb_weight_change_insert(new_antisegment, orbital, color_temp);

  //Leo: compute the dissipation weight change
  double dissipation_weight_change = 1.0; //TODO
  //double dissipation_weight_change=ohmic_config.dissipation_weight_change(new_antisegment, orbital, true, local_config);
  
  //Leo: compute the number of segments and antisegments of the new configuration
  //std::vector<int> n_segments_temp = local_config.get_new_n_segments_insert_antisegment(new_antisegment, orbital);

  //compute the proposal probability ratio
  //Leo: the old algorithm must be modified when n_env>1
  double permutation_factor = beta * t_next_segment_end;
  //double permutation_factor=beta*t_next_segment_end/(n_segments_temp[color_temp+n_env]);

  //perform metropolis
  double weight_change = eta * local_weight_change * permutation_factor * dissipation_weight_change;

  if(std::abs(weight_change)>random())
  {
    nacc[10]++;
    has_worm = true;
    //std::cout<<cred<<"accepting insert antisegment."<<cblack<<std::endl;
    if(weight_change < 0) sign*=-1.;
    local_config.insert_antisegment(new_antisegment, orbital);
    local_config.set_worm(orbital, t_end, t_start);
    //local_config.set_n_segments(orbital, n_segments_temp); //Leo: update the number of segments and antisegments 
    //hyb_config.insert_antisegment(new_antisegment, orbital, color_temp);
    dissipation_weight_ratio = 1.0/dissipation_weight_change; //Leo: keep the weight change (of removal!) for measure_G //TODO

    //Leo: record the updated color and set color_updated to true
    //color = color_temp;
    //updated_colors[color]++;
    //ncolor[color]++;
    //ncolor_diff[color]++;  //Leo: + for insertion, - for removal 
    //color_updated = true;
    //std::cout<<cred<<"done accepting insert antisegment."<<cblack<<std::endl;

    /* Leo Fang: for test purpose, print out the segment map */
    if(VERY_VERBOSE && sweeps<=debug_number) 
    {
        debug_output(10, local_weight_change, 1.0, dissipation_weight_change, permutation_factor); 
    }
  }
}


void hybridization::remove_worm_antisegment_update(int orbital)
{
  nprop[12]++; N_G++;

  //cannot remove worm if there's none
  if(!has_worm) return;

  //if worm's head and tail are adjacent, get the segment index of the worm, otherwise leave
  int k = local_config.order(orbital);
  int segment_nr = 0;
  if(!local_config.is_worm_antisegment(segment_nr, orbital)) return;

  //try to merge segment k and segment k+1
  segment segment_earlier = local_config.get_segment(segment_nr, orbital);
  segment segment_later   = local_config.get_segment(segment_nr==k-1?0:segment_nr+1, orbital);

  //check we really get the worm //TODO: this is redundant if everything works correctly
  assert(segment_earlier.c_end_ == segment_later.c_start_);
  assert(segment_earlier.c_end_ == WORM_COLOR);

  //compute local weight of the antisegment. note that time direction here has to be forward
  segment segment_forward(segment_earlier.t_end_, segment_later.t_start_, WORM_COLOR, WORM_COLOR);
  double local_weight_change=1./local_config.local_weight_change(segment_forward, orbital, true);
  
  //compute hybridization weight change
  //Leo: need to check!
  segment antisegment(segment_later.t_start_, segment_earlier.t_end_, WORM_COLOR, WORM_COLOR);
  //double hybridization_weight_change=1.0/hyb_config.hyb_weight_change_remove(antisegment, orbital, color_temp);

  //Leo: compute the dissipation weight change
  double dissipation_weight_change = 1.0; //TODO
  //double dissipation_weight_change=1.0/ohmic_config.dissipation_weight_change(antisegment, orbital, false, local_config);

  //Leo: get the current number of segments and antisegments
  //std::vector<int> n_segments = local_config.get_n_segments(orbital); 
  
  //compute the proposal probability ratio
  double t_next_segment_end = local_config.order(orbital)==1?beta:segment_later.t_end_-segment_earlier.t_end_; 
  if(t_next_segment_end<0.) t_next_segment_end += beta;
  //Leo: the old algorithm must be modified when n_env>1
  double permutation_factor = 1.0/(beta*t_next_segment_end);
  //double permutation_factor=n_segments[color_temp+n_env]/(beta*t_next_segment_end);

  //perform metropolis;
  double weight_change = local_weight_change * permutation_factor * dissipation_weight_change / eta;

  if(std::abs(weight_change)>random())
  {
    nacc[12]++;
    has_worm = false;
    if(weight_change < 0) sign*=-1.;
    local_config.remove_antisegment(antisegment, orbital);
    // //Leo: compute the number of segments and antisegments of the new configuration
    // std::vector<int> n_segments_temp = local_config.rebuild_n_segments(orbital);
    // local_config.set_n_segments(orbital, n_segments_temp); //Leo: update the number of segments and antisegments 
    //hyb_config.remove_antisegment(antisegment, orbital, color_temp);
    dissipation_weight_ratio = dissipation_weight_change; //Leo: keep the weight change for measure_G //TODO

    //Leo: record the updated color and set color_updated to true
    //color = color_temp;
    //updated_colors[color]++;
    //ncolor[color]++;
    //ncolor_diff[color]--;  //Leo: + for insertion, - for removal 
    //color_updated = true;
    //std::cout<<cred<<"done accepting remove antisegment."<<cblack<<std::endl;

    /* Leo Fang: for test purpose, print out the segment map */
    if(VERY_VERBOSE && sweeps<=debug_number) 
    {
        debug_output(12, local_weight_change, 1.0, dissipation_weight_change, permutation_factor); 
    }
  }
}


//this update attemps to swap the worm head with any d^dagger attached to a hybridization line
void hybridization::worm_creep_update()
{
  nprop[17]++; N_G++;

  //cannot move the worm if there's none
  if(!has_worm) return;

  //the worm cannot creep to anywhere if there is no other d^dagger
  int orbital = local_config.get_worm_orbital();
  int k = local_config.order(orbital); //number of segments in the worm orbital (including the worm!)
  if(k==1) return; 

  //catch the worm
  std::set<segment>::iterator old_worm_head = local_config.get_worm_head();
  double t_start = old_worm_head->t_start_;

  //find the new head for the worm
  int segment_nr = (int)(random()*k);
  std::set<segment>::iterator new_worm_head = local_config.get_segment_iterator(segment_nr, orbital);

  if(new_worm_head->c_start_ == WORM_COLOR) return; //caught the worm...
  
  //compute hybridization weight change
  double worm_tail = local_config.get_worm_tail();
  double hybridization_weight_change = hyb_config.hyb_weight_change_worm_creep(*new_worm_head, t_start, worm_tail, orbital);

  //Leo: compute the dissipation weight change //TODO
  double dissipation_weight_change = 1.0;
  //double dissipation_weight_change=ohmic_config.dissipation_weight_change(new_segment, orbital, true, local_config);
  
  //perform metropolis
  double weight_change = hybridization_weight_change * dissipation_weight_change;

  if(std::abs(weight_change)>random())
  {
    nacc[17]++;
    if(weight_change < 0) sign*=-1.;
    if(hybridization_weight_change < 0) std::cout << "***** need to flip sign *****" << std::endl;
    hyb_config.worm_creep(*new_worm_head, t_start, worm_tail, orbital);
    local_config.set_worm(orbital, new_worm_head->t_start_);
    std::swap(new_worm_head->c_start_, old_worm_head->c_start_);
    dissipation_weight_ratio = 1.0/dissipation_weight_change; //Leo: keep the weight change (of removal!) for measure_G //TODO

    //Leo: record the updated color and set color_updated to true
    //color = color_temp;
    //updated_colors[color]++;
    //ncolor[color]++;
    //ncolor_diff[color]++;  //Leo: + for insertion, - for removal 
    ////color_updated = true;

    //Leo: for test purpose, print out a lot of things...
    if(VERY_VERBOSE && sweeps<=debug_number) 
    { 
       debug_output(17, 1.0, hybridization_weight_change, dissipation_weight_change, 1.0); 
    }
  }
}


//Leo: experimental local color swapping in a randomly chosen orbital
//Notes:
//1. currently only two colors are supported, so this part should be extended in the future //TODO
//2. local weight will not be affected by this update
//3. color_tar is to be replaced by color_des
void hybridization::color_swap_update(int orbital)
{
  nprop[8]++; N_Z++;
  if(has_worm) return;

  if( local_config.order(orbital)==0 ) return; //no segment for swapping
  
  //Leo: choose the color in which we do the update
  //Now only two colors (red/1 and blue/0) are considered, but it can be easily changed
  size_t color_tar = (int)(random()*n_env);
  size_t color_des = (color_tar ? 0 : 1); //TODO: change this line for more colors

  int color_tar_mat_size = hyb_config.color_matrix_size(orbital, color_tar);
  int color_des_mat_size = hyb_config.color_matrix_size(orbital, color_des);
  if(color_tar_mat_size==0) return; //target color does not exist

  //compute hybridization weight change
  double times[2] = {random(), random()}; //pass two random numbers to get times = {t_start, t_end}
  double hybridization_weight_change = hyb_config.hyb_weight_change_swap(orbital, color_tar, color_des, times);

  //compute the dissipation weight change //TODO
  double dissipation_weight_change = 1.;//ohmic_config.color_flip_weight_change(orbital, color_1, color_2, local_config);

  //Leo: test!
  //double permutation_factor = (color_tar_mat_size*color_tar_mat_size) / ((color_tar_mat_size+1.0)*(color_tar_mat_size+1.0));
  double permutation_factor = (color_tar_mat_size*color_tar_mat_size) / ((color_des_mat_size+1.0)*(color_des_mat_size+1.0));
  //double permutation_factor = (color_tar_mat_size) / (color_tar_mat_size+1.0);
  
  //perform metropolis
  double weight_change = hybridization_weight_change * dissipation_weight_change * permutation_factor;

  //std::cout << "cdagger time: " << times[0] << std::endl;
  //std::cout << "c time: " << times[1] << std::endl;
  //std::cout << local_config << std::endl;
  //std::cout << hyb_config << std::endl;

  if(std::abs(weight_change)>random())
  {
    nacc[8]++;
    if(weight_change < 0) sign*=-1.;
    local_config.swap_color(orbital, color_tar, color_des, times); 
    hyb_config.swap_color(orbital, color_tar, color_des, times);
    dissipation_weight_ratio = dissipation_weight_change; //TODO: check if necessary

    //Leo: record the updated color 
    ncolor[color_tar]++;
    ncolor[color_des]++;
    ncolor_diff[color_tar]--;  //Leo: + for insertion, - for removal 
    ncolor_diff[color_des]++;  //Leo: + for insertion, - for removal 
    updated_colors[color_tar]++; 
    updated_colors[color_des]++; 

    //Leo: for test purpose, print out a lot of things...
    if(VERY_VERBOSE && sweeps<=debug_number) 
    { 
       debug_output(8, 1, hybridization_weight_change, dissipation_weight_change, permutation_factor); 
    }
  }
}
