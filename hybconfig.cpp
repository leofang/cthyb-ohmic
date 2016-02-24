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

#include"hybconfig.hpp"

hybridization_configuration::hybridization_configuration(const alps::params &p):
//  Delta((int)(p["N_ENV"]), p),
  n_env_(p.defined("N_ENV")?(int)p["N_ENV"]:1),
  hybmat_((int)(p["N_ORBITALS"]), std::vector<hybmatrix>(n_env_, p))
{
  //n_env_ = (int)p["N_ENV"]|1;
  //if(p.defined("VERY_VERBOSE") && p["VERY_VERBOSE"].cast<bool>()==true)
  //{
  //   std::cout << "hybridization_configuration is initializing..." << std::endl;
  //   std::cout << "N_ENV = " << n_env_ << std::endl;
  //}
  //hybmat_((int)(p["N_ORBITALS"]), std::vector<hybmatrix>(n_env_, p));
  std::cout << "Initialize Delta for " << n_env_ << " color(s)..." << std::endl << std::endl;
  initialize_Delta(p);
}


// initialize Delta from files
void hybridization_configuration::initialize_Delta(const alps::params &p)
{
//    std::cout << "n_env_ = " << n_env_ << std::endl;
//    std::cout << "Initialize Delta..." << std::endl;
    // read multiple hybridization files upon request
    if(p.defined("N_ENV") && n_env_>1) //N_ENV>=2
    {
       // detect if hybridization input files are given
       for(int i=0; i< n_env_; i++)
       {
         std::stringstream temp_stream; // stringstream used for the conversion
         temp_stream << i;              // add the value of i to the characters in the stream
         std::string temp_filename = "DELTA" + temp_stream.str();

         if(!p.defined(temp_filename))
              throw  std::runtime_error("Parameter " + temp_filename + \
              " (filename for hybridization function of color " + temp_stream.str() + ") is not provided. Abort!");
       }
       //prepare parameter set for each color by pretending other colors do not exist
       //and call "DELTAi" as "DELTA"
       for(int i=0; i< n_env_; i++)
       {
         alps::params temp_params = p; // temp_params is a "buffer"
         for(int j=0; j< n_env_; j++)
         {
             std::stringstream temp_stream;
             temp_stream << j;
             std::string temp_string = "DELTA" + temp_stream.str();
             if (j==i) { temp_params["DELTA"] = p[temp_string].cast<std::string>(); }
             temp_params.erase(temp_string);
         }
         temp_params.erase("N_ENV"); //Just in case something is wrong when reading hybridization file...
         
	 //Leo: depreciated
         //std::cout << "the parameter set has " << temp_params.size() << " parameters." << std::endl;
         //std::cout << temp_params << std::endl;
         std::cout << "reading the delta file " << temp_params["DELTA"].cast<std::string>() << " ..." << std::endl;

         //use temp_params to initialize hybfun constructor
         Delta.push_back(temp_params);
       }
       //check the size of Delta
       if (Delta.size() != n_env_)
         throw std::runtime_error("Delta should have size N_ENV but it does not! Abort!");
    }
    else if ( (p.defined("N_ENV") && n_env_==1) || !p.defined("N_ENV") ) //N_ENV=1
    {
        alps::params temp_params = p;
        if(!p.defined("DELTA") && !p.defined("DELTA0"))
          throw  std::runtime_error("Parameter DELTA or DELTA0 (filename for hybridization function) is not provided. Abort!");
        if(!p.defined("DELTA")) // read from DELTA0
        {
           temp_params["DELTA"]=temp_params["DELTA0"].cast<std::string>();
           temp_params.erase("DELTA0");
           if(p.defined("N_ENV")) { temp_params.erase("N_ENV"); } //Just in case
        }
 
        //Leo: depreciated
        //std::cout << "the parameter set has " << temp_params.size() << " parameters." << std::endl;
        //std::cout << temp_params << std::endl;
        std::cout << "reading the delta file " << temp_params["DELTA"].cast<std::string>() << " ..." << std::endl;

        //use temp_params to initialize hybfun constructor
        Delta.push_back(temp_params);
    }
    else //N_ENV is something nonsense
    {
        throw  std::invalid_argument("N_ENV is invalid, abort."); //TODO: move this to sanity_check
    }
}


void hybridization_configuration::dump() 
{
    for (int i=0; i<hybmat_.size(); i++) //N_orbital
    {
        for (int j=0; j<hybmat_[i].size(); j++) //Leo: N_ENV
	   std::cerr << "Weight for orbital " << i << " and reservoir " << j << " : " << hybmat_[i][j].full_weight() << std::endl;
    }
}


void hybridization_configuration::rebuild() 
{
  for (int i=0; i<hybmat_.size(); i++) //N_orbital
  {
     for (int j=0; j<hybmat_[i].size(); j++) //Leo: N_ENV
     {
    	hybmat_[i][j].rebuild_hyb_matrix(i, Delta[j]);
     }
  }
}


//Leo: for debug purpose
void hybridization_configuration::rebuild_ordered() 
{
  for (int i=0; i<hybmat_.size(); i++) //N_orbital
  {
     for (int j=0; j<hybmat_[i].size(); j++) //Leo: N_ENV
     {
    	hybmat_[i][j].rebuild_ordered_hyb_matrix(i, Delta[j]);
     }
  }
}


void hybridization_configuration::rebuild(int orbital) 
{
   for (int j=0; j<hybmat_[orbital].size(); j++) //Leo: N_ENV
       hybmat_[orbital][j].rebuild_hyb_matrix(orbital, Delta[j]);
}


void hybridization_configuration::rebuild(std::vector<int> orbital) 
{
  for (int i=0;i<orbital.size();i++)
  {
     for (int j=0; j<hybmat_[orbital[i]].size(); j++) //Leo: N_ENV
        hybmat_[orbital[i]][j].rebuild_hyb_matrix(orbital[i], Delta[j]);
  }
}


//Leo: since the size of reservoir matrices coupled to the same orbital are correlated (n_L+n_R=n),
//     there should be a consistency check after each update to gaurantee we're doing it correctly
int hybridization_configuration::total_color_matrix_size(int orbital) const
{
   int total_size=0;
   for (int i=0; i<n_env_; i++)
   {
       total_size+=hybmat_[orbital][i].size();
   }
   return total_size;
}


//Leo: when N_ENV>1, it is possible to have the sign problem due to operator time ordering,
//     therefore, this function returns the overall sign by multiplying those of each matrix.
//     This sign should agree with that of hybridization_configuration::full_weight!
//     TOTALLY EXPERIMENTAL! BE CAREFUL!
int hybridization_configuration::overall_color_matrix_sign() const
{
   int sign=1;
   for (int i=0; i<hybmat_.size(); i++) //N_orbital
   {
      for (int j=0; j<hybmat_[i].size(); j++) //N_ENV
      {
     	sign*=hybmat_[i][j].sign();
      }
      //sign*=(total_color_matrix_size(i)%2?-1:1);
   }
   return sign;
}


void hybridization_configuration::haunt_missing_sign(int orbital)
{
    int orbital_order = total_color_matrix_size(orbital);
    if(orbital_order==0) return;

    std::vector<int> matrix_size(n_env_, 0);
    for(int i=0; i<n_env_; i++)
       matrix_size[i] = hybmat_[orbital][i].size();
  
    //select the case in which n is even, and both n_R and n_L are odd
    if( !(orbital_order%2) && (matrix_size[0]%2) )
        std::cout << "SPECIAL CASE: n is even, and both n_R and n_L are odd!" << std::endl;
}




double hybridization_configuration::hyb_weight_change_insert(const segment &new_segment, int orbital, size_t color)
{
//  if(new_segment.c_start_ == new_segment.c_end_) //Leo: the colors of both ends of a (anti-)segment should be the same!
//  {int color = new_segment.c_start_; }
//  else 
//  {throw std::runtime_error("hybridization_configuration: the colors to be inserted are not the same!")}

  return hybmat_[orbital][color].hyb_weight_change_insert(new_segment, orbital, Delta[color]); //hand this off to the determinant matrix
}


void hybridization_configuration::insert_segment(const segment &new_segment, int orbital, size_t color)
{
  //std::cout<<clmagenta<<"before insert recompute "<<cblack<<std::endl;
  //if(hybmat_[orbital].size()>0) hybmat_[orbital].rebuild_hyb_matrix(orbital, Delta);
  //std::cout<<clmagenta<<"done before insert recompute "<<cblack<<std::endl;
  hybmat_[orbital][color].insert_segment(new_segment, orbital); //hand this off to the determinant matrix
  //std::cout<<clmagenta<<"after insert recompute "<<cblack<<std::endl;
  //hybmat_[orbital].rebuild_hyb_matrix(orbital, Delta);
  //std::cout<<clmagenta<<"done after insert recompute "<<cblack<<std::endl;
}


double hybridization_configuration::hyb_weight_change_remove(const segment &new_segment, int orbital, size_t color)
{
  return hybmat_[orbital][color].hyb_weight_change_remove(new_segment, orbital, Delta[color]); //hand this off to the determinant matrix
}


void hybridization_configuration::remove_segment(const segment &new_segment, int orbital, size_t color)
{
  //std::cout<<clmagenta<<"before remove recompute "<<cblack<<std::endl;
  //hybmat_[orbital].rebuild_hyb_matrix(orbital, Delta);
  //std::cout<<clmagenta<<"done before remove recompute "<<cblack<<std::endl;
  hybmat_[orbital][color].remove_segment(new_segment, orbital); //hand this off to the determinant matrix
  //std::cout<<clmagenta<<"after remove recompute "<<cblack<<std::endl;
  //if(hybmat_[orbital].size()>0) hybmat_[orbital].rebuild_hyb_matrix(orbital, Delta);
  //std::cout<<clmagenta<<"done after remove recompute "<<cblack<<std::endl;
}


void hybridization_configuration::remove_antisegment(const segment &new_antisegment, int orbital, size_t color)
{
  hybmat_[orbital][color].remove_segment(new_antisegment, orbital); //hand this off to the determinant matrix
  //std::cout<<clmagenta<<"after as remove recompute "<<cblack<<std::endl;
  if(hybmat_[orbital][color].size()>0) hybmat_[orbital][color].rebuild_hyb_matrix(orbital, Delta[color]);
  //std::cout<<clmagenta<<"done after as remove recompute "<<cblack<<std::endl;
}


void hybridization_configuration::insert_antisegment(const segment &new_antisegment, int orbital, size_t color)
{
  hybmat_[orbital][color].insert_segment(new_antisegment, orbital); //hand this off to the determinant matrix
  //std::cout<<clmagenta<<"after as insert recompute "<<cblack<<std::endl;
  //hybmat_[orbital].rebuild_hyb_matrix(orbital, Delta);
  //std::cout<<clmagenta<<"done after as insert recompute "<<cblack<<std::endl;
}


void hybridization_configuration::measure_G(std::vector<std::vector<double> > &G, std::vector<std::vector<double> > &F, const std::vector<std::map<double,double> > &F_prefactor, double sign, double dissipation_weight_ratio) const
{
  for(std::size_t orbital=0; orbital<hybmat_.size(); ++orbital)
  {
     //Leo: not sure if the sign here is correct for N_ENV=2, need to check!
     for(int color=0; color<hybmat_[orbital].size(); color++) //Leo: not sure if color here works...
     {
    //   int color = 1; //(n_env_==1?0:1);
    //   if(color!=0 && overall_color_matrix_sign()<0)  sign*=-1; //Leo: test!
    //   if(overall_color_matrix_sign()<0)
    //        hybmat_[orbital][color].measure_G(G[orbital], F[orbital], F_prefactor[orbital],\
    //                                           sign*hybmat_[orbital][color].sign(), dissipation_weight_ratio);
    //   else
            hybmat_[orbital][color].measure_G(G[orbital], F[orbital], F_prefactor[orbital],\
                                               sign, dissipation_weight_ratio);
                                         // (sign<0?sign+color:sign), dissipation_weight_ratio);
     }
  }
}


void hybridization_configuration::measure_Gw(std::vector<std::vector<double> > &Gwr, std::vector<std::vector<double> > &Gwi, std::vector<std::vector<double> > &Fwr, std::vector<std::vector<double> > &Fwi, const std::vector<std::map<double,double> > &F_prefactor, double sign) const
{
  for(std::size_t orbital=0;orbital<hybmat_.size();++orbital)
  { //Leo: not sure if color here works...
    for(size_t color=0; color<hybmat_[orbital].size(); color++) //Leo: not sure if color here works...
    {
       hybmat_[orbital][color].measure_Gw(Gwr[orbital], Gwi[orbital], Fwr[orbital], Fwi[orbital], F_prefactor[orbital], sign);
    }
  }
}


void hybridization_configuration::measure_G2w(std::vector<std::vector<std::complex<double> > >&G2w, std::vector<std::vector<std::complex<double> > > &F2w, int N_w2, int N_w_aux, const std::vector<std::map<double,double> > &F_prefactor) const
{
  for(std::size_t orbital=0;orbital<hybmat_.size();++orbital)
  {
    for(size_t color=0; color<hybmat_[orbital].size(); color++) //Leo: not sure if color here works...
    {
       hybmat_[orbital][color].measure_G2w(G2w[orbital], F2w[orbital], N_w2, N_w_aux, F_prefactor[orbital]);
    }
  }
}


void hybridization_configuration::measure_Gl(std::vector<std::vector<double> > &Gl, std::vector<std::vector<double> > &Fl, const std::vector<std::map<double,double> > &F_prefactor, double sign) const
{
  for(std::size_t orbital=0;orbital<hybmat_.size();++orbital)
  {
     for(size_t color=0; color<hybmat_[orbital].size(); color++) //Leo: not sure if color here works...
     {
        hybmat_[orbital][color].measure_Gl(Gl[orbital], Fl[orbital], F_prefactor[orbital], sign);
     }
  }
}


double hybridization_configuration::full_weight() const
{
  double weight=1.;
  for(std::size_t orbital=0;orbital<hybmat_.size();++orbital)
  {
    for(size_t color=0; color<hybmat_[orbital].size(); color++) //Leo: not sure if color here works...
    {
       weight*=hybmat_[orbital][color].full_weight();
    }
  }
  return weight;
}


std::ostream &operator<<(std::ostream &os, const hybridization_configuration &hyb_config)
{
  for(std::size_t i=0;i<hyb_config.hybmat_.size();++i)
  {
    for(size_t j=0; j<hyb_config.hybmat_[i].size(); j++) //Leo: not sure if j=color here works...
    {
        //Leo: get rid of text colors
	//os<<cblue<<"------- "<<"orbital: "<<i<< ", color: " << j << " ------"<<cblack<<std::endl;
        os << "------- orbital: " << i << ", color: " << j << " -------" << std::endl;
        os << hyb_config.hybmat_[i][j];
        if(i!=hyb_config.hybmat_.size()-1)  os << std::endl;
    }
  }
  return os;
}
