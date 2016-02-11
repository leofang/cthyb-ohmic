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
#ifndef OHMIC_CONFIG_HPP
#define OHMIC_CONFIG_HPP

#include <alps/ngs/params.hpp>
#include "hybint.hpp"
#include "hybsegment.hpp"
#include "hyblocal.hpp"
#include <set>
#include <vector>

//This is a friend class of local_configuration that handles the ohmic
//dissipative environment. It reads the segment information and then 
//compute the phase correlators which show up in the dissipation weight. 
//
//Note that unlike the local and hyb classes, this class stores *nothing* 
//about the dissipative environment, since everything (mainly weight 
//change) can be computed on the fly. The purpose of this class is just
//to make the code more structured and readable. In particular, one can
//think of this class as not existing at all when Dissipation=0 is set.

class local_configuration; //forward declaration

class dissipation_configuration{
public:
  dissipation_configuration(const alps::params &p);
  //friend std::ostream &operator<<(std::ostream &os, const dissipation_configuration &ohmic_conf);

  double dissipation_weight_change(const segment &seg, int orbital, bool insert, const local_configuration &local_config) const;
  inline double phase_correlator_J(double tau) const;

private:
  int n_env_; //number of reservoirs
  double beta_;
  int n_orbitals_;
  bool dissipation_; //turn on dissipation or not
  double r_;         //dissipation strength
  double C0_;        //capacitance ratio of the 0-th capacitor to total capacitance
  double wc_;        //cutoff frequency
  double kappa_;     //kappa=1/(beta*wc)
  double gamma_;     //Euler constant
  std::vector< std::vector<double> > dissipation_coeff_;
};

//std::ostream &operator<<(std::ostream &os, const dissipation_configuration &ohmic_conf);



#ifndef COLORS
#define COLORS
#define cblack "\033[22;30m"
#define cred "\033[22;31m"
#define cgreen "\033[22;32m"
#define cbrown "\033[22;33m"
#define cblue "\033[22;34m"
#define cmagenta "\033[22;35m"
#define ccyan "\033[22;36m"
#define cgray "\033[22;37m"
#define cdgray "\033[01;30m"
#define clred "\033[01;31m"
#define clgreen "\033[01;32m"
#define clyellow "\033[01;33m"
#define clblue "\033[01;34m"
#define clmagenta "\033[01;35m"
#define clcyan "\033[01;36m"
#define cwhite "\033[01;37m"
#endif

#endif
