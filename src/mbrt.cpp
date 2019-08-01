//     mbrt.cpp: Mean tree BART model class methods.
//     Copyright (C) 2012-2016 Matthew T. Pratola, Robert E. McCulloch and Hugh A. Chipman
//
//     This file is part of BART.
//
//     BART is free software: you can redistribute it and/or modify
//     it under the terms of the GNU Affero General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
//
//     BART is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU Affero General Public License for more details.
//
//     You should have received a copy of the GNU Affero General Public License
//     along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//     Author contact information
//     Matthew T. Pratola: mpratola@gmail.com
//     Robert E. McCulloch: robert.e.mculloch@gmail.com
//     Hugh A. Chipman: hughchipman@gmail.com


#include "mbrt.h"
//#include "brtfuns.h"
#include <iostream>
#include <map>
#include <vector>

#ifndef NotInR
#include <Rcpp.h>
#define COUT Rcpp::Rcout 
#else
#define COUT std::cout  << "std::cout "
#endif

using std::cout;
using std::endl;


//--------------------------------------------------
//a single iteration of the MCMC for brt model
void mbrt::draw(rn& gen)
{
   //All the usual steps
   brt::draw(gen);

   // Update the in-sample predicted vector
   setf();

   // Update the in-sample residual vector
   setr();

}
//--------------------------------------------------
//draw theta for a single bottom node for the brt model
double mbrt::drawnodetheta(sinfo& si, rn& gen)
{
   msinfo& msi=static_cast<msinfo&>(si);
   double muhat = msi.sumwy/msi.sumw;
   double a = 1.0/(ci.tau*ci.tau);
   return (msi.sumw*muhat)/(a+msi.sumw) + gen.normal()/sqrt(a+msi.sumw);
}
//--------------------------------------------------
//lm: log of integrated likelihood, depends on prior and suff stats
double mbrt::lm(sinfo& si)
{
   msinfo& msi=static_cast<msinfo&>(si);
   double t2 =ci.tau*ci.tau;
   double k = msi.sumw*t2+1;
   return -.5*log(k)+.5*msi.sumwy*msi.sumwy*t2/k;
}
//--------------------------------------------------
//Add in an observation, this has to be changed for every model.
//Note that this may well depend on information in brt with our leading example
//being double *sigma in cinfo for the case of e~N(0,sigma_i^2).
// Note that we are using the training data and the brt object knows the training data
//     so all we need to specify is the row of the data (argument size_t i).
void mbrt::add_observation_to_suff(diterator& diter, sinfo& si)
{
   msinfo& msi=static_cast<msinfo&>(si);
   double w;
   w=1.0/(ci.sigma[*diter]*ci.sigma[*diter]);
   msi.n+=1;
   msi.sumw+=w;
   msi.sumwy+=w*diter.gety();
}
//--------------------------------------------------
//pr for brt
void mbrt::pr()
{
   COUT << "***** mbrt object:\n";
   COUT << "Conditioning info:" << endl;
   COUT << "   mean:   tau=" << ci.tau << endl;
   if(!ci.sigma)
     COUT << "         sigma=[]" << endl;
   else
     COUT << "         sigma=[" << ci.sigma[0] << ",...," << ci.sigma[di->n-1] << "]" << endl;
   brt::pr();
}
