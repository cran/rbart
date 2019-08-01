//     ambrt.h: Additive mean BART model class definition.
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


#ifndef GUARD_ambrt_h
#define GUARD_ambrt_h

#include "tree.h"
#include "treefuns.h"
#include "dinfo.h"
#include "mbrt.h"


class ambrt : public mbrt 
{
public:
   //--------------------
   //classes
   // tprior and mcmcinfo are same as in brt
   // cinfo same as in mbrt
   //--------------------
   //constructors/destructors
   ambrt(): mbrt(),m(200),mb(m),notjmus(m),divec(m) {}
   ambrt(size_t im): mbrt(),m(im),mb(m),notjmus(m),divec(m) {}
   virtual ~ambrt() {
      if(!notjmus.empty()) {
         for(size_t j=0;j<m;j++) notjmus[j].clear();
         notjmus.clear();
         for(size_t j=0;j<m;j++) delete divec[j];
      }
   }

   //--------------------
   //methods
   void draw(rn& gen);
   void adapt();
   void setci(double tau, double* sigma) { ci.tau=tau; ci.sigma=sigma; for(size_t j=0;j<m;j++) mb[j].setci(tau,sigma); }
   void settc(int tc) { this->tc = tc; for(size_t j=0;j<m;j++) mb[j].settc(tc); }
   void setxi(xinfo *xi) { this->xi=xi; for(size_t j=0;j<m;j++) mb[j].setxi(xi); }
   void setdata(dinfo *di);
   void settp(double alpha, double beta) { tp.alpha=alpha;tp.beta=beta; for(size_t j=0;j<m;j++) mb[j].settp(alpha,beta); }
   tree::tree_p gettree(size_t i) { return &mb[i].t; }
   void setmi(double pbd, double pb, size_t minperbot, bool dopert, double pertalpha, double pchgv, std::vector<std::vector<double> >* chgv)
             { mi.pbd=pbd; mi.pb=pb; mi.minperbot=minperbot; mi.dopert=dopert;
               mi.pertalpha=pertalpha; mi.pchgv=pchgv; mi.corv=chgv; 
               for(size_t j=0;j<m;j++) mb[j].setmi(pbd,pb,minperbot,dopert,pertalpha,pchgv,chgv); }
   void setstats(bool dostats) { mi.dostats=dostats; for(size_t j=0;j<m;j++) mb[j].setstats(dostats); if(dostats) mi.varcount=new unsigned int[xi->size()]; }
   void pr();
   // drawnodetheta, lm, add_observation_to_suff and newsinfo/newsinfovec unused here.

   //--------------------
   //data
   //--------------------------------------------------
   //stuff that maybe should be protected
   // std::vector<mbrt> mb;  // the vector of individual mu trees for sum representation
protected:
   //--------------------
   //model information
   size_t m;  //number of trees in sum representation
   std::vector<mbrt> mb;  // the vector of individual mu trees for sum representation
   //--------------------
   //data
   std::vector<std::vector<double> > notjmus;
   std::vector<dinfo*> divec;
   //--------------------
   //mcmc info
   //--------------------
   //methods
   virtual void local_setf(diterator& diter);  //set the vector of predicted values
   virtual void local_setr(diterator& diter);  //set the vector of residuals
   virtual void local_predict(diterator& diter); // predict y at the (npred x p) settings *di.x
   virtual void local_savetree(size_t iter, int beg, int end, std::vector<int>& nn, std::vector<std::vector<int> >& id, std::vector<std::vector<int> >& v,
                  std::vector<std::vector<int> >& c, std::vector<std::vector<double> >& theta);
   virtual void local_loadtree(size_t iter, int beg, int end, std::vector<int>& nn, std::vector<std::vector<int> >& id, std::vector<std::vector<int> >& v,
                  std::vector<std::vector<int> >& c, std::vector<std::vector<double> >& theta);

};


#endif
