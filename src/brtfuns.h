//     brtfuns.h: Base BART model class help functions header file.
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


#ifndef GUARD_brtfuns_h
#define GUARD_brtfuns_h

#include <iostream>
#include "tree.h"
#include "treefuns.h"
#include "brt.h"

using std::cout;
using std::endl;

//--------------------------------------------------
//make xinfo = cutpoints
void makexinfo(size_t p, size_t n, double *x, xinfo& xi, size_t nc);
//--------------------------------------------------
//compute prob of a birth, goodbots will contain all the good bottom nodes
double getpb(tree& t, xinfo& xi, double pipb, tree::npv& goodbots);
//--------------------------------------------------
//bprop: function to generate birth proposal
void bprop(tree& x, xinfo& xi, brt::tprior& tp, double pb, tree::npv& goodbots, double& PBx, tree::tree_p& nx, size_t& v, size_t& c, double& pr, rn& gen);
//--------------------------------------------------
// death proposal
void dprop(tree& x, xinfo& xi, brt::tprior& tp, double pb, tree::npv& goodbots, double& PBx, tree::tree_p& nx, double& pr, rn& gen);
//--------------------------------------------------
//get prob a node grows, 0 if no good vars, else alpha/(1+d)^beta
double pgrow(tree::tree_p n, xinfo& xi, brt::tprior& tp);
//--------------------------------------------------
//calculate beginning and end points of data vector to be accessed in parallel computations
void calcbegend(int n, int my_rank, int thread_count, int* beg, int* end);

//--------------------------------------------------
// Functions to support change-of-variable proposal
//--------------------------------------------------
// update the correlation matrix for chgv move taking into account that not all
// variables may be eligible at pertnode.
void updatecormat(tree::tree_p pertnode, xinfo& xi, std::vector<std::vector<double> >& chgv);
//--------------------------------------------------
// renormalize the correlation matrix so that the probability of row sums to 1.
void normchgvrow(size_t row, std::vector<std::vector<double> >& chgv);
//--------------------------------------------------
// randomly choose a new variable to transition to from oldv
size_t getchgv(size_t oldv, std::vector<std::vector<double> >& chgv, rn& gen);


//--------------------------------------------------
// Functions to support rotate proposal
//--------------------------------------------------
//setup the initial right rotation
void rotright(tree::tree_p n);
//--------------------------------------------------
//setup the initial left rotation
void rotleft(tree::tree_p n);
//--------------------------------------------------
//eliminate immediate dead ends from the rotate
void reduceleft(tree::tree_p n, size_t v, size_t c);
//--------------------------------------------------
//eliminate immediate dead ends from the rotate
void reduceright(tree::tree_p n, size_t v, size_t c);
//--------------------------------------------------
//split tree along variable v at cutpoint c retaining only 
//part of the tree that is ``left'' of this v,c rule
void splitleft(tree::tree_p t, size_t v, size_t c);
//--------------------------------------------------
//split tree along variable v at cutpoint c retaining only 
//part of the tree that is ``right'' of this v,c rule
void splitright(tree::tree_p t, size_t v, size_t c);
//--------------------------------------------------
//does an actual merge (randomly chosen) 
bool merge(tree::tree_p tl, tree::tree_p tr, tree::tree_p t, size_t v, size_t c, rn& gen);
//--------------------------------------------------
// only to get nways, not to actually do the merge.
bool mergecount(tree::tree_p tl, tree::tree_p tr, size_t v, size_t c, int* nways);
//--------------------------------------------------
// End of functions to support rotate proposal
//--------------------------------------------------


#endif
