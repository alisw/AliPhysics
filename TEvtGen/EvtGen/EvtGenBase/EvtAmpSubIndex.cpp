//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2002      Caltech
//
// Module: EvtAmpSubIndex.cc
//
// Description: Class to manipulate the amplitudes in the decays.
//
// Modification history:
//
//    RYD     Nov 22, 2002         Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtAmpIndex.hh"
#include "EvtGenBase/EvtAmpSubIndex.hh"
#include <vector>
using std::vector;


EvtAmpSubIndex::EvtAmpSubIndex(EvtAmpIndex* ind,std::vector<int> sub):
  _ind(ind),
  _sub(sub),
  _size(sub.size()),
  _nstate(sub.size())
{
  int i;
  
  for(i=0;i<_size;i++) {
    if (i==0){
      _nstate[i]=1;
    }
    else{
      _nstate[i]=_nstate[i-1]*_ind->_ind[sub[i-1]];
    }
  }
}


int EvtAmpSubIndex::index(){

  int i;
  int ind=0;

  for(i=0;i<_size;i++) {
    ind+=_ind->_state[_ind->_ind[i]]*_nstate[i];
  }

  return ind;

}











