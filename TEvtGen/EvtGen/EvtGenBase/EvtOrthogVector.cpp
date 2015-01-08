//--------------------------------------------------------------------------
// 
// Environment: 
// This software is part of the EvtGen package developed jointly 
// for the BaBar and CLEO collaborations.  If you use all or part 
// of it, please give an appropriate acknowledgement.
// 
// Copyright Information: See EvtGen/COPYRIGHT 
// Copyright (C) 2000 Caltech, LLNL
// 
// Module: EvtGen/EvtOrthogVector.hh
// 
// Description:
// 
// Modification history: 
//
// Lange August 11, 2000 Created
//
//------------------------------------------------------------------------
#include "EvtGenBase/EvtPatches.hh"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "EvtGenBase/EvtOrthogVector.hh"
using std::fstream;

EvtOrthogVector::EvtOrthogVector(int n, std::vector<double> *vectors){

  _dimen=n;
  _holder.resize(n);

  std::vector<int> temp;
  
  int i;
  for (i=0;i<n;i++) {
    _orthogVector.push_back(0.);
    temp.push_back(i);
  }

  findOrthog(_dimen,temp, vectors);

}

EvtOrthogVector::~EvtOrthogVector(){
}

void EvtOrthogVector::findOrthog(int dim, std::vector<int> invect, 
			    std::vector<double> *vectors) {


  if ( dim==2 ) {
    _holder[0]=invect[0];
    _holder[1]=invect[1];
    int sign=findEvenOddSwaps();
    {
      double addition=1;
      int i;
      for (i=1; i<_dimen; i++){
	addition*=vectors[i-1][_holder[i]];
      }
      addition*=sign;
      _orthogVector[_holder[0]]+=addition;
    }
    
    _holder[0]=invect[1];
    _holder[1]=invect[0];
    
    {
      double addition=1;
      int i;
      for (i=1; i<_dimen; i++){
	addition*=vectors[i-1][_holder[i]];
      }
      addition*=sign;
      _orthogVector[_holder[0]]-=addition;
    }
    
    return;
  }
  else{
    std::vector<int> temp((2*dim));

    int i;
    for (i=0; i<dim; i++) temp[i]=invect[i];
    for (i=0; i<dim; i++) temp[i+dim]=invect[i];

    for (i=0; i<dim; i++) {
      _holder[dim-1]=temp[dim-1+i];
      std::vector<int> tempDim((dim-1));

      int j;
      for (j=0; j<(dim-1); j++) tempDim[j]=temp[j+i];
      findOrthog(dim-1, tempDim, vectors); 
    }
  }
 
  return;
}

int EvtOrthogVector::findEvenOddSwaps() {

  std::vector<int> temp(_dimen);

  int i,j,nSwap;
  for (i=0; i<_dimen; i++) temp[i]=_holder[i];

  nSwap=0;
  for (i=0; i<(_dimen-1); i++) {
    for (j=i+1; j<_dimen; j++) {

      if ( temp[i]>temp[j] ) {
	int duh=temp[j];
	temp[j]=temp[i];
	temp[i]=duh;
	nSwap+=1;
      }
    }
  }
  nSwap-= (nSwap/2)*2;

  if ( nSwap ) return -1;
  
  return 1;

}
