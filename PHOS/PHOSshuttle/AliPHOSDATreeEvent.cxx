/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/* $Id$ */

// --
// --
// Implementation for TTree output in PHOS DA
// for calibrating energy by pi0 and MIP.
// --
// -- Author: Hisayuki Torii (Hiroshima Univ.)
// --

#include <iostream>
#include "AliPHOSDATreeEvent.h"
#include "AliPHOSDATreeCluster.h"
#include "AliPHOSDATreeDigit.h"

ClassImp(AliPHOSDATreeEvent)
//------------------------------------------------------------------------
AliPHOSDATreeEvent::AliPHOSDATreeEvent(const AliPHOSDATreeEvent& evt)
: fTime(0),
  fNDigits(0),
  fDigits(0),
  fNClusters(0),
  fClusters(0)
{
  // Copy Constructor

  fNDigits = evt.fNDigits;
  if( fNDigits > 0 ){
    fDigits = new AliPHOSDATreeDigit[fNDigits];
    int ndigits = fNDigits;
    while( ndigits-- ){
      fDigits[ndigits] = evt.fDigits[ndigits];
    }
  } else {
    fDigits = 0;
  }
  //
  fNClusters = evt.fNClusters;
  if( fNClusters > 0 ){
    fClusters = new AliPHOSDATreeCluster[fNClusters];
    int nclusters = fNClusters;
    while( nclusters-- ){
      fClusters[nclusters] = evt.fClusters[nclusters];
    }
  } else {
    fClusters = 0;
  }
}

AliPHOSDATreeEvent::~AliPHOSDATreeEvent(){
  delete[] fDigits;
  delete[] fClusters;
};

//------------------------------------------------------------------------
AliPHOSDATreeEvent& AliPHOSDATreeEvent::operator=(const AliPHOSDATreeEvent& evt){
  // Copy Operator

  if (this != &evt) {
    if( fNDigits > 0 ) delete[] fDigits;
    fNDigits = evt.fNDigits;
    if( fNDigits > 0 ){
      fDigits = new AliPHOSDATreeDigit[fNDigits];
      int ndigits = fNDigits;
      while( ndigits-- ){
	fDigits[ndigits] = evt.fDigits[ndigits];
      }
    } else {
      fDigits = 0;
    }
    //
    if( fNClusters > 0 ) delete[] fClusters;
    fNClusters = evt.fNClusters;
    if( fNClusters > 0 ){
      fClusters = new AliPHOSDATreeCluster[fNClusters];
      int nclusters = fNClusters;
      while( nclusters-- ){
	fClusters[nclusters] = evt.fClusters[nclusters];
      }
    } else {
      fClusters = 0;
    }
  }
  return *this;
}
//------------------------------------------------------------------------
bool AliPHOSDATreeEvent::Fill(float energy,int row,int col){
  // Fill new digit information

  AliPHOSDATreeDigit digit(energy,row,col);
  return Fill(digit);
}
//------------------------------------------------------------------------
bool AliPHOSDATreeEvent::Fill(AliPHOSDATreeDigit& digit){
  // Fill new digit information

  AliPHOSDATreeDigit* newDigits = new AliPHOSDATreeDigit[fNDigits+1];
  int ndigits = fNDigits;
  newDigits[ndigits] = digit;
  while( ndigits-- ){
    newDigits[ndigits] = fDigits[ndigits];
  }
  if( fNDigits>0 ) delete[] fDigits;
  fDigits = newDigits;
  fNDigits++;
  return true;
}
//------------------------------------------------------------------------
bool AliPHOSDATreeEvent::Clusterize(AliPHOSDATreeDigit& digit){
  //Input digit information into clustering

  bool status = false;
  int nclusters = fNClusters;
  while( nclusters-- && !status ){
    if( fClusters[nclusters].IsNeighbor(digit) ){
      status = fClusters[nclusters].Append(digit);
    }
  }
  //std::cout<<" DEBUGDEBUG:: status = "<<status<<" fNClusters="<<fNClusters<<std::endl;
  if( !status ){
    AliPHOSDATreeCluster* newClusters = new AliPHOSDATreeCluster[fNClusters+1];
    nclusters = fNClusters;
    status = newClusters[nclusters].Append(digit);
    while( nclusters-- ){
      newClusters[nclusters] = fClusters[nclusters];
    }
    if( fNClusters>0 ) delete[] fClusters;
    fClusters = newClusters;
    fNClusters++;
  }
  return status;
}
//------------------------------------------------------------------------
bool AliPHOSDATreeEvent::ExecuteClustering(){
  // Run clustering algorithm.
  // With the filled digit information
  //

  // Loop over all digit information
  if( fNDigits <= 0 ) return true;
  bool status = true;
  int ndigits = fNDigits;
  while( ndigits-- && status ){
    status &= Clusterize(fDigits[ndigits]);
  }
  if(! status ) return status;

  // Cluster merging
  int nclusters0 = fNClusters;
  while( nclusters0-- ){
    AliPHOSDATreeCluster& cluster0 = fClusters[nclusters0];
    int nclusters1 = nclusters0;
    while( nclusters1-- ){
      AliPHOSDATreeCluster& cluster1 = fClusters[nclusters1];
      if( cluster0.IsNeighbor(cluster1) ){
	cluster0.Append(cluster1);
	cluster1.Reset();
      }
    }
  }
  int nclustersnonzero = 0;
  int nclusters = fNClusters;
  while( nclusters-- ){
    fClusters[nclusters].CalculateProperty();
    if( fClusters[nclusters].GetEnergy()>0 ) nclustersnonzero++;
  }
  //std::cout<<" nclustersnonzero = "<<nclustersnonzero<<std::endl;
  if( nclustersnonzero > 0 ){
    AliPHOSDATreeCluster* newClusters = new AliPHOSDATreeCluster[nclustersnonzero];
    nclusters = fNClusters;
    fNClusters = nclustersnonzero;
    while( nclusters-- ){
      if( fClusters[nclusters].GetEnergy()>0 ){
	newClusters[--nclustersnonzero] = fClusters[nclusters];
	//std::cout<<"     : "<<nclusters<<" --> "<<nclustersnonzero<<std::endl;
      }
    }
    if( fNClusters>0 ) delete[] fClusters;
    fClusters = newClusters;
  } else {
    if( fNClusters>0 ) delete[] fClusters;
    fClusters = 0;
    fNClusters = 0;
  }
  return true;
}
//------------------------------------------------------------------------
void AliPHOSDATreeEvent::Print(Option_t *opt) const
{
  // Print Out
  
  char* when = ctime(&fTime);
  std::cout<<" AliPHOSDATreeEvent:: fNDigits = "<<fNDigits<<" fNClusters="<<fNClusters<<" : "<<when;
  int ndigits = fNDigits;
  while( ndigits-- ){
    std::cout<<" =>["<<ndigits<<"] : ";
    fDigits[ndigits].Print(opt);
  }
  int nclusters = fNClusters;
  while( nclusters-- ){
    std::cout<<" =>["<<nclusters<<"] : ";
    fClusters[nclusters].Print(opt);
  }
}
//------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& out, const AliPHOSDATreeEvent& event){
  // Print Out

  char* when = ctime(&(event.fTime));
  out<<" AliPHOSDATreeEvent:: fNClusters="<<event.fNClusters<<" : "<<when;
  int ndigits = event.fNDigits;
  while( ndigits-- ){
    out<<" =>["<<ndigits<<"] : "<<event.fDigits[ndigits];
    if( ndigits!=0 ) out<<std::endl;
  }
  int nclusters = event.fNClusters;
  while( nclusters-- ){
    out<<" =>["<<nclusters<<"] : "<<event.fClusters[nclusters];
    if( nclusters!=0 ) out<<std::endl;
  }
  return out;
}
//------------------------------------------------------------------------
