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

#include <math.h>
#include <Rtypes.h>
#include <iostream>
#include "AliPHOSDATreeDigit.h"
#include "AliPHOSDATreeCluster.h"
ClassImp(AliPHOSDATreeCluster)
//------------------------------------------------------------------------
AliPHOSDATreeCluster::AliPHOSDATreeCluster(const AliPHOSDATreeCluster& cluster):
//fEnergy(cluster.fEnergy),fX(cluster.fX),fY(cluster.fY),fZ(cluster.fZ),fNDigits(cluster.fNDigits){
fEnergy(cluster.fEnergy),fRow(cluster.fRow),fCol(cluster.fCol),fNDigits(cluster.fNDigits),fDigits(0){
  // Copy Constructor

  if( fNDigits > 0 ){
    fDigits = new AliPHOSDATreeDigit[fNDigits];
    int ndigits = fNDigits;
    while( ndigits-- ){
      fDigits[ndigits] = cluster.fDigits[ndigits];
    }
  } else {
    fDigits = 0;
  }
}
//------------------------------------------------------------------------
AliPHOSDATreeCluster& AliPHOSDATreeCluster::operator=(const AliPHOSDATreeCluster& cluster){
  // Copy Operator

  if (this != &cluster) {
    if( fNDigits> 0 ) delete[] fDigits;
    fEnergy = cluster.fEnergy;
    fNDigits = cluster.fNDigits;
    fRow = cluster.fRow;
    fCol = cluster.fCol;
    //fX = cluster.fX;
    //fY = cluster.fY;
    //fZ = cluster.fZ;
    if( fNDigits > 0 ){
      fDigits = new AliPHOSDATreeDigit[fNDigits];
      int ndigits = fNDigits;
      while( ndigits-- ){
	fDigits[ndigits] = cluster.fDigits[ndigits];
      }
    } else {
      fDigits = 0;
    }
  }
  return *this;
}
//------------------------------------------------------------------------
void AliPHOSDATreeCluster::Print(Option_t *opt) const
{
  // Print out
  std::cout<<" AliPHOSDATreeCluster:: Energy="<<fEnergy<<" NDigits="<<fNDigits
	   <<" (row,col)=("<<fRow<<","<<fCol<<")"<<std::endl;
    //<<" (x,y,z)=("<<fX<<","<<fY<<","<<fZ<<")"<<std::endl;
  int ndigits = fNDigits;
  while( ndigits-- ){
    std::cout<<"   -->["<<ndigits<<"] : ";
    fDigits[ndigits].Print(opt);
  }
}
//------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& out, const AliPHOSDATreeCluster& cluster){
  // Print out
  out<<" AliPHOSDATreeCluster:: Energy="<<cluster.fEnergy<<" NDigits="<<cluster.fNDigits
     <<" (row,col)=("<<cluster.fRow<<","<<cluster.fCol<<")"<<std::endl;
  int ndigits = cluster.fNDigits;
  while( ndigits-- ){
    out<<"   -->["<<ndigits<<"] : "<<cluster.fDigits[ndigits];
    if( ndigits!=0 ) out<<std::endl;
  }
  return out;
}
//------------------------------------------------------------------------
void AliPHOSDATreeCluster::Reset(){
  // Reset information

  fEnergy = 0;
  fRow = -100;
  fCol = -100;
  //fX = 0;
  //fY = 0;
  //fZ = 0;
  if( fNDigits> 0 ) delete[] fDigits;
  fNDigits = 0;
}
//------------------------------------------------------------------------
bool AliPHOSDATreeCluster::Append(AliPHOSDATreeDigit& digit){
  // Add digit information and sum all energy
  //
  if(! digit.IsValid() ){
    std::cout<<" AliPHOSDATreeCluster::Append():: Error!! Digit is not valid.."<<std::endl;
    return false;
  }
  AliPHOSDATreeDigit* newfDigits = new AliPHOSDATreeDigit[fNDigits+1];
  bool bsearching = true;
  int ndigit = fNDigits;
  while( ndigit-- ){
    if( fDigits[ndigit].GetAbsId() == digit.GetAbsId() ){
      std::cout<<" AliPHOSDATreeCluster::Append():: Error!! The channel already exist."<<std::endl;
      std::cout<<" Add "<<digit<<std::endl;
      std::cout<<" into *this"<<*this<<std::endl;
      delete[] newfDigits;
      return false;
    }
    if( fDigits[ndigit].GetEnergy() < digit.GetEnergy() ){
      newfDigits[ndigit+1] = fDigits[ndigit];
    } else {
      if( bsearching ) {
	bsearching = false;
	newfDigits[ndigit+1] = digit;
      }
      newfDigits[ndigit] = fDigits[ndigit];
    }
  }
  if( bsearching ) newfDigits[0] = digit;
  if( fNDigits>0 ) delete[] fDigits;
  fNDigits++;
  fDigits = newfDigits;
  fEnergy += digit.GetEnergy();
  return true;
}
//------------------------------------------------------------------------
bool AliPHOSDATreeCluster::Append(AliPHOSDATreeCluster& cluster){
  // Add another cluster information and sum all energy
  //
  AliPHOSDATreeDigit* newfDigits = new AliPHOSDATreeDigit[fNDigits+cluster.fNDigits];
  int ndigits1 = fNDigits;
  int ndigits2 = cluster.fNDigits;
  int ndigitsall = ndigits1 + ndigits2;
  while( ndigitsall-- ){
    //std::cout<<" ------ ndigits1:"<<ndigits1<<" ndigits2:"<<ndigits2<<std::endl;
    if( ndigits1 && ndigits2 ){
      if( fDigits[ndigits1-1].GetEnergy() < cluster.fDigits[ndigits2-1].GetEnergy() ){
	newfDigits[ndigitsall] = fDigits[--ndigits1];
      } else {
	newfDigits[ndigitsall]= cluster.fDigits[--ndigits2];
      }
    } else if ( ndigits1 && ndigits2==0 ){
      newfDigits[ndigitsall] = fDigits[--ndigits1];
    } else if ( ndigits2 && ndigits1==0 ){
      newfDigits[ndigitsall]= cluster.fDigits[--ndigits2];
    } else {
      std::cout<<" AliPHOSDATreeCluster::Append() Something wrong.. "<<std::endl;
      delete newfDigits;
      return false;
    }
  }
  if(fNDigits>0) delete[] fDigits;
  fDigits = newfDigits;
  fNDigits += cluster.fNDigits;
  fEnergy += cluster.GetEnergy();
  return true;
    
}
//------------------------------------------------------------------------
bool AliPHOSDATreeCluster::IsNeighbor(const AliPHOSDATreeDigit& digit) const{
  // Check wether the given digit is neighboring to this cluster.
  // Return true if yes.

  bool status = false;
  int ndigits = fNDigits;
  while( ndigits-- && !status ){
    status = digit.IsNeighbor(fDigits[ndigits]);
  }
  return status;
}
//------------------------------------------------------------------------
bool AliPHOSDATreeCluster::IsNeighbor(const AliPHOSDATreeCluster& cluster) const{
  // Check wether the given cluster is neighboring to this cluster.
  // Return true if yes.

  bool status = false;
  int ndigits = fNDigits;
  while( ndigits-- && !status ){
    status = cluster.IsNeighbor(fDigits[ndigits]);
  }
  return status;
}
//------------------------------------------------------------------------
bool AliPHOSDATreeCluster::CalculateProperty(){
  // Calculate the hit position
  // (calculation of dispersion is not valid)
  
  fCol = 0;
  fRow = 0;
  float totweight = 0;
  float weight;
  int ndigits = fNDigits;
  while( ndigits-- ){
    weight = log(fDigits[ndigits].GetEnergy()/fEnergy) + 4.5; //4.5 is for PHOS
    //std::cout<<" AliPHOSDATreeCluster::CalculateProperty() DEBUG: ndigits="<<ndigits<<" weight="<<weight<<std::endl;
    if( weight > 0 ){
      totweight += weight;
      fRow += fDigits[ndigits].GetRow() * weight;
      fCol += fDigits[ndigits].GetCol() * weight;
    }
  }
  //std::cout<<" AliPHOSDATreeCluster::CalculateProperty() DEBUG: totweight="<<totweight<<std::endl;
  if( totweight > 0 ){
    fRow /= totweight;
    fCol /= totweight;
  } else {
    fRow = 0;
    fCol = 0;
  }
  /*
  float disp = 0;
  if( totweight > ){
    ndigits = fNDigits;
    while( ndigits-- ){
      weight = log(fDigits[ndigits].GetEnergy()/fEnergy) + 4.5; //4.5 is for PHOS
      disp += weight * ( (fDigits[ndigits].GetRow()-fRow)*(fDigits[ndigits].GetRow()-fRow) +
			 (fDigits[ndigits].GetCol()-fCol)*(fDigits[ndigits].GetCol()-fCol) );
    }
    disp /= totweight;
  }
  */
  return true;
}
//------------------------------------------------------------------------
