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

#include "AliRICHRawCluster.h"
#include <iostream.h>


ClassImp(AliRICHRawCluster)
//__________________________________________________________________________________________________

AliRICHRawCluster :: AliRICHRawCluster() 
{//default ctor
  fTracks[0]=fTracks[1]=fTracks[2]=-1; 
  fX=fY=0;
  fQ=fMultiplicity=0;
  for (int k=0;k<50;k++) {
    fIndexMap[k]=-1;
    fContMap[k]=0;
    fPhysicsMap[k]=-1;
    fCtype=-1;
  }
  fNcluster[0]=fNcluster[1]=-1;
}//ctor
//__________________________________________________________________________________________________
Int_t AliRICHRawCluster::Compare(const TObject *obj) const
{//Compare two clusters
  AliRICHRawCluster *raw=(AliRICHRawCluster *)obj;
  Float_t y=fY;
  Float_t yo=raw->fY;
  if(y>yo)      return  1;
  else if(y<yo) return -1;
  else          return  0;
}//Compare()
//__________________________________________________________________________________________________
Int_t AliRICHRawCluster::PhysicsContribution()
{//Type of physics processes
  Int_t iPhys=0;
  Int_t iBg=0;
  Int_t iMixed=0;
  for (Int_t i=0; i<fMultiplicity; i++){
    if(fPhysicsMap[i]==2) iPhys++;
    if(fPhysicsMap[i]==1) iMixed++;
    if(fPhysicsMap[i]==0) iBg++;
  }
  if(iMixed==0 && iBg==0)                         return 2;
  else if((iPhys != 0 && iBg !=0) || iMixed != 0) return 1;
  else	                                          return 0;
}//PhysicsContribution
//__________________________________________________________________________________________________
void AliRICHRawCluster::Print(Option_t*)const
{
  Info("","X=%7.2f, Y=%7.2f, Qdc=%4i, Peak=%4i, Multip=%2i,      T0=%5i T1=%5i T2=%5i",
           fX,      fY,      fQ,      fPeakSignal,fMultiplicity, fTracks[0], fTracks[1], fTracks[2]);
  
  for(int i=0;i<fMultiplicity;i++)
    cout<<"D"<<i<<"="<<fIndexMap[i]<<" C="<<fContMap[i]<<endl;
}//void AliRICHRawCluster::Print(Option_t *option)const
