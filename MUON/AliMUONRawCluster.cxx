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

/*
$Log$
Revision 1.1.2.1  2000/06/09 22:04:50  morsch
Was before in DataStructures.cxx

*/

#include "AliMUONRawCluster.h"
#include <TObjArray.h>
#include <TArrayF.h>

ClassImp(AliMUONRawCluster);


AliMUONRawCluster::AliMUONRawCluster() {
// Constructor
    fTracks[0]=fTracks[1]=fTracks[2]=-1; 
    for (int j=0;j<2;j++) {
	fQ[j]=0;
	fX[j]=0;
	fY[j]=0;
	fMultiplicity[j]=0;
	fPeakSignal[j]=-1;
	fChi2[j]=-1;
	
	for (int k=0;k<50;k++) {
	    fIndexMap[k][j]=-1;
	    fOffsetMap[k][j]=0;
	    fContMap[k][j]=0;
	    fPhysicsMap[k]=-1;
	}
    }
    fNcluster[0]=fNcluster[1]=-1;
}

Int_t AliMUONRawCluster::Compare(TObject *obj)
{
  /*
         AliMUONRawCluster *raw=(AliMUONRawCluster *)obj;
	 Float_t r=GetRadius();
         Float_t ro=raw->GetRadius();
         if (r>ro) return 1;
         else if (r<ro) return -1;
         else return 0;
  */
         AliMUONRawCluster *raw=(AliMUONRawCluster *)obj;
	 Float_t y=fY[0];
         Float_t yo=raw->fY[0];
         if (y>yo) return 1;
         else if (y<yo) return -1;
         else return 0;

}

Int_t AliMUONRawCluster::
BinarySearch(Float_t y, TArrayF coord, Int_t from, Int_t upto)
{
   // Find object using a binary search. Array must first have been sorted.
   // Search can be limited by setting upto to desired index.

   Int_t low=from, high=upto-1, half;
   while(high-low>1) {
        half=(high+low)/2;
        if(y>coord[half]) low=half;
        else high=half;
   }
   return low;
}

void AliMUONRawCluster::SortMin(Int_t *idx,Float_t *xdarray,Float_t *xarray,Float_t *yarray,Float_t *qarray, Int_t ntr)
{
  //
  // Get the 3 closest points(cog) one can find on the second cathode 
  // starting from a given cog on first cathode
  //
  
  //
  //  Loop over deltax, only 3 times
  //
  
    Float_t xmin;
    Int_t jmin;
    Int_t id[3] = {-2,-2,-2};
    Float_t jx[3] = {0.,0.,0.};
    Float_t jy[3] = {0.,0.,0.};
    Float_t jq[3] = {0.,0.,0.};
    Int_t jid[3] = {-2,-2,-2};
    Int_t i,j,imax;
  
    if (ntr<3) imax=ntr;
    else imax=3;
    for(i=0;i<imax;i++){
        xmin=1001.;
        jmin=0;
    
        for(j=0;j<ntr;j++){
            if ((i == 1 && j == id[i-1]) 
	          ||(i == 2 && (j == id[i-1] || j == id[i-2]))) continue;
           if (TMath::Abs(xdarray[j]) < xmin) {
	      xmin = TMath::Abs(xdarray[j]);
	      jmin=j;
           }       
        } // j
        if (xmin != 1001.) {    
           id[i]=jmin;
           jx[i]=xarray[jmin]; 
           jy[i]=yarray[jmin]; 
           jq[i]=qarray[jmin]; 
           jid[i]=idx[jmin];
        } 
    
    }  // i
  
    for (i=0;i<3;i++){
        if (jid[i] == -2) {
            xarray[i]=1001.;
            yarray[i]=1001.;
            qarray[i]=1001.;
            idx[i]=-1;
        } else {
            xarray[i]=jx[i];
            yarray[i]=jy[i];
            qarray[i]=jq[i];
            idx[i]=jid[i];
        }
    }

}


Int_t AliMUONRawCluster::PhysicsContribution()
{
// Evaluate physics contribution to cluster
  Int_t iPhys=0;
  Int_t iBg=0;
  Int_t iMixed=0;
  for (Int_t i=0; i<fMultiplicity[0]; i++) {
    if (fPhysicsMap[i]==2) iPhys++;
    if (fPhysicsMap[i]==1) iMixed++;
    if (fPhysicsMap[i]==0) iBg++;
  }
  if (iMixed==0 && iBg==0) {
    return 2;
  } else if ((iPhys != 0 && iBg !=0) || iMixed != 0) {
    return 1;
  } else {
    return 0;
  }
}
