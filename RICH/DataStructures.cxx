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
*/


#include "DataStructures.h"
#include "AliRun.h"
#include <TObjArray.h>
#include <TParticle.h>

ClassImp(AliRICHPadHit)

ClassImp(AliRICHHit)
    
//___________________________________________
AliRICHHit::AliRICHHit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
    AliHit(shunt, track)
{
    fChamber=vol[0];
    fParticle=hits[0];
    fX=hits[1];
    fY=hits[2];
    fZ=hits[3];
    fTheta=hits[4];
    fPhi=hits[5];
    fTlength=hits[6];
    fEloss=hits[7];
    fPHfirst=(Int_t) hits[8];
    fPHlast=(Int_t) hits[9];
    fLoss=hits[13];
    fMomX=hits[14];
    fMomY=hits[15];
    fMomZ=hits[16];
}

ClassImp(AliRICHCerenkov)
//___________________________________________
AliRICHCerenkov::AliRICHCerenkov(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
    AliHit(shunt, track)
{
    fChamber=vol[0];
    fX=hits[1];
    fY=hits[2];
    fZ=hits[3];
    fTheta=hits[4];
    fPhi=hits[5];
    fEloss=hits[7];
    fPHfirst=(Int_t) hits[8];
    fPHlast=(Int_t) hits[9];
    fCMother=Int_t(hits[10]);
    fIndex = hits[11];
    fProduction = hits[12];  
    fLoss=hits[13];
    fMomX=hits[14];
    fMomY=hits[15];
    fMomZ=hits[16];
}

//___________________________________________
AliRICHPadHit::AliRICHPadHit(Int_t *clhits)
{
    fHitNumber=clhits[0];
    fCathode=clhits[1];
    fQ=clhits[2];
    fPadX=clhits[3];
    fPadY=clhits[4];
    fQpad=clhits[5];
    fRSec=clhits[6];
}

ClassImp(AliRICHDigit)
//_____________________________________________________________________________
AliRICHDigit::AliRICHDigit(Int_t *digits)
{
    //
    // Creates a RICH digit object to be updated
    //
    fPadX        = digits[0];
    fPadY        = digits[1];
    fSignal      = digits[2];
    
}
//_____________________________________________________________________________
AliRICHDigit::AliRICHDigit(Int_t *tracks, Int_t *charges, Int_t *digits)
{
    //
    // Creates a RICH digit object
    //
    fPadX        = digits[0];
    fPadY        = digits[1];
    fSignal      = digits[2];
    for(Int_t i=0; i<100; i++) {
	fTcharges[i]  = charges[i];
	fTracks[i]    = tracks[i];
    }
}

ClassImp(AliRICHTransientDigit)
    
//____________________________________________________________________________
AliRICHTransientDigit::AliRICHTransientDigit(Int_t ich, Int_t *digits): 
    AliRICHDigit(digits)
{
    //
    // Creates a RICH digit list object
    //
    
    fChamber = ich;
    fTrackList   = new TObjArray;
    
}
//_____________________________________________________________________________


ClassImp(AliRICHRawCluster)
Int_t AliRICHRawCluster::Compare(TObject *obj)
{
         AliRICHRawCluster *raw=(AliRICHRawCluster *)obj;
	 Float_t y=fY;
         Float_t yo=raw->fY;
         if (y>yo) return 1;
         else if (y<yo) return -1;
         else return 0;

}

Int_t AliRICHRawCluster::
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

void AliRICHRawCluster::SortMin(Int_t *idx,Float_t *xdarray,Float_t *xarray,Float_t *yarray,Float_t *qarray, Int_t ntr)
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


Int_t AliRICHRawCluster::PhysicsContribution()
{
    Int_t iPhys=0;
    Int_t iBg=0;
    Int_t iMixed=0;
    for (Int_t i=0; i<fMultiplicity; i++) {
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
//_____________________________________________________________________________


ClassImp(AliRICHRecHit)


//_____________________________________________________________________________
AliRICHRecHit::AliRICHRecHit(Int_t id, Float_t *rechit)
{
    //
    // Creates a RICH rec. hit object
    //
    Theta        = rechit[0];
    Phi          = rechit[1];
    Omega        = rechit[2];
    fX           = rechit[3];
    fY           = rechit[4];
}
