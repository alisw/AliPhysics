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
Revision 1.2  2000/06/28 12:19:18  morsch
More consequent seperation of global input data services (AliMUONClusterInput singleton) and the
cluster and hit reconstruction algorithms in AliMUONClusterFinderVS.
AliMUONClusterFinderVS becomes the base class for clustering and hit reconstruction.
It requires two cathode planes. Small modifications in the code will make it usable for
one cathode plane and, hence, more general (for test beam data).
AliMUONClusterFinder is now obsolete.

Revision 1.1  2000/06/28 08:06:10  morsch
Avoid global variables in AliMUONClusterFinderVS by seperating the input data for the fit from the
algorithmic part of the class. Input data resides inside the AliMUONClusterInput singleton.
It also naturally takes care of the TMinuit instance.

*/
#include "AliRun.h"
#include "AliMUON.h"
#include "AliMUONChamber.h"
#include "AliMUONClusterInput.h"
#include "AliMUONSegmentation.h"
#include "AliMUONResponse.h"
#include "AliMUONRawCluster.h"
#include "AliMUONDigit.h"

#include <TClonesArray.h>
#include <TMinuit.h>

ClassImp(AliMUONClusterInput)

AliMUONClusterInput* AliMUONClusterInput::fgClusterInput = 0; 
TMinuit* AliMUONClusterInput::fgMinuit = 0; 

AliMUONClusterInput* AliMUONClusterInput::Instance()
{
// return pointer to the singleton instance
    if (fgClusterInput == 0) {
	fgClusterInput = new AliMUONClusterInput();
	fgMinuit = new TMinuit(5);
    }
    
    return fgClusterInput;
}

void AliMUONClusterInput::SetDigits(Int_t chamber, TClonesArray* dig1, TClonesArray* dig2)
{
// Set pointer to digits with corresponding segmentations and responses (two cathode planes)
    fChamber=chamber;
    fDigits[0]=dig1;
    fDigits[1]=dig2; 
    fNDigits[0]=dig1->GetEntriesFast();
    fNDigits[1]=dig2->GetEntriesFast();
    
    AliMUON *pMUON;
    AliMUONChamber* iChamber;

    pMUON = (AliMUON*) gAlice->GetModule("MUON");
    iChamber =  &(pMUON->Chamber(chamber));

    fSegmentation[0]=iChamber->SegmentationModel(1);
    fSegmentation[1]=iChamber->SegmentationModel(2);
    fResponse=iChamber->ResponseModel();
    fNseg = 2;
}

void AliMUONClusterInput::SetDigits(Int_t chamber, TClonesArray* dig)
{
// Set pointer to digits with corresponding segmentations and responses (one cathode plane)
    fDigits[0]=dig;
    AliMUON *pMUON;
    AliMUONChamber* iChamber;

    pMUON = (AliMUON*) gAlice->GetModule("MUON");
    iChamber =  &(pMUON->Chamber(chamber));

    fSegmentation[0]=iChamber->SegmentationModel(1);
    fResponse=iChamber->ResponseModel();
    fNseg=1;
}

void  AliMUONClusterInput::SetCluster(AliMUONRawCluster* cluster)
{
// Set the current cluster
    printf("\n %p \n", cluster);
    fCluster=cluster;
    Float_t qtot;
    Int_t   i, cath, ix, iy;
    AliMUONDigit* digit;
    fNmul[0]=cluster->fMultiplicity[0];
    fNmul[1]=cluster->fMultiplicity[1];
    printf("\n %p %p ", fDigits[0], fDigits[1]);
    
    for (cath=0; cath<2; cath++) {
	qtot=0;
	for (i=0; i<fNmul[cath]; i++) {
	    // pointer to digit
	    digit =(AliMUONDigit*)
		(fDigits[cath]->UncheckedAt(cluster->fIndexMap[i][cath]));
	    // pad coordinates
	    ix = digit->fPadX;
	    iy = digit->fPadY;
	    // pad charge
	    fCharge[i][cath] = digit->fSignal;
	    // pad centre coordinates
//	    fSegmentation[cath]->GetPadCxy(ix, iy, x, y);
            // globals kUsed in fitting functions
	    fix[i][cath]=ix;
	    fiy[i][cath]=iy;
	    // total charge per cluster
	    qtot+=fCharge[i][cath];
	} // loop over cluster digits
	fQtot[cath]=qtot;
	fChargeTot[cath]=Int_t(qtot);  
    }  // loop over cathodes
}



Float_t AliMUONClusterInput::DiscrChargeS1(Int_t i,Double_t *par) 
{
// par[0]    x-position of cluster
// par[1]    y-position of cluster

   fSegmentation[0]->SetPad(fix[i][0], fiy[i][0]);
//  First Cluster
   fSegmentation[0]->SetHit(par[0],par[1],0);
   Float_t q1=fResponse->IntXY(fSegmentation[0]);
    
   Float_t value = fQtot[0]*q1;
   return value;
}

Float_t AliMUONClusterInput::DiscrChargeCombiS1(Int_t i,Double_t *par, Int_t cath) 
{
// par[0]    x-position of cluster
// par[1]    y-position of cluster

   fSegmentation[cath]->SetPad(fix[i][cath], fiy[i][cath]);
//  First Cluster
   fSegmentation[cath]->SetHit(par[0],par[1],0);
   Float_t q1=fResponse->IntXY(fSegmentation[cath]);
    
   Float_t value = fQtot[cath]*q1;
   return value;
}


Float_t AliMUONClusterInput::DiscrChargeS2(Int_t i,Double_t *par) 
{
// par[0]    x-position of first  cluster
// par[1]    y-position of first  cluster
// par[2]    x-position of second cluster
// par[3]    y-position of second cluster
// par[4]    charge fraction of first  cluster
// 1-par[4]  charge fraction of second cluster

   fSegmentation[0]->SetPad(fix[i][0], fiy[i][0]);
//  First Cluster
   fSegmentation[0]->SetHit(par[0],par[1],0);
   Float_t q1=fResponse->IntXY(fSegmentation[0]);
    
//  Second Cluster
   fSegmentation[0]->SetHit(par[2],par[3],0);
   Float_t q2=fResponse->IntXY(fSegmentation[0]);
    
   Float_t value = fQtot[0]*(par[4]*q1+(1.-par[4])*q2);
   return value;
}

Float_t AliMUONClusterInput::DiscrChargeCombiS2(Int_t i,Double_t *par, Int_t cath) 
{
// par[0]    x-position of first  cluster
// par[1]    y-position of first  cluster
// par[2]    x-position of second cluster
// par[3]    y-position of second cluster
// par[4]    charge fraction of first  cluster
// 1-par[4]  charge fraction of second cluster

   fSegmentation[cath]->SetPad(fix[i][cath], fiy[i][cath]);
//  First Cluster
   fSegmentation[cath]->SetHit(par[0],par[1],0);
   Float_t q1=fResponse->IntXY(fSegmentation[cath]);
    
//  Second Cluster
   fSegmentation[cath]->SetHit(par[2],par[3],0);
   Float_t q2=fResponse->IntXY(fSegmentation[cath]);
   Float_t value;
   if (cath==0) {
       value = fQtot[0]*(par[4]*q1+(1.-par[4])*q2);
   } else {
       value = fQtot[1]*(par[5]*q1+(1.-par[5])*q2);
   }
   return value;
}

void AliMUONClusterInput::Streamer(TBuffer &R__b) {} 


