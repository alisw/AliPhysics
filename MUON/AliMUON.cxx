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


///////////////////////////////////////////////
//  Manager and hits classes for set:MUON     //
////////////////////////////////////////////////

#include "Riostream.h"

#include <AliPDG.h>
#include <TBRIK.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGeometry.h>
#include <TMinuit.h>
#include <TNode.h> 
#include <TNtuple.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TObjectTable.h>
#include <TPad.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TRandom.h> 
#include <TRotMatrix.h>
#include <TTUBE.h>
#include <TTUBE.h>
#include <TTree.h> 
#include <TVector.h>
#include <TVirtualMC.h>

#include "AliConst.h" 
#include "AliHeader.h"
#include "AliHitMap.h"
#include "AliLoader.h"
#include "AliMUONLoader.h"
#include "AliMUON.h"
#include "AliMUONChamberTrigger.h"
#include "AliMUONClusterFinderVS.h"
#include "AliMUONClusterInput.h"
#include "AliMUONConstants.h"
#include "AliMUONDigit.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONHit.h"
#include "AliMUONHitMapA1.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONMerger.h"	
#include "AliMUONPadHit.h"
#include "AliMUONRawCluster.h"
#include "AliMUONTransientDigit.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONTriggerDecision.h"
#include "AliRun.h"	


// Defaults parameters for Z positions of chambers
// taken from values for "stations" in AliMUON::AliMUON
//     const Float_t zch[7]={528, 690., 975., 1249., 1449., 1610, 1710.};
// and from array "dstation" in AliMUONv1::CreateGeometry
//          Float_t dstation[5]={20., 20., 20, 20., 20.};
//     for tracking chambers,
//          according to (Z1 = zch - dstation) and  (Z2 = zch + dstation)
//          for the first and second chambers in the station, respectively,
// and from "DTPLANES" in AliMUONv1::CreateGeometry
//           const Float_t DTPLANES = 15.;
//     for trigger chambers,
//          according to (Z1 = zch) and  (Z2 = zch + DTPLANES)
//          for the first and second chambers in the station, respectively

ClassImp(AliMUON)
//___________________________________________
AliMUON::AliMUON()
{
// Default Constructor
//
    fNCh             = 0;
    fNTrackingCh     = 0;
    fIshunt          = 0;
    fPadHits         = 0;
    fNPadHits        = 0;
    fChambers        = 0;
    fDchambers       = 0;
    fTriggerCircuits = 0;  
    fNdch            = 0;
    fRawClusters     = 0;
    fNrawch          = 0;
    fGlobalTrigger   = 0;
    fNLocalTrigger   = 0;
    fLocalTrigger    = 0;
    fNLocalTrigger   = 0;
    fAccMin          = 0.;
    fAccMax          = 0.;   
    fAccCut          = kFALSE;
    fMerger          = 0;
    fFileName        = 0;
}
 
//___________________________________________
AliMUON::AliMUON(const char *name, const char *title)
       : AliDetector(name,title)
{
//Begin_Html
/*
<img src="gif/alimuon.gif">
*/
//End_Html

   fHits     = new TClonesArray("AliMUONHit",1000);
   gAlice->AddHitList(fHits);
   fPadHits = new TClonesArray("AliMUONPadHit",10000);
   fNPadHits  =  0;
   fIshunt     =  0;

   fNCh             = AliMUONConstants::NCh(); 
   fNTrackingCh     = AliMUONConstants::NTrackingCh();
   fNdch            = new Int_t[fNCh];

   fDchambers = new TObjArray(AliMUONConstants::NCh());

   Int_t i;
   
   for (i=0; i<AliMUONConstants::NCh() ;i++) {
       fDchambers->AddAt(new TClonesArray("AliMUONDigit",10000),i); 
       fNdch[i]=0;
   }

   fNrawch      = new Int_t[fNTrackingCh];

   fRawClusters = new TObjArray(AliMUONConstants::NTrackingCh());

   for (i=0; i<AliMUONConstants::NTrackingCh();i++) {
       fRawClusters->AddAt(new TClonesArray("AliMUONRawCluster",10000),i); 
       fNrawch[i]=0;
   }

   fGlobalTrigger = new TClonesArray("AliMUONGlobalTrigger",1);    
   fNGlobalTrigger = 0;
   fLocalTrigger  = new TClonesArray("AliMUONLocalTrigger",234);    
   fNLocalTrigger = 0;

   SetMarkerColor(kRed);
//
//
//
//

    Int_t ch;

    fChambers = new TObjArray(AliMUONConstants::NCh());

    // Loop over stations
    for (Int_t st = 0; st < AliMUONConstants::NCh() / 2; st++) {
      // Loop over 2 chambers in the station
	for (Int_t stCH = 0; stCH < 2; stCH++) {
//
//    
//    Default Parameters for Muon Tracking Stations


	    ch = 2 * st + stCH;
//
	    if (ch < AliMUONConstants::NTrackingCh()) {
	      fChambers->AddAt(new AliMUONChamber(ch),ch);
	    } else {
	      fChambers->AddAt(new AliMUONChamberTrigger(ch),ch);
	    }
	    
        //PH	    AliMUONChamber* chamber = (AliMUONChamber*) (*fChambers)[ch];
	    AliMUONChamber* chamber = (AliMUONChamber*) fChambers->At(ch);
	    
	    chamber->SetGid(0);
	    // Default values for Z of chambers
	    chamber->SetZ(AliMUONConstants::DefaultChamberZ(ch));
//
	    chamber->InitGeo(AliMUONConstants::DefaultChamberZ(ch));
//          Set chamber inner and outer radius to default
	    chamber->SetRInner(AliMUONConstants::Dmin(st)/2);
	    chamber->SetROuter(AliMUONConstants::Dmax(st)/2);
//
	} // Chamber stCH (0, 1) in 
    }     // Station st (0...)
//    fChambers->SetLast(AliMUONConstants::NCh());
    fMaxStepGas=0.01; 
    fMaxStepAlu=0.1; 
    fMaxDestepGas=-1;
    fMaxDestepAlu=-1;
//
   fMaxIterPad   = 0;
   fCurIterPad   = 0;
//
   fAccMin          = 0.;
   fAccMax          = 0.;   
   fAccCut          = kFALSE;

   // cp new design of AliMUONTriggerDecision
   fTriggerCircuits = new TObjArray(AliMUONConstants::NTriggerCircuit());
   for (Int_t circ=0; circ<AliMUONConstants::NTriggerCircuit(); circ++) {
     fTriggerCircuits->AddAt(new AliMUONTriggerCircuit(),circ);     

   }
     fMerger = 0;
}
 
//___________________________________________
AliMUON::AliMUON(const AliMUON& rMUON):AliDetector(rMUON)
{
// Dummy copy constructor
    ;
    
}

AliMUON::~AliMUON()
{
// Destructor
    if(fDebug) printf("%s: Calling AliMUON destructor !!!\n",ClassName());
    
    fIshunt  = 0;
 
    // Delete TObjArrays
 
    if (fChambers){
      fChambers->Delete();
      delete fChambers;
    }
 
    if (fTriggerCircuits){
      fTriggerCircuits->Delete();
      delete fTriggerCircuits;
    }
 
    if (fDchambers){
      fDchambers->Delete();
      delete fDchambers;
    }
 
    if (fRawClusters){
      fRawClusters->Delete();
      delete fRawClusters;
    }

    if (fNrawch) delete [] fNrawch;
 
    // Delete TClonesArrays
 
    if (fPadHits){
      fPadHits->Delete();
      delete fPadHits;
    }
 
    if (fGlobalTrigger){
      fGlobalTrigger->Delete();
      delete fGlobalTrigger;
    }
    fNGlobalTrigger = 0;
    
    if (fLocalTrigger){
      fLocalTrigger->Delete();
      delete fLocalTrigger;
    }
    fNLocalTrigger = 0;

    if (fHits) {
      fHits->Delete();
      delete fHits;
    }

    if (fMerger) delete fMerger;
    if (fNdch) delete [] fNdch;

}
 
//___________________________________________
void AliMUON::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliMUONHit(fIshunt,track,vol,hits);
}
//___________________________________________
void AliMUON::AddHit(Int_t fIshunt, Int_t track, Int_t iChamber, 
	      Int_t idpart, Float_t X, Float_t Y, Float_t Z, 
	      Float_t tof, Float_t momentum, Float_t theta, 
	      Float_t phi, Float_t length, Float_t destep)
{
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliMUONHit(fIshunt, track, iChamber, 
	       idpart, X, Y, Z, 
	       tof, momentum, theta, 
	       phi, length, destep);
}
//___________________________________________
void AliMUON::AddPadHit(Int_t *clhits)  // To be removed
{
   TClonesArray &lclusters = *fPadHits;
   new(lclusters[fNPadHits++]) AliMUONPadHit(clhits);
}
//_____________________________________________________________________________
void AliMUON::AddDigits(Int_t id, Int_t *tracks, Int_t *charges, Int_t *digits)
{
    //
    // Add a MUON digit to the list
    //

  //PH    TClonesArray &ldigits = * ((TClonesArray*)(*fDchambers)[id]);
    TClonesArray &ldigits = * ( (TClonesArray*) fDchambers->At(id) );
    new(ldigits[fNdch[id]++]) AliMUONDigit(tracks,charges,digits);
}

//_____________________________________________________________________________
void AliMUON::AddRawCluster(Int_t id, const AliMUONRawCluster& c)
{
    //
    // Add a MUON digit to the list
    //

  //PH    TClonesArray &lrawcl = *((TClonesArray*)(*fRawClusters)[id]);
    TClonesArray &lrawcl = *((TClonesArray*)fRawClusters->At(id));
    new(lrawcl[fNrawch[id]++]) AliMUONRawCluster(c);
}

//___________________________________________
void AliMUON::AddGlobalTrigger(Int_t *singlePlus, Int_t *singleMinus,
			       Int_t *singleUndef,
			       Int_t *pairUnlike, Int_t *pairLike)
{
// add a MUON Global Trigger to the list (only one GlobalTrigger per event !)
  TClonesArray &globalTrigger = *fGlobalTrigger;
  new(globalTrigger[fNGlobalTrigger++]) 
    AliMUONGlobalTrigger(singlePlus, singleMinus,  singleUndef, pairUnlike, 
			 pairLike);
}
//___________________________________________
void AliMUON::AddLocalTrigger(Int_t *localtr)
{
// add a MUON Local Trigger to the list
  TClonesArray &localTrigger = *fLocalTrigger;
  new(localTrigger[fNLocalTrigger++]) AliMUONLocalTrigger(localtr);
}

//___________________________________________
void AliMUON::BuildGeometry()
{
// Geometry for event display
  for (Int_t i=0; i<7; i++) {
    for (Int_t j=0; j<2; j++) {
      Int_t id=2*i+j+1;
      this->Chamber(id-1).SegmentationModel(1)->Draw("eventdisplay");
    }
  }
}

//___________________________________________
Int_t AliMUON::DistancetoPrimitive(Int_t , Int_t )
{
   return 9999;
}

//___________________________________________
void AliMUON::MakeBranch(Option_t* option)
{
    //
    // Create Tree branches for the MUON.
    //
    const Int_t kBufferSize = 4000;
    char branchname[30];
    sprintf(branchname,"%sCluster",GetName());
    
    
    const char *cD = strstr(option,"D");
    const char *cR = strstr(option,"R");
    const char *cH = strstr(option,"H");

    if (TreeH() && cH) 
     {
      if (fPadHits == 0x0) fPadHits = new TClonesArray("AliMUONPadHit",10000);
      MakeBranchInTree(TreeH(), branchname, &fPadHits, kBufferSize, 0);
      if (fHits == 0x0) fHits     = new TClonesArray("AliMUONHit",1000);
     }
    //it must be under fHits creation
    AliDetector::MakeBranch(option);
    
    if (cD  && fLoader->TreeD()) {
      //
      // one branch for digits per chamber
      // 
      Int_t i;
      if (fDchambers  == 0x0) 
        {
          fDchambers = new TObjArray(AliMUONConstants::NCh());
          for (Int_t i=0; i<AliMUONConstants::NCh() ;i++) {
              fDchambers->AddAt(new TClonesArray("AliMUONDigit",10000),i); 
          }
        }
    
      for (i=0; i<AliMUONConstants::NCh() ;i++) 
       {
        sprintf(branchname,"%sDigits%d",GetName(),i+1);	
        MakeBranchInTree(fLoader->TreeD(), branchname, &((*fDchambers)[i]), kBufferSize, 0);
        printf("Making Branch %s for digits in chamber %d\n",branchname,i+1);
       }
    }
    
    if (cR  && fLoader->TreeR()) {
      //     
      // one branch for raw clusters per chamber
      //  
      printf("Make Branch - TreeR address %p\n",fLoader->TreeR());
      
      Int_t i;
      if (fRawClusters == 0x0)
      {
        fRawClusters = new TObjArray(AliMUONConstants::NTrackingCh());
        for (Int_t i=0; i<AliMUONConstants::NTrackingCh();i++) {
            fRawClusters->AddAt(new TClonesArray("AliMUONRawCluster",10000),i); 
        }
      }

      for (i=0; i<AliMUONConstants::NTrackingCh() ;i++) 
       {
         sprintf(branchname,"%sRawClusters%d",GetName(),i+1);	
         MakeBranchInTree(fLoader->TreeR(), branchname, &((*fRawClusters)[i]), kBufferSize, 0);
         printf("Making Branch %s for raw clusters in chamber %d\n",branchname,i+1);
      }
      //
      // one branch for global trigger
      //
      sprintf(branchname,"%sGlobalTrigger",GetName());
      
      if (fGlobalTrigger == 0x0) {
        fGlobalTrigger = new TClonesArray("AliMUONGlobalTrigger",1); 
      }
      MakeBranchInTree(fLoader->TreeR(), branchname, &fGlobalTrigger, kBufferSize, 0);
      printf("Making Branch %s for Global Trigger\n",branchname);
      //
      // one branch for local trigger
      //  
      sprintf(branchname,"%sLocalTrigger",GetName());
      
      if (fLocalTrigger == 0x0) {
        fLocalTrigger  = new TClonesArray("AliMUONLocalTrigger",234);
      }
      
      MakeBranchInTree(fLoader->TreeR(), branchname, &fLocalTrigger, kBufferSize, 0);
      printf("Making Branch %s for Local Trigger\n",branchname);
   }
}

//___________________________________________
void AliMUON::SetTreeAddress()
{
  // Set branch address for the Hits and Digits Tree.
  char branchname[30];

  TBranch *branch;
  TTree *treeH = fLoader->TreeH();
  TTree *treeD = fLoader->TreeD();
  TTree *treeR = fLoader->TreeR();

  if (treeH) {
    if (fPadHits == 0x0) fPadHits = new TClonesArray("AliMUONPadHit",10000);
    if (fPadHits) {
      branch = treeH->GetBranch("MUONCluster");
      if (branch) branch->SetAddress(&fPadHits);
    }
    if (fHits == 0x0) fHits     = new TClonesArray("AliMUONHit",1000);
  }
  //it must be under fHits creation
  AliDetector::SetTreeAddress();

  if (treeD) {
      if (fDchambers == 0x0) 
        {
          fDchambers = new TObjArray(AliMUONConstants::NCh());
          for (Int_t i=0; i<AliMUONConstants::NCh() ;i++) {
              fDchambers->AddAt(new TClonesArray("AliMUONDigit",10000),i); 
          }
        }
      for (int i=0; i<AliMUONConstants::NCh(); i++) {
	  sprintf(branchname,"%sDigits%d",GetName(),i+1);
                        
                      if (fDchambers) {
	      branch = treeD->GetBranch(branchname);
	      if (branch) branch->SetAddress(&((*fDchambers)[i]));
	  }
      }
  }

  // printf("SetTreeAddress --- treeR address  %p \n",treeR);

  if (treeR) {
      if (fRawClusters == 0x0)
      {
        fRawClusters = new TObjArray(AliMUONConstants::NTrackingCh());
        for (Int_t i=0; i<AliMUONConstants::NTrackingCh();i++) {
            fRawClusters->AddAt(new TClonesArray("AliMUONRawCluster",10000),i); 
        }
      }
      
      for (int i=0; i<AliMUONConstants::NTrackingCh(); i++) {
	  sprintf(branchname,"%sRawClusters%d",GetName(),i+1);
	  if (fRawClusters) {
	      branch = treeR->GetBranch(branchname);
	      if (branch) branch->SetAddress(&((*fRawClusters)[i]));
	  }
      }

      if (fLocalTrigger == 0x0) {
        fLocalTrigger  = new TClonesArray("AliMUONLocalTrigger",234);
      }

      if (fLocalTrigger) {
	branch = treeR->GetBranch("MUONLocalTrigger");
	if (branch) branch->SetAddress(&fLocalTrigger);
      }
      
      if (fGlobalTrigger == 0x0) {
        fGlobalTrigger = new TClonesArray("AliMUONGlobalTrigger",1); 
      }
      
      if (fGlobalTrigger) {
	branch = treeR->GetBranch("MUONGlobalTrigger");
	if (branch) branch->SetAddress(&fGlobalTrigger);
      }
  }
}
//___________________________________________
void AliMUON::ResetHits()
{
  // Reset number of clusters and the cluster array for this detector
  AliDetector::ResetHits();
  fNPadHits = 0;
  if (fPadHits) fPadHits->Clear();
}

//____________________________________________
void AliMUON::ResetDigits()
{
    //
    // Reset number of digits and the digits array for this detector
    //
    if (fDchambers == 0x0) return;
    for ( int i=0;i<AliMUONConstants::NCh();i++ ) {
      //PH	if ((*fDchambers)[i])    ((TClonesArray*)(*fDchambers)[i])->Clear();
	if ((*fDchambers)[i])    ((TClonesArray*)fDchambers->At(i))->Clear();
	if (fNdch)  fNdch[i]=0;
    }
}
//____________________________________________
void AliMUON::ResetRawClusters()
{
    //
    // Reset number of raw clusters and the raw clust array for this detector
    //
    for ( int i=0;i<AliMUONConstants::NTrackingCh();i++ ) {
      //PH	if ((*fRawClusters)[i])    ((TClonesArray*)(*fRawClusters)[i])->Clear();
	if ((*fRawClusters)[i])    ((TClonesArray*)fRawClusters->At(i))->Clear();
	if (fNrawch)  fNrawch[i]=0;
    }
}

//____________________________________________
void AliMUON::ResetTrigger()
{
  //  Reset Local and Global Trigger 
  fNGlobalTrigger = 0;
  if (fGlobalTrigger) fGlobalTrigger->Clear();
  fNLocalTrigger = 0;
  if (fLocalTrigger) fLocalTrigger->Clear();
}

//____________________________________________
void AliMUON::SetPadSize(Int_t id, Int_t isec, Float_t p1, Float_t p2)
{
// Set the pad size for chamber id and cathode isec
    Int_t i=2*(id-1);
    //PH    ((AliMUONChamber*) (*fChambers)[i])  ->SetPadSize(isec,p1,p2);
    //PH    ((AliMUONChamber*) (*fChambers)[i+1])->SetPadSize(isec,p1,p2);
    ((AliMUONChamber*) fChambers->At(i))  ->SetPadSize(isec,p1,p2);
    ((AliMUONChamber*) fChambers->At(i+1))->SetPadSize(isec,p1,p2);
}

//___________________________________________
void AliMUON::SetChambersZ(const Float_t *Z)
{
  // Set Z values for all chambers (tracking and trigger)
  // from the array pointed to by "Z"
    for (Int_t ch = 0; ch < AliMUONConstants::NCh(); ch++)
      //PH	((AliMUONChamber*) ((*fChambers)[ch]))->SetZ(Z[ch]);
	((AliMUONChamber*) fChambers->At(ch))->SetZ(Z[ch]);
    return;
}

//___________________________________________
void AliMUON::SetChambersZToDefault()
{
  // Set Z values for all chambers (tracking and trigger)
  // to default values
  SetChambersZ(AliMUONConstants::DefaultChamberZ());
  return;
}

//___________________________________________
void AliMUON::SetChargeSlope(Int_t id, Float_t p1)
{
// Set the inverse charge slope for chamber id
    Int_t i=2*(id-1);
    //PH    ((AliMUONChamber*) (*fChambers)[i])->SetChargeSlope(p1);
    //PH    ((AliMUONChamber*) (*fChambers)[i+1])->SetChargeSlope(p1);
    ((AliMUONChamber*) fChambers->At(i))->SetChargeSlope(p1);
    ((AliMUONChamber*) fChambers->At(i+1))->SetChargeSlope(p1);
}

//___________________________________________
void AliMUON::SetChargeSpread(Int_t id, Float_t p1, Float_t p2)
{
// Set sigma of charge spread for chamber id
    Int_t i=2*(id-1);
    //PH    ((AliMUONChamber*) fChambers->At(i))->SetChargeSpread(p1,p2);
    //PH    ((AliMUONChamber*) fChambers->Ati+1])->SetChargeSpread(p1,p2);
    ((AliMUONChamber*) fChambers->At(i))->SetChargeSpread(p1,p2);
    ((AliMUONChamber*) fChambers->At(i+1))->SetChargeSpread(p1,p2);
}

//___________________________________________
void AliMUON::SetSigmaIntegration(Int_t id, Float_t p1)
{
// Set integration limits for charge spread
    Int_t i=2*(id-1);
    //PH    ((AliMUONChamber*) (*fChambers)[i])->SetSigmaIntegration(p1);
    //PH    ((AliMUONChamber*) (*fChambers)[i+1])->SetSigmaIntegration(p1);
    ((AliMUONChamber*) fChambers->At(i))->SetSigmaIntegration(p1);
    ((AliMUONChamber*) fChambers->At(i+1))->SetSigmaIntegration(p1);
}

//___________________________________________
void AliMUON::SetMaxAdc(Int_t id, Int_t p1)
{
// Set maximum number for ADCcounts (saturation)
    Int_t i=2*(id-1);
    //PH    ((AliMUONChamber*) (*fChambers)[i])->SetMaxAdc(p1);
    //PH    ((AliMUONChamber*) (*fChambers)[i+1])->SetMaxAdc(p1);
    ((AliMUONChamber*) fChambers->At(i))->SetMaxAdc(p1);
    ((AliMUONChamber*) fChambers->At(i+1))->SetMaxAdc(p1);
}

//___________________________________________
void AliMUON::SetMaxStepGas(Float_t p1)
{
// Set stepsize in gas
     fMaxStepGas=p1;
}

//___________________________________________
void AliMUON::SetMaxStepAlu(Float_t p1)
{
// Set step size in Alu
    fMaxStepAlu=p1;
}

//___________________________________________
void AliMUON::SetMaxDestepGas(Float_t p1)
{
// Set maximum step size in Gas
    fMaxDestepGas=p1;
}

//___________________________________________
void AliMUON::SetMaxDestepAlu(Float_t p1)
{
// Set maximum step size in Alu
    fMaxDestepAlu=p1;
}
//___________________________________________
void AliMUON::SetAcceptance(Bool_t acc, Float_t angmin, Float_t angmax)
{
// Set acceptance cuts 
   fAccCut=acc;
   fAccMin=angmin*TMath::Pi()/180;
   fAccMax=angmax*TMath::Pi()/180;
   Int_t ch;
   if (acc) {
       for (Int_t st = 0; st < AliMUONConstants::NCh() / 2; st++) {
	   // Loop over 2 chambers in the station
	   for (Int_t stCH = 0; stCH < 2; stCH++) {
	       ch = 2 * st + stCH;
//         Set chamber inner and outer radius according to acceptance cuts
	       Chamber(ch).SetRInner(AliMUONConstants::DefaultChamberZ(ch)*TMath::Tan(fAccMin));
	       Chamber(ch).SetROuter(AliMUONConstants::DefaultChamberZ(ch)*TMath::Tan(fAccMax));
	   } // chamber loop
       } // station loop
   }
}
//___________________________________________
void   AliMUON::SetSegmentationModel(Int_t id, Int_t isec, AliSegmentation *segmentation)
{
// Set the segmentation for chamber id cathode isec
  //PH    ((AliMUONChamber*) (*fChambers)[id])->SetSegmentationModel(isec, segmentation);
    ((AliMUONChamber*) fChambers->At(id))->SetSegmentationModel(isec, segmentation);

}
//___________________________________________
void   AliMUON::SetResponseModel(Int_t id, AliMUONResponse *response)
{
// Set the response for chamber id
  //PH    ((AliMUONChamber*) (*fChambers)[id])->SetResponseModel(response);
    ((AliMUONChamber*) fChambers->At(id))->SetResponseModel(response);
}

void   AliMUON::SetReconstructionModel(Int_t id, AliMUONClusterFinderVS *reconst)
{
// Set ClusterFinder for chamber id
  //PH    ((AliMUONChamber*) (*fChambers)[id])->SetReconstructionModel(reconst);
    ((AliMUONChamber*) fChambers->At(id))->SetReconstructionModel(reconst);
}

void   AliMUON::SetNsec(Int_t id, Int_t nsec)
{
// Set number of segmented cathods for chamber id
  //PH    ((AliMUONChamber*) (*fChambers)[id])->SetNsec(nsec);
    ((AliMUONChamber*) fChambers->At(id))->SetNsec(nsec);
}

//___________________________________________
void AliMUON::SDigits2Digits()
{

// write TreeD here 

    if (!fMerger) {
      if (gAlice->GetDebug()>0) {
	cerr<<"AliMUON::SDigits2Digits: create default AliMUONMerger "<<endl;
	cerr<<" no merging, just digitization of 1 event will be done"<<endl;
      }
      fMerger = new AliMUONMerger();
    }
    fMerger->Init();
    fMerger->Digitise();
    char hname[30];
    //    sprintf(hname,"TreeD%d",fLoader->GetHeader()->GetEvent());
    fLoader->TreeD()->Write(hname,TObject::kOverwrite);
    fLoader->TreeD()->Reset();
}


//__________________________________________________________________________
AliLoader* AliMUON::MakeLoader(const char* topfoldername)
{ 
//builds standard getter (AliLoader type)
//if detector wants to use castomized getter, it must overload this method

 if (GetDebug())
   Info("MakeLoader",
        "Creating standard getter for detector %s. Top folder is %s.",
         GetName(),topfoldername);
     
 fLoader = new AliMUONLoader(GetName(),topfoldername);
 return fLoader;
}
//__________________________________________________________________________
// To be removed
void AliMUON::MakePadHits(Float_t xhit,Float_t yhit, Float_t zhit,
			  Float_t eloss, Float_t tof,  Int_t idvol)
{
//
//  Calls the charge disintegration method of the current chamber and adds
//  the simulated cluster to the root treee 
//
    Int_t clhits[7];
    Float_t newclust[6][500];
    Int_t nnew;
    
    
//
//  Integrated pulse height on chamber

    
    clhits[0]=fNhits+1;
//
//
//    if (idvol == 6) printf("\n ->Disintegration %f %f %f", xhit, yhit, eloss );
    

    //PH    ((AliMUONChamber*) (*fChambers)[idvol])
    ((AliMUONChamber*) fChambers->At(idvol))
	->DisIntegration(eloss, tof, xhit, yhit, zhit, nnew, newclust);
    Int_t ic=0;
//    if (idvol == 6) printf("\n nnew  %d \n", nnew);
//
//  Add new clusters
    for (Int_t i=0; i<nnew; i++) {
	if (Int_t(newclust[3][i]) > 0) {
	    ic++;
// Cathode plane
	    clhits[1] = Int_t(newclust[5][i]);
//  Cluster Charge
	    clhits[2] = Int_t(newclust[0][i]);
//  Pad: ix
	    clhits[3] = Int_t(newclust[1][i]);
//  Pad: iy 
	    clhits[4] = Int_t(newclust[2][i]);
//  Pad: charge
	    clhits[5] = Int_t(newclust[3][i]);
//  Pad: chamber sector
	    clhits[6] = Int_t(newclust[4][i]);
	    
	    AddPadHit(clhits);
	}
    }
}

//___________________________________________
void AliMUON::Trigger(Int_t nev){
// call the Trigger Algorithm and fill TreeR

  Int_t singlePlus[3]  = {0,0,0}; 
  Int_t singleMinus[3] = {0,0,0}; 
  Int_t singleUndef[3] = {0,0,0};
  Int_t pairUnlike[3]  = {0,0,0}; 
  Int_t pairLike[3]    = {0,0,0};

  ResetTrigger();
  AliMUONTriggerDecision* decision= new AliMUONTriggerDecision(1);
  decision->Trigger();   
  decision->GetGlobalTrigger(singlePlus, singleMinus, singleUndef,
			     pairUnlike, pairLike);

// add a local trigger in the list 
  AddGlobalTrigger(singlePlus, singleMinus, singleUndef, pairUnlike, pairLike);
  Int_t i;
  
  for (Int_t icirc=0; icirc<AliMUONConstants::NTriggerCircuit(); icirc++) { 
      if(decision->GetITrigger(icirc)==1) {
	  Int_t localtr[7]={0,0,0,0,0,0,0};      
	  Int_t loLpt[2]={0,0}; Int_t loHpt[2]={0,0}; Int_t loApt[2]={0,0};
	  decision->GetLutOutput(icirc, loLpt, loHpt, loApt);
	  localtr[0] = icirc;
	  localtr[1] = decision->GetStripX11(icirc);
	  localtr[2] = decision->GetDev(icirc);
	  localtr[3] = decision->GetStripY11(icirc);
	  for (i=0; i<2; i++) {    // convert the Lut output in 1 digit 
	      localtr[4] = localtr[4]+Int_t(loLpt[i]*TMath::Power(2,i));
	      localtr[5] = localtr[5]+Int_t(loHpt[i]*TMath::Power(2,i));
	      localtr[6] = localtr[6]+Int_t(loApt[i]*TMath::Power(2,i));
	  }
	  AddLocalTrigger(localtr);  // add a local trigger in the list
      }
  }

  delete decision;

  fLoader->TreeR()->Fill();
//  char hname[30];
//  sprintf(hname,"TreeR%d",nev);
//  fLoader->TreeR()->Write(hname,TObject::kOverwrite);
//  fLoader->TreeR()->Reset();
  fLoader->WriteRecPoints("OVERWRITE");
  ResetTrigger();
  
  printf("\n End of trigger for event %d", nev);
}


//____________________________________________
void AliMUON::Digits2Reco()
{
  FindClusters();
  Int_t nev = gAlice->GetHeader()->GetEvent();
  fLoader->TreeR()->Fill();
  // char hname[30];
  // sprintf(hname,"TreeR%d", nev);
  //fLoader->TreeR()->Write(hname);
  //fLoader->TreeR()->Reset();
  fLoader->WriteRecPoints("OVERWRITE");
  ResetRawClusters();        
  printf("\n End of cluster finding for event %d", nev);
}

void AliMUON::FindClusters()
{
//
//  Perform cluster finding
//
    TClonesArray *dig1, *dig2;
    Int_t ndig, k;
    dig1 = new TClonesArray("AliMUONDigit",1000);
    dig2 = new TClonesArray("AliMUONDigit",1000);
    AliMUONDigit *digit;
// Loop on chambers and on cathode planes
//
    ResetRawClusters();        
    TClonesArray * muonDigits;

    for (Int_t ich = 0; ich < 10; ich++) {
      //PH	AliMUONChamber* iChamber = (AliMUONChamber*) (*fChambers)[ich];
	AliMUONChamber* iChamber = (AliMUONChamber*) fChambers->At(ich);
	AliMUONClusterFinderVS* rec = iChamber->ReconstructionModel();
    
	ResetDigits();
	fLoader->TreeD()->GetEvent(0);
	//TClonesArray *
	muonDigits = (TClonesArray *) Dchambers()->At(ich);
	ndig=muonDigits->GetEntriesFast();
	printf("\n 1 Found %d digits in %p chamber %d", ndig, muonDigits,ich);
	TClonesArray &lhits1 = *dig1;
	Int_t n = 0;
	for (k = 0; k < ndig; k++) {
	    digit = (AliMUONDigit*) muonDigits->UncheckedAt(k);
	    if (rec->TestTrack(digit->Track(0)))
		new(lhits1[n++]) AliMUONDigit(*digit);
	}
	ResetDigits();
	fLoader->TreeD()->GetEvent(1);
	//muonDigits  = this->DigitsAddress(ich);
	muonDigits = (TClonesArray *) Dchambers()->At(ich);
	ndig=muonDigits->GetEntriesFast();
	printf("\n 2 Found %d digits in %p %d", ndig, muonDigits, ich);
	TClonesArray &lhits2 = *dig2;
	n=0;
	
	for (k=0; k<ndig; k++) {
	    digit= (AliMUONDigit*) muonDigits->UncheckedAt(k);
	    if (rec->TestTrack(digit->Track(0)))
	    new(lhits2[n++]) AliMUONDigit(*digit);
	}

	if (rec) {	 
	    AliMUONClusterInput::Instance()->SetDigits(ich, dig1, dig2);
	    rec->FindRawClusters();
	}
	dig1->Delete();
	dig2->Delete();
    } // for ich
    delete dig1;
    delete dig2;
}
 
#ifdef never
void AliMUON::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliMUON.
      AliMUONChamber        *iChamber;
      AliMUONTriggerCircuit *iTriggerCircuit;
      AliSegmentation       *segmentation;
      AliMUONResponse       *response;
      TClonesArray          *digitsaddress;
      TClonesArray          *rawcladdress;
      Int_t i;
      if (R__b.IsReading()) {
	  Version_t R__v = R__b.ReadVersion(); if (R__v) { }
	  AliDetector::Streamer(R__b);
	  R__b >> fNPadHits;
	  R__b >> fPadHits; // diff
	  R__b >> fNLocalTrigger;       
	  R__b >> fLocalTrigger;       
	  R__b >> fNGlobalTrigger;       
	  R__b >> fGlobalTrigger;   
	  R__b >> fDchambers;
	  R__b >> fRawClusters;
	  R__b.ReadArray(fNdch);
	  R__b.ReadArray(fNrawch);
	  R__b >> fAccCut;
	  R__b >> fAccMin;
	  R__b >> fAccMax; 
	  R__b >> fChambers;
	  R__b >> fTriggerCircuits;
	  for (i =0; i<AliMUONConstants::NTriggerCircuit(); i++) {
	      iTriggerCircuit=(AliMUONTriggerCircuit*) (*fTriggerCircuits)[i];
	      iTriggerCircuit->Streamer(R__b);
	  }
// Stream chamber related information
	  for (i =0; i<AliMUONConstants::NCh(); i++) {
	      iChamber=(AliMUONChamber*) (*fChambers)[i];
	      iChamber->Streamer(R__b);
	      if (iChamber->Nsec()==1) {
		  segmentation=iChamber->SegmentationModel(1);
		  if (segmentation)
		      segmentation->Streamer(R__b);
	      } else {
		  segmentation=iChamber->SegmentationModel(1);
		  if (segmentation)
		      segmentation->Streamer(R__b);
		  if (segmentation)
		      segmentation=iChamber->SegmentationModel(2);
		  segmentation->Streamer(R__b);
	      }
	      response=iChamber->ResponseModel();
	      if (response)
		  response->Streamer(R__b);	  
	      digitsaddress=(TClonesArray*) (*fDchambers)[i];
	      digitsaddress->Streamer(R__b);
	      if (i < AliMUONConstants::NTrackingCh()) {
		  rawcladdress=(TClonesArray*) (*fRawClusters)[i];
		  rawcladdress->Streamer(R__b);
	      }
	  }
	  
      } else {
	  R__b.WriteVersion(AliMUON::IsA());
	  AliDetector::Streamer(R__b);
	  R__b << fNPadHits;
	  R__b << fPadHits; // diff
	  R__b << fNLocalTrigger;       
	  R__b << fLocalTrigger;       
	  R__b << fNGlobalTrigger;       
	  R__b << fGlobalTrigger; 
	  R__b << fDchambers;
	  R__b << fRawClusters;
	  R__b.WriteArray(fNdch, AliMUONConstants::NCh());
	  R__b.WriteArray(fNrawch, AliMUONConstants::NTrackingCh());
	  
	  R__b << fAccCut;
	  R__b << fAccMin;
	  R__b << fAccMax; 
	  
	  R__b << fChambers;
	  R__b << fTriggerCircuits;
	  for (i =0; i<AliMUONConstants::NTriggerCircuit(); i++) {
	      iTriggerCircuit=(AliMUONTriggerCircuit*) (*fTriggerCircuits)[i];
	      iTriggerCircuit->Streamer(R__b);
	  }
	  for (i =0; i<AliMUONConstants::NCh(); i++) {
	      iChamber=(AliMUONChamber*) (*fChambers)[i];
	      iChamber->Streamer(R__b);
	      if (iChamber->Nsec()==1) {
		  segmentation=iChamber->SegmentationModel(1);
		  if (segmentation)
		      segmentation->Streamer(R__b);
	      } else {
		  segmentation=iChamber->SegmentationModel(1);
		  if (segmentation)
		      segmentation->Streamer(R__b);
		  segmentation=iChamber->SegmentationModel(2);
		  if (segmentation)
		      segmentation->Streamer(R__b);
	      }
	      response=iChamber->ResponseModel();
	      if (response)
		  response->Streamer(R__b);
	      digitsaddress=(TClonesArray*) (*fDchambers)[i];
	      digitsaddress->Streamer(R__b);
	      if (i < AliMUONConstants::NTrackingCh()) {
		  rawcladdress=(TClonesArray*) (*fRawClusters)[i];
		  rawcladdress->Streamer(R__b);
	      }
	  }
      }
}
#endif

AliMUONPadHit* AliMUON::FirstPad(AliMUONHit*  hit, TClonesArray *clusters) 
{
// to be removed
    // Initialise the pad iterator
    // Return the address of the first padhit for hit
    TClonesArray *theClusters = clusters;
    Int_t nclust = theClusters->GetEntriesFast();
    if (nclust && hit->PHlast() > 0) {
	AliMUON::fMaxIterPad=hit->PHlast();
	AliMUON::fCurIterPad=hit->PHfirst();
	return (AliMUONPadHit*) clusters->UncheckedAt(AliMUON::fCurIterPad-1);
    } else {
	return 0;
    }
}

AliMUONPadHit* AliMUON::NextPad(TClonesArray *clusters) 
{
  // To be removed
// Get next pad (in iterator) 
//
    AliMUON::fCurIterPad++;
    if (AliMUON::fCurIterPad <= AliMUON::fMaxIterPad) {
	return (AliMUONPadHit*) clusters->UncheckedAt(AliMUON::fCurIterPad-1);
    } else {
	return 0;
    }
}


AliMUONRawCluster *AliMUON::RawCluster(Int_t ichamber, Int_t icathod, Int_t icluster)
{
//
//  Return rawcluster (icluster) for chamber ichamber and cathode icathod
//  Obsolete ??
    TClonesArray *muonRawCluster  = RawClustAddress(ichamber);
    ResetRawClusters();
    TTree *treeR = fLoader->TreeR();
    Int_t nent=(Int_t)treeR->GetEntries();
    treeR->GetEvent(nent-2+icathod-1);
    //treeR->GetEvent(icathod);
    //Int_t nrawcl = (Int_t)muonRawCluster->GetEntriesFast();

    AliMUONRawCluster * mRaw = (AliMUONRawCluster*)muonRawCluster->UncheckedAt(icluster);
    //printf("RawCluster _ nent nrawcl icluster mRaw %d %d %d%p\n",nent,nrawcl,icluster,mRaw);
    
    return  mRaw;
}
 
void   AliMUON::SetMerger(AliMUONMerger* merger)
{
// Set pointer to merger 
    fMerger = merger;
}

AliMUONMerger*  AliMUON::Merger()
{
// Return pointer to merger
    return fMerger;
}



AliMUON& AliMUON::operator = (const AliMUON& /*rhs*/)
{
// copy operator
// dummy version
    return *this;
}

////////////////////////////////////////////////////////////////////////
void AliMUON::MakeBranchInTreeD(TTree *treeD, const char *file)
{
    //
    // Create TreeD branches for the MUON.
    //

  const Int_t kBufferSize = 4000;
  char branchname[30];
    
  if (fDchambers  == 0x0)   {
    fDchambers = new TObjArray(AliMUONConstants::NCh());
    for (Int_t i=0; i<AliMUONConstants::NCh() ;i++) {
      fDchambers->AddAt(new TClonesArray("AliMUONDigit",10000),i); 
    }
  }
  //
  // one branch for digits per chamber
  // 
  for (Int_t i=0; i<AliMUONConstants::NCh() ;i++) {
    sprintf(branchname,"%sDigits%d",GetName(),i+1);	
    if (fDchambers && treeD) {
      MakeBranchInTree(treeD, 
		       branchname, &((*fDchambers)[i]), kBufferSize, file);
//      printf("Making Branch %s for digits in chamber %d\n",branchname,i+1);
    }
  }
}

//___________________________________________

