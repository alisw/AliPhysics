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
//__________________________________________________________________
AliMUON::AliMUON()
{
// Default Constructor
//
    fNCh             = 0;
    fNTrackingCh     = 0;
    fIshunt          = 0;
    fChambers        = 0;
    fTriggerCircuits = 0;
    fAccMin          = 0.;
    fAccMax          = 0.;   
    fAccCut          = kFALSE;
    fMerger          = 0;
    fFileName        = 0;
    fMUONData        = 0;
}
//__________________________________________________________________
AliMUON::AliMUON(const char *name, const char *title)
  : AliDetector(name,title)
{
//Begin_Html
/*
<img src="gif/alimuon.gif">
*/
//End_Html
  fMUONData  = 0x0;
  fIshunt     =  0;

  fNCh             = AliMUONConstants::NCh(); 
  fNTrackingCh     = AliMUONConstants::NTrackingCh();

  SetMarkerColor(kRed);//
//
// Creating List of Chambers
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
	if (ch < AliMUONConstants::NTrackingCh()) {
	  fChambers->AddAt(new AliMUONChamber(ch),ch);
	} else {
	  fChambers->AddAt(new AliMUONChamberTrigger(ch),ch);
	}
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
    
    fMaxStepGas=0.01; 
    fMaxStepAlu=0.1; 
    fMaxDestepGas=-1;
    fMaxDestepAlu=-1;
    
    fMaxIterPad   = 0;
    fCurIterPad   = 0;
    
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
//____________________________________________________________________
AliMUON::AliMUON(const AliMUON& rMUON):AliDetector(rMUON)
{
// Dummy copy constructor
    ;
    
}
//____________________________________________________________________
AliMUON::~AliMUON()
{
// Destructor
  if(fDebug) printf("%s: Calling AliMUON destructor !!!\n",ClassName());
  fIshunt  = 0;
  if (fMerger) delete fMerger;
}
//____________________________________________________________________
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
//___________________________________________________________________
Int_t AliMUON::DistancetoPrimitive(Int_t , Int_t )
{
  return 9999;
}
//__________________________________________________________________
void  AliMUON::SetTreeAddress()
{
  GetMUONData()->SetLoader(fLoader); 
  GetMUONData()->SetTreeAddress("H,D,RC");
}

//____________________________________________________________________
void AliMUON::SetPadSize(Int_t id, Int_t isec, Float_t p1, Float_t p2)
{
// Set the pad size for chamber id and cathode isec
    Int_t i=2*(id-1);
    ((AliMUONChamber*) fChambers->At(i))  ->SetPadSize(isec,p1,p2);
    ((AliMUONChamber*) fChambers->At(i+1))->SetPadSize(isec,p1,p2);
}

//___________________________________________
void AliMUON::SetChambersZ(const Float_t *Z)
{
  // Set Z values for all chambers (tracking and trigger)
  // from the array pointed to by "Z"
    for (Int_t ch = 0; ch < AliMUONConstants::NCh(); ch++)
	((AliMUONChamber*) fChambers->At(ch))->SetZ(Z[ch]);
    return;
}
//_________________________________________________________________
void AliMUON::SetChambersZToDefault()
{
  // Set Z values for all chambers (tracking and trigger)
  // to default values
  SetChambersZ(AliMUONConstants::DefaultChamberZ());
  return;
}
//_________________________________________________________________
void AliMUON::SetChargeSlope(Int_t id, Float_t p1)
{
// Set the inverse charge slope for chamber id
    Int_t i=2*(id-1);    //PH    ((AliMUONChamber*) (*fChambers)[i])->SetSigmaIntegration(p1);
    //PH    ((AliMUONChamber*) (*fChambers)[i+1])->SetSigmaIntegration(p1);
    ((AliMUONChamber*) fChambers->At(i))->SetChargeSlope(p1);
    ((AliMUONChamber*) fChambers->At(i+1))->SetChargeSlope(p1);
}
//__________________________________________________________________
void AliMUON::SetChargeSpread(Int_t id, Float_t p1, Float_t p2)
{
// Set sigma of charge spread for chamber id
    Int_t i=2*(id-1);
    ((AliMUONChamber*) fChambers->At(i))->SetChargeSpread(p1,p2);
    ((AliMUONChamber*) fChambers->At(i+1))->SetChargeSpread(p1,p2);
}
//___________________________________________________________________
void AliMUON::SetSigmaIntegration(Int_t id, Float_t p1)
{
// Set integration limits for charge spread
    Int_t i=2*(id-1);
    ((AliMUONChamber*) fChambers->At(i))->SetSigmaIntegration(p1);
    ((AliMUONChamber*) fChambers->At(i+1))->SetSigmaIntegration(p1);
}

//__________________________________________________________________
void AliMUON::SetMaxAdc(Int_t id, Int_t p1)
{
// Set maximum number for ADCcounts (saturation)
    Int_t i=2*(id-1);
    ((AliMUONChamber*) fChambers->At(i))->SetMaxAdc(p1);
    ((AliMUONChamber*) fChambers->At(i+1))->SetMaxAdc(p1);
}
//__________________________________________________________________
void AliMUON::SetMaxStepGas(Float_t p1)
{
// Set stepsize in gas
  fMaxStepGas=p1;
}
//__________________________________________________________________
void AliMUON::SetMaxStepAlu(Float_t p1)
{
// Set step size in Alu
    fMaxStepAlu=p1;
}
//__________________________________________________________________
void AliMUON::SetMaxDestepGas(Float_t p1)
{
// Set maximum step size in Gas
    fMaxDestepGas=p1;
}
//__________________________________________________________________
void AliMUON::SetMaxDestepAlu(Float_t p1)
{
// Set maximum step size in Alu
  fMaxDestepAlu=p1;
}
//___________________________________________________________________
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
//____________________________________________________________________
void   AliMUON::SetSegmentationModel(Int_t id, Int_t isec, AliSegmentation *segmentation)
{
// Set the segmentation for chamber id cathode isec
    ((AliMUONChamber*) fChambers->At(id))->SetSegmentationModel(isec, segmentation);

}
//____________________________________________________________________
void   AliMUON::SetResponseModel(Int_t id, AliMUONResponse *response)
{
// Set the response for chamber id
    ((AliMUONChamber*) fChambers->At(id))->SetResponseModel(response);
}
//____________________________________________________________________
void   AliMUON::SetReconstructionModel(Int_t id, AliMUONClusterFinderVS *reconst)
{
// Set ClusterFinder for chamber id
    ((AliMUONChamber*) fChambers->At(id))->SetReconstructionModel(reconst);
}
//____________________________________________________________________
void   AliMUON::SetNsec(Int_t id, Int_t nsec)
{
// Set number of segmented cathods for chamber id
    ((AliMUONChamber*) fChambers->At(id))->SetNsec(nsec);
}
//_____________________________________________________________________
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

//_______________________________________________________________________
AliLoader* AliMUON::MakeLoader(const char* topfoldername)
{ 
//builds standard getter (AliLoader type)
//if detector wants to use castomized getter, it must overload this method

 if (GetDebug())
   Info("MakeLoader",
        "Creating standard getter for detector %s. Top folder is %s.",
         GetName(),topfoldername);
     
 fLoader   = new AliLoader(GetName(),topfoldername);
 fMUONData = new AliMUONData(fLoader,GetName(),GetName()); 
 return fLoader;
}

//_______________________________________________________________________
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
  GetMUONData()->AddGlobalTrigger(singlePlus, singleMinus, singleUndef, pairUnlike, pairLike);
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
      GetMUONData()->AddLocalTrigger(localtr);  // add a local trigger in the list
    }
  }
  
  delete decision;

  //  fLoader->TreeR()->Fill();
  GetMUONData()->Fill("GLT"); //Filling Global and Local Trigger GLT
  //  char hname[30];
  //  sprintf(hname,"TreeR%d",nev);
  //  fLoader->TreeR()->Write(hname,TObject::kOverwrite);
    //  fLoader->TreeR()->Reset();
  fLoader->WriteRecPoints("OVERWRITE");
  
  printf("\n End of trigger for event %d", nev);
}

//____________________________________________________________________
void AliMUON::Digits2Reco()
{
  FindClusters();
  Int_t nev = gAlice->GetHeader()->GetEvent();
  GetMUONData()->Fill("RC"); //Filling Reconstructed Cluster
  fLoader->WriteRecPoints("OVERWRITE");
  GetMUONData()->ResetRawClusters();        
  Info("Digits2Reco","End of cluster finding for event %d", nev);
}
//____________________________________________________________________
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
	GetMUONData()->GetCathode(0);
	//TClonesArray *
	muonDigits = GetMUONData()->Digits(ich); 
	ndig=muonDigits->GetEntriesFast();
	printf("\n 1 Found %d digits in %p chamber %d", ndig, muonDigits,ich);
	TClonesArray &lhits1 = *dig1;
	Int_t n = 0;
	for (k = 0; k < ndig; k++) {
	    digit = (AliMUONDigit*) muonDigits->UncheckedAt(k);
	    if (rec->TestTrack(digit->Track(0)))
		new(lhits1[n++]) AliMUONDigit(*digit);
	}
	GetMUONData()->ResetDigits();
	GetMUONData()->GetCathode(1);
	muonDigits =  GetMUONData()->Digits(ich);  
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
//______________________________________________________________________
#ifdef never
void AliMUON::Streamer(TBuffer &R__b)_
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
//_______________________________________________________________________
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
//_______________________________________________________________________
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
//_______________________________________________________________________

AliMUONRawCluster *AliMUON::RawCluster(Int_t ichamber, Int_t icathod, Int_t icluster)
{
//
//  Return rawcluster (icluster) for chamber ichamber and cathode icathod
//  Obsolete ??
    TClonesArray *muonRawCluster  = GetMUONData()->RawClusters(ichamber);
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
//________________________________________________________________________
void   AliMUON::SetMerger(AliMUONMerger* merger)
{
// Set pointer to merger 
    fMerger = merger;
}
//________________________________________________________________________
AliMUONMerger*  AliMUON::Merger()
{
// Return pointer to merger
    return fMerger;
}
//________________________________________________________________________
AliMUON& AliMUON::operator = (const AliMUON& /*rhs*/)
{
// copy operator
// dummy version
    return *this;
}

