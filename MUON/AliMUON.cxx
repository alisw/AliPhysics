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

//#include "AliHeader.h"
#include "AliLoader.h"
#include "AliRunDigitizer.h"
#include "AliMC.h"
#include "AliRun.h"	
#include "AliMUON.h"
#include "AliMUONChamberTrigger.h"
#include "AliMUONConstants.h"
#include "AliMUONHit.h"
#include "AliMUONMerger.h"	
#include "AliMUONPadHit.h"
#include "AliMUONRawCluster.h"
#include "AliMUONTransientDigit.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONGeometryBuilder.h"
#include "AliMUONVGeometryBuilder.h"	
#include "AliMUONDigitizerv2.h"
#include "AliMUONSDigitizerv1.h"
#include "AliMUONRawData.h"

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
  : AliDetector(),
    fNCh(0),
    fNTrackingCh(0),
    fMUONData(0),
    fSplitLevel(0),
    fChambers(0),
    fTriggerCircuits(0),
    fGeometryBuilder(0),
    fAccCut(kFALSE),
    fAccMin(0.),
    fAccMax(0.),   
    fMaxStepGas(0.),
    fMaxStepAlu(0.),
    fMaxDestepGas(0.),
    fMaxDestepAlu(0.),
    fMaxIterPad(0),
    fCurIterPad(0),
    fMerger(0)
{
// Default Constructor
//
    fIshunt          = 0;
}

//__________________________________________________________________
AliMUON::AliMUON(const char *name, const char *title)
  : AliDetector(name,title),
    fNCh(AliMUONConstants::NCh()),
    fNTrackingCh(AliMUONConstants::NTrackingCh()),
    fMUONData(0),
    fSplitLevel(0),
    fChambers(0),
    fTriggerCircuits(0),
    fGeometryBuilder(0),
    fAccCut(kFALSE),
    fAccMin(0.),
    fAccMax(0.),   
    fMaxStepGas(0.1),
    fMaxStepAlu(0.1),
    fMaxDestepGas(-1), // Negatives values are ignored by geant3 CONS200 
    fMaxDestepAlu(-1), // in the calculation of the tracking parameters
    fMaxIterPad(0),
    fCurIterPad(0),
    fMerger(0)
{

  fIshunt =  0;

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
	//chamber->SetGid(0);
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
    
    // cp new design of AliMUONTriggerDecision
    fTriggerCircuits = new TObjArray(AliMUONConstants::NTriggerCircuit());
    for (Int_t circ=0; circ<AliMUONConstants::NTriggerCircuit(); circ++) {
      fTriggerCircuits->AddAt(new AliMUONTriggerCircuit(),circ);          
    }
    
    // Geometry builder
    fGeometryBuilder = new AliMUONGeometryBuilder(this);
}

//____________________________________________________________________
AliMUON::AliMUON(const AliMUON& rMUON)
 : AliDetector(rMUON)
{
// Protected copy constructor

  Fatal("AliMUONMergerModule", "Not implemented.");
}

//____________________________________________________________________
AliMUON::~AliMUON()
{
// Destructor
  if(fDebug) printf("%s: Calling AliMUON destructor !!!\n",ClassName());
  fIshunt  = 0;
  if (fMerger) delete fMerger;

  if (fChambers){
    fChambers->Delete();
    delete fChambers;
  }
  if (fTriggerCircuits){
    fTriggerCircuits->Delete();
    delete fTriggerCircuits;
  }
  delete fMUONData;
  delete fGeometryBuilder;
}

//________________________________________________________________________
AliMUON& AliMUON::operator = (const AliMUON& rhs)
{
// Protected assignement operator

  if (this == &rhs) return *this;

  Fatal("operator=", "Not implemented.");
    
  return *this;  
}

//_____________________________________________________________________________
void AliMUON::AddGeometryBuilder(AliMUONVGeometryBuilder* geomBuilder)
{
// Adds the geometry builder to the list
// ---

  fGeometryBuilder->AddBuilder(geomBuilder);
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

//__________________________________________________________________
void  AliMUON::SetTreeAddress()
{
  GetMUONData()->SetLoader(fLoader); 
  //  GetMUONData()->MakeBranch("D,S,RC");
  //  GetMUONData()->SetTreeAddress("H,D,S,RC");
  GetMUONData()->SetTreeAddress("H");
  if (fHits !=  GetMUONData()->Hits())  {
    if ( gAlice->GetMCApp() )
      if ( gAlice->GetMCApp()->GetHitLists() ) {
	fHits = GetMUONData()->Hits();
	gAlice->GetMCApp()->AddHitList(fHits); // For purifyKine, only necessary when Hit list is created in AliMUONData
      }  
  }
  fHits = GetMUONData()->Hits(); // Added by Ivana to use the methods FisrtHit, NextHit of AliDetector    
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
	Chamber(ch).SetRInner(TMath::Abs(AliMUONConstants::DefaultChamberZ(ch)*TMath::Tan(fAccMin)));
	Chamber(ch).SetROuter(TMath::Abs(AliMUONConstants::DefaultChamberZ(ch)*TMath::Tan(fAccMax)));
      } // chamber loop
    } // station loop
  }
}

//____________________________________________________________________
Float_t  AliMUON::GetMaxStepGas() const
{
// Return stepsize in gas
  
  return fMaxStepGas;
}  

//____________________________________________________________________
Float_t  AliMUON::GetMaxStepAlu() const
{
// Return step size in Alu
  
  return fMaxStepAlu;
}
  
//____________________________________________________________________
Float_t  AliMUON::GetMaxDestepGas() const
{
// Return maximum step size in Gas
  
  return fMaxDestepGas;
}
  
//____________________________________________________________________
Float_t  AliMUON::GetMaxDestepAlu() const
{
// Return maximum step size in Gas
  
  return fMaxDestepAlu;
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
void   AliMUON::SetNsec(Int_t id, Int_t nsec)
{
// Set number of segmented cathods for chamber id
    ((AliMUONChamber*) fChambers->At(id))->SetNsec(nsec);
}
//____________________________________________________________________
AliDigitizer* AliMUON::CreateDigitizer(AliRunDigitizer* manager) const
{
  return new AliMUONDigitizerv2(manager);
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

//_____________________________________________________________________
void AliMUON::Hits2SDigits()
{
  // Adaption of AliMUONSDigitizerv1 to be excuted by the AliSimulation framework
  AliRunLoader* runLoader = fLoader->GetRunLoader();
  AliRunDigitizer   * manager = new AliRunDigitizer(1,1);
  manager->SetInputStream(0,runLoader->GetFileName(),AliConfig::GetDefaultEventFolderName());
  AliMUONDigitizer * dMUON   = new AliMUONSDigitizerv1(manager);
  fLoader->LoadHits("READ");
  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++) {
    runLoader->GetEvent(iEvent);
    dMUON->Exec("");
  }
  fLoader->UnloadHits();
}
//_____________________________________________________________________
void AliMUON::Digits2Raw()
{
  // convert digits of the current event to raw data
  AliMUONRawData* rawData;

  rawData = new AliMUONRawData(fLoader);
  if (!rawData->WriteRawData()) Info("MUON","pb writting raw data");
  delete rawData;
  return;
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
 fMUONData->SetSplitLevel(fSplitLevel);
 return fLoader;
}
//_______________________________________________________________________
AliMUONPadHit* AliMUON::FirstPad(AliMUONHit*  hit, TClonesArray *clusters) 
{
// to be removed
    // Initialise the pad iterator
    // Return the address of the first padhit for hit
    TClonesArray *theClusters = clusters;
    Int_t nclust = theClusters->GetEntriesFast();
    if (nclust && hit->PHlast() > 0) {
	fMaxIterPad=hit->PHlast();
	fCurIterPad=hit->PHfirst();
	return (AliMUONPadHit*) clusters->UncheckedAt(fCurIterPad-1);
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
    fCurIterPad++;
    if (fCurIterPad <= fMaxIterPad) {
	return (AliMUONPadHit*) clusters->UncheckedAt(fCurIterPad-1);
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
/* PH Commented out waiting for correct implementation
//________________________________________________________________________
void AliMUON::RemapTrackHitIDs(Int_t* map)
{
// Remaps the track numbers in the hits arrays, so that they correspond
// to the entry indices in the Kine tree.
// The correspondance is not direct. To get the real index into the Kine tree
// compute the particle index as follows:
//
//   num_primaries = AliStack::GetNprimary();
//   num_tracks = AliStack::GetNtracks();
//   track = AliMUONHit::Track()
//
//   if (track < num_primaries)
//       particleindex = track + num_tracks - num_primaries;
//   else
//       particleindex = track - num_primaries;
	
	// Remap the track numbers based on the specified map.
	AliMUONData* data = GetMUONData();
	TClonesArray* hits = data->Hits();
	for (Int_t i = 0; i < hits->GetEntriesFast(); i++)
	{
		AliMUONHit* hit = static_cast<AliMUONHit*>( hits->At(i) );
		hit->SetTrack( map[hit->Track()] );
	};
};
*/
