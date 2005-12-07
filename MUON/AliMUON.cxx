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
#include "AliMUONRawCluster.h"
#include "AliMUONTransientDigit.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONGeometry.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryBuilder.h"
#include "AliMUONCommonGeometryBuilder.h"
#include "AliMUONVGeometryBuilder.h"	
#include "AliMUONGeometrySegmentation.h"
#include "AliMUONDigitizerv2.h"
#include "AliMUONSDigitizerv1.h"
#include "AliMUONRawData.h"
#include "AliMUONSegmentation.h"
#include "AliLog.h"

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
    fSegmentation(0),
    fAccCut(kFALSE),
    fAccMin(0.),
    fAccMax(0.),   
    fMaxStepGas(0.),
    fMaxStepAlu(0.),
    fMaxDestepGas(0.),
    fMaxDestepAlu(0.),
    fMaxIterPad(0),
    fCurIterPad(0)
{
// Default Constructor
//
	AliDebug(1,Form("default (empty) ctor this=%p",this));
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
    fSegmentation(0),
    fAccCut(kFALSE),
    fAccMin(0.),
    fAccMax(0.),   
    fMaxStepGas(0.1),
    fMaxStepAlu(0.1),
    fMaxDestepGas(-1), // Negatives values are ignored by geant3 CONS200 
    fMaxDestepAlu(-1), // in the calculation of the tracking parameters
    fMaxIterPad(0),
    fCurIterPad(0)
{
	AliDebug(1,Form("ctor this=%p",this));
  fIshunt =  0;

  SetMarkerColor(kRed);//
    
  // Geometry builder
  fGeometryBuilder = new AliMUONGeometryBuilder(this);
  
  // Common geometry definitions
  fGeometryBuilder
    ->AddBuilder(new AliMUONCommonGeometryBuilder(this));

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
	  fChambers->AddAt(new AliMUONChamberTrigger(ch, GetGeometryTransformer()),ch);
	}
      } // Chamber stCH (0, 1) in 
    }     // Station st (0...)
    
    // cp new design of AliMUONTriggerDecision
    fTriggerCircuits = new TObjArray(AliMUONConstants::NTriggerCircuit());
    for (Int_t circ=0; circ<AliMUONConstants::NTriggerCircuit(); circ++) {
      fTriggerCircuits->AddAt(new AliMUONTriggerCircuit(),circ);          
    }
}

//____________________________________________________________________
AliMUON::AliMUON(const AliMUON& rMUON)
 : AliDetector(rMUON)
{
// Protected copy constructor

  AliFatal("Not implemented.");
}

//____________________________________________________________________
AliMUON::~AliMUON()
{
// Destructor
  AliDebug(1,Form("dtor this=%p",this));
  fIshunt  = 0;

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
  delete fSegmentation;
}

//________________________________________________________________________
AliMUON& AliMUON::operator = (const AliMUON& rhs)
{
// Protected assignement operator

  if (this == &rhs) return *this;

  AliFatal("Not implemented.");
    
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


//     for (Int_t i = 0; i < AliMUONConstants::NCh(); i++)     
//       this->Chamber(i).SegmentationModel2(1)->Draw("eventdisplay");// to be check !
     
  
}

//____________________________________________________________________
const AliMUONGeometry*  AliMUON::GetGeometry() const
{
// Return geometry parametrisation

  if ( !fGeometryBuilder) {
    AliWarningStream() << "GeometryBuilder not defined." << std::endl;
    return 0;
  }
  
  return fGeometryBuilder->GetGeometry();
}   

//____________________________________________________________________
const AliMUONGeometryTransformer*  AliMUON::GetGeometryTransformer() const
{
// Return geometry parametrisation

  const AliMUONGeometry* kGeometry = GetGeometry();
  
  if ( !kGeometry) return 0;

  return kGeometry->GetTransformer();
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
 void  AliMUON::SetAlign(Bool_t align)
{
 // Sets option for alignement to geometry builder
 
   fGeometryBuilder->SetAlign(align);
}   

//____________________________________________________________________
 void  AliMUON::SetAlign(const TString& fileName, Bool_t align)
{
 // Sets option for alignement to geometry builder
 
   fGeometryBuilder->SetAlign(fileName, align);
}   

//____________________________________________________________________
void   AliMUON::SetResponseModel(Int_t id, AliMUONResponse *response)
{
// Set the response for chamber id
    ((AliMUONChamber*) fChambers->At(id))->SetResponseModel(response);
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
  if (!rawData->Digits2Raw()) AliInfo("pb writting raw data");
  delete rawData;
  return;
}
//_______________________________________________________________________
AliLoader* AliMUON::MakeLoader(const char* topfoldername)
{ 
//builds standard getter (AliLoader type)
//if detector wants to use castomized getter, it must overload this method

 
 AliDebug(1,Form("Creating standard getter for detector %s. Top folder is %s.",
         GetName(),topfoldername));
 fLoader   = new AliLoader(GetName(),topfoldername);
 fMUONData = new AliMUONData(fLoader,GetName(),GetName()); 
 fMUONData->SetSplitLevel(fSplitLevel);
 return fLoader;
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

