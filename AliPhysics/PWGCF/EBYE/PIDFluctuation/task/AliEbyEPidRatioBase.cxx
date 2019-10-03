/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: ALICE Offline.                                                 *
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

//=========================================================================//
//             AliEbyE Analysis for Particle Ratio Fluctuation             //
//                   Deepika Rathee  | Satyajit Jena                       //
//                   drathee@cern.ch | sjena@cern.ch                       //
//                  Date: Wed Jul  9 18:38:30 CEST 2014                    //
//          New approch to find particle ratio to reduce memory            //
//                             (Test Only)                                 //
//        Copied from NetParticle Classes
//        Origin: Authors: Jochen Thaeder <jochen@thaeder.de>
//                         Michael Weber <m.weber@cern.ch>
//=========================================================================//

#include "TMath.h"
#include "TAxis.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliESDtrackCuts.h"


#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliEbyEPidRatioBase.h"

using namespace std;
ClassImp(AliEbyEPidRatioBase)
//________________________________________________________________________
AliEbyEPidRatioBase::AliEbyEPidRatioBase() :
  TNamed(),
  fHelper(NULL),
  
  fESD(NULL), 
  fESDTrackCuts(NULL),
  fAOD(NULL),
  fArrayMC(NULL),
  fAODtrackCutBit(1024),
  fIsMC(kFALSE),
  fMCEvent(NULL),
  fStack(NULL),
  
  fCentralityBin(-1.),
  fNTracks(0), fIsRatio(kFALSE), fIsPtBin(kFALSE), fIsDetectorWise(kFALSE) {
  // Constructor   

  AliLog::SetClassDebugLevel("AliEbyEPidRatioBase",10);
}

//________________________________________________________________________
AliEbyEPidRatioBase::AliEbyEPidRatioBase(const Char_t* name, const Char_t* title) :
  TNamed(name, title),
  fHelper(NULL),
 
  fESD(NULL), 
  fESDTrackCuts(NULL),
  fAOD(NULL),
  fArrayMC(NULL),
  fAODtrackCutBit(1024),
  fIsMC(kFALSE),
  
  
  fMCEvent(NULL),
  fStack(NULL),
  
  fCentralityBin(-1.),
  fNTracks(0), fIsRatio(kFALSE),fIsPtBin(kFALSE), fIsDetectorWise(kFALSE){
  // Constructor   

  AliLog::SetClassDebugLevel("AliEbyEPidRatioBase",10);
}

//________________________________________________________________________
AliEbyEPidRatioBase::~AliEbyEPidRatioBase() {
  // Destructor
}

//________________________________________________________________________
void AliEbyEPidRatioBase::Initialize(AliEbyEPidRatioHelper* helper, AliESDtrackCuts* cuts) {
  fHelper           = helper;
  fESDTrackCuts     = (cuts) ? cuts : helper->GetESDTrackCuts();
  fIsMC             = helper->GetIsMC();
  fIsRatio          = helper->GetIsRatio();
  fIsPtBin          = helper->GetIsPtBin();
  fIsDetectorWise   = helper->GetDetWise();
  fAODtrackCutBit   = helper->GetAODtrackCutBit();
  Init();
  CreateHistograms();
 
  Float_t ptRange[2];
  fESDTrackCuts->GetPtRange(ptRange[0],ptRange[1]);
  Printf(">>>> Pt Initialisation:  [%f,%f]",ptRange[0],ptRange[1]);
  fESDTrackCuts->GetEtaRange(ptRange[0],ptRange[1]);
  Printf(">>>> Eta Initialisation: [%f,%f]",ptRange[0],ptRange[1]);
 
  return;
}

//________________________________________________________________________
Int_t AliEbyEPidRatioBase::SetupEvent() {
  // -- Setup event


  ResetEvent();

  // -- Get ESD objects
  if (dynamic_cast<AliESDInputHandler*>(fHelper->GetInputEventHandler())) {
    fESD         = dynamic_cast<AliESDEvent*>(fHelper->GetInputEventHandler()->GetEvent());
    fNTracks     = fESD->GetNumberOfTracks();
  }

  // -- Get AOD objects
  else if (dynamic_cast<AliAODInputHandler*>(fHelper->GetInputEventHandler())) {
    fAOD         = dynamic_cast<AliAODEvent*>(fHelper->GetInputEventHandler()->GetEvent());
    fNTracks     = fAOD->GetNumberOfTracks();
    
    if (fIsMC) {
      fArrayMC = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!fArrayMC)
	AliFatal("No array of MC particles found !!!"); // MW  no AliFatal use return values
    }
  }

  if (fIsMC) {
    fMCEvent     = fHelper->GetMCEvent();
    if (fMCEvent)
      fStack     = fMCEvent->Stack();
  }

  fCentralityBin = fHelper->GetCentralityBin();

  return Setup();
}

//________________________________________________________________________
void AliEbyEPidRatioBase::ResetEvent() {
  // -- Reset ESD Event
  fESD       = NULL;
  // -- Reset AOD Event
  fAOD       = NULL;
  // -- Reset MC Event
  if (fIsMC)
    fMCEvent = NULL;
  // -- Reset in class
  Reset();
  return;
}


