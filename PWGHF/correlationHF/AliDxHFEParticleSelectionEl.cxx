// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE Project            * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  Sedat Altinpinar <Sedat.Altinpinar@cern.ch>           *
//*                  Hege Erdal       <hege.erdal@gmail.com>               *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliDxHFEParticleSelectionEl.cxx
/// @author Sedat Altinpinar, Hege Erdal, Matthias Richter
/// @date   2012-03-19
/// @brief  D0 selection for D0-HFE correlation
///
#include "AliDxHFEParticleSelectionEl.h"
#include "AliVParticle.h"
#include "AliVEvent.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliHFEcontainer.h"
#include "AliHFEpid.h"
#include "AliHFEpidBase.h"
#include "AliHFEtools.h"
#include "AliHFEcuts.h"
#include "AliAODTrack.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"
#include "AliCFManager.h"
#include "THnSparse.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TObjArray.h"
#include <iostream>
#include <cerrno>
#include <memory>

using namespace std;

/// ROOT macro for the implementation of ROOT specific class methods
ClassImp(AliDxHFEParticleSelectionEl)

AliDxHFEParticleSelectionEl::AliDxHFEParticleSelectionEl(const char* opt)
  : AliDxHFEParticleSelection("electron", opt)
  , fPID(NULL)
  , fElectronProperties(NULL)
  , fWhichCut(NULL)  
  , fCuts(NULL)
  , fCFM(NULL)
{
  // constructor
  // 
  // 
  // 
  // 
}

AliDxHFEParticleSelectionEl::~AliDxHFEParticleSelectionEl()
{
  // destructor
  if (fElectronProperties) {
    delete fElectronProperties;
    fElectronProperties=NULL;
  }
  if(fCFM){
    delete fCFM;
    fCFM=NULL;
  }

  // NOTE: external objects fPID and fCuts are not deleted here
  fPID=NULL;
  fCuts=NULL;
}
 
int AliDxHFEParticleSelectionEl::InitControlObjects()
{
  /// init control and monitoring objects
  AliInfo("Electron THnSparse");
  const int thnSize = 3;
  const double Pi=TMath::Pi();
  TString name;// ="e information";
  //     		       0    1       2
  // 	 	               Pt   Phi    Eta
  int    thnBins[thnSize] = { 1000,  200, 500};
  double thnMin [thnSize] = {    0,    0, -1.};
  double thnMax [thnSize] = {  100, 2*Pi,  1.};

  name.Form("%s info", GetName());
  std::auto_ptr<THnSparseF> ElectronProperties(new THnSparseF(name, name, thnSize, thnBins, thnMin, thnMax));

  if (ElectronProperties.get()==NULL) {
    return -ENOMEM;
  }
  int axis=0;
  ElectronProperties->GetAxis(axis++)->SetTitle("Pt");
  ElectronProperties->GetAxis(axis++)->SetTitle("Phi");
  ElectronProperties->GetAxis(axis++)->SetTitle("Eta"); 

  fElectronProperties=ElectronProperties.release();

  AddControlObject(fElectronProperties);

  fWhichCut= new TH1F("fWhichCut","effective cut for a rejected particle",6,-0.5,5.5);
  AddControlObject(fWhichCut);

  //--------Initialize correction Framework and Cuts
  // Consider moving this, either to separate function or
  // add a set function for AliCFManager
  // Do we need this? Can we just call AliHFEcuts::CheckParticleCuts
  AliInfo("Setting up CFM");
  fCFM = new AliCFManager;
  // the setup of cut objects is done in AliHFEcuts::Initialize
  // the ids used by this class must be the same, the code might be broken if
  // the sequence in AliHFEcuts::Initialize is changed
  const Int_t kNcutSteps = AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kNcutStepsRecTrack + AliHFEcuts::kNcutStepsDETrack;
  // reset pointers in the CF manager
  fCFM->SetNStepParticle(kNcutSteps);
  for(Int_t istep = 0; istep < kNcutSteps; istep++) {
    fCFM->SetParticleCutsList(istep, NULL);
  }
  if(!fCuts) {
    AliWarning("Cuts not available. Default cuts will be used");
    fCuts = new AliHFEcuts;
    fCuts->CreateStandardCuts();
  }
  // TODO: error handling?
  fCuts->Initialize(fCFM);

  return 0;
}

int AliDxHFEParticleSelectionEl::HistogramParticleProperties(AliVParticle* p, int selectionCode)
{
  /// histogram particle properties
  if (!p) return -EINVAL;
  //if (!fControlObjects) return 0;
  if(selectionCode==0) return  0;

  AliAODTrack *track=(AliAODTrack*)p;
  Double_t eProperties[]={track->Pt(),track->Phi(),track->Eta()};
  if(fElectronProperties) fElectronProperties->Fill(eProperties);
  
  return 0;
}

int AliDxHFEParticleSelectionEl::IsSelected(AliVParticle* pEl, const AliVEvent*)
{
  /// select El candidates
  AliAODTrack *track=(AliAODTrack*)pEl;
 

  //--------track cut selection-----------------------
  //Using AliHFECuts:
  // RecKine: ITSTPC cuts  
  if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track)){
    fWhichCut->Fill(0);
    return 0;
  }
  
  // RecPrim
  if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) {
    fWhichCut->Fill(1);
    return 0;
  }
  
  // HFEcuts: ITS layers cuts
  if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track)) {
    fWhichCut->Fill(2);
    return 0;
  }
  
  // HFE cuts: TOF PID and mismatch flag
  if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTOF, track)) {
    fWhichCut->Fill(3);
    return 0;
  }
  
  // HFE cuts: TPC PID cleanup
  if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, track)){
    fWhichCut->Fill(4);
    return 0;
  } 
 
  // HFEcuts: Nb of tracklets TRD0
  //if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTRD, track)) continue;


  //--------PID selection-----------------------
  AliHFEpidObject hfetrack;
  hfetrack.SetAnalysisType(AliHFEpidObject::kAODanalysis);
  hfetrack.SetRecTrack(track);

  // TODO: configurable colliding system
  //if(IsPbPb()) hfetrack.SetPbPb();
  hfetrack.SetPP();

  if(fPID && fPID->IsSelected(&hfetrack)) {
    AliDebug(3,"Inside FilldPhi, electron is selected");

    return 1;
  }
  else{
    fWhichCut->Fill(5);
    return 0;
  }
}

void AliDxHFEParticleSelectionEl::SetCuts(TObject* cuts, int level)
{
  /// set cut objects
  if (level==kCutHFE) {
    fCuts=dynamic_cast<AliHFEcuts*>(cuts);
    if (!fCuts && cuts) {
      AliError(Form("Cut object is not of required type AliHFEcuts but %s", cuts->ClassName()));
    }
    return;
  }

  if (level==kCutPID) {
    fPID=dynamic_cast<AliHFEpid*>(cuts);
    if (!fPID && cuts) {
      AliError(Form("cuts object is not of required type AliHFEpid but %s", cuts->ClassName()));
    }
    return;
  }
}

//________________________________________________________________________
Bool_t AliDxHFEParticleSelectionEl::ProcessCutStep(Int_t cutStep, AliVParticle *track)
{
  // Check single track cuts for a given cut step
  const Int_t kMCOffset = AliHFEcuts::kNcutStepsMCTrack;
  if(!fCFM->CheckParticleCuts(cutStep + kMCOffset, track)) return kFALSE;
  //if(!fCuts->CheckParticleCuts(cutStep + kMCOffset, track)) return kFALSE;
  return kTRUE;
}
