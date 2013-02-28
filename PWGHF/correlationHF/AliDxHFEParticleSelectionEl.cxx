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
#include "TH2F.h"
#include "TAxis.h"
#include "TList.h"
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
  , fPIDTOF(NULL)
  , fElectronProperties(NULL)
  , fHistoList(NULL)
  , fCutPidList(NULL)
  , fPIDResponse(NULL)
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
  if(fHistoList){
    delete fHistoList;
    fHistoList=NULL;
  }
  if(fCFM){
    delete fCFM;
    fCFM=NULL;
  }
  if(fPIDResponse){
    delete fPIDResponse;
    fPIDResponse=NULL;
  }
  // NOTE: external objects fPID and fCuts are not deleted here
  fPID=NULL;
  fPIDTOF=NULL;
  fCuts=NULL;
}

const char* AliDxHFEParticleSelectionEl::fgkCutBinNames[]={
  "kRecKineITSTPC",
  "kRecPrim",
  "kHFEcuts",
  "kHFEcutsTOFPID",
  "kHFEcutsTPCPID",
  "kPIDTOF",
  "kPIDTPCTOF",
  "Selected e"

};

int AliDxHFEParticleSelectionEl::Init()
{
  //
  // init function
  // 
  int iResult=0;

  // Implicit call to InitControlObjects() before setting up CFM and fCuts
  // (if not there)
  iResult=AliDxHFEParticleSelection::Init();
  if (iResult<0) return iResult;

  //--------Initialize correction Framework and Cuts-------------------------
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

THnSparse* AliDxHFEParticleSelectionEl::DefineTHnSparse()
{
  //
  // Defines the THnSparse. For now, only calls CreatControlTHnSparse

  // here is the only place to change the dimension
  const int thnSize = 3;
  InitTHnSparseArray(thnSize);
  const double Pi=TMath::Pi();
  TString name;
  name.Form("%s info", GetName());

  //                                 0    1       2
  // 	 	                     Pt   Phi    Eta
  int         thnBins [thnSize] = { 1000,  200, 500};
  double      thnMin  [thnSize] = {    0,    0, -1.};
  double      thnMax  [thnSize] = {  100, 2*Pi,  1.};
  const char* thnNames[thnSize] = { "Pt","Phi","Eta"};

  return CreateControlTHnSparse(name,thnSize,thnBins,thnMin,thnMax,thnNames);
}

int AliDxHFEParticleSelectionEl::InitControlObjects()
{
  /// init control and monitoring objects
  AliInfo("Electron THnSparse");

  fElectronProperties=DefineTHnSparse();
  AddControlObject(fElectronProperties);

  double dEdxBins[6]={1000,0.,10.,200,0.,200.};
  double nSigBins[6]={1000,0.,10.,200,-10.,10.};

  fHistoList=new TList;
  fHistoList->SetName("HFE Histograms");
  fHistoList->SetOwner();

  // Histogram storing which cuts have been applied to the tracks
  fHistoList->Add(CreateControlHistogram("fWhichCut","effective cut for a rejected particle", kNCutLabels, fgkCutBinNames));

  // dEdx plots, TPC signal vs momentum
  fHistoList->Add(CreateControl2DHistogram("fdEdx", "dEdx before cuts", dEdxBins,"momentum (GeV/c)","dE/dx in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fdEdxCut", "dEdx after cuts",dEdxBins,"momentum (GeV/c)","dE/dx in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fdEdxPidTOF", "dEdx after TOF pid",dEdxBins,"momentum (GeV/c)","dE/dx in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fdEdxPid", "dEdx after pid",dEdxBins,"momentum (GeV/c)","dE/dx in TPC (a.u.)"));

  // nSigmaTPC vs momentum
  fHistoList->Add(CreateControl2DHistogram("fnSigTPC", "nSigTPC before cuts",nSigBins,"momentum (GeV/c)","nSigma in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fnSigTPCCut", "nSigmaTPC after cuts",nSigBins,"momentum (GeV/c)","nSigma in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fnSigTPCPidTOF", "nSigmaTPC after TOF PID", nSigBins,"momentum (GeV/c)","nSigma in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fnSigTPCPid", "nSigmaTPC after PID", nSigBins,"momentum (GeV/c)","nSigma in TPC (a.u.)"));

  // nSigmaTOF vs momentum
  fHistoList->Add(CreateControl2DHistogram("fnSigTOF", "nSigmaTOF before cuts",nSigBins,"momentum (GeV/c)","nSigma in TOF (a.u.)"));  
  fHistoList->Add(CreateControl2DHistogram("fnSigTOFCut", "nSigmaTOF after cuts",nSigBins,"momentum (GeV/c)","nSigma in TOF (a.u.)"));  
  fHistoList->Add(CreateControl2DHistogram("fnSigTOFPidTOF", "nSigmaTOF after TOF PID", nSigBins,"momentum (GeV/c)","nSigma in TOF (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fnSigTOFPid", "nSigmaTOF after PID", nSigBins,"momentum (GeV/c)","nSigma in TOF (a.u.)"));
 
  AddControlObject(fHistoList);

  return AliDxHFEParticleSelection::InitControlObjects();
}

int AliDxHFEParticleSelectionEl::HistogramParticleProperties(AliVParticle* p, int selectionCode)
{
  /// histogram particle properties
  if (!p) return -EINVAL;
  //if (!fControlObjects) return 0;
  if(selectionCode==0) return  0;
  
  if(fElectronProperties && ParticleProperties()) {
    FillParticleProperties(p, ParticleProperties(), GetDimTHnSparse());
    fElectronProperties->Fill(ParticleProperties());
  }

  return 0;
}

int AliDxHFEParticleSelectionEl::FillParticleProperties(AliVParticle* p, Double_t* data, int dimension) const
{
  // fill the data array from the particle data
  if (!data) return -EINVAL;
  AliAODTrack *track=(AliAODTrack*)p;
  if (!track) return -ENODATA;
  int i=0;
  if (dimension!=GetDimTHnSparse()) {
    // TODO: think about filling only the available data and throwing a warning
    return -ENOSPC;
  }
  data[i++]=track->Pt();
  data[i++]=track->Phi();
  data[i++]=track->Eta();
  return i;
}


int AliDxHFEParticleSelectionEl::IsSelected(AliVParticle* pEl, const AliVEvent* pEvent)
{
  /// select El candidates
  // TODO: How to handle MC? would be too much duplicate code if copy entire IsSelected. 

  AliAODTrack *track=(AliAODTrack*)pEl;
  fCFM->SetRecEventInfo(pEvent);
  
  ((TH2D*)fHistoList->FindObject("fdEdx"))->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
  ((TH2D*)fHistoList->FindObject("fnSigTPC"))->Fill(track->GetTPCmomentum(), fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
  ((TH2D*)fHistoList->FindObject("fnSigTOF"))->Fill(track->P(), fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron));
  
  //--------track cut selection-----------------------
  //Using AliHFECuts:
  // RecKine: ITSTPC cuts  
  if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track)){
    ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kRecKineITSTPC);
    return 0;
  }
  
  // RecPrim
  if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) {
    ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kRecPrim);
    return 0;
  }
  
  // HFEcuts: ITS layers cuts
  if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track)) {
    ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kHFEcutsITS);
    return 0;
  }
  
  // HFE cuts: TOF PID and mismatch flag
  if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTOF, track)) {
    ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kHFEcutsTOF);
    return 0;
  }
  
  // HFE cuts: TPC PID cleanup
  if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, track)){
    ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kHFEcutsTPC);
    return 0;
  } 
 
  // HFEcuts: Nb of tracklets TRD0
  //if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTRD, track)) continue;
  ((TH2D*)fHistoList->FindObject("fdEdxCut"))->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
  ((TH2D*)fHistoList->FindObject("fnSigTPCCut"))->Fill(track->GetTPCmomentum(), fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
  ((TH2D*)fHistoList->FindObject("fnSigTOFCut"))->Fill(track->P(), fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron));


  //--------PID selection-----------------------
  AliHFEpidObject hfetrack;
  hfetrack.SetAnalysisType(AliHFEpidObject::kAODanalysis);
  hfetrack.SetRecTrack(track);

  // TODO: configurable colliding system
  //if(IsPbPb()) hfetrack.SetPbPb();
  hfetrack.SetPP();

  // TODO: Put this into a while-loop instead, looping over the number of pid objects in the cut-list?
  // This needs a bit of thinking and finetuning (wrt histogramming)
  if(fPIDTOF && fPIDTOF->IsSelected(&hfetrack)) {
    ((TH2D*)fHistoList->FindObject("fdEdxPidTOF"))->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
    ((TH2D*)fHistoList->FindObject("fnSigTPCPidTOF"))->Fill(track->GetTPCmomentum(), fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
    ((TH2D*)fHistoList->FindObject("fnSigTOFPidTOF"))->Fill(track->P(), fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron));
  }
  else{
    ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kPIDTOF);
    return 0;
  }

  //Combined tof & tpc pid
  if(fPID && fPID->IsSelected(&hfetrack)) {
    AliDebug(3,"Inside FilldPhi, electron is selected");
    ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kSelected);
    ((TH2D*)fHistoList->FindObject("fdEdxPid"))->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
    ((TH2D*)fHistoList->FindObject("fnSigTPCPid"))->Fill(track->GetTPCmomentum(), fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
    ((TH2D*)fHistoList->FindObject("fnSigTOFPid"))->Fill(track->P(), fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron));
    return 1;
  }
  else{
    ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kPID);
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

  if (level==kCutPIDTOF) {
    fPIDTOF=dynamic_cast<AliHFEpid*>(cuts);
    if (!fPIDTOF && cuts) {
      AliError(Form("cuts object is not of required type AliHFEpid but %s", cuts->ClassName()));
    }
    return;
  }
  if(level==kCutList){
    fCutPidList=dynamic_cast<TList*>(cuts);
    if (!fCutPidList && cuts) {
      AliError(Form("cuts object is not of required type TList but %s", cuts->ClassName()));
    }
    else{
      // TODO: Could be done more elegantly, at the moment requires that the cut and pid objects are in 
      // a specific order..
      TObject *obj=NULL;
      int iii=0;
      TIter next(fCutPidList);
      while((obj = next())){
	iii++;
	if(iii==1) {
	  fCuts=dynamic_cast<AliHFEcuts*>(obj);
	  if (!fCuts) 
	    AliError(Form("Cut object is not of required type AliHFEcuts but %s", cuts->ClassName()));
	}
	if(iii==2){
	  fPID=dynamic_cast<AliHFEpid*>(obj);
	  if (!fPID) 
	    AliError(Form("cuts object is not of required type AliHFEpid but %s", cuts->ClassName()));
	}
	if(iii==3){ 
	  fPIDTOF=dynamic_cast<AliHFEpid*>(obj);
	  if (!fPIDTOF) 
	    AliError(Form("cuts object is not of required type AliHFEpid but %s", cuts->ClassName()));
	}
      }
      
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
