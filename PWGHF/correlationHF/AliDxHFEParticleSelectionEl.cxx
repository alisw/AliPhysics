//  $Id$

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
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include <iostream>
#include <cerrno>
#include <memory>

using namespace std;

/// ROOT macro for the implementation of ROOT specific class methods
ClassImp(AliDxHFEParticleSelectionEl)

AliDxHFEParticleSelectionEl::AliDxHFEParticleSelectionEl(const char* opt)
  : AliDxHFEParticleSelection("electron", opt)
  , fPIDTOFTPC(NULL)
  , fPIDTOF(NULL)
  , fElectronProperties(NULL)
  , fHistoList(NULL)
  , fCutPidList(NULL)
  , fPIDResponse(NULL)
  , fCuts(NULL)
  , fCFM(NULL)
  , fFinalCutStep(kPIDTOFTPC)
  , fInvMassLow(0.1)
  , fUseInvMassCut(kFALSE)
{
  // constructor
  // 
  // 
  // 
  // 
  ParseArguments(opt);
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
  fPIDTOFTPC=NULL;
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
  // TODO: Add also invariant mass? here and in correlation, to do cut afterwards..

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

  fHistoList=new TList;
  fHistoList->SetName("HFE Histograms");
  fHistoList->SetOwner();

  // Histogram storing which cuts have been applied to the tracks
  fHistoList->Add(CreateControlHistogram("fWhichCut","effective cut for a rejected particle", kNCutLabels, fgkCutBinNames));
  fHistoList->Add(CreateControlHistogram("fTPCnClAOD","Nr TPC clusters/track for all tracks", 160, 0., 159.));
  fHistoList->Add(CreateControlHistogram("fTPCnClSingleTrackCuts","Nr TPC clusters/track After SingleTrackCuts", 160, 0., 159.));
  fHistoList->Add(CreateControlHistogram("fTPCnClTPCTOFPID","Nr TPC clusters/track After TPCTOF PID", 160, 0., 159.));

  double dEdxBins[6]={1000,0.,10.,200,0.,200.};
  double nSigBins[6]={1000,0.,10.,200,-10.,10.};

  // dEdx plots, TPC signal vs momentum
  fHistoList->Add(CreateControl2DHistogram("fdEdx", "dEdx before cuts", dEdxBins,"momentum (GeV/c)","dE/dx in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fdEdxCut", "dEdx after cuts",dEdxBins,"momentum (GeV/c)","dE/dx in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fdEdxPidTOF", "dEdx after TOF pid",dEdxBins,"momentum (GeV/c)","dE/dx in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fdEdxPid", "dEdx after pid",dEdxBins,"momentum (GeV/c)","dE/dx in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fdEdxIM", "dEdx after Inv Mass",dEdxBins,"momentum (GeV/c)","dE/dx in TPC (a.u.)"));

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

  // Invariant mass LS and ULS without cut
  fHistoList->Add(CreateControlHistogram("fInvMassLS", "Invariant mass LS", 1000, 0., 0.5));
  fHistoList->Add(CreateControlHistogram("fInvMassULS", "Invariant mass ULS", 1000, 0., 0.5));

  // Invariant mass LS and ULS without cut, extended x-axis
  fHistoList->Add(CreateControlHistogram("fInvMassLSExt", "Invariant mass LS", 1000, 0., 5.));
  fHistoList->Add(CreateControlHistogram("fInvMassULSExt", "Invariant mass ULS", 1000, 0., 5.));

  // Invariant mass LS and ULS with cut
  fHistoList->Add(CreateControlHistogram("fInvMassLScut", "Invariant mass LS (cut)", 1000, 0., 0.5));
  fHistoList->Add(CreateControlHistogram("fInvMassULScut", "Invariant mass ULS (cut)", 1000, 0., 0.5));

 
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


TObjArray* AliDxHFEParticleSelectionEl::Select(const AliVEvent* pEvent)
{
  /// create selection from 'Tracks' member of the event,
  /// array contains only pointers but does not own the objects
  /// object array needs to be deleted by caller
  if (!pEvent) return NULL;

  TList* selectedTracks=new TList;

  TObjArray* finalTracks=new TObjArray;
  if (!finalTracks) return NULL;
  finalTracks->SetOwner(kFALSE); // creating new track objects below

  int nofTracks=pEvent->GetNumberOfTracks();
  int selectionCode=0;
  for (int itrack=0; itrack<nofTracks; itrack++) {
    selectionCode=0;
    AliVParticle* track=pEvent->GetTrack(itrack);
    selectionCode=IsSelected(track,pEvent);
    // If fUseInvMassCut is set to true, only call HistogramParticleProperties
    // (here) for those who didn't pass single track cuts and pid. 
    // if fUseInvMassCut is not set, then call the function for all tracks
    // TODO: Change this when have added invariant mass i thnsparse
    if(!fUseInvMassCut || selectionCode==0)
      HistogramParticleProperties(track, selectionCode);
    if (selectionCode==0) continue;
    selectedTracks->Add(track);
    //if fUseInvMassCut, fill the array below
    if(!fUseInvMassCut) finalTracks->Add(CreateParticle(track));
  }

  // Calculating invariant mass for electron pairs with same selection criteria
  // TODO: Could be removed if it is verified that looser cuts is better, but keep
  // for now
  if(selectedTracks->GetSize()>0)
    {
      Bool_t* selTrackIndex=new Bool_t[selectedTracks->GetSize()];
      InvMassFilter(selectedTracks, selTrackIndex);

      //Only remove electrons if fUseInvMassCut is set
      if(fUseInvMassCut){
	for(Int_t k=0; k<selectedTracks->GetSize(); k++)
	  {
	    selectionCode=0; //Reset selectionCode here to be based on invariant mass cut 
	    AliAODTrack *trackIM=(AliAODTrack*)selectedTracks->At(k);
	    if(selTrackIndex[k]==kTRUE)
	      {
		((TH2D*)fHistoList->FindObject("fdEdxIM"))->Fill(trackIM->GetTPCmomentum(), trackIM->GetTPCsignal());
		finalTracks->Add(CreateParticle(trackIM));
		selectionCode=1;
	      }
	    HistogramParticleProperties(trackIM, selectionCode);
	  }
      }
      
    }
  return finalTracks;
}

int AliDxHFEParticleSelectionEl::IsSelected(AliVParticle* pEl, const AliVEvent* pEvent)
{
  /// select El candidates
  // TODO: How to handle MC? would be too much duplicate code if copy entire IsSelected. 
  if(!pEvent){
    AliError("No event information");
    return 0;
  }

  AliAODTrack *track=(AliAODTrack*)pEl;
  fCFM->SetRecEventInfo(pEvent);
  
  ((TH2D*)fHistoList->FindObject("fdEdx"))->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
  ((TH2D*)fHistoList->FindObject("fnSigTPC"))->Fill(track->GetTPCmomentum(), fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
  ((TH2D*)fHistoList->FindObject("fnSigTOF"))->Fill(track->P(), fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron));
  ((TH1D*)fHistoList->FindObject("fTPCnClAOD"))->Fill(track->GetTPCNcls());

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

  ((TH2D*)fHistoList->FindObject("fdEdxCut"))->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
  ((TH2D*)fHistoList->FindObject("fnSigTPCCut"))->Fill(track->GetTPCmomentum(), fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
  ((TH2D*)fHistoList->FindObject("fnSigTOFCut"))->Fill(track->P(), fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron));
  ((TH1D*)fHistoList->FindObject("fTPCnClSingleTrackCuts"))->Fill(track->GetTPCNcls());

  if(fFinalCutStep==kHFEcutsTPC) {AliDebug(2,"Returns after track cuts"); return 1;}


  //--------PID selection-----------------------
  AliHFEpidObject hfetrack;
  hfetrack.SetAnalysisType(AliHFEpidObject::kAODanalysis);
  hfetrack.SetRecTrack(track);

  // TODO: configurable colliding system
  if(GetSystem()==1) hfetrack.SetPbPb();
  else  hfetrack.SetPP();

  // TODO: Put this into a while-loop instead, looping over the number of pid objects in the cut-list?
  // This needs a bit of thinking and finetuning (wrt histogramming)
  if(fPIDTOF && fPIDTOF->IsSelected(&hfetrack)) {
    ((TH2D*)fHistoList->FindObject("fdEdxPidTOF"))->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
    ((TH2D*)fHistoList->FindObject("fnSigTPCPidTOF"))->Fill(track->GetTPCmomentum(), fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
    ((TH2D*)fHistoList->FindObject("fnSigTOFPidTOF"))->Fill(track->P(), fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron));

    if(fFinalCutStep==kPIDTOF) {AliDebug(2,"Returns at PIDTOF"); return 1;}
  }
  else{
    ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kPIDTOF);
  }

  //Combined tof & tpc pid
  if(fPIDTOFTPC && fPIDTOFTPC->IsSelected(&hfetrack)) {
    AliDebug(3,"Inside FilldPhi, electron is selected");
    ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kSelected);
    ((TH2D*)fHistoList->FindObject("fdEdxPid"))->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
    ((TH2D*)fHistoList->FindObject("fnSigTPCPid"))->Fill(track->GetTPCmomentum(), fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
    ((TH2D*)fHistoList->FindObject("fnSigTOFPid"))->Fill(track->P(), fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron));
    ((TH1D*)fHistoList->FindObject("fTPCnClTPCTOFPID"))->Fill(track->GetTPCNcls());

    return 1;
  }
  else{
    ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kPIDTOFTPC);
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

  if (level==kCutPIDTOFTPC) {
    fPIDTOFTPC=dynamic_cast<AliHFEpid*>(cuts);
    if (!fPIDTOFTPC && cuts) {
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
	  fPIDTOFTPC=dynamic_cast<AliHFEpid*>(obj);
	  if (!fPIDTOFTPC) 
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

int AliDxHFEParticleSelectionEl::ParseArguments(const char* arguments)
{
  // parse arguments and set internal flags
  TString strArguments(arguments);
  auto_ptr<TObjArray> tokens(strArguments.Tokenize(" "));
  if (!tokens.get()) return 0;

  AliInfo(strArguments);
  TIter next(tokens.get());
  TObject* token;
  while ((token=next())) {
    TString argument=token->GetName();
    if(argument.BeginsWith("elreco=")){
      argument.ReplaceAll("elreco=", "");
      int selectionStep=kPIDTOFTPC; //Default
      if(argument.CompareTo("aftertrackcuts")==0) selectionStep=kHFEcutsTPC;
      if(argument.CompareTo("aftertofpid")==0) selectionStep=kPIDTOF;
      if(argument.CompareTo("afterfullpid")==0) selectionStep=kPIDTOFTPC;
      SetFinalCutStep(selectionStep); 
      continue;   
    }
    if(argument.BeginsWith("useinvmasscut")){
      fUseInvMassCut=kTRUE;
      AliInfo("Using Invariant mass cut");
      continue;   
    }
    if(argument.BeginsWith("invmasscut=")){
      argument.ReplaceAll("invmasscut=", "");
      fInvMassLow=argument.Atoi();
      fUseInvMassCut=kTRUE;
      continue;   
    }

    // forwarding of single argument works, unless key-option pairs separated
    // by blanks are introduced
    AliDxHFEParticleSelection::ParseArguments(argument);
  }
  
  return 0;
}


void AliDxHFEParticleSelectionEl::InvMassFilter(TList *elList, Bool_t *selIndx)
{

  //Function for getting invariant mass of electron pairs
   
  for(int i=0; i<((elList->GetSize())-1); i++)
    {
      AliAODTrack *trackAsso=(AliAODTrack*)elList->At(i);
      for(int j=i+1; j<elList->GetSize(); j++)
	{
	  
	  AliAODTrack *track=(AliAODTrack*)elList->At(j);
	  if(trackAsso && track)
	    {
	      
	      Double_t mass=-999., width = -999;
	      Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
	      
	      Int_t chargeAsso = trackAsso->Charge();
	      Int_t charge = track->Charge();
	      
	      Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;
	      if(charge>0) fPDGe1 = -11;
	      if(chargeAsso>0) fPDGe2 = -11;
	      
	      if(charge == chargeAsso) fFlagLS = kTRUE;
	      if(charge != chargeAsso) fFlagULS = kTRUE;
	      
	      AliKFParticle ge1 = AliKFParticle(*track, fPDGe1);
	      AliKFParticle ge2 = AliKFParticle(*trackAsso, fPDGe2);
	      AliKFParticle recg(ge1, ge2);
	      
	      if(recg.GetNDF()<1) continue;
	      Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
	      if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;
	      
	      recg.GetMass(mass,width);

	      if(fFlagLS) {
		if(mass<0.5)
		  ((TH1D*)fHistoList->FindObject("fInvMassLS"))->Fill(mass);
		((TH1D*)fHistoList->FindObject("fInvMassLSExt"))->Fill(mass);
		if(mass>fInvMassLow)
		  {
		    ((TH1D*)fHistoList->FindObject("fInvMassLScut"))->Fill(mass);
		    selIndx[i]=kTRUE;
		    selIndx[j]=kTRUE;
		  }
	      }
	      if(fFlagULS) {
		if(mass<0.5)
		  ((TH1D*)fHistoList->FindObject("fInvMassULS"))->Fill(mass);
		((TH1D*)fHistoList->FindObject("fInvMassULSExt"))->Fill(mass);
		if(mass>fInvMassLow)
		  {
		    ((TH1D*)fHistoList->FindObject("fInvMassULScut"))->Fill(mass);
		    selIndx[i]=kTRUE;
		    selIndx[j]=kTRUE;
		  }

	      }
	      
	    }
	}
    }
}
