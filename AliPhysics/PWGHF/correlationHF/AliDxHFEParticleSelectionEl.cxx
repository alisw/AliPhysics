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
#include "AliSelectNonHFE.h"
#include "AliReducedParticle.h"
#include "AliESDtrackCuts.h"
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
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"
#include "AliExternalTrackParam.h"
#include "AliAODVertex.h"
#include "AliHFEextraCuts.h"
#include "AliCFManager.h"
#include "THnSparse.h"
#include "AliLog.h"
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
  , fPIDTPC(NULL)
  , fPIDTPCEMCAL(NULL)
  , fElectronProperties(NULL)
  , fHistoList(NULL)
  , fCutPidList(NULL)
  , fPIDResponse(NULL)
  , fCuts(NULL)
  , fSelNHFE(NULL)
  , fTrackCuts(NULL)
  , fCFM(NULL)
  , fFinalCutStep(kINVMASS)
  , fInvMassLow(0.1)
  , fUseInvMassCut(kNoInvMass)
  , fSystem(0)
  , fTrackNum(0)
  , fImpactParamCutRadial(3)
  , fEtaCut(0.8)
  , fSurvivedCutStep(kNotSelected)
  , fStoreCutStepInfo(kFALSE)
  , fSetFilterBit(kTRUE)
  , fBit(0)
  , fMaxPtCombinedPID(999)
  , fMaxPTOFWhenPresent(-1)
  , fUseEMCAL(kFALSE)
  , fMinM20(0.03)
  , fMaxM20(0.3)
  , fMinM02(0.03)
  , fMaxM02(0.5)
  , fDispersion(1)
  , fEovPMin(0.8)
  , fEovPMax(1.2)
  , fUseTOFonlyWhenPresent(false)
  , fStopAfterFilterBit(false)
  , fCutLS(kFALSE)
  , fTPCnSigmaPRej(3)
  , fTPCnSigmaPiRej(3.5)
  , fPRejPMin(1.3)
  , fPRejPMax(2.6)
  , fPiRejPMin(1.0)
  , fPiRejPMax(999)
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
  fPIDTPC=NULL;
  fPIDTPCEMCAL=NULL;
  fCuts=NULL;
  if(fSelNHFE){
    delete fSelNHFE;
    fSelNHFE=NULL;
  }
  if(fTrackCuts){
    delete fTrackCuts;
    fTrackCuts=NULL;
  }
}

const char* AliDxHFEParticleSelectionEl::fgkCutBinNames[]={
  "kRecKineITSTPC",
  "kRecPrim",
  "kHFEcuts",
  "kHFEcutsTOFPID",
  "kHFEcutsTPCPID",
  "kPIDTOF",
  "kPIDTPC",
  "kPIDTPCTOF",
  "kRejPi",
  "kRejProton",
  "kINVMASS",
  "Selected e"
};

int AliDxHFEParticleSelectionEl::Init()
{
  //
  // init function
  // 
  int iResult=0;
  // TODO: think about initializing fSelNHFE and fTrackCuts from the macro
  //Initialization of invariant mass cut function (AliSelectNonHFE)
  fSelNHFE= new AliSelectNonHFE("IM","IM");
  // Cuts for associated track in Inv Mass cut
  fTrackCuts = new AliESDtrackCuts();
  fTrackCuts->SetAcceptKinkDaughters(kFALSE);
  fTrackCuts->SetRequireTPCRefit(kTRUE);
  fTrackCuts->SetEtaRange(-0.9,0.9);
  fTrackCuts->SetRequireSigmaToVertex(kTRUE);
  fTrackCuts->SetMaxChi2PerClusterTPC(4.0);
  fTrackCuts->SetMinNClustersTPC(80);
  fTrackCuts->SetPtRange(0.3,1e10);
  
  fSelNHFE->SetTrackCuts(-3, 3, fTrackCuts);
  fSelNHFE->SetInvariantMassCut(fInvMassLow);
  fSelNHFE->SetAlgorithm("KF");
  fSelNHFE->SetAODanalysis(kTRUE);
  
  
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
  fCuts->SetPIDResponse(fPIDResponse);
  fCuts->Initialize(fCFM);

  if(fUseInvMassCut>kNoInvMass)  AliDebug(2,Form("Setting up with invariant mass cut of %f",fInvMassLow));
  if(fStoreCutStepInfo) AliDebug(2,"will store cut step info");
  AliDebug(2,Form("Stopping selection after step: %d", fFinalCutStep));
  
  return 0;
  
}

THnSparse* AliDxHFEParticleSelectionEl::DefineTHnSparse()
{
  
  //
  // Defines the THnSparse. For now, only calls CreatControlTHnSparse
  // TODO: Add also invariant mass? here and in correlation, to do cut afterwards..
  
  // here is the only place to change the dimension
  const double Pi=TMath::Pi();
  TString name;
  THnSparse* thn=NULL;
  name.Form("%s info", GetName());
  
  if(fStoreCutStepInfo){
    const int thnSizeExt =4;
    
    InitTHnSparseArray(thnSizeExt);
    
    // TODO: Redo binning of distributions more?
    //     		               0    1      2   
    // 	 	                       Pt   Phi   Eta 
    int    thnBinsExt[thnSizeExt] = { 50,  100, 100,   kNCutLabels-1};
    double thnMinExt [thnSizeExt] = {   0,    0, -1.,   kRecKineITSTPC-0.5 };
    double thnMaxExt [thnSizeExt] = {  10, 2*Pi,  1.,   kSelected-0.5};
    const char* thnNamesExt[thnSizeExt]={
      "Pt",
      "Phi",
      "Eta", 
      "Last survived cut step"
    };
    thn=(THnSparse*)CreateControlTHnSparse(name,thnSizeExt,thnBinsExt,thnMinExt,thnMaxExt,thnNamesExt);
  }
  else{
    
    const int thnSize =3;
    InitTHnSparseArray(thnSize);
    
    // TODO: Redo binning of distributions more?
    //     		       0    1      2    
    // 	 	               Pt   Phi   Eta   
    int    thnBins[thnSize] = { 50,  100, 100 };
    double thnMin [thnSize] = {   0,    0, -1. };
    double thnMax [thnSize] = {  10, 2*Pi,  1. };
    const char* thnNames[thnSize]={
      "Pt",
      "Phi",
      "Eta", 
    };
    thn=(THnSparse*)CreateControlTHnSparse(name,thnSize,thnBins,thnMin,thnMax,thnNames);
    
  }
  return thn;
}

int AliDxHFEParticleSelectionEl::InitControlObjects()
{
  /// init control and monitoring objects
  AliInfo("Electron THnSparse");
  
  fElectronProperties=DefineTHnSparse();
  fElectronProperties->Sumw2();
  AddControlObject(fElectronProperties);
  
  fHistoList=new TList;
  fHistoList->SetName("HFE Histograms");
  fHistoList->SetOwner();
  
  // Histogram storing which cuts have been applied to the tracks
  // TODO: revise the names of the histograms, starting with 'f'
  // confuses the reader to think of class members
  fHistoList->Add(CreateControlHistogram("fWhichCut","effective cut for a rejected particle", kNCutLabels, fgkCutBinNames));
  fHistoList->Add(CreateControlHistogram("fTPCnClAOD","Nr TPC clusters/track for all tracks", 160, 0., 159.));
  fHistoList->Add(CreateControlHistogram("fTPCnClSingleTrackCuts","Nr TPC clusters/track After SingleTrackCuts", 160, 0., 159.));
  fHistoList->Add(CreateControlHistogram("fTPCnClTPCTOFPID","Nr TPC clusters/track After TPCTOF PID", 160, 0., 159.));
  
  double dEdxBins[6]={1000,0.,10.,200,0.,200.};
  double nSigBins[6]={1000,0.,10.,200,-10.,10.};
  double eovpBins[6]={100 ,0.,10.,350,-15.,20.};
  double eovppTBins[6]={500,0.,20.,300,0.,5.};
  
  // dEdx plots, TPC signal vs momentum
  fHistoList->Add(CreateControl2DHistogram("fdEdx", "dEdx before cuts", dEdxBins,"momentum (GeV/c)","dE/dx in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fdEdxCut", "dEdx after cuts",dEdxBins,"momentum (GeV/c)","dE/dx in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fdEdxPidTOF", "dEdx after TOF pid",dEdxBins,"momentum (GeV/c)","dE/dx in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fdEdxPidTPC", "dEdx after TPC pid",dEdxBins,"momentum (GeV/c)","dE/dx in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fdEdxPidEMCAL", "dEdx after EMCAL pid",dEdxBins,"momentum (GeV/c)","dE/dx in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fdEdxPid", "dEdx after pid",dEdxBins,"momentum (GeV/c)","dE/dx in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fdEdxRejPi", "dEdx after pion rejection",dEdxBins,"momentum (GeV/c)","dE/dx in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fdEdxRejProton", "dEdx after proton rejection",dEdxBins,"momentum (GeV/c)","dE/dx in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fdEdxIM", "dEdx after Inv Mass",dEdxBins,"momentum (GeV/c)","dE/dx in TPC (a.u.)"));
  
  // nSigmaTPC vs momentum
  fHistoList->Add(CreateControl2DHistogram("fnSigTPC", "nSigTPC before cuts",nSigBins,"momentum (GeV/c)","nSigma in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fnSigTPCCut", "nSigmaTPC after cuts",nSigBins,"momentum (GeV/c)","nSigma in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fnSigTPCPidTPC", "nSigmaTPC after TPC PID", nSigBins,"momentum (GeV/c)","nSigma in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fnSigTPCPidEMCAL", "nSigmaTPC after EMCAL PID", nSigBins,"momentum (GeV/c)","nSigma in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fnSigTPCPidTOF", "nSigmaTPC after TOF PID", nSigBins,"momentum (GeV/c)","nSigma in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fnSigTPCPid", "nSigmaTPC after PID", nSigBins,"momentum (GeV/c)","nSigma in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fnSigTPCRejPi", "nSigmaTPC after pion rejection", nSigBins,"momentum (GeV/c)","nSigma in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fnSigTPCRejProton", "nSigmaTPC after proton rejection", nSigBins,"momentum (GeV/c)","nSigma in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fnSigTPCIM", "nSigmaTPC after Inv Mass", nSigBins,"momentum (GeV/c)","nSigma in TPC (a.u.)"));
  
  // nSigmaTOF vs momentum
  fHistoList->Add(CreateControl2DHistogram("fnSigTOF", "nSigmaTOF before cuts",nSigBins,"momentum (GeV/c)","nSigma in TOF (a.u.)"));  
  fHistoList->Add(CreateControl2DHistogram("fnSigTOFCut", "nSigmaTOF after cuts",nSigBins,"momentum (GeV/c)","nSigma in TOF (a.u.)"));  
  fHistoList->Add(CreateControl2DHistogram("fnSigTOFPidTOF", "nSigmaTOF after TOF PID", nSigBins,"momentum (GeV/c)","nSigma in TOF (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fnSigTOFPidTPC", "nSigmaTOF after TPC PID", nSigBins,"momentum (GeV/c)","nSigma in TOF (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fnSigTOFPidEMCAL", "nSigmaTOF after EMCAL PID", nSigBins,"momentum (GeV/c)","nSigma in TOF (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fnSigTOFPid", "nSigmaTOF after PID", nSigBins,"momentum (GeV/c)","nSigma in TOF (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fnSigTOFRejPi", "nSigmaTOF after pion rejection", nSigBins,"momentum (GeV/c)","nSigma in TOF (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fnSigTOFRejProton", "nSigmaTOF after proton rejection", nSigBins,"momentum (GeV/c)","nSigma in TOF (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fnSigTOFIM", "nSigmaTOF after Inv Mass", nSigBins,"momentum (GeV/c)","nSigma in TOF (a.u.)"));
  
  // E/p vs TPC nSig
  fHistoList->Add(CreateControl2DHistogram("feopCut", "E/p after single track cuts", eovpBins,"E/p","nSigma in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("feopEMCAL", "E/p after EMCAL PID", eovpBins,"E/p","nSigma in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("feopTPC", "E/p after TPC PID", eovpBins,"E/p","nSigma in TPC (a.u.)"));
  
  // E/p vs pT
  fHistoList->Add(CreateControl2DHistogram("feoppTCut", "E/p vs pT after single track cuts", eovppTBins,"momentum (GeV/c)","E/p"));
  fHistoList->Add(CreateControl2DHistogram("feoppTEMCAL", "E/p vs pT after EMCAL PID", eovppTBins,"momentum (GeV/c)","E/p"));
  fHistoList->Add(CreateControl2DHistogram("feoppTTPC", "E/p vs pT after TPC PID", eovppTBins,"momentum (GeV/c)","E/p"));
  
  // Invariant mass LS and ULS without cut
  fHistoList->Add(CreateControlHistogram("fInvMassLS", "Invariant mass LS", 1000, 0., 0.5));
  fHistoList->Add(CreateControlHistogram("fInvMassULS", "Invariant mass ULS", 1000, 0., 0.5));
  
  // Invariant mass LS and ULS without cut, extended x-axis
  fHistoList->Add(CreateControlHistogram("fInvMass2SelLS", "Invariant mass two selected particles LS", 1000, 0., 5.));
  fHistoList->Add(CreateControlHistogram("fInvMass2SelULS", "Invariant mass two selected particles ULS", 1000, 0., 5.));
  
  // Invariant mass LS and ULS with cut
  // TODO: Remove if remove first try at invariant mass
  fHistoList->Add(CreateControlHistogram("fInvMass2SelLScut", "Invariant mass two selected particles LS (cut)", 1000, 0., 0.5));
  fHistoList->Add(CreateControlHistogram("fInvMass2SelULScut", "Invariant mass two selected particles ULS (cut)", 1000, 0., 0.5));

  fHistoList->Add(CreateControlHistogram("fPtULScut", "Pt particle cut from ULS ", 100, 0., 10.));
  fHistoList->Add(CreateControlHistogram("fPtLScut", "Pt particle cut from LS ", 100, 0., 10.));

  fSelNHFE->SetHistMass((TH1F*)fHistoList->FindObject("fInvMassULS"));
  fSelNHFE->SetHistMassBack((TH1F*)fHistoList->FindObject("fInvMassLS"));
  
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
  memset(data, 0, dimension*sizeof(data[0]));
  data[i++]=track->Pt();
  data[i++]=track->Phi();
  data[i++]=track->Eta();
  if (i<dimension) {
    if(fStoreCutStepInfo) data[i]=GetLastSurvivedCutsStep();
    i++; // take out of conditionals to be save
  }
  
  return i;
}


TObjArray* AliDxHFEParticleSelectionEl::Select(const AliVEvent* pEvent)
{
  /// create selection from 'Tracks' member of the event,
  /// array contains only pointers but does not own the objects
  /// object array needs to be deleted by caller
  if (!pEvent) return NULL;
  
  TList* selectedTracks=new TList;
  fSelNHFE->SetPIDresponse(fPIDResponse);
  // TODO: move out to analysis task
  fPIDTPC->SetPIDResponse(fPIDResponse);
  //Previously only fPIDTPC had fPIDResponse set.
  fPIDTOFTPC->SetPIDResponse(fPIDResponse); 
  fPIDTOF->SetPIDResponse(fPIDResponse);
  TObjArray* finalTracks=new TObjArray;
  if (!finalTracks) return NULL;
  finalTracks->SetOwner(kFALSE); // creating new track objects below
  TObjArray* finalTracksIM=new TObjArray;
  if (!finalTracksIM) return NULL;
  finalTracksIM->SetOwner(kFALSE); // creating new track objects below
  
  int nofTracks=pEvent->GetNumberOfTracks();
  fTrackNum=0;
  int selectionCode=0;
  //  int origin
  for (int itrack=0; itrack<nofTracks; itrack++) {
    selectionCode=0;
    AliVParticle* track=pEvent->GetTrack(itrack);
    selectionCode=IsSelected(track,pEvent);
    fTrackNum++; 
    // If fUseInvMassCut is set to true, only call HistogramParticleProperties
    // (here) for those who didn't pass single track cuts and pid. 
    // if fUseInvMassCut is not set, then call the function for all tracks
    // TODO: Change this when have added invariant mass i thnsparse
    //CHANGED FOR OLD INV MASS
    if (selectionCode==0) continue;
    AliReducedParticle *trackReduced2 =(AliReducedParticle*)CreateParticle(track);
    finalTracks->Add(trackReduced2);
    if(fUseInvMassCut!=kInvMassTwoSelected || selectionCode==0)
      HistogramParticleProperties(trackReduced2, selectionCode);
    selectedTracks->Add(track);
  }

  // Calculating invariant mass for electron pairs with same selection criteria
  // TODO: Remove
  if(selectedTracks->GetSize()>0)
    {
      Bool_t* selTrackIndex=new Bool_t[selectedTracks->GetSize()];
      if (selTrackIndex) {
      InvMassFilter(selectedTracks, selTrackIndex);

      //Only remove electrons if fUseInvMassCut is set
      if(fUseInvMassCut==kInvMassTwoSelected){
	for(Int_t k=0; k<selectedTracks->GetSize(); k++)
	  {
	    selectionCode=0; //Reset selectionCode here to be based on invariant mass cut 
	    //On the basis that selectedTracks and finalTracks contain same particles
	    AliAODTrack *trackIM=(AliAODTrack*)selectedTracks->At(k);
	    AliReducedParticle *trackReduced=(AliReducedParticle*)finalTracks->At(k);
	    if(selTrackIndex[k]==kTRUE)
	      {
		((TH2D*)fHistoList->FindObject("fdEdxIM"))->Fill(trackIM->GetTPCmomentum(), trackIM->GetTPCsignal());
		finalTracksIM->Add(trackReduced);
		selectionCode=1;
	      }
	    HistogramParticleProperties(trackReduced, selectionCode);
	  }
      }
      delete [] selTrackIndex;
      }
      
    }

  HistogramEventProperties(AliDxHFEParticleSelection::kHistoNrTracksPrEvent,finalTracks->GetEntries());
  if(fUseInvMassCut!=kInvMassTwoSelected) 
    return finalTracks; 
  else
    return finalTracksIM;

}

int AliDxHFEParticleSelectionEl::IsSelected(AliVParticle* pEl, const AliVEvent* pEvent)
{
  /// select El candidates
  if(!pEvent){
    AliError("No event information");
    return 0;
  }  
  fSurvivedCutStep=kNotSelected;
  
  AliAODTrack *track=(AliAODTrack*)pEl;
  fCFM->SetRecEventInfo(pEvent);
  Double_t pPart = track->P();
  
  ((TH2D*)fHistoList->FindObject("fdEdx"))->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
  ((TH2D*)fHistoList->FindObject("fnSigTPC"))->Fill(track->GetTPCmomentum(), fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
  ((TH2D*)fHistoList->FindObject("fnSigTOF"))->Fill(pPart, fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron));
  ((TH1D*)fHistoList->FindObject("fTPCnClAOD"))->Fill(track->GetTPCNcls());
  Double_t fClsE = -999, p = -999, fEovP=-999, pt = -999;
  
  if(fFinalCutStep==kNoCuts){
    Float_t radial=999;
    Bool_t acceptTrack=kFALSE;
    
    //Cut on Eta
    if(TMath::Abs(track->Eta()) < fEtaCut) {
      acceptTrack=kTRUE;
    }
    
    const Double_t kBeampiperadius=fImpactParamCutRadial;
    Double_t dcaD[2]={-999.,-999.}, covD[3]={-999.,-999.,-999.};
    
    AliAODEvent *aodevent = (AliAODEvent*)(pEvent);
    if(!aodevent) {
      AliDebug(1, "No aod event available\n");
      return 0;
    }
    
    AliAODVertex *vtxAODSkip  = aodevent->GetPrimaryVertex();
    if(!vtxAODSkip) return 0;
    AliExternalTrackParam etp; etp.CopyFromVTrack(track);
    if(etp.PropagateToDCA(vtxAODSkip, aodevent->GetMagneticField(), kBeampiperadius, dcaD, covD)) {
      radial = dcaD[0];
      //z = dcaD[1];
    }
    
    // Cut on z  
    if(acceptTrack && TMath::Abs(radial) < fImpactParamCutRadial){
      acceptTrack=kTRUE;
    }
    
    if(acceptTrack)
      fSurvivedCutStep=kNoCuts;
    
    if(!acceptTrack) return 0;
    
    // Should only return at this point if you only want to store all tracks
    AliDebug(2,"Returns after kNoCuts"); 
    ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kSelected); 
    return 1;
  }
  
  //--------track cut selection-----------------------
  // Filter Bit:
  if(fSetFilterBit){
    if (!track->TestFilterMask(BIT(fBit))){
      AliDebug(2,Form("cut due to filter bit %d",fBit));
      return 0;
    }
  }
  if(fStopAfterFilterBit) return 1;
  
  //Using AliHFECuts:
  // RecKine: ITSTPC cuts  
  if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track)){
    AliDebug(4,"Cut: kStepRecKineITSTPC");
    ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kRecKineITSTPC);
    return 0;
  }
  if(fStoreCutStepInfo) fSurvivedCutStep=kRecKineITSTPC;
  if(fFinalCutStep==kRecKineITSTPC) {AliDebug(2,"Returns after kStepRecKineITSTPC"); ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kSelected); return 1;}
  
  // RecPrim
  if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) {
    AliDebug(4,"Cut: kStepRecPrim");
    ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kRecPrim);
    if(!fStoreCutStepInfo) return 0;
    else return 1; //return 1 because it passed cuts above, but not this (no need to go further)
  }
  if(fStoreCutStepInfo) fSurvivedCutStep=kRecPrim;
  if(fFinalCutStep==kRecPrim) {AliDebug(2,"Returns after kRecPrim"); ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kSelected); return 1;}
  
  // HFEcuts: ITS layers cuts
  if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track)) {
    AliDebug(4,"Cut: kStepHFEcutsITS");
    ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kHFEcutsITS);
    if(!fStoreCutStepInfo) return 0;
    else return 1; //return 1 because it passed cuts above, but not this (no need to go further)
    
  }
  if(fStoreCutStepInfo) fSurvivedCutStep=kHFEcutsITS;
  if(fFinalCutStep==kHFEcutsITS) {AliDebug(2,"Returns after kHFEcutsITS "); ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kSelected); return 1;}

  
  // if fMaxPtCombinedPID is set to lower than upper Ptlimit (10GeV/c), will separate
  // PID into two regions: below fMaxptCombinedPID - both TPC and TOF, above only TPC
  Bool_t useCombinedPID=kTRUE;
  Bool_t okTOF = CheckTOFPIDStatus(track);
  if(okTOF){
    if((pPart<fMaxPTOFWhenPresent && fMaxPTOFWhenPresent>-1) || pPart<fMaxPtCombinedPID){ //If the momentum is within either of the two regions, use TOF PID
      useCombinedPID=kTRUE;
    } else {//If not within the range where tof is to be used
      useCombinedPID=kFALSE; //Do not use combined pid (tof + tpc)
    }
  }

  if(useCombinedPID) AliDebug(2,Form("P: %f, use CombinedPID (fMaxPtCombinedPID= %f)",pPart,fMaxPtCombinedPID));
  else AliDebug(2,Form("P: %f, use only TPC PID (fMaxPtCombinedPID= %f)",pPart,fMaxPtCombinedPID));

  if(!fUseEMCAL){
    // HFE cuts: TOF PID and mismatch flag
    if(useCombinedPID && !okTOF){//Using TOF PID, but okTOF fails
      if(pPart<fMaxPtCombinedPID){//TOF forced to be used
	AliDebug(4,"Cut: kStepHFEcutsTOF");
	((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kHFEcutsTOF);
	if(!fStoreCutStepInfo){ return 0;}
      	else {return 1;} //return 1 because it passed cuts above, but not this (no need to go further)
      }else{//Above the pPart region where tof is forced, turn off combined PID
	AliDebug(2,"Don't use forced TOF anymore, no TOF mathing of the track");
	useCombinedPID=kFALSE; 
      }
    }
    if(fStoreCutStepInfo) fSurvivedCutStep=kHFEcutsTOF;
    if(fFinalCutStep==kHFEcutsTOF) {AliDebug(2,"Returns after kHFEcutsTOF"); ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kSelected); return 1;}
  }
  
  // HFE cuts: TPC PID cleanup
  if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, track)){
    AliDebug(4,"Cut: kStepHFEcutsTPC");
    ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kHFEcutsTPC);
    if(!fStoreCutStepInfo) return 0;
    else return 1; //return 1 because it passed cuts above, but not this (no need to go further)
  } 
  if(fStoreCutStepInfo) fSurvivedCutStep=kHFEcutsTPC;
  
  ((TH2D*)fHistoList->FindObject("fdEdxCut"))->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
  ((TH2D*)fHistoList->FindObject("fnSigTPCCut"))->Fill(track->GetTPCmomentum(), fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
  ((TH2D*)fHistoList->FindObject("fnSigTOFCut"))->Fill(pPart, fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron));
  ((TH1D*)fHistoList->FindObject("fTPCnClSingleTrackCuts"))->Fill(track->GetTPCNcls());
  
  if(fFinalCutStep==kHFEcutsTPC) {AliDebug(2,"Returns after track cuts"); ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kSelected); return 1;}
  
  
  //--------PID selection-----------------------
  AliHFEpidObject hfetrack;
  hfetrack.SetAnalysisType(AliHFEpidObject::kAODanalysis);
  hfetrack.SetRecTrack(track);
  
  // TODO: configurable colliding system
  // TODO: Check problem with PbPb, for now workaround with setting multiplicity to -1
  if(GetSystem()==1) {
    hfetrack.SetPbPb();
    hfetrack.SetCentrality(-1); 
  }
  else  hfetrack.SetPP();
  //  hfetrack.SetMulitplicity(ncontribVtx);
  
  // TODO: Put this into a while-loop instead, looping over the number of pid objects in the cut-list?
  // This needs a bit of thinking and finetuning (wrt histogramming)

  //---------EMCAL PID----------//
  if(fUseEMCAL)
    {
      pt = track->Pt();      
      // Track extrapolation to EMCAL
      Int_t fClsId = track->GetEMCALcluster();
      if(fClsId <0) return 0;
      AliVCluster *cluster = pEvent->GetCaloCluster(fClsId);
      if(!cluster->IsEMCAL()) return 0;
      if(TMath::Abs(cluster->GetTrackDx())>0.05 || TMath::Abs(cluster->GetTrackDz())>0.05) return 0;    
      
      fClsE = cluster->E();
      fEovP = fClsE/pPart;
      
      //Electron id with shower shape  
      
      if(cluster->GetM20()< fMinM20 || cluster->GetM20()> fMaxM20 || 
	 cluster->GetM02()< fMinM02 || cluster->GetM02()> fMaxM02 || 
	 cluster->GetDispersion()> fDispersion) return 0;
      
      //Electron id with E/p
      ((TH2D*)fHistoList->FindObject("feopCut"))->Fill(fEovP, fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
      ((TH2D*)fHistoList->FindObject("feoppTCut"))->Fill(pt, fEovP);//track->GetTPCmomentum(), fEovP);
      if(fEovP < fEovPMin || fEovP > fEovPMax) return 0;
      ((TH2D*)fHistoList->FindObject("fdEdxPidEMCAL"))->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
      ((TH2D*)fHistoList->FindObject("fnSigTPCPidEMCAL"))->Fill(track->GetTPCmomentum(), fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
      ((TH2D*)fHistoList->FindObject("fnSigTOFPidEMCAL"))->Fill(pPart, fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron));
      ((TH2D*)fHistoList->FindObject("feopEMCAL"))->Fill(fEovP, fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
      ((TH2D*)fHistoList->FindObject("feoppTEMCAL"))->Fill(pt, fEovP);//track->GetTPCmomentum(), fEovP);
    }
  
  // TPC PID
  if(fPIDTPC && fPIDTPC->IsSelected(&hfetrack)) {
    ((TH2D*)fHistoList->FindObject("fdEdxPidTPC"))->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
    ((TH2D*)fHistoList->FindObject("fnSigTPCPidTPC"))->Fill(track->GetTPCmomentum(), fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
    ((TH2D*)fHistoList->FindObject("fnSigTOFPidTPC"))->Fill(pPart, fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron));
    if(fStoreCutStepInfo){fSurvivedCutStep=kPIDTPC; }
    if(fFinalCutStep==kPIDTPC) {AliDebug(2,"Returns at PIDTPC"); ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kSelected); return 1;}
  }
  else{
    ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kPIDTPC);
    if(!useCombinedPID) {return 0; }// if only use combined PID, return 0 here (not selected by TPC PID)
    if(fUseEMCAL) {return 0; }// if use EMCAL and TPC return here
  }
  if(fUseEMCAL)
    { ((TH2D*)fHistoList->FindObject("feopTPC"))->Fill(fEovP, fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
      ((TH2D*)fHistoList->FindObject("feoppTTPC"))->Fill(pt, fEovP);//track->GetTPCmomentum(), fEovP);
    }
  if(fFinalCutStep==kPIDTPC) {AliDebug(2,"Returns at PIDTPC"); return 0;}
  

  if(!fUseEMCAL && useCombinedPID){

    //TOF PID
    if(fPIDTOF && fPIDTOF->IsSelected(&hfetrack)) {
      ((TH2D*)fHistoList->FindObject("fdEdxPidTOF"))->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
      ((TH2D*)fHistoList->FindObject("fnSigTPCPidTOF"))->Fill(track->GetTPCmomentum(), fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
      ((TH2D*)fHistoList->FindObject("fnSigTOFPidTOF"))->Fill(pPart, fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron));
      if(fStoreCutStepInfo){fSurvivedCutStep=kPIDTOF; }
      if(fFinalCutStep==kPIDTOF) {AliDebug(2,"Returns at PIDTOF");((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kSelected); return 1;}
    }
    else{
      ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kPIDTOF);
    }
    if(fFinalCutStep==kPIDTOF) {AliDebug(2,"Returns at PIDTOF"); return 0;}

    //Combined TOF & TPC PID
    if(fPIDTOFTPC && fPIDTOFTPC->IsSelected(&hfetrack)) {
      AliDebug(3,"Inside FilldPhi, electron is selected");
      if(fStoreCutStepInfo){
	fSurvivedCutStep=kPIDTOFTPC;
      }
    }
    else{
      ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kPIDTOFTPC);
      AliDebug(4,"Cut: kTPCTOFPID");
      if(!fStoreCutStepInfo) { return 0;}
      else if(fStoreCutStepInfo){  return 1; }//return 1 because it passed cuts above, but not this (no need to go further)
    }
  }

  // Filling histograms with particles passing PID criteria
  // (Filled here due to the option of separating regions for TPC+TOF and TPC)
  ((TH2D*)fHistoList->FindObject("fdEdxPid"))->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
  ((TH2D*)fHistoList->FindObject("fnSigTPCPid"))->Fill(track->GetTPCmomentum(), fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
  if(!fUseEMCAL) ((TH2D*)fHistoList->FindObject("fnSigTOFPid"))->Fill(pPart, fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron));
  ((TH1D*)fHistoList->FindObject("fTPCnClTPCTOFPID"))->Fill(track->GetTPCNcls());
  if(fFinalCutStep==kPIDTOFTPC) {AliDebug(2,"Returns at PIDTOFTPC"); return 1;}

  // [FIXME]: Under development
  // Rejection of pions and protons based on 3. nsigma cut in tpc
  Double_t nsigmaTPCpi=fPIDResponse->NumberOfSigmasTPC((AliVParticle*)track,(AliPID::EParticleType)AliPID::kPion);
  Double_t nsigmaTPCp=fPIDResponse->NumberOfSigmasTPC((AliVParticle*)track,(AliPID::EParticleType)AliPID::kProton);

  //Pion rejection
  if(fTPCnSigmaPiRej>0){ //Pion rejection on by default, but can be turned off by setting fTPCnSigmaPiRej to -1
    if(TMath::Abs(nsigmaTPCpi)<fTPCnSigmaPiRej && pPart>fPiRejPMin && pPart<fPiRejPMax){ //Only applied above between given P() limits 
      // Track within pion identification range. Reject
      ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kRejPi);
      return 0; 
    }else{
      ((TH2D*)fHistoList->FindObject("fdEdxRejPi"))->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
      ((TH2D*)fHistoList->FindObject("fnSigTPCRejPi"))->Fill(track->GetTPCmomentum(), fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
      if(!fUseEMCAL) ((TH2D*)fHistoList->FindObject("fnSigTOFRejPi"))->Fill(pPart, fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron));
    }
  }  
  if(fFinalCutStep==kRejPi) {AliDebug(2,"Returns after pion rejection"); return 1;}  
  //Proton rejection
  if(fTPCnSigmaPRej>0){ //Proton rejection on by default, but can be turned off by setting fTPCnSigmaPRej to -1
    if(TMath::Abs(nsigmaTPCp)<fTPCnSigmaPRej && pPart>fPRejPMin && pPart<fPRejPMax){ //Only applied above 1GeV
      // Track within proton identification range. Reject
      ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kRejProton);
      return 0;
    }else{
      ((TH2D*)fHistoList->FindObject("fdEdxRejProton"))->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
      ((TH2D*)fHistoList->FindObject("fnSigTPCRejProton"))->Fill(track->GetTPCmomentum(), fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
      if(!fUseEMCAL) ((TH2D*)fHistoList->FindObject("fnSigTOFRejProton"))->Fill(pPart, fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron));
    }
  }
  if(fFinalCutStep==kRejProton) {AliDebug(2,"Returns after proton rejection"); return 1;}  

  //Removing electrons based on invariant mass method
  if(fUseInvMassCut==kInvMassSingleSelected)
    {
      AliDebug(4,"invmass check");
      bool isLS=kFALSE;
      bool isULS=kFALSE;
      fSelNHFE->FindNonHFE(fTrackNum, track, const_cast<AliVEvent*>(pEvent));
      if(fSelNHFE->IsULS())
	{
	  //Not selected
	  isULS=kTRUE;
	  ((TH1F*)fHistoList->FindObject("fPtULScut"))->Fill(track->Pt());
	  ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kINVMASS);
	  AliDebug(4,"Cut: Invmass");
	  if(!fStoreCutStepInfo) return 0;
	  
	}
      if(fSelNHFE->IsLS()){
	((TH1F*)fHistoList->FindObject("fPtLScut"))->Fill(track->Pt());
	if(fCutLS){
	  //Not selected
	  isLS=kTRUE;
	  ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kINVMASS);
	  AliDebug(4,"Cut: Invmass");
	  if(!fStoreCutStepInfo) return 0;
	}
      }
      
      if(fStoreCutStepInfo && !(isLS || isULS)) fSurvivedCutStep=kINVMASS;
      
    }
  ((TH2D*)fHistoList->FindObject("fdEdxIM"))->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
  ((TH2D*)fHistoList->FindObject("fnSigTPCIM"))->Fill(track->GetTPCmomentum(), fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
  if(!fUseEMCAL) ((TH2D*)fHistoList->FindObject("fnSigTOFIM"))->Fill(pPart, fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron));
  ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kSelected);
  if(fFinalCutStep==kINVMASS) {AliDebug(2,"Returns after invariant mass"); return 1;}  
  return 1;
  
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
	    AliError(Form("Cut object is not of required type AliHFEcuts but %s", obj->ClassName()));
	}
	if(iii==2){
	  fPIDTOFTPC=dynamic_cast<AliHFEpid*>(obj);
	  if (!fPIDTOFTPC) 
	    AliError(Form("(TOFTPC) cuts object is not of required type AliHFEpid but %s", obj->ClassName()));
	}
	if(iii==3){ 
	  fPIDTOF=dynamic_cast<AliHFEpid*>(obj);
	  if (!fPIDTOF) 
	    AliError(Form("(TOF) cuts object is not of required type AliHFEpid but %s", obj->ClassName()));
	}
	if(iii==4){
	  fPIDTPC=dynamic_cast<AliHFEpid*>(obj);
	  if (!fPIDTPC) 
	    AliError(Form("(TPC) cuts object is not of required type AliHFEpid but %s", obj->ClassName()));
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
  unique_ptr<TObjArray> tokens(strArguments.Tokenize(" "));
  if (!tokens.get()) return 0;
  
  AliInfo(strArguments);
  TIter next(tokens.get());
  TObject* token;
  while ((token=next())) {
    TString argument=token->GetName();
    if(argument.BeginsWith("elreco=")){
      argument.ReplaceAll("elreco=", "");
      int selectionStep=kPIDTOFTPC; //Default
      if(argument.CompareTo("alltracks")==0) selectionStep=kNoCuts;      
      else if(argument.CompareTo("afterreckineitstpc")==0) selectionStep=kRecKineITSTPC;
      else if(argument.CompareTo("afterrecprim")==0) selectionStep=kRecPrim;
      else if(argument.CompareTo("afterhfeits")==0) selectionStep=kHFEcutsITS;
      else if(argument.CompareTo("afterhfetof")==0) selectionStep=kHFEcutsTOF;
      else if(argument.CompareTo("afterhfetpc")==0) selectionStep=kHFEcutsTPC;
      else if(argument.CompareTo("aftertrackcuts")==0) selectionStep=kHFEcutsTPC;
      else if(argument.CompareTo("aftertofpid")==0) selectionStep=kPIDTOF;
      else if(argument.CompareTo("aftertpcpid")==0) selectionStep=kPIDTPC;
      else if(argument.CompareTo("afterfullpid")==0) selectionStep=kPIDTOFTPC;
      else AliFatal(Form("unknown argument '%s'", argument.Data()));
      SetFinalCutStep(selectionStep); 
      continue;   
    }
    if(argument.BeginsWith("twoselectedinvmasscut")){
      fUseInvMassCut=kInvMassTwoSelected;
      AliInfo("Using Invariant mass cut for two selected particles");
      continue;   
    }
    if(argument.BeginsWith("useinvmasscut")){
      fUseInvMassCut=kInvMassSingleSelected;
      AliInfo("Using Invariant mass cut for single selected particle and looser cuts on partner");
      continue;   
    }
    if(argument.BeginsWith("EMCALPID")){
      fUseEMCAL=kTRUE;
      AliInfo("Using EMCAL PID");
      continue;   
    }
    if(argument.BeginsWith("maxPtCombinedPID=")){
      argument.ReplaceAll("maxPtCombinedPID=", "");
      fMaxPtCombinedPID=argument.Atof();
      AliInfo(Form("Using only TPC PID over %f GeV/c",fMaxPtCombinedPID));
      continue;
    }
    if(argument.BeginsWith("maxPTOFWhenPresent=")){
      argument.ReplaceAll("maxPTOFWhenPresent=", "");
      fMaxPTOFWhenPresent=argument.Atof();
      AliInfo(Form("TOF PID used when present up to %f GeV/c",fMaxPTOFWhenPresent));
      fUseTOFonlyWhenPresent=kTRUE;
      AliInfo("Using TOF when present");
      continue;
    }
    //v--- [FIXME] Should no longer be used
    if(argument.BeginsWith("onlyTOFwhenpresent") || argument.BeginsWith("TOFwhenpresent") || argument.BeginsWith("TOFwhenthere")){
      fUseTOFonlyWhenPresent=kTRUE;
      AliInfo("Only using TOF when present");
      continue;
    }
    if(argument.BeginsWith("StopFB") || argument.BeginsWith("stopafterfilterbit") || argument.BeginsWith("stopfb")){
      fStopAfterFilterBit=kTRUE;
      AliInfo("Only using TOF when present");
      continue;
    }
    if(argument.BeginsWith("notusefilterbit")){
      fSetFilterBit=kFALSE;
      AliInfo("Switching off use of Filter Bit");
      continue;   
    }
    if(argument.BeginsWith("filterbit=")){
      argument.ReplaceAll("filterbit=", "");
      fBit=argument.Atoi();
      AliInfo(Form("Using filter bit: %d",fBit));
      fSetFilterBit=kTRUE;
      continue;   
    }
    if(argument.BeginsWith("invmasscut=")){
      argument.ReplaceAll("invmasscut=", "");
      fInvMassLow=argument.Atof();
      AliInfo(Form("Using invariant mass cut: %f",fInvMassLow));
      //      fUseInvMassCut=kInvMassSingleSelected;
      continue;   
    }
    if(argument.BeginsWith("impactparamcut=")){
      argument.ReplaceAll("impactparamcut=", "");
      fImpactParamCutRadial=argument.Atof();
      AliInfo(Form("Using impact parameter cut: %f",fImpactParamCutRadial));
      continue;   
    }
    if(argument.BeginsWith("etacut=")){
      argument.ReplaceAll("etacut=", "");
      fEtaCut=argument.Atof();
      AliInfo(Form("Using Eta cut: %f",fEtaCut));
      continue;   
    }
    if(argument.BeginsWith("M20min=")){
      argument.ReplaceAll("M20min=", "");
      fMinM20=argument.Atof();
      AliInfo(Form("Using M20min cut: %f",fMinM20));
      continue;   
    }
    if(argument.BeginsWith("M20max=")){
      argument.ReplaceAll("M20max=", "");
      fMaxM20=argument.Atof();
      AliInfo(Form("Using M20max cut: %f",fMaxM20));
      continue;   
    }
    
    if(argument.BeginsWith("M02min=")){
      argument.ReplaceAll("M02min=", "");
      fMinM02=argument.Atof();
      AliInfo(Form("Using M02min cut: %f",fMinM02));
      continue;   
    }
    if(argument.BeginsWith("M02max=")){
      argument.ReplaceAll("M02max=", "");
      fMaxM02=argument.Atof();
      AliInfo(Form("Using M02max cut: %f",fMaxM02));
      continue;   
    }
    
    if(argument.BeginsWith("Dispersion=")){
      argument.ReplaceAll("Dispersion=", "");
      fDispersion=argument.Atof();
      AliInfo(Form("Using Dispersion cut: %f",fDispersion));
      continue;   
    }
    
    if(argument.BeginsWith("EoPmin=")){
      argument.ReplaceAll("EoPmin=", "");
      fEovPMin=argument.Atof();
      AliInfo(Form("Using E/p min cut: %f",fEovPMin));
      continue;   
    }
    if(argument.BeginsWith("EoPmax=")){
      argument.ReplaceAll("EoPmax=", "");
      fEovPMax=argument.Atof();
      AliInfo(Form("Using E/p max cut: %f",fEovPMax));
      continue;   
    }
    if(argument.BeginsWith("storelastcutstep")){
      AliInfo("Stores the last cut step");
      SetStoreLastCutStep(kTRUE);
      continue;
    }
    if(argument.BeginsWith("cutLS")){
      AliInfo("Cut also on LS distribution of electron pairs");
      fCutLS=kTRUE;
      continue;
    }
    if(argument.BeginsWith("RejP=")){ //Disable proton rejection by setting to -1
      argument.ReplaceAll("RejP=", "");
      fTPCnSigmaPRej=argument.Atof();
      AliInfo(Form("Using proton rejection in TPC, nsigma+/-: %f",fTPCnSigmaPRej));
      continue;   
    }
    if(argument.BeginsWith("Rejpi=")){ //Disable pion rejection by setting to -1
      argument.ReplaceAll("Rejpi=", "");
      fTPCnSigmaPiRej=argument.Atof();
      AliInfo(Form("Using pion rejection in TPC, nsigma+/-: %f",fTPCnSigmaPiRej));
      continue;   
    }
    if(argument.BeginsWith("minPRejP=")){
      argument.ReplaceAll("minPRejP=", "");
      fPRejPMin=argument.Atof();
      AliInfo(Form("Lower p limit for use of Proton rejection: %f GeV/c",fPRejPMin));
      continue;
    }
    if(argument.BeginsWith("maxPRejP=")){
      argument.ReplaceAll("maxPRejP=", "");
      fPRejPMax=argument.Atof();
      AliInfo(Form("Upper p limit for use of Proton rejection: %f GeV/c",fPRejPMax));
      continue;
    }
    if(argument.BeginsWith("minPRejpi=")){
      argument.ReplaceAll("minPRejpi=", "");
      fPiRejPMin=argument.Atof();
      AliInfo(Form("Lower p limit for use of Pion rejection: %f GeV/c",fPiRejPMin));
      continue;
    }
    if(argument.BeginsWith("maxPRejpi=")){
      argument.ReplaceAll("maxPRejpi=", "");
      fPiRejPMax=argument.Atof();
      AliInfo(Form("Upper p limit for use of Pion rejection: %f GeV/c",fPiRejPMax));
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
		//		if(mass<0.5)
		// ((TH1D*)fHistoList->FindObject("fInvMass2SelLScut"))->Fill(mass);
		((TH1D*)fHistoList->FindObject("fInvMass2SelLS"))->Fill(mass);
		if(mass>fInvMassLow)
		  {
		    ((TH1D*)fHistoList->FindObject("fInvMass2SelLScut"))->Fill(mass);
		    selIndx[i]=kTRUE;
		    selIndx[j]=kTRUE;
		  }
	      }
	      if(fFlagULS) {
		//if(mass<0.5)
		//  ((TH1D*)fHistoList->FindObject("fInvMass2SelULScut"))->Fill(mass);
		((TH1D*)fHistoList->FindObject("fInvMass2SelULS"))->Fill(mass);
		if(mass>fInvMassLow)
		  {
		    ((TH1D*)fHistoList->FindObject("fInvMass2SelULScut"))->Fill(mass);
		    selIndx[i]=kTRUE;
		    selIndx[j]=kTRUE;
		  }

	      }
	      
	    }
	}
    }
}

//___________________________________________________________________________
Bool_t AliDxHFEParticleSelectionEl::CheckTOFPIDStatus(AliAODTrack *track) const{
  // Borrowed from AliSingleTrackEffCuts.cxx
  // Check TOC PID status
  //

  if ((track->GetStatus()&AliESDtrack::kTOFout)==0)   return kFALSE;
  if ((track->GetStatus()&AliESDtrack::kTIME)==0)     return kFALSE;
  if ((track->GetStatus()&AliESDtrack::kTOFpid)==0)   return kFALSE;
  if ((track->GetStatus()&AliESDtrack::kTOFmismatch)!=0) return kFALSE;
  return kTRUE;
}
