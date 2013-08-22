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
  , fFinalCutStep(kPIDTOFTPC)
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
  , fUseEMCAL(kFALSE)
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
  fCuts->Initialize(fCFM);

  //Setting up TPC PID
  //Add settings for asymmetric cut on nSigma TPC
  //TODO: have this completely set up from addtask 
  const int paramSize=4;
  Double_t params[paramSize];
  memset(params, 0, sizeof(Double_t)*paramSize);
  params[0]=-1.;
  fPIDTPC = new AliHFEpid("hfePidTPC");
  if(!fPIDTPC->GetNumberOfPIDdetectors()) { 
    fPIDTPC->AddDetector("TPC",1);
  }
  fPIDTPC->ConfigureTPCdefaultCut(NULL, params, 3.);
  fPIDTPC->InitializePID();

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
    int    thnBinsExt[thnSizeExt] = { 100,  100, 100,   kNCutLabels-1};
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
    int    thnBins[thnSize] = { 100,  100, 100 };
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

  // dEdx plots, TPC signal vs momentum
  fHistoList->Add(CreateControl2DHistogram("fdEdx", "dEdx before cuts", dEdxBins,"momentum (GeV/c)","dE/dx in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fdEdxCut", "dEdx after cuts",dEdxBins,"momentum (GeV/c)","dE/dx in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fdEdxPidTOF", "dEdx after TOF pid",dEdxBins,"momentum (GeV/c)","dE/dx in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fdEdxPidTPC", "dEdx after TPC pid",dEdxBins,"momentum (GeV/c)","dE/dx in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fdEdxPidEMCAL", "dEdx after EMCAL pid",dEdxBins,"momentum (GeV/c)","dE/dx in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fdEdxPid", "dEdx after pid",dEdxBins,"momentum (GeV/c)","dE/dx in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fdEdxIM", "dEdx after Inv Mass",dEdxBins,"momentum (GeV/c)","dE/dx in TPC (a.u.)"));

  // nSigmaTPC vs momentum
  fHistoList->Add(CreateControl2DHistogram("fnSigTPC", "nSigTPC before cuts",nSigBins,"momentum (GeV/c)","nSigma in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fnSigTPCCut", "nSigmaTPC after cuts",nSigBins,"momentum (GeV/c)","nSigma in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fnSigTPCPidTPC", "nSigmaTPC after TPC PID", nSigBins,"momentum (GeV/c)","nSigma in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fnSigTPCPidEMCAL", "nSigmaTPC after EMCAL PID", nSigBins,"momentum (GeV/c)","nSigma in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fnSigTPCPidTOF", "nSigmaTPC after TOF PID", nSigBins,"momentum (GeV/c)","nSigma in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fnSigTPCPid", "nSigmaTPC after PID", nSigBins,"momentum (GeV/c)","nSigma in TPC (a.u.)"));

  // nSigmaTOF vs momentum
  fHistoList->Add(CreateControl2DHistogram("fnSigTOF", "nSigmaTOF before cuts",nSigBins,"momentum (GeV/c)","nSigma in TOF (a.u.)"));  
  fHistoList->Add(CreateControl2DHistogram("fnSigTOFCut", "nSigmaTOF after cuts",nSigBins,"momentum (GeV/c)","nSigma in TOF (a.u.)"));  
  fHistoList->Add(CreateControl2DHistogram("fnSigTOFPidTOF", "nSigmaTOF after TOF PID", nSigBins,"momentum (GeV/c)","nSigma in TOF (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fnSigTOFPidTPC", "nSigmaTOF after TPC PID", nSigBins,"momentum (GeV/c)","nSigma in TOF (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fnSigTOFPidEMCAL", "nSigmaTOF after EMCAL PID", nSigBins,"momentum (GeV/c)","nSigma in TOF (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("fnSigTOFPid", "nSigmaTOF after PID", nSigBins,"momentum (GeV/c)","nSigma in TOF (a.u.)"));

  // E/p
  fHistoList->Add(CreateControl2DHistogram("feopEMCAL", "E/p after EMCAL PID", eovpBins,"E/p","nSigma in TPC (a.u.)"));
  fHistoList->Add(CreateControl2DHistogram("feopTPC", "E/p after TPC PID", eovpBins,"E/p","nSigma in TPC (a.u.)"));

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
  fPIDTPC->SetPIDResponse(fPIDResponse);
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
  // TODO: Could be removed if it is verified that looser cuts is better, but keep
  // for now
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
  if(!pEvent){
    AliError("No event information");
    return 0;
  }  
  fSurvivedCutStep=kNotSelected;

  AliAODTrack *track=(AliAODTrack*)pEl;
  //  AliVTrack *trackv = dynamic_cast<AliVTrack*>(track);
  fCFM->SetRecEventInfo(pEvent);
  
  ((TH2D*)fHistoList->FindObject("fdEdx"))->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
  ((TH2D*)fHistoList->FindObject("fnSigTPC"))->Fill(track->GetTPCmomentum(), fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
  ((TH2D*)fHistoList->FindObject("fnSigTOF"))->Fill(track->P(), fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron));
  ((TH1D*)fHistoList->FindObject("fTPCnClAOD"))->Fill(track->GetTPCNcls());
  Double_t fClsE = -999, p = -999, fEovP=-999, pt = -999, dEdx=-999, fTPCnSigma=0;

  if(fFinalCutStep==kNoCuts){
    Float_t radial=999;
    //Float_t z=999;
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

  // if fMaxPtCombinedPID is set to lower than upper Ptlimit (10GeV/c), will separate
  // PID into two regions: below fMaxptCombinedPID - both TPC and TOF, above only TPC
  Bool_t useCombinedPID=kTRUE;
  if(track->Pt() > fMaxPtCombinedPID) useCombinedPID=kFALSE;
  if(useCombinedPID) AliDebug(2,Form("Pt: %f, use CombinedPID (fMaxPtCombinedPID= %f)",track->Pt(),fMaxPtCombinedPID));
  else AliDebug(2,Form("Pt: %f, use only TPC PID (fMaxPtCombinedPID= %f)",track->Pt(),fMaxPtCombinedPID));
  
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

  if(!fUseEMCAL){
    // HFE cuts: TOF PID and mismatch flag
    if(useCombinedPID && !ProcessCutStep(AliHFEcuts::kStepHFEcutsTOF, track)) {
      if(!useCombinedPID) cout << "should not be here "<< track->Pt() << endl;
      AliDebug(4,"Cut: kStepHFEcutsTOF");
      ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kHFEcutsTOF);
      if(!fStoreCutStepInfo) return 0;
      else return 1; //return 1 because it passed cuts above, but not this (no need to go further)
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
  ((TH2D*)fHistoList->FindObject("fnSigTOFCut"))->Fill(track->P(), fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron));
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

  if(fUseEMCAL)
    {
      // Development of EMCAL PID //

      pt = track->Pt();
      p = track->P();
      //    dEdx = track->GetTPCsignal();
      //fTPCnSigma = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron);//fPIDTPCEMCAL->GetPIDResponse() ? fPIDTPCEMCAL->GetPIDResponse()->NumberOfSigmasTPC(track,AliPID::kElectron) : 1000; [FIX] check if this change is ok
      //if(pt<2){
      //Cutreason: Out of bounds, pt
      //return 0;
      //}
      // Crude cut on TPC nSigma, -3,3
      //     if(fTPCnSigma < -3 || fTPCnSigma > 3){
      //Cutreason: Out of bounds, TPC nSigma
      // return 0;
      //}
      //eta cut (-0.7,0.7)
      //    if(track->Eta() < -0.7 || track->Eta() > 0.7){
      //Cutreason: Out of bounds, eta
      //     return 0;
      //   }
      // printf("EMCAL4\n");
      // Track extrapolation to EMCAL
      Int_t fClsId = track->GetEMCALcluster();
      //      printf("%d\n",fClsId);
      if(fClsId <0) return 0;
      AliVCluster *cluster = pEvent->GetCaloCluster(fClsId);
      if(!cluster->IsEMCAL()) return 0;
      //      printf("EMCAL5\n");
      if(TMath::Abs(cluster->GetTrackDx())>0.05 || TMath::Abs(cluster->GetTrackDz())>0.05) return 0;    
      //fdEdxBef->Fill(p,dEdx);
      //fTPCnsigma->Fill(p,fTPCnSigma);
      
      //     fTrkpt->Fill(pt);
      fClsE = cluster->E();
      fEovP = fClsE/p;
            
      //Electron id with TPC
      //    if(fTPCnSigma < fTPCnsigEleMin || fTPCnSigma > fTPCnsigEleMax) continue;
      //    fEovPWoSS->Fill(pt,fEovP);
      //    fElecPhiTPCEovP->Fill(track->Phi());
      
      //Electron id with shower shape  
      /*[FIX] Look into this later    
	if(cluster->GetM20()< fM20CutMin || cluster->GetM20()> fM20CutMax || cluster->GetM02()< fM02CutMin || cluster->GetM02()> fM02CutMax || cluster->GetDispersion()> fDispCutMax) continue;
	fEovPWSS->Fill(pt,fEovP);
      */
      //Electron id with E/p
      Double_t fEovPMin=0.8;
      Double_t fEovPMax=1.2;
      if(fEovP < fEovPMin || fEovP > fEovPMax) return 0;
      //      printf("Track selected by EMCAL only\n");
      //    fTrkEovPAft->Fill(pt,fEovP);
      //    fElecPhi->Fill(track->Phi());
      //    fElecPhiPt->Fill(track->Phi(),track->Pt());
      //    if (track->Eta() >0 && track->Eta() <0.7) fElecPhiTPChalf->Fill(track->Phi());
      ((TH2D*)fHistoList->FindObject("fdEdxPidEMCAL"))->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
      ((TH2D*)fHistoList->FindObject("fnSigTPCPidEMCAL"))->Fill(track->GetTPCmomentum(), fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
      ((TH2D*)fHistoList->FindObject("fnSigTOFPidEMCAL"))->Fill(track->P(), fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron));
      ((TH2D*)fHistoList->FindObject("feopEMCAL"))->Fill(fEovP, fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
      // /Development of EMCAL PID //
    }
  
  // TODO: Put this into a while-loop instead, looping over the number of pid objects in the cut-list?
  // This needs a bit of thinking and finetuning (wrt histogramming)
  if(!fUseEMCAL){
    if(fPIDTOF && fPIDTOF->IsSelected(&hfetrack)) {
      ((TH2D*)fHistoList->FindObject("fdEdxPidTOF"))->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
      ((TH2D*)fHistoList->FindObject("fnSigTPCPidTOF"))->Fill(track->GetTPCmomentum(), fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
      ((TH2D*)fHistoList->FindObject("fnSigTOFPidTOF"))->Fill(track->P(), fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron));
      if(fStoreCutStepInfo){fSurvivedCutStep=kPIDTOF; }
      if(fFinalCutStep==kPIDTOF) {AliDebug(2,"Returns at PIDTOF");((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kSelected); return 1;}
    }
    else{
      ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kPIDTOF);
    }
    
    if(fFinalCutStep==kPIDTOF) {AliDebug(2,"Returns at PIDTOF"); return 0;}
  }
  //if(useCombinedPID) AliInfo(Form("Pt: %f, use CombinedPID (fMaxPtCombinedPID= %f)",track->Pt(),fMaxPtCombinedPID));
  //else AliInfo(Form("Pt: %f, use only TPC PID (fMaxPtCombinedPID= %f)",track->Pt(),fMaxPtCombinedPID));

  if(fPIDTPC && fPIDTPC->IsSelected(&hfetrack)) {
    ((TH2D*)fHistoList->FindObject("fdEdxPidTPC"))->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
    ((TH2D*)fHistoList->FindObject("fnSigTPCPidTPC"))->Fill(track->GetTPCmomentum(), fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
    ((TH2D*)fHistoList->FindObject("fnSigTOFPidTPC"))->Fill(track->P(), fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron));
    if(fStoreCutStepInfo){fSurvivedCutStep=kPIDTPC; }
    if(fFinalCutStep==kPIDTPC) {AliDebug(2,"Returns at PIDTPC"); ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kSelected); return 1;}
    //if(!useCombinedPID)  cout << "Will be selected" << endl;
  }
  else{
    ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kPIDTPC);
    if(!useCombinedPID) {return 0; }// if only use combined PID, return 0 here (not selected by TPC PID)
  }
  if(fUseEMCAL) ((TH2D*)fHistoList->FindObject("feopTPC"))->Fill(fEovP, fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
  if(fFinalCutStep==kPIDTPC) {AliDebug(2,"Returns at PIDTPC"); return 0;}

  //Combined tof & tpc pid
  if(!fUseEMCAL){
  if(fPIDTOFTPC && fPIDTOFTPC->IsSelected(&hfetrack)) {
    AliDebug(3,"Inside FilldPhi, electron is selected");
    if(fStoreCutStepInfo){
      fSurvivedCutStep=kPIDTOFTPC;
    }
  }
  else{
    ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kPIDTOFTPC);
    AliDebug(4,"Cut: kTPCTOFPID");
    //Should do nothing if not use combinedPID
    //    if(!useCombinedPID) cout << "HERE" << endl;
    if(!fStoreCutStepInfo && useCombinedPID) { return 0;}
    else if(fStoreCutStepInfo){  return 1; }//return 1 because it passed cuts above, but not this (no need to go further)
  }
  //if(!useCombinedPID) cout << "HERE" << endl;
  }
  // Filling histograms with particles passing PID criteria
  // (Filled here due to the option of separating regions for TPC+TOF and TPC)
  ((TH2D*)fHistoList->FindObject("fdEdxPid"))->Fill(track->GetTPCmomentum(), track->GetTPCsignal());
  ((TH2D*)fHistoList->FindObject("fnSigTPCPid"))->Fill(track->GetTPCmomentum(), fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron));
  if(!fUseEMCAL) ((TH2D*)fHistoList->FindObject("fnSigTOFPid"))->Fill(track->P(), fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron));
  ((TH1D*)fHistoList->FindObject("fTPCnClTPCTOFPID"))->Fill(track->GetTPCNcls());

  //if(fStoreCut

  //Removing electrons based on invariant mass method
  if(fUseInvMassCut==kInvMassSingleSelected)
    {
      AliDebug(4,"invmass check");
      fSelNHFE->FindNonHFE(fTrackNum, track, const_cast<AliVEvent*>(pEvent));
      if(fSelNHFE->IsLS() || fSelNHFE->IsULS())
	{
	  //Not selected
	  ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kINVMASS);
	  AliDebug(4,"Cut: Invmass");
	  if(!fStoreCutStepInfo) return 0;
	}
      else
	if(fStoreCutStepInfo) fSurvivedCutStep=kINVMASS;

    }
  ((TH1D*)fHistoList->FindObject("fWhichCut"))->Fill(kSelected);
  
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
	    AliError(Form("Cut object is not of required type AliHFEcuts but %s", cuts->ClassName()));
	}
	if(iii==2){
	  fPIDTOFTPC=dynamic_cast<AliHFEpid*>(obj);
	  if (!fPIDTOFTPC) 
	    AliError(Form("(TOFTPC) cuts object is not of required type AliHFEpid but %s", cuts->ClassName()));
	}
	if(iii==3){ 
	  fPIDTOF=dynamic_cast<AliHFEpid*>(obj);
	  if (!fPIDTOF) 
	    AliError(Form("(TOF) cuts object is not of required type AliHFEpid but %s", cuts->ClassName()));
	}
	/*if(iii=4){
	  fPIDTPC=dynamic_cast<AliHFEpid*>(obj);
	  if (!fPIDTPC) 
	    AliError(Form("(TPC) cuts object is not of required type AliHFEpid but %s", cuts->ClassName()));
	    }*/
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
      fMaxPtCombinedPID=argument.Atoi();
      AliInfo(Form("Using only TPC PID over %f GeV/c",fMaxPtCombinedPID));
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
    if(argument.BeginsWith("storelastcutstep")){
      AliInfo("Stores the last cut step");
      SetStoreLastCutStep(kTRUE);
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
