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

/// @file   AliDxHFECorrelation.cxx
/// @author Sedat Altinpinar, Hege Erdal, Matthias Richter
/// @date   2012-04-25
/// @brief  Worker class for D0-HF electron correlation
///

#include "AliDxHFECorrelation.h"
#include "AliVParticle.h"
#include "AliLog.h"
//#include "AliAnalysisCuts.h"         // required dependency libANALYSISalice.so
//#include "AliFlowTrackSimple.h"      // required dependency libPWGflowBase.so
//#include "AliFlowCandidateTrack.h"   // required dependency libPWGflowTasks.so
//#include "AliCFContainer.h"          // required dependency libCORRFW.so
#include "TObjArray.h"
#include "AliHFCorrelator.h"
#include "AliHFAssociatedTrackCuts.h"
#include "AliVertexingHFUtils.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliRDHFCutsD0toKpi.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "AliReducedParticle.h"
#include "AliDxHFEParticleSelection.h"
#include <iostream>
#include <cerrno>
#include <memory>

using namespace std;

ClassImp(AliDxHFECorrelation)

AliDxHFECorrelation::AliDxHFECorrelation(const char* name)
  : TNamed(name?name:"AliDxHFECorrelation", "")
  , fHistograms(NULL)  
  , fControlObjects(NULL)
  , fCorrProperties(NULL)
  , fhEventControlCorr(NULL)
  , fCuts(NULL)
  , fCutsD0(NULL)
  , fUseMC(kFALSE)
  , fCorrelator(NULL)
  , fUseEventMixing(kFALSE)
  , fSystem(0)
  , fMinPhi(-TMath::Pi()/2)
  , fMaxPhi(3*TMath::Pi()/2)
  , fDeltaPhi(0)
  , fDeltaEta(0)
  , fDimThn(0)
  , fCorrArray(NULL)
  , fEventType(0)
  , fTriggerParticleType(kD)
  , fUseTrackEfficiency(kFALSE)
  , fUseD0Efficiency(kFALSE)
  , fRunMode(kReducedMode)
  , fUseCentrality(0)
  , fPoolBin(-1)
{
  // default constructor
  // 
  //

}

const char* AliDxHFECorrelation::fgkEventControlBinNames[]={
  "nEventsAll",
  "nEventsSelected",
  "nEventsTriggerd",
  "nEventsCorrelated"
};

AliDxHFECorrelation::~AliDxHFECorrelation()
{
  // destructor
  //
  //
  if (fHistograms) delete fHistograms;
  fHistograms=NULL;

  // NOTE: fControlObjects owns the object, and they are deleted in the
  // destructor of TList
  if (fControlObjects) delete fControlObjects;
  fControlObjects=NULL;
  if (fCorrProperties) delete fCorrProperties;
  fCorrProperties=NULL;
  fhEventControlCorr=NULL;
  if(fCorrelator) delete fCorrelator;
  fCorrelator=NULL;
  if(fCorrArray) delete fCorrArray;
  fCorrArray=NULL;

  // NOTE: the external object is deleted elsewhere
  fCuts=NULL;
  fCutsD0=NULL;
}

int AliDxHFECorrelation::Init(const char* arguments)
{
  //
  // Will initialize thnsparse, histogram and AliHFCorrelator
  //
  AliInfo("Initializing correlation objects");
  ParseArguments(arguments);

  //----------------------------------------------
  // Setting up THnSparse 
  fCorrProperties=DefineTHnSparse();
  fCorrProperties->Sumw2();
  AddControlObject(fCorrProperties);

  //----------------------------------------------
  // Histogram for storing event information

  TString histoname="";
  if(fTriggerParticleType==kElectron)
    histoname="hEventControlHFExDCorr";
  else
    histoname="hEventControlDxHFECorr";
  std::unique_ptr<TH1D> hEventControl(new TH1D(histoname.Data(), histoname.Data(), 10, 0, 10));
  if (!hEventControl.get()) {
    return -ENOMEM;
  }
  int iLabel=0;
  for (iLabel=0; iLabel<kNEventControlLabels; iLabel++)
    hEventControl->GetXaxis()->SetBinLabel(iLabel+1, fgkEventControlBinNames[iLabel]);

  fhEventControlCorr=hEventControl.release();
  AddControlObject(fhEventControlCorr);

  //----------------------------------------------
  // AliHFCorrelator for Event Mixing and correlation
  // 
  // fCuts is the hadron cut object, fSystem to switch between pp or PbPb
  AliHFAssociatedTrackCuts* cuts=dynamic_cast<AliHFAssociatedTrackCuts*>(fCuts);
  if (!cuts) {
    if (fCuts)
      AliError(Form("cuts object of wrong type %s, required AliHFAssociatedTrackCuts", fCuts->ClassName()));
    else
      AliError("mandatory cuts object missing");
    return -EINVAL;
  }
  if (cuts->GetNCentPoolBins()==0 || cuts->GetCentPoolBins()==NULL ||
      cuts->GetNZvtxPoolBins()==0 || cuts->GetZvtxPoolBins()==NULL) {
    // the bin information is used further downstream so it
    // needs to be available in order to continue
    AliError(Form("inavlid object %s: bin configuration is mandatory", cuts->GetName()));
    cuts->Dump();
    return -EINVAL;
  }
  cuts->Print("");
  // Fetching out the RDHF-cut objects for D0 to store in output stream
  TObject *obj2=NULL;
  TIter next2(fCutsD0);
  obj2 = next2();
  AliRDHFCutsD0toKpi* copyfCuts=new AliRDHFCutsD0toKpi(dynamic_cast<AliRDHFCutsD0toKpi&>(*obj2));
  if(!copyfCuts){
    AliFatal(Form("cut object is of incorrect type %s, expecting AliRDHFCutsD0toKpi", obj2->ClassName()));
  }
  fCorrelator = new AliHFCorrelator("Correlator", cuts, fUseCentrality, copyfCuts); 
  fCorrelator->SetDeltaPhiInterval(fMinPhi,fMaxPhi); //Correct Phi Interval
  fCorrelator->SetEventMixing(fUseEventMixing);      // mixing Off/On 
  fCorrelator->SetAssociatedParticleType(AliHFCorrelator::kElectron);
  // 0: don't calculate d0; 1: return d0; 2: return d0/d0err
  fCorrelator->SetApplyDisplacementCut(kFALSE); 
  fCorrelator->SetUseMC(fUseMC);
  fCorrelator->SetUseReco(kTRUE); // Reco/MCTruth
  Bool_t pooldef = fCorrelator->DefineEventPool();
	
  if(!pooldef) AliInfo("Warning:: Event pool not defined properly");


  // ============================= EVENT MIXING CHECKS ======================================
  // TODO: Not sure if all 4 histos are needed. Keep for now..	
  // TODO: Set them up more nicely
  Double_t MaxNofEvents = cuts->GetMaxNEventsInPool();
  Double_t MinNofTracks = cuts->GetMinNTracksInPool();
  Int_t NofCentBins = cuts->GetNCentPoolBins();
  const Double_t * CentBins = cuts->GetCentPoolBins();
  const Double_t defaultCentBins[] = {0,100};
  if (NofCentBins==0 || CentBins==NULL) {
    NofCentBins=1; // note: array dimension minus one, because of bin is bound by upper and lower
    CentBins=defaultCentBins;
  }
  Int_t NofZVrtxBins = cuts->GetNZvtxPoolBins();
  const Double_t *ZVrtxBins = cuts->GetZvtxPoolBins();
  const Double_t defaultZVrtxBins[] = {-10,10};
  if (NofZVrtxBins==0 || ZVrtxBins==NULL) {
    NofZVrtxBins=1; // note: array dimension minus one, because of bin is bound by upper and lower
    ZVrtxBins=defaultZVrtxBins;
  }

  Int_t nofEventPropBins =0;

  if(fUseCentrality) nofEventPropBins = 100; // PbPb centrality
  if(!fUseCentrality) nofEventPropBins = NofCentBins; // pp multiplicity

  Double_t minvalue = CentBins[0];
  Double_t maxvalue = CentBins[NofCentBins];
  Double_t Zminvalue = ZVrtxBins[0];
  Double_t Zmaxvalue = ZVrtxBins[NofZVrtxBins];

  const Double_t Nevents[]={0,2*MaxNofEvents/10,4*MaxNofEvents/10,6*MaxNofEvents/10,8*MaxNofEvents/10,MaxNofEvents};
  const Double_t * events = Nevents;

  TH3D * EventsPerPoolBin = new TH3D("EventsPerPoolBin","Number of events in bin pool",NofCentBins,CentBins,NofZVrtxBins,ZVrtxBins,5,events);
  EventsPerPoolBin->GetXaxis()->SetTitle("Centrality/multiplicity ");
  EventsPerPoolBin->GetYaxis()->SetTitle("Z vertex [cm]");
  EventsPerPoolBin->GetZaxis()->SetTitle("Number of events in pool bin");
  if(fUseEventMixing) AddControlObject(EventsPerPoolBin);

  Double_t MaxNofTracks = (MaxNofEvents+1)*MinNofTracks;
  Double_t Diff = MaxNofTracks-MinNofTracks;

  Double_t Ntracks[]={MinNofTracks,MinNofTracks+Diff/5,MinNofTracks+2*Diff/5,MinNofTracks+3*Diff/5,MinNofTracks+4*Diff/5,MaxNofTracks};
  Double_t  * trackN = Ntracks;

  TH3D * NofTracksPerPoolBin = new TH3D("NofTracksPerPoolBin","Number of tracks in bin pool",NofCentBins,CentBins,NofZVrtxBins,ZVrtxBins,5,trackN);
  NofTracksPerPoolBin->GetXaxis()->SetTitle("Centrality/multiplicity ");
  NofTracksPerPoolBin->GetYaxis()->SetTitle("Z vertex [cm]");
  NofTracksPerPoolBin->GetZaxis()->SetTitle("Number of tracks per bin");

  if(fUseEventMixing) AddControlObject(NofTracksPerPoolBin);

  TH2D * NofPoolBinCalls = new TH2D("NofPoolBinCalls","Number of tracks in bin pool",NofCentBins,CentBins,NofZVrtxBins,ZVrtxBins);
  NofPoolBinCalls->GetXaxis()->SetTitle("Centrality/multiplicity ");
  NofPoolBinCalls->GetYaxis()->SetTitle("Z vertex [cm]");
  if(fUseEventMixing) AddControlObject(NofPoolBinCalls);

  TH2D * EventProps = new TH2D("EventProps","Number of tracks in bin pool",nofEventPropBins,minvalue,maxvalue,100,Zminvalue,Zmaxvalue);
  EventProps->GetXaxis()->SetTitle("Centrality/multiplicity ");
  EventProps->GetYaxis()->SetTitle("Z vertex [cm]");
  if(fUseEventMixing) AddControlObject(EventProps);

  TH2D * MultvsVtx = new TH2D("MultvsVtx","Mult vs Vtx",501,0,500,200,Zminvalue,Zmaxvalue);
  MultvsVtx->GetXaxis()->SetTitle("Centrality/multiplicity ");
  MultvsVtx->GetYaxis()->SetTitle("Z vertex [cm]");
  AddControlObject(MultvsVtx);
  return 0;
}

int AliDxHFECorrelation::ParseArguments(const char* arguments)
{
  // parse arguments and set internal flags
  TString strArguments(arguments);
  unique_ptr<TObjArray> tokens(strArguments.Tokenize(" "));
  if (!tokens.get()) return -ENOMEM;

  TIter next(tokens.get());
  TObject* token;
  while ((token=next())) {
    TString argument=token->GetName();
   
    if (argument.BeginsWith("event-mixing")) {
      fUseEventMixing=true;
      continue;
    }
      
    if (argument.BeginsWith("mc") ||
	argument.BeginsWith("use-mc")) {
      fUseMC=true;
      continue;
    }
    if (argument.BeginsWith("useTrackEff")){
      fUseTrackEfficiency=true;
      AliInfo("Applying track efficiency");
      continue;
    }
    if (argument.BeginsWith("reducedMode") || argument.BeginsWith("reducedmode")){
      fRunMode=kReducedMode;
      AliInfo("Running in Reduced mode");
      continue;
    }
    if (argument.BeginsWith("FullMode") || argument.BeginsWith("fullmode")){
      fRunMode=kFullMode;
      AliInfo("Running in Full mode");
      continue;
    }
    if (argument.BeginsWith("useD0Eff")){
      fUseD0Efficiency=true;
      AliInfo("Applying Correction for D0 efficiency");
      continue;
    }
    if (argument.BeginsWith("system=")) {
      argument.ReplaceAll("system=", "");
      if (argument.CompareTo("pp")==0) {fSystem=0; fUseCentrality=0;}
      else if (argument.CompareTo("Pb-Pb")==0) {fSystem=1; fUseCentrality=1;}
      else if (argument.CompareTo("p-Pb")==0) {fSystem=2; fUseCentrality=0;}
      else {
	AliWarning(Form("can not set collision system, unknown parameter '%s'", argument.Data()));
	// TODO: check what makes sense
	fSystem=0;
      }
      continue;
    }
    if (argument.BeginsWith("minphi=")) {
      argument.ReplaceAll("minphi=", "");
      fMinPhi=argument.Atof();
      AliInfo(Form("changing fMinPhi: %f ",fMinPhi));
      continue;
    }
    if (argument.BeginsWith("maxphi=")) {
      argument.ReplaceAll("maxphi=", "");
      if (argument.BeginsWith("2pi")) {
	fMaxPhi=2*TMath::Pi();
      }
      else
	fMaxPhi=argument.Atof();
      AliInfo(Form("changing fMaxPhi: %f ",fMaxPhi));
      continue;
    }
    if (argument.BeginsWith("trigger=")){
      argument.ReplaceAll("trigger=", "");
      if (argument.CompareTo("D")==0) { fTriggerParticleType=kD; AliInfo("Trigger on D"); }
      if (argument.CompareTo("D0")==0) { fTriggerParticleType=kD; AliInfo("Trigger on D");}
      else if (argument.CompareTo("electron")==0){ fTriggerParticleType=kElectron; AliInfo("trigger on electron");}
      continue;
    }  
    if(argument.BeginsWith("runmode=")){
      argument.ReplaceAll("runmode=","");
      cout << argument.Data() << endl;
      if (argument.CompareTo("full")==0) {AliInfo("Run in Full mode"); SetRunFullMode(kTRUE);}
      if (argument.CompareTo("reduced")==0) {AliInfo("Run in Reduced mode"); SetRunFullMode(kFALSE);}
      continue;
    }
    AliWarning(Form("unknown argument '%s'", argument.Data()));
      
  }

  return 0;
}

THnSparse* AliDxHFECorrelation::DefineTHnSparse()
{
  //
  //Defines the THnSparse. For now, only calls CreateControlTHnSparse
  if(RunFullMode())
    AliDebug(1, "Creating Full Corr THnSparse");
  else
    AliDebug(1, "Creating Reduced Corr THnSparse");
  // here is the only place to change the dimension
  static const int sizeEventdphi = 7;  
  static const int sizeEventdphipPb = 6;
  static const int sizeEventdphiReduced = 5;  
  InitTHnSparseArray(sizeEventdphi);
  const double pi=TMath::Pi();
  THnSparse* thn=NULL;

  TString name;
  name.Form("%s info", GetName());

  if(fRunMode==kFullMode){

    // 			                        0           1       2      3         4     5      6   
    // 			                      D0invmass   PtD0    PhiD0  PtbinD0    Pte   dphi   deta   
    int         binsEventdphi[sizeEventdphi] = {   200,      500,   100,    21,     100,   64,    100};
    double      minEventdphi [sizeEventdphi] = { 1.5648,      0,      0,     0,       0,  fMinPhi, -2};
    double      maxEventdphi [sizeEventdphi] = { 2.1648,     50,   2*pi,    20,      10,  fMaxPhi,  2};
    const char* nameEventdphi[sizeEventdphi] = {
      "D0InvMass",
      "PtD0",
      "PhiD0",
      "PtBinD0",
      "PtEl",
      "#Delta#Phi", 
      "#Delta#eta"
    };
    thn=(THnSparse*)CreateControlTHnSparse(name,sizeEventdphi,binsEventdphi,minEventdphi,maxEventdphi,nameEventdphi);
  }
  else if(2==fSystem){ //Reduced bins for p-Pb
    // 			                                   0          1      2       3      4      5    
    // 			                                D0invmass   PtD0     Pte    dphi   deta  poolbin
    int         binsEventdphipPb[sizeEventdphipPb] = {   150,      28,   50,   16,          20,     6   };
    double      minEventdphipPb [sizeEventdphipPb] = { 1.5848,      2,      0,  fMinPhi,    -2,    -0.5  };
    double      maxEventdphipPb [sizeEventdphipPb] = { 2.1848,      16,    10,  fMaxPhi,     2,    5.5  }; 
    const char* nameEventdphipPb[sizeEventdphipPb] = {
      "D0InvMass",
      "PtD0",
      "PtEl",
      "#Delta#Phi", 
      "#Delta#eta",
      "Poolbin"
    };
    thn=(THnSparse*)CreateControlTHnSparse(name,sizeEventdphipPb,binsEventdphipPb,minEventdphipPb,maxEventdphipPb,nameEventdphipPb);

  }
  else{
    // 			                                   0          1      2       3      4 
    // 			                                D0invmass   PtD0     Pte   dphi   deta   
    int         binsEventdphiRed[sizeEventdphiReduced] = {   200,      500,   100,    64,    100};
    double      minEventdphiRed [sizeEventdphiReduced] = { 1.5648,      0,      0,  fMinPhi,  -2};
    double      maxEventdphiRed [sizeEventdphiReduced] = { 2.1648,      50,    10,  fMaxPhi,   2};
    const char* nameEventdphiRed[sizeEventdphiReduced] = {
      "D0InvMass",
      "PtD0",
      "PtEl",
      "#Delta#Phi", 
      "#Delta#eta"
    };
    thn=(THnSparse*)CreateControlTHnSparse(name,sizeEventdphiReduced,binsEventdphiRed,minEventdphiRed,maxEventdphiRed,nameEventdphiRed);

  }
  return thn;

}

THnSparse* AliDxHFECorrelation::CreateControlTHnSparse(const char* name,
							     int thnSize,
							     int* thnBins,
							     double* thnMin,
							     double* thnMax,
							     const char** binLabels) const
{
  //
  // Creates THnSparse.
  //

  AliInfo("Setting up THnSparse");

  std::unique_ptr<THnSparseD> th(new THnSparseD(name, name, thnSize, thnBins, thnMin, thnMax));
  if (th.get()==NULL) {
    return NULL;
  }
  for (int iLabel=0; iLabel<thnSize; iLabel++) {
    th->GetAxis(iLabel)->SetTitle(binLabels[iLabel]);    
   
  }
  return th.release();

}

int AliDxHFECorrelation::AddControlObject(TObject* pObj)
{
  AliInfo("Adding object");
  /// add control object to list, the base class becomes owner of the object
  if (!pObj) return -EINVAL;
  if (!fControlObjects) {
    fControlObjects=new TList;
    if (!fControlObjects) return -ENOMEM;
    fControlObjects->SetOwner();
  }
  if (fControlObjects->FindObject(pObj->GetName())) {
    AliError(Form("ignoring duplicate object '%s' of type %s", pObj->GetName(), pObj->ClassName()));
    return -EEXIST;
  }
  fControlObjects->Add(pObj);
  return 0;
}

int AliDxHFECorrelation::HistogramEventProperties(int bin)
{
  /// histogram event properties
  if (!fhEventControlCorr) return 0;
  fhEventControlCorr->Fill(bin);

  return 0;
}

int AliDxHFECorrelation::Fill(const TObjArray* triggerCandidates, const TObjArray* associatedTracks, const AliVEvent* pEvent)
{
  //
  // will use AliHFCorrelator to process D0-electron pair and then fill THnSparse.
  //
  if (!triggerCandidates || !associatedTracks) return -EINVAL;
  if (!fControlObjects) {
    Init();
  }
  if (!fControlObjects) {
    AliError("Initialisation failed, can not fill THnSparse");
    return -ENODEV;
  }
  // set the event to be processed
  // TODO: change the correlator class to take the const pointer
  if (!fCorrelator) {
    AliError("working class instance fCorrelator missing");
    return -ENODEV;
  }
  AliHFAssociatedTrackCuts* cuts=dynamic_cast<AliHFAssociatedTrackCuts*>(fCuts);
  if (!cuts) {
    if (fCuts)
      AliError(Form("cuts object of wrong type %s, required AliHFAssociatedTrackCuts", fCuts->ClassName()));
    else
      AliError("mandatory cuts object missing");
    return -EINVAL;
  }
  AliAODEvent *AOD= (AliAODEvent*)(pEvent);
  Double_t multEv = (Double_t)(AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(AOD,-1.,1.));
  AliAODVertex *vtx = AOD->GetPrimaryVertex();
  Double_t zvertex = vtx->GetZ();

  fCorrelator->SetAODEvent(dynamic_cast<AliAODEvent*>(const_cast<AliVEvent*>(pEvent))); 

  Bool_t correlatorON = fCorrelator->Initialize(); //define the pool for mixing
  if(!correlatorON) {
    AliError("AliHFCorrelator didn't initialize the pool correctly or processed a bad event");
    return 1;
  }

  TIter itrigger(triggerCandidates);
  TObject* otrigger=NULL;
  int ctrigger=-1;

  // For the moment this is very specific to D0-electron correlation. Should be 
  // changed to be less specific. 
  while ((otrigger=itrigger())!=NULL) {
    // loop over trigger D0 particle
    ctrigger++;
    AliReducedParticle* ptrigger=dynamic_cast<AliReducedParticle*>(otrigger);
    if (!ptrigger)  continue;

    if(AliDxHFECorrelation::GetTriggerParticleType()==kD){
      if(! TestParticle(ptrigger,kD)) continue;
    }
    else if(AliDxHFECorrelation::GetTriggerParticleType()==kElectron)
      if(! TestParticle(ptrigger,kElectron)) continue;

    Double_t phiTrigger = ptrigger->Phi();
    Double_t ptTrigger = ptrigger->Pt();
    Double_t etaTrigger = ptrigger->Eta();

    // set the phi of the D meson in the correct range
    // TODO: Is this correct to do this??
    phiTrigger = fCorrelator->SetCorrectPhiRange(phiTrigger);
    // pass to the object the necessary trigger part parameters
    fCorrelator->SetTriggerParticleProperties(ptTrigger,phiTrigger,etaTrigger); 

    Bool_t execPool = fCorrelator->ProcessEventPool();
    if(fUseEventMixing && !execPool) {
      AliDebug(1,"Mixed event analysis: pool is not ready");
      continue;
    }
    Int_t NofEventsinPool = 1;
    if(fUseEventMixing) NofEventsinPool = fCorrelator->GetNofEventsInPool();
		
    // loop on events in the pool; if it is SE analysis, stops at one
    for (Int_t jMix =0; jMix < NofEventsinPool; jMix++){
      Bool_t analyzetracks = fCorrelator->ProcessAssociatedTracks(jMix, associatedTracks); //[FIXME] This line causes a segfault from AliHFCorrelator in some cases. Needs to be understood
			
      if(!analyzetracks) {
	AliError("AliHFCorrelator::Cannot process the track array");
	continue;
      }

      Int_t NofTracks = fCorrelator->GetNofTracks();

      // looping on track candidates
      for(Int_t iTrack = 0; iTrack<NofTracks; iTrack++){ 
	fPoolBin=-1;
	Bool_t runcorrelation = fCorrelator->Correlate(iTrack);
	if(!runcorrelation) continue;
			
	fDeltaPhi = fCorrelator->GetDeltaPhi();
	fDeltaEta = fCorrelator->GetDeltaEta();
	
	AliReducedParticle *assoc = fCorrelator->GetAssociatedParticle();
	//	Double_t efficiency = assoc->GetWeight();
	if(!assoc) continue;

	// Test associated particle, which is electron if trigger is D, 
	// or D if trigger is electron
	if(AliDxHFECorrelation::GetTriggerParticleType()==kD){
	  if(! TestParticle(assoc,kElectron)) continue;
	}
	else if(AliDxHFECorrelation::GetTriggerParticleType()==kElectron)
	  if(! TestParticle(assoc,kD)) continue;
	Int_t multiplicity = -1;
	Double_t MultipOrCent = -1;
	AliCentrality *centralityObj = 0;
	
	// get the pool for event mixing
	if(!fUseCentrality){ // pp or pPb
	  //[FIXME][Old method] multiplicity = AOD->GetNumberOfTracks();
	  //  fCorrelatorTr = new AliHFCorrelator("CorrelatorTr",fCutsTracks,fSys,fCutsD0);//fSys=0 use multiplicity, =1 use centrality
	  //MultipOrCent = multiplicity; // convert from Int_t to Double_t
	  MultipOrCent=fCorrelator->GetCentrality();
	}
	if(fUseCentrality){ // PbPb		
	  centralityObj = ((AliVAODHeader*)AOD->GetHeader())->GetCentralityP();
	  MultipOrCent = centralityObj->GetCentralityPercentileUnchecked("V0M");
	  //AliInfo(Form("Centrality is %f", MultipOrCent));
	}
	fPoolBin=cuts->GetPoolBin(MultipOrCent, zvertex);
	
	Double_t weight =1.;
	if(fUseTrackEfficiency){
	  if(AliDxHFECorrelation::GetTriggerParticleType()==kElectron)
	    weight=cuts->GetTrackWeight(ptrigger->Pt(),ptrigger->Eta(),zvertex);
	  else
	    weight=cuts->GetTrackWeight(assoc->Pt(),assoc->Eta(),zvertex);
	  AliDebug(2,Form("Vertex: %f  weight: %f ",zvertex, weight));
	}
	if(fUseD0Efficiency){
	  Double_t D0eff=1;
	  if(AliDxHFECorrelation::GetTriggerParticleType()==kD)
 	    D0eff=GetD0Eff(ptrigger, multEv);
	  else
	    D0eff=GetD0Eff(assoc, multEv);
	  weight=weight*D0eff;
	  AliDebug(2,Form("D0eff: %f, combined efficiency: %f",D0eff, weight));
	}

	FillParticleProperties(ptrigger,assoc,ParticleProperties(),GetDimTHnSparse());
	if(weight!=0){ fCorrProperties->Fill(ParticleProperties(),1./weight);} //Due to possibility of empty regions in eff-map, weights may come out as 0, which crashes the code here.
      } // loop over electron tracks in event
    } // loop over events in pool
    ((TH2D*)fControlObjects->FindObject("MultvsVtx"))->Fill(multEv,zvertex); // event properties
  } // loop over trigger particle

  Bool_t updated = fCorrelator->PoolUpdate(associatedTracks);
  if(fUseEventMixing){
    if(!updated) AliDebug(1,"Pool was not updated");
    else {
      EventMixingChecks(pEvent);
      AliDebug(1,"Pool was updated");
    }
  }
  return 0;
}

double AliDxHFECorrelation::GetD0Eff(AliVParticle* tr, Double_t evMult){

  AliReducedParticle *track=(AliReducedParticle*)tr;
  if (!track) return -ENODATA;
  Double_t pt=track->Pt();

  AliHFAssociatedTrackCuts* cuts=dynamic_cast<AliHFAssociatedTrackCuts*>(fCuts);
  if (!cuts) {
    if (fCuts)
      AliError(Form("cuts object of wrong type %s, required AliHFAssociatedTrackCuts", fCuts->ClassName()));
    else
      AliError("mandatory cuts object missing");
    return -EINVAL;
  }

  return cuts->GetTrigWeight(pt,evMult);

}


int AliDxHFECorrelation::FillParticleProperties(AliVParticle* tr, AliVParticle *as, Double_t* data, int dimension) const
{
  // fill the data array from the particle data
  if (!data) return -EINVAL;
  AliReducedParticle *ptrigger=(AliReducedParticle*)tr;
  AliReducedParticle *assoc=(AliReducedParticle*)as;
  if (!ptrigger || !assoc) return -ENODATA;
  int i=0;
  if (dimension!=GetDimTHnSparse()) {
    // TODO: think about filling only the available data and throwing a warning
    return -ENOSPC;
  }
  if(fTriggerParticleType==kD){
    data[i++]=ptrigger->GetInvMass();
    data[i++]=ptrigger->Pt();
    if(RunFullMode()) data[i++]=ptrigger->Phi();
    if(RunFullMode()) data[i++]=ptrigger->GetPtBin(); 
    data[i++]=assoc->Pt();
  } 
  else{
    data[i++]=assoc->GetInvMass();
    data[i++]=assoc->Pt();
    if(RunFullMode()) data[i++]=assoc->Phi();
    if(RunFullMode()) data[i++]=assoc->GetPtBin(); 
    data[i++]=ptrigger->Pt();
  }
  data[i++]=GetDeltaPhi();
  data[i++]=GetDeltaEta();
  if(fSystem==2) data[i++]=fPoolBin;

  return i;
}

void AliDxHFECorrelation::Clear(Option_t * /*option*/)
{
  /// overloaded from TObject: cleanup

  // nothing to be done so far
  return TObject::Clear();
}

void AliDxHFECorrelation::Print(Option_t */*option*/) const
{
  /// overloaded from TObject: print info
  cout << "====================================================================" << endl;
  TObject::Print();
  if (fHistograms) {
    fHistograms->Print();
  }
}

void AliDxHFECorrelation::Draw(Option_t */*option*/)
{
  /// overloaded from TObject: draw histograms
}

TObject* AliDxHFECorrelation::FindObject(const char *name) const
{
  /// overloaded from TObject: find object by name
  if (fControlObjects) {
    return fControlObjects->FindObject(name);
  }
  return NULL;
}

TObject* AliDxHFECorrelation::FindObject(const TObject *obj) const
{
  /// overloaded from TObject: find object by pointer
  if (fControlObjects) {
    return fControlObjects->FindObject(obj);
  }
  return NULL;
}

void AliDxHFECorrelation::SaveAs(const char *filename, Option_t */*option*/) const
{
  /// overloaded from TObject: save to file
  std::unique_ptr<TFile> output(TFile::Open(filename, "RECREATE"));
  if (!output.get() || output->IsZombie()) {
    AliError(Form("can not open file %s from writing", filename));
    return;
  }
  output->cd();
  if (fControlObjects) fControlObjects->Write();
  output->Close();
}

AliDxHFECorrelation& AliDxHFECorrelation::operator+=(const AliDxHFECorrelation& other)
{
  /// add histograms from another instance
  // TODO - need to change this to ThnSparse?
  if (!fHistograms || !other.fHistograms) return *this;
  
  for (int i=0; i<kNofHistograms; i++) {
    if (fHistograms->At(i)==NULL || other.fHistograms->At(i)==NULL) continue;
    TH1* target=reinterpret_cast<TH1*>(fHistograms->At(i));
    TH1* source=reinterpret_cast<TH1*>(other.fHistograms->At(i));
    if (!target || !source) continue;
    TString name(fHistograms->At(i)->GetName());
    if (name.CompareTo(target->GetName())!=0) {
      AliWarning(Form("skipping incompatible objects at position %d: %s vs %s", i, source->GetName(), target->GetName()));
      continue;
    }
    if (source->IsA()!=target->IsA()) {
      AliWarning(Form("skipping incompatible classes at position %d: %s vs %s", i, source->ClassName(), target->ClassName()));
      continue;
    }
    target->Add(source);
  }
  return *this;
}


//____________________________  Run checks on event mixing ___________________________________________________
void AliDxHFECorrelation::EventMixingChecks(const AliVEvent* pEvent){
	
  AliAODEvent *AOD= (AliAODEvent*)(pEvent);
  AliCentrality *centralityObj = 0;
  Int_t multiplicity = -1;
  Double_t MultipOrCent = -1;
	
  // get the pool for event mixing
  if(!fUseCentrality){ // pp
    /* Old method
    multiplicity = AOD->GetNumberOfTracks();
    MultipOrCent = multiplicity; // convert from Int_t to Double_t */
    MultipOrCent=fCorrelator->GetCentrality();
  }
  if(fUseCentrality){ // PbPb		
    centralityObj = ((AliVAODHeader*)AOD->GetHeader())->GetCentralityP();
    MultipOrCent = centralityObj->GetCentralityPercentileUnchecked("V0M");
    AliInfo(Form("Centrality is %f", MultipOrCent));
  }
	
  AliAODVertex *vtx = AOD->GetPrimaryVertex();
  Double_t zvertex = vtx->GetZ(); // zvertex

  AliEventPool *pool = fCorrelator->GetPool();

  ((TH2D*)fControlObjects->FindObject("NofPoolBinCalls"))->Fill(MultipOrCent,zvertex); // number of calls of pool
  ((TH2D*)fControlObjects->FindObject("EventProps"))->Fill(MultipOrCent,zvertex); // event properties
  ((TH3D*)fControlObjects->FindObject("EventsPerPoolBin"))->Fill(MultipOrCent,zvertex,pool->NTracksInPool()); // number of events in the pool
  ((TH3D*)fControlObjects->FindObject("NofTracksPerPoolBin"))->Fill(MultipOrCent,zvertex,pool->GetCurrentNEvents()); // number of calls of pool
}
	
Bool_t AliDxHFECorrelation::TestParticle(AliVParticle* /*p*/, Int_t /*id*/){

  // Testing particle - for now mainly needed in MC
  return kTRUE;
}
