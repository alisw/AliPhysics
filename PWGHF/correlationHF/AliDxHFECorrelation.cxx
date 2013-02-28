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
#include "AliAODEvent.h"
#include "AliAODVertex.h"
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
  fCorrProperties=NULL;
  fhEventControlCorr=NULL;
  if(fCorrelator) delete fCorrelator;
  fCorrelator=NULL;
  if(fCorrArray) delete fCorrArray;
  fCorrArray=NULL;

  // NOTE: the external object is deleted elsewhere
  fCuts=NULL;
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
  AddControlObject(fCorrProperties);

  //----------------------------------------------
  // Histogram for storing event information

  TString histoname="";
  if(fTriggerParticleType==kElectron)
    histoname="hEventControlHFExDCorr";
  else
    histoname="hEventControlDxHFECorr";
  std::auto_ptr<TH1D> hEventControl(new TH1D(histoname.Data(), histoname.Data(), 10, 0, 10));
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
  fCorrelator = new AliHFCorrelator("Correlator", cuts, fSystem); 
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
  Int_t MaxNofEvents = cuts->GetMaxNEventsInPool();
  Int_t MinNofTracks = cuts->GetMinNTracksInPool();
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

  if(fSystem) nofEventPropBins = 100; // PbPb centrality
  if(!fSystem) nofEventPropBins = NofCentBins; // pp multiplicity

  Double_t minvalue = CentBins[0];
  Double_t maxvalue = CentBins[NofCentBins];
  Double_t Zminvalue = ZVrtxBins[0];
  Double_t Zmaxvalue = ZVrtxBins[NofCentBins];

  const Double_t Nevents[]={0,2*MaxNofEvents/10,4*MaxNofEvents/10,6*MaxNofEvents/10,8*MaxNofEvents/10,MaxNofEvents};
  const Double_t * events = Nevents;

  TH3D * EventsPerPoolBin = new TH3D("EventsPerPoolBin","Number of events in bin pool",NofCentBins,CentBins,NofZVrtxBins,ZVrtxBins,5,events);
  EventsPerPoolBin->GetXaxis()->SetTitle("Centrality/multiplicity ");
  EventsPerPoolBin->GetYaxis()->SetTitle("Z vertex [cm]");
  EventsPerPoolBin->GetZaxis()->SetTitle("Number of events in pool bin");
  if(fUseEventMixing) AddControlObject(EventsPerPoolBin);

  Int_t MaxNofTracks = (MaxNofEvents+1)*MinNofTracks;
  Int_t Diff = MaxNofTracks-MinNofTracks;

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

  return 0;
}

int AliDxHFECorrelation::ParseArguments(const char* arguments)
{
  // parse arguments and set internal flags
  TString strArguments(arguments);
  auto_ptr<TObjArray> tokens(strArguments.Tokenize(" "));
  if (!tokens.get()) return -ENOMEM;

  TIter next(tokens.get());
  TObject* token;
  while ((token=next())) {
    TString argument=token->GetName();
   
    if (argument.BeginsWith("event-mixing")) {
      fUseEventMixing=true;
      continue;
    }
      
    if (argument.BeginsWith("use-mc")) {
      fUseMC=true;
      continue;
    }
    if (argument.BeginsWith("system=")) {
      argument.ReplaceAll("system=", "");
      if (argument.CompareTo("pp")==0) fSystem=0;
      else if (argument.CompareTo("Pb-Pb")==0) fSystem=1;
      else {
	AliWarning(Form("can not set collision system, unknown parameter '%s'", argument.Data()));
	// TODO: check what makes sense
	fSystem=0;
      }
      continue;
    }
    if (argument.BeginsWith("trigger=")){
      argument.ReplaceAll("trigger=", "");
      if (argument.CompareTo("D")==0) { fTriggerParticleType=kD; AliInfo("Trigger on D"); }
      if (argument.CompareTo("D0")==0) { fTriggerParticleType=kD; AliInfo("Trigger on D");}
      else if (argument.CompareTo("electron")==0){ fTriggerParticleType=kElectron; AliInfo("trigger on electron");}
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
  AliDebug(1, "Creating Corr THnSparse");
  // here is the only place to change the dimension
  static const int sizeEventdphi = 7;  
  InitTHnSparseArray(sizeEventdphi);
  const double pi=TMath::Pi();

  //TODO: add phi for electron??
  // 			                        0           1       2      3         4     5      6   
  // 			                      D0invmass   PtD0    PhiD0  PtbinD0    Pte   dphi   deta   
  int         binsEventdphi[sizeEventdphi] = {   200,      1000,   100,    21,     1000,   100,    100};
  double      minEventdphi [sizeEventdphi] = { 1.5648,      0,       0,     0,       0,  fMinPhi, -2};
  double      maxEventdphi [sizeEventdphi] = { 2.1648,     100,   2*pi,    20,      100, fMaxPhi, 2};
  const char* nameEventdphi[sizeEventdphi] = {
    "D0InvMass",
    "PtD0",
    "PhiD0",
    "PtBinD0",
    "PtEl",
    "#Delta#Phi", 
    "#Delta#eta"
  };

  TString name;
  name.Form("%s info", GetName());


  return CreateControlTHnSparse(name,sizeEventdphi,binsEventdphi,minEventdphi,maxEventdphi,nameEventdphi);

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

  std::auto_ptr<THnSparseD> th(new THnSparseD(name, name, thnSize, thnBins, thnMin, thnMax));
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
      Bool_t analyzetracks = fCorrelator->ProcessAssociatedTracks(jMix, associatedTracks);
			
      if(!analyzetracks) {
	AliError("AliHFCorrelator::Cannot process the track array");
	continue;
      }

      Int_t NofTracks = fCorrelator->GetNofTracks();

      // looping on track candidates
      for(Int_t iTrack = 0; iTrack<NofTracks; iTrack++){ 
	Bool_t runcorrelation = fCorrelator->Correlate(iTrack);
	if(!runcorrelation) continue;
			
	fDeltaPhi = fCorrelator->GetDeltaPhi();
	fDeltaEta = fCorrelator->GetDeltaEta();
	
	AliReducedParticle *assoc = fCorrelator->GetAssociatedParticle();

	FillParticleProperties(ptrigger,assoc,ParticleProperties(),GetDimTHnSparse());
	fCorrProperties->Fill(ParticleProperties());

      } // loop over electron tracks in event
    } // loop over events in pool
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
    data[i++]=ptrigger->Phi();
    data[i++]=ptrigger->GetPtBin(); 
    data[i++]=assoc->Pt();
  } 
  else{
    data[i++]=assoc->GetInvMass();
    data[i++]=assoc->Pt();
    data[i++]=assoc->Phi();
    data[i++]=assoc->GetPtBin(); 
    data[i++]=ptrigger->Pt();
  }
  data[i++]=GetDeltaPhi();
  data[i++]=GetDeltaEta();

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
  std::auto_ptr<TFile> output(TFile::Open(filename, "RECREATE"));
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
  if(!fSystem){ // pp
    multiplicity = AOD->GetNTracks();
    MultipOrCent = multiplicity; // convert from Int_t to Double_t
  }
  if(fSystem){ // PbPb		
    centralityObj = AOD->GetHeader()->GetCentralityP();
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
	
