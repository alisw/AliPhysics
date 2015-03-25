

/**************************************************************************
 * Copyright(c) 1998-2012, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id: AliAnalysisTaskSEDplusCorrelations.cxx 58883 2012-10-02 09:41:01Z prino $ */
/*  AliAnalysisTask for HF(DPlus:3Prongs)-Hadron/Kaon/K0 azimuthal correlations
    By: Jitendra Kumar
*/

#include <Riostream.h>
#include <TClonesArray.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TChain.h>
#include <TString.h>
#include <THnSparse.h>
#include <TDatabasePDG.h>

#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliAnalysisManager.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODPidHF.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSEDplusCorrelations.h"
#include "AliNormalizationCounter.h"
#include "AliVParticle.h"
#include "AliHFAssociatedTrackCuts.h"
#include "AliReducedParticle.h"
#include "AliHFCorrelator.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliCentrality.h"
#include "AliMultiplicity.h"
#include "AliAODTracklets.h"
#include "AliVertexingHFUtils.h"


using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskSEDplusCorrelations)

//____________________| Default Constructor
AliAnalysisTaskSEDplusCorrelations::AliAnalysisTaskSEDplusCorrelations():
AliAnalysisTaskSE(),
  fSelect(0),
  fOutput(0x0),
  fOutputCorr(0X0),
  fReadMC(kFALSE),
  fReco(kFALSE),
  fMontecarlo(kFALSE),
  fMCGenEvType(kFALSE),
  fMixing(kFALSE),
  farrayMC(0x0),
  fSystem(kFALSE),
  fUseBit(kTRUE),
  fTCconfig(kFALSE),
  fHistNEvents(0),
  fCounter(0x0),
  fDplusCuts(0),
  fAssoCuts(0),
  fCorrelator(0x0),
  fNPtBins(0),
  fBinWidth(0),
  fCentrOrMult(-99),
  fMultiplicity(1),
  fEffTrack(kFALSE),
  fEffDplus(kFALSE),
  fCentralityEstimator(0),
  fEvalCentrality(kFALSE),
  fMinCentrality(0),
  fMaxCentrality(100)
{
  //  Default constructor
}


//____________________| Specific Constructor
AliAnalysisTaskSEDplusCorrelations::AliAnalysisTaskSEDplusCorrelations(const Char_t* name,AliRDHFCutsDplustoKpipi* DplusCuts, AliHFAssociatedTrackCuts *AsscCuts) :
  AliAnalysisTaskSE(name),
  fSelect(0),
  fOutput(0x0),
  fOutputCorr(0X0),
  fReadMC(kFALSE),
  fReco(kFALSE),
  fMontecarlo(kFALSE),
  fMCGenEvType(kFALSE),
  fMixing(kFALSE),
  farrayMC(0x0),
  fSystem(kFALSE),
  fUseBit(kTRUE),
  fTCconfig(kFALSE),
  fHistNEvents(0),
  fCounter(0x0),
  fDplusCuts(0),
  fAssoCuts(AsscCuts),
  fCorrelator(0x0),
  fNPtBins(0),
  fBinWidth(0.002),
  fCentrOrMult(-99),
  fMultiplicity(1),
  fEffTrack(kFALSE),
  fEffDplus(kFALSE),
  fCentralityEstimator(0),
  fEvalCentrality(kFALSE),
  fMinCentrality(0),
  fMaxCentrality(100)
{
  
  Info("AliAnalysisTaskSEDplusCorrelations","Calling Constructor");
  fNPtBins=DplusCuts->GetNPtBins();
  fDplusCuts=DplusCuts;
  
  DefineInput(0, TChain::Class());
  DefineOutput(1,TList::Class()); // Basic output slot (more needed)
  DefineOutput(2,TList::Class()); // Correlations form Data and MC
  DefineOutput(3,AliRDHFCutsDplustoKpipi::Class());  //D meson cuts
  DefineOutput(4,AliHFAssociatedTrackCuts::Class());  // Associated tracks cuts
  DefineOutput(5,AliNormalizationCounter::Class()); //Norm 
  
}


//____________________| Soruce Operator
AliAnalysisTaskSEDplusCorrelations::AliAnalysisTaskSEDplusCorrelations(const AliAnalysisTaskSEDplusCorrelations &source):
  AliAnalysisTaskSE(source),
  fSelect(source.fSelect),
  fOutput(source.fOutput),
  fOutputCorr(source.fOutputCorr),
  fReadMC(source.fReadMC),
  fReco(source.fReadMC),
  fMontecarlo(source.fMontecarlo),
  fMCGenEvType(source.fMCGenEvType),
  fMixing(source.fMixing),
  farrayMC(source.farrayMC),
  fSystem(source.fMixing),
  fUseBit(source.fUseBit),
  fTCconfig(source.fTCconfig),
  fHistNEvents(source.fHistNEvents),
  fCounter(source.fCounter),
  fDplusCuts(source.fDplusCuts),
  fAssoCuts(source.fAssoCuts),
  fCorrelator(source.fCorrelator),
  fNPtBins(source.fNPtBins),
  fBinWidth(source.fBinWidth),
  fCentrOrMult(source.fCentrOrMult),
  fMultiplicity(source.fMultiplicity),
  fEffTrack(source.fEffTrack),
  fEffDplus(source.fEffDplus),
  fCentralityEstimator(source.fCentralityEstimator),
  fEvalCentrality(source.fEvalCentrality),
  fMinCentrality(source.fMinCentrality),
  fMaxCentrality(source.fMaxCentrality)
{

}


//____________________| Destructor
AliAnalysisTaskSEDplusCorrelations::~AliAnalysisTaskSEDplusCorrelations()
{
  Info("AliAnalysisTaskSEDplusCorrelations","Calling Destructor");
  
  if(fOutputCorr) {delete fOutputCorr; fOutputCorr = 0;}
  if(fOutput) {delete fOutput; fOutput = 0;}
  if(farrayMC) {delete farrayMC; farrayMC = 0;}
  if(fHistNEvents) {delete fHistNEvents; fHistNEvents = 0;}
  if(fCounter) {delete fCounter; fCounter = 0;}
  if(fDplusCuts) {delete fDplusCuts; fDplusCuts = 0;}
  if(fAssoCuts) {delete fAssoCuts; fAssoCuts=0;}
  if(fCorrelator) {delete fCorrelator; fCorrelator = 0;}
  
}

//____________________| Assignment Operator
AliAnalysisTaskSEDplusCorrelations& AliAnalysisTaskSEDplusCorrelations::operator=(const AliAnalysisTaskSEDplusCorrelations& orig)
{  
  if (&orig == this) return *this; //if address is the same (same object), returns itself
  AliAnalysisTaskSE::operator=(orig); //Uses the AliAnalysisTaskSE operator to assign the inherited part of the class
  fSelect = orig.fSelect;
  fOutput = orig.fOutput;
  fOutputCorr = orig.fOutputCorr;
  fReadMC = orig.fReadMC;
  fReco = orig.fReco;
  fMontecarlo = orig.fMontecarlo;
  fMCGenEvType = orig.fMCGenEvType;
  fMixing = orig.fMixing;
  farrayMC = orig.farrayMC;
  fSystem = orig.fSystem;
  fUseBit = orig.fUseBit;
  fTCconfig=orig.fTCconfig;
  fHistNEvents = orig.fHistNEvents;
  fCounter = orig.fCounter;
  fDplusCuts = orig.fDplusCuts;
  fAssoCuts = orig.fAssoCuts;
  fCorrelator = orig.fCorrelator;
  fNPtBins = orig.fNPtBins;
  fBinWidth = orig.fBinWidth;
  fCentrOrMult = orig.fCentrOrMult;
  fMultiplicity = orig.fMultiplicity;
  fEffTrack = orig.fEffTrack;
  fEffDplus = orig.fEffDplus;
  fCentralityEstimator = orig.fCentralityEstimator;
  fEvalCentrality=orig.fEvalCentrality;
  fMinCentrality=orig.fMinCentrality;
  fMaxCentrality=orig.fMaxCentrality;
  return *this; //returns pointer of the class

}


//____________________| Inilizations
void AliAnalysisTaskSEDplusCorrelations::Init()
{
  
  Info("AliAnalysisTaskSEDplusCorrelations","Calling Inilizations");

  
  //Copy of cuts objects
  AliRDHFCutsDplustoKpipi* DplusCutsCopy = new AliRDHFCutsDplustoKpipi(*fDplusCuts);
  if(!DplusCutsCopy)return;
  const char* name=GetOutputSlot(3)->GetContainer()->GetName();
  DplusCutsCopy->SetName(name);

  
  if(fEvalCentrality){
    Bool_t isCentEstimatorOk = kFALSE;
    if(fCentralityEstimator==fDplusCuts->GetUseCentrality())isCentEstimatorOk = kTRUE; //
    if(!isCentEstimatorOk) {
      AliFatal(Form("Chosen centrality estimator is conflict with D cuts file"));
      return;
    }
  }

  
  //Posting cuts for D+ and Associated tracks
  PostData(3,DplusCutsCopy); 
  PostData(4,fAssoCuts); 
 
  return;
}



//____________________| UserCreateOutputObjects
void AliAnalysisTaskSEDplusCorrelations::UserCreateOutputObjects()
{
  
  Info("AliAnalysisTaskSEDplusCorrelations","Creating UserCreateOutputObjects");
  
  // Category one: Basic O/P
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("BasicHistograms");
  
  fOutputCorr = new TList();
  fOutputCorr->SetOwner();
  fOutputCorr->SetName("CorrelationsHistolist");
  
  fHistNEvents = new TH1F("fHistNEvents", "number of events ", 11, -0.5 , 9.5);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"nEvents analyzed");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"Rejected due to Physics Sel");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"Within choosen centrality");
  fHistNEvents->GetXaxis()->SetBinLabel(4,"Event after correlator interface");
  fHistNEvents->GetXaxis()->SetBinLabel(5,"Events after bad vertex rejection");
  fHistNEvents->GetXaxis()->SetBinLabel(6,"Total Dplus candidate");
  fHistNEvents->GetXaxis()->SetBinLabel(7,"Dplus cand not passing DplusCuts");
  fHistNEvents->GetXaxis()->SetBinLabel(8,"Dplus after passing DplusCuts");
  fHistNEvents->GetXaxis()->SetBinLabel(9,"Dplus after topological cuts");
  fHistNEvents->GetXaxis()->SetBinLabel(10,"Dplus after Vertex cleaning and TC");
  fHistNEvents->GetXaxis()->SetBinLabel(11,"Dplus after fiducial accptnc: finally accpeted");
  fHistNEvents->GetXaxis()->SetNdivisions(1,kFALSE);
  fHistNEvents->SetMinimum(0);
  fHistNEvents->Sumw2();
  fOutput->Add(fHistNEvents);
  
  HistoNomenclature();// Hist Nomenclature function
  
  TString normName="NormalizationCounter"; 
  AliAnalysisDataContainer *cont = GetOutputSlot(5)->GetContainer();
  if(cont)normName=(TString)cont->GetName();
  fCounter = new AliNormalizationCounter(normName.Data());
  fCounter->Init();
  
  //Setting up from Correlator  
  Double_t Pi = TMath::Pi();
  fCorrelator = new AliHFCorrelator("Correlator",fAssoCuts,fSystem,fDplusCuts);
  fCorrelator->SetDeltaPhiInterval(-0.5*Pi, 1.5*Pi); 
  fCorrelator->SetEventMixing(fMixing); 
  fCorrelator->SetAssociatedParticleType(fSelect);
  // fCorrelator->SetApplyDisplacementCut(fDisplacement); //set kFALSE/kTRUE for using the displacement cut
  fCorrelator->SetUseReco(fReco);
  fCorrelator->SetUseMC(fMontecarlo);
  
  Bool_t pooldef = fCorrelator->DefineEventPool();
  if(!pooldef) AliInfo("Warning:: Event pool not defined properly");
  
  // Post Data
  PostData(1,fOutput);
  PostData(2,fOutputCorr);
  PostData(5,fCounter);  
  
}



//____________________| UserExec 
void AliAnalysisTaskSEDplusCorrelations::UserExec(Option_t *) {
  
  Info("AliAnalysisTaskSEDplusCorrelations","Running ... UserExec");
  
  AliAODEvent *aod = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!aod) return;
  
  TClonesArray *array3Prong = 0;
  array3Prong=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
  if(!aod && AODEvent() && IsStandardAOD()) {
    
    aod = dynamic_cast<AliAODEvent*> (AODEvent());
    AliAODHandler* aodHandler = (AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    
    if(aodHandler->GetExtensions()) { 
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      array3Prong=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");   
    }
  }
  else if(!aod || !array3Prong){
    printf("AliAnalysisTaskSEDplusCorrelationselation::UserExec: AOD  Charm3Prong branch not found!\n");
    return;
  }
  
 
  if(!aod->GetPrimaryVertex()||TMath::Abs(aod->GetMagneticField())<0.001) return;
  
  fCounter->StoreEvent(aod,fDplusCuts,fMontecarlo);
  fHistNEvents->Fill(0);

  
  //MC Kinemactics
  if(fMontecarlo){
    
    farrayMC =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!farrayMC){
      printf("AliAnalysisTaskSEDplusCorrelations::UserExec: MC particles branch not found!\n");
      return;
    }
    fCorrelator->SetMCArray(farrayMC);
    
    AliAODMCHeader  *mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskSEDplusCorrelations::UserExec: MC header branch not found!\n");
      return;
    }

    Double_t zVtxMC = mcHeader->GetVtxZ();
    if(TMath::Abs(zVtxMC)>10) return;
    if(aod->GetTriggerMask()==0 && (aod->GetRunNumber()>=195344 && aod->GetRunNumber()<=195677)) return;
    
    if(fMCGenEvType){
      Bool_t isMCeventgood = kFALSE;
      Int_t      eventType = mcHeader->GetEventType();
      Int_t      NMCevents = fAssoCuts->GetNofMCEventType();
      for(Int_t k=0; k<NMCevents; k++){
	Int_t * MCEventType = fAssoCuts->GetMCEventType();
	if(eventType == MCEventType[k]) isMCeventgood= kTRUE;
	// ((TH1D*)fOutputBasic->FindObject("EventTypeMC"))->Fill(eventType);
      }
      if(NMCevents && !isMCeventgood){
	if(fDebug>2) std::cout << "The MC event " << eventType << " not interesting for this analysis: skipping" << std::endl;
	return;
      }
    }
  }
  
  
  //Event Selection
  if(!fReadMC && !fDplusCuts->IsEventSelected(aod))return;
  else if(fReadMC){
    AliAODVertex *vertex = (AliAODVertex*)aod->GetPrimaryVertex();
    if(!vertex)return;
    Double_t zVtxMCreco = vertex->GetZ();
    if(TMath::Abs(zVtxMCreco)>10) return;
  }
  fHistNEvents->Fill(1);
  
    
  if(fEvalCentrality){
    fMinCentrality = fDplusCuts->GetMinCentrality(); 
    fMaxCentrality = fDplusCuts->GetMaxCentrality(); 
    fCentrOrMult = fDplusCuts->GetCentrality(aod);
  }  
  else if(!fEvalCentrality){
    Double_t count = -1, mineta = -1.0, maxeta = 1.0;
    AliAODTracklets* tracklets = aod->GetTracklets();
    if(!tracklets)return;
    Int_t nTr=tracklets->GetNumberOfTracklets();
    for(Int_t iTr=0; iTr<nTr; iTr++){
      Double_t theta=tracklets->GetTheta(iTr);
      Double_t eta=-TMath::Log(TMath::Tan(theta/2.));
      if(eta>mineta && eta<maxeta) count++;
    }
    //if(!fEvalCentrality)fCentrOrMult = (Double_t)(AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aod,mineta,maxeta));
    fMultiplicity = count;
    fCentrOrMult = fMultiplicity;
  }
  fHistNEvents->Fill(2);
  //Printf("Value of fCentrOrMult ----> %0.01f",fCentrOrMult);
  
  
  //HFCorrelator interface
  fCorrelator->SetPidAssociated();
  fCorrelator->SetAODEvent(aod);
  Bool_t correlatorON = fCorrelator->Initialize(); //
  if(!correlatorON) {
    AliInfo("AliHFCorrelator didn't initialize the pool correctly or processed a bad event");
    return;
  }
  fHistNEvents->Fill(3); // count event after correlator interface
  
  
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  Bool_t isGoodVtx=kFALSE;
  TString primTitle = vtx1->GetTitle();
  if(primTitle.Contains("VertexerTracks") && vtx1->GetNContributors()>0) {
    isGoodVtx=kTRUE;
    fHistNEvents->Fill(4);
  }       
  if(!isGoodVtx) return; 
    
  //D+ Particle loop and Correlation	  
  Int_t nDplusKpipi = array3Prong->GetEntriesFast();
  Int_t       labDp = -1;
  Bool_t    isDplus = kFALSE;    
  Int_t pdgDgDplustoKpipi[3]={321,211,211};
  
  
  //printf("Number of D+->Kpipi: %d and of tracks: %d\n",nDplusKpipi,aod->GetNumberOfTracks());    
  for (Int_t iDplusKpipi = 0; iDplusKpipi < nDplusKpipi; iDplusKpipi++){
    
    // D+ Primary Vertex Setting
    AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(iDplusKpipi); // D+
    if(!d)continue;
    fHistNEvents->Fill(5);
    
    if(d->Pt()<2.) continue; 
    
    if(fUseBit && !d->HasSelectionBit(AliRDHFCuts::kDplusCuts)){
      fHistNEvents->Fill(6);
      continue;
    }	 
    fHistNEvents->Fill(7);
    
    //Tight cuts
    Int_t passTightCuts = fDplusCuts->IsSelected(d,AliRDHFCuts::kAll,aod);
    
    //Topological cuts
    Int_t passTopologicalCuts=fDplusCuts->GetIsSelectedCuts();
    if(!passTopologicalCuts) continue;
    fHistNEvents->Fill(8);
    
    //Int_t passPIDCuts=fDplusCuts->GetIsSelectedPID();
    //if(!passPIDCuts) continue; //not needed
    
    if(fTCconfig) {if(!passTightCuts) continue;}	

    
    Bool_t unsetvtx=kFALSE;
    if(!d->GetOwnPrimaryVtx()){
      d->SetOwnPrimaryVtx(vtx1);
      unsetvtx=kTRUE;
    }
    
    
    Bool_t recVtx=kFALSE;
    AliAODVertex *origownvtx=0x0;
    if(fDplusCuts->GetIsPrimaryWithoutDaughters()){
      if(d->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*d->GetOwnPrimaryVtx());	
      if(fDplusCuts->RecalcOwnPrimaryVtx(d,aod))recVtx=kTRUE;
      else fDplusCuts->CleanOwnPrimaryVtx(d,aod,origownvtx);
    }
    fHistNEvents->Fill(9);    
    
    
    Double_t ptDplusCand = d->Pt();
    Double_t rapDplusCand=d->YDplus();    
    Int_t iPtDplus = fDplusCuts->PtBin(ptDplusCand);
    if(iPtDplus<0)continue;
    
    // D+ Fiducial Acceptance in pt and eta
    if(!fDplusCuts->IsInFiducialAcceptance(ptDplusCand,rapDplusCand)) continue;
    fHistNEvents->Fill(10); 
    
    if(fMontecarlo){
      labDp = d->MatchToMC(411,farrayMC,3,pdgDgDplustoKpipi);
      if(labDp>=0){
	  isDplus = kTRUE;
      }
    }
    
    //Computing Correlations: Below Funtions 
    HadronCorrelations(d,farrayMC,isDplus);
    
    //Vertexing cleaning... 
    if(recVtx)fDplusCuts->CleanOwnPrimaryVtx(d,aod,origownvtx);	 
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
    
  }// loop over D+ 
  
  
  //EVT MIXING PART
  if(fMixing ) {
    Bool_t updated = fCorrelator->PoolUpdate();
    if(!updated) AliInfo("Pool was not updated");
  }

  
  // Posting Slots
  PostData(1,fOutput);
  PostData(2,fOutputCorr);
  PostData(5,fCounter); 
  return;
  
} 



//____________________| Calculating Correlation
void AliAnalysisTaskSEDplusCorrelations::HadronCorrelations(AliAODRecoDecayHF3Prong* d,  TClonesArray *arrayMC, Bool_t isDplus) {
  
  //Switch for String for Hist name 
  TString parttype = "", datatype = "";
  if(fSelect==1)parttype += "Hdron";
  else if(fSelect==2)parttype += "Kaons";
  else if(fSelect==3)parttype += "kZero";
  
  if(fReco){
    if(!fReadMC)datatype += "Data";
    else if(fReadMC)datatype += "Reco";
  }

  
  // D+ Quantities
  Double_t mDplus       = d->InvMassDplus();
  Double_t ptDplusCand  = d->Pt();
  Double_t etaDplusCand = d->YDplus();
  Double_t phiDplusCand = d->Phi();

  
  Double_t effDplus =1.0;
  if(fEffDplus)effDplus = fAssoCuts->GetTrigWeight(ptDplusCand,fMultiplicity); //Dplus efficiency
  if(effDplus<1.0e-9) effDplus=1; // case of 0 bin content
  
  Int_t iPtBin = fDplusCuts->PtBin(ptDplusCand); // Pt bins
  if(iPtBin<0) return ;

  
  //Correlation with SE or ME
  Double_t nSparceCorrDplusinfo[1] = {mDplus};
  
  if(!fMixing){
    if(fReco){ // data/MC-Reco
      ((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseM_%s_%s_Bin_%d", parttype.Data(), datatype.Data(),iPtBin)))->Fill(nSparceCorrDplusinfo,1.0/effDplus);
    }
    if(fMontecarlo){
      ((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseM_%s_%s_Bin_%d", parttype.Data(), datatype.Data(),iPtBin)))->Fill(nSparceCorrDplusinfo);
    }
  }
  
  
  //Dplus parameter resetting and storage
  Double_t phiDplus = fCorrelator->SetCorrectPhiRange(phiDplusCand);
  fCorrelator->SetTriggerParticleProperties(ptDplusCand,phiDplus,etaDplusCand);   
  
  //Pool Setting and Pool events
  Bool_t execPool = fCorrelator->ProcessEventPool();
  if(fMixing && !execPool) {
    AliInfo("Mixed event analysis: pool is not ready");
  }
  
  Int_t NofEventsinPool = 1; // SE
  if(fMixing) {
    NofEventsinPool = fCorrelator->GetNofEventsInPool();
  }
    
  //D+ Daughters ids
  Int_t nIDs[3] = {-9999999};
  nIDs[0] = ((AliAODTrack*)d->GetDaughter(0))->GetID();
  nIDs[1] = ((AliAODTrack*)d->GetDaughter(1))->GetID();
  nIDs[2] = ((AliAODTrack*)d->GetDaughter(2))->GetID();
    
  //cout << "-----" <<phiDplusCand << ",  "<<ptDplusCand << endl;
  //Correlation with associated tracks 
  for (Int_t jMix =0; jMix < NofEventsinPool; jMix++){
    
    Bool_t analyzetracks = fCorrelator->ProcessAssociatedTracks(jMix);
    if(!analyzetracks) {
      AliInfo("AliHFCorrelator::Cannot process the track array");
      continue;
    }
    
    // Leading Particle varibales
    Double_t refpt = 0, effLead = 1, TotaleffLead =1;
    Double_t ptleadHadron=0;
    Double_t DeltaphiLead=0, DeltaetaLead=0;
    
    //Correlation with associated tracks 
    for (Int_t iTrack = 0; iTrack<fCorrelator->GetNofTracks(); iTrack++) {
      
      // Correlation using HFCorrelator
      Bool_t runcorrelation = fCorrelator->Correlate(iTrack);
      if(!runcorrelation) continue;
      
      AliReducedParticle* redpart = fCorrelator->GetAssociatedParticle();
      if(!redpart)continue;
      
      Double_t phiHad=redpart->Phi();
      Double_t ptHad=redpart->Pt();
      Int_t label = redpart->GetLabel();
      Int_t trackid = redpart->GetID();

      phiHad = fCorrelator->SetCorrectPhiRange(phiHad);
      
      Double_t effTr =1;
      if(fEffTrack)effTr = redpart->GetWeight(); //track efficiency
      if(effTr<1.e-9) effTr=1;    
      Double_t effTotal = effTr*effDplus;
      
      if(trackid < 0) continue;      
      if (!fMixing)if( trackid == nIDs[0] || trackid == nIDs[1] || trackid == nIDs[2]) continue;
      
      // leading particle correlation
      if (ptHad > refpt) {
	refpt = ptHad;
	DeltaphiLead  = fCorrelator->GetDeltaPhi();
	DeltaetaLead  = fCorrelator->GetDeltaEta();
	ptleadHadron  = ptHad;
	effLead = redpart->GetWeight();
	TotaleffLead = effLead*effDplus;
      }
      
      //filling correlations in ThnSparce format
      if(fReco){
	Bool_t* partSource = NULL;
	CorrelationNSparsePlots(d, redpart, iPtBin, partSource, 1.0/effTotal);
	delete[] partSource;
      }
      
      if(fMontecarlo && isDplus){
	Bool_t* partSource = fAssoCuts->IsMCpartFromHF(label,arrayMC); // check source for 1/2/3
	CorrelationNSparsePlots(d, redpart, iPtBin, partSource, 1.0/effTotal);
	delete[] partSource;
      } // MC Gen Origin ends
    }//track loop end here
    
    //Leading particle correlations
    Double_t nSparceCorrLeadPart[4] = {mDplus, DeltaphiLead, DeltaetaLead, ptleadHadron};
    if(!fMixing) {	
      if(fReco)((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseL_%s_%s_Bin_%d", parttype.Data(), datatype.Data(),iPtBin)))->Fill(nSparceCorrLeadPart);
      if(fMontecarlo) ((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseL_%s_%s_Bin_%d", parttype.Data(), datatype.Data(),iPtBin)))->Fill(nSparceCorrLeadPart); 
    }
    
    if(fMixing) {
      ((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseL_%s_%s_Bin_%d_evMix", parttype.Data(), datatype.Data(),iPtBin)))->Fill(nSparceCorrLeadPart);
    }   
  }  // event loop.  
  
  
}


//____________________|Correlation with ThnSparse Histograms
void AliAnalysisTaskSEDplusCorrelations::CorrelationNSparsePlots(AliAODRecoDecayHF3Prong *d, AliReducedParticle* track, Int_t iPtBin, Bool_t *origDplus, Double_t weightEff) {
  
  iPtBin = fDplusCuts->PtBin(d->Pt());
  Double_t mDplus = d->InvMassDplus();
  Double_t deltaPhi = fCorrelator->GetDeltaPhi();
  Double_t deltaEta = fCorrelator->GetDeltaEta();
  Double_t ptTrack = track->Pt(); // replacing by multiplicity for the time being
  
  TString partype = "", datatype = "";
  if(fSelect==1)partype += "Hdron";
  else if(fSelect==2)partype += "Kaons";
  else if(fSelect==3)partype += "kZero";
  
  
  if(fReco){
    if(!fReadMC)datatype += "Data";
    else if(fReadMC)datatype += "Reco";
  }
  
  if(fMontecarlo){
    origDplus = new Bool_t[4];
    if(origDplus[0] && origDplus[2]) datatype += "Frmc"; // is from charm ->D
    else if(origDplus[1] && origDplus[2])datatype += "Frmb"; // is from beauty ->D
    else if(origDplus[1] && origDplus[3])datatype += "FrmB"; // is from beauty ->B
    else   datatype += "Trth";
  }

  
  //pT bins setting
  Double_t ptLim_Sparse = 0.0;
  if(!fMixing) {        
    if(fReco)ptLim_Sparse = ((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseC_%s_%s_Bin_%d",partype.Data(), datatype.Data(),iPtBin)))->GetAxis(3)->GetXmax();
    if(fMontecarlo)ptLim_Sparse = ((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseC_%s_%s_Bin_%d",partype.Data(), datatype.Data(),iPtBin)))->GetAxis(3)->GetXmax();
  }
  else if(fMixing) {
    ptLim_Sparse=((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseC_%s_%s_Bin_%d_evMix",partype.Data(), datatype.Data(), iPtBin )))->GetAxis(3)->GetXmax();
  }
  if(ptTrack > ptLim_Sparse) ptTrack = ptLim_Sparse-0.01; //filling all above pT in last bin 

  

  //Correlations
  Double_t nSparceCorr[4] = {mDplus, deltaPhi, deltaEta, ptTrack};
  if(!fMixing) {
    if(fReco)((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseC_%s_%s_Bin_%d",partype.Data(), datatype.Data(),iPtBin)))->Fill(nSparceCorr,weightEff);
    if(fMontecarlo)((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseC_%s_%s_Bin_%d",partype.Data(), datatype.Data(),iPtBin)))->Fill(nSparceCorr,weightEff);
  }
  else if(fMixing) {
    ((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseC_%s_%s_Bin_%d_evMix",partype.Data(), datatype.Data(), iPtBin )))->Fill(nSparceCorr);
  }              
  
  delete[] origDplus; 
  return;	        
  
}


//____________________|  Histograms Nomenclature
void AliAnalysisTaskSEDplusCorrelations::HistoNomenclature() {
  
  Double_t Pi = TMath::Pi();
  
  //D+ Candidate Vars: InvMass
  Float_t         range  = 0.200;
  Float_t fLowmasslimit  = 1.865 - range;
  Float_t  fUpmasslimit  = 1.865 + range;
  Float_t         width  = fBinWidth;       
  Int_t           nbins  = (Int_t)((fUpmasslimit-fLowmasslimit)/width+0.5);
  Int_t     missingbins  = 4 - nbins%4;
  nbins  =  nbins + missingbins; 
  width  = (fUpmasslimit-fLowmasslimit)/nbins; // new width
  fBinWidth=width;
  
  
  //1. Invariant Mass 
  Int_t     nBinsDinfo[1] = {nbins};
  Double_t binMinDinfo[1] = {fLowmasslimit};
  Double_t binMaxDinfo[1] = {fUpmasslimit};
  
  //2. SEorME Correlations Particle Vars: InvMass, DeltaPhi, DeltaEta, Centrality
  Int_t           nVarsBins[4] = {nbins,              32,      16,    6   };
  Double_t nMinimumEdgeBins[4] = {fLowmasslimit,   -Pi/2,    -1.6,    0.00};  //is the minimum for all the bins
  Double_t nMaximumEdgeBins[4] = {fUpmasslimit,   3*Pi/2,     1.6,    3.00};  //is the maximum for all the bins
  
  
  TString nameMasterStting   = "";
  TString namePlotThnDplus   = "hnSparse" ;
  TString namePlotThnDplusMC = "hnSparse" ;
  TString namePlotThn        = "hnSparse" ;
  TString namePlotThnL       = "hnSparse" ;
  TString namePlotThnMC      = "hnSparse" ;
  TString namePlotThnLMC     = "hnSparse" ;
  TString namePlotThnMCc     = "hnSparse" ;
  TString namePlotThnMCb     = "hnSparse" ;
  TString namePlotThnMCB     = "hnSparse" ;
  
  if(fSelect==1)nameMasterStting += "_Hdron";
  else  if(fSelect==2)nameMasterStting += "_Kaons";
  else  if(fSelect==3)nameMasterStting += "_kZero";
  else  nameMasterStting = "_Nulls";
  
  if(fReco){
    if(!fReadMC)nameMasterStting +="_Data";
    if(fReadMC)nameMasterStting +="_Reco";
    namePlotThnDplus += Form("M%s_Bin_", nameMasterStting.Data());
    namePlotThn      += Form("C%s_Bin_", nameMasterStting.Data());
    namePlotThnL     += Form("L%s_Bin_", nameMasterStting.Data());
  }

  if(fMontecarlo){  
    namePlotThnDplusMC  +=  Form("M%s_Trth_Bin_", nameMasterStting.Data());
    namePlotThnMCc      +=  Form("C%s_Frmc_Bin_", nameMasterStting.Data());
    namePlotThnMCb      +=  Form("C%s_Frmb_Bin_", nameMasterStting.Data());
    namePlotThnMCB      +=  Form("C%s_FrmB_Bin_", nameMasterStting.Data());
    namePlotThnLMC      +=  Form("L%s_Trth_Bin_", nameMasterStting.Data());
  }
  
  
  for(Int_t i=0;  i<fNPtBins;  i++){
    
    namePlotThnDplus.Resize(25);
    namePlotThnDplusMC.Resize(25);
    namePlotThn.Resize(25);
    namePlotThnL.Resize(25);
    namePlotThnMC.Resize(25);
    namePlotThnLMC.Resize(25);
    namePlotThnMCc.Resize(25);
    namePlotThnMCb.Resize(25);
    namePlotThnMCB.Resize(25);
    
    namePlotThnDplus+=i;
    namePlotThnDplusMC+=i;
    namePlotThn+=i;
    namePlotThnL+=i;
    namePlotThnMC+=i;
    namePlotThnLMC+=i;
    namePlotThnMCc+=i;
    namePlotThnMCb+=i;
    namePlotThnMCB+=i;
    
    
    // Correlations in ThnSparce
    if(!fMixing) {
      if(fReco){
	// Histograms for D+ Mass
	THnSparseI *DplusHisto = new THnSparseI(namePlotThnDplus.Data(), "Reco_Dplus_Mass_Histo; D^{+} Inv Mass;",1,nBinsDinfo,binMinDinfo,binMaxDinfo);
	DplusHisto->Sumw2();
	fOutputCorr->Add(DplusHisto);
        
	//D-h
	THnSparseI *hPhiCorr = new THnSparseI(namePlotThn.Data(), "Reco_SE_Corr_ThnSparse; D^{+} Inv Mass; #Delta#phi; #Delta#eta; track p_T;",4,nVarsBins,nMinimumEdgeBins,nMaximumEdgeBins);
	hPhiCorr->Sumw2();
	fOutputCorr->Add(hPhiCorr);
        
	//D-Leading Particles
	THnSparseI *hCorrLead = new THnSparseI(namePlotThnL.Data(), "Reco_SE_LeadingCorr_ThnSparse; D^{+} Inv Mass; #Delta#phi; #Delta#eta; track p_T;", 4, nVarsBins,nMinimumEdgeBins,nMaximumEdgeBins);
	hCorrLead->Sumw2();
	fOutputCorr->Add(hCorrLead);
      }
      
      if(fMontecarlo){
	// Histograms for D+ Candidates
	THnSparseI *DplusHistoMC = new THnSparseI(namePlotThnDplusMC.Data(), "MCGen_Dplus_Mass_Histo; D^{+} Inv Mass;",1,nBinsDinfo,binMinDinfo,binMaxDinfo);
	DplusHistoMC->Sumw2();
	fOutputCorr->Add(DplusHistoMC);
        
	THnSparseI *hPhiCorrMCc = new THnSparseI(namePlotThnMCc.Data(), "MCGen_SE_CorrCharm_ThnSparse; D^{+} Inv Mass; #Delta#phi; #Delta#eta; track p_T;",4,nVarsBins,nMinimumEdgeBins,nMaximumEdgeBins);
	hPhiCorrMCc->Sumw2();
	fOutputCorr->Add(hPhiCorrMCc);
        
	THnSparseI *hPhiCorrMCb = new THnSparseI(namePlotThnMCb.Data(), "MCGen_SE_CorrBeauty_ThnSparse; D^{+} Inv Mass; #Delta#phi; #Delta#eta; track p_T;",4,nVarsBins,nMinimumEdgeBins,nMaximumEdgeBins);
	hPhiCorrMCb->Sumw2();
	fOutputCorr->Add(hPhiCorrMCb);
        
	THnSparseI *hPhiCorrMCB = new THnSparseI(namePlotThnMCB.Data(), "MCGen_SE_CorrBMeson_ThnSparse; D^{+} Inv Mass; #Delta#phi; #Delta#eta; track p_T;",4,nVarsBins,nMinimumEdgeBins,nMaximumEdgeBins);
	hPhiCorrMCB->Sumw2();
	fOutputCorr->Add(hPhiCorrMCB);
        
	THnSparseI *hCorrLeadMC = new THnSparseI(namePlotThnLMC.Data(), "MCGen_SE_LeadingCorr_ThnSparse; D^{+} Inv Mass; #Delta#phi; #Delta#eta; track p_T;",4,nVarsBins,nMinimumEdgeBins,nMaximumEdgeBins);
	hCorrLeadMC->Sumw2();
	fOutputCorr->Add(hCorrLeadMC);
      }
    }

    
    
    if(fMixing) {  // Histogram for Event Mixing	  
      if(fReco){
	namePlotThn+="_evMix";
	THnSparseI *hPhiCorrMixing = new THnSparseI(namePlotThn.Data(), "Reco_SE_Corr_ThnSparse_EvMix; D^{+} Inv. Mass; #Delta#phi; #Delta#eta; track p_T;",4,nVarsBins,nMinimumEdgeBins,nMaximumEdgeBins);
	hPhiCorrMixing->Sumw2();
	fOutputCorr->Add(hPhiCorrMixing);
        
	// Correlations 2D for Leading Particles
	namePlotThnL+="_evMix";
	THnSparseI *hCorrLead = new THnSparseI(namePlotThnL.Data(), "Reco_SE_LeadingCorr_ThnSparse_EvMix; D^{+} Inv Mass; #Delta#phi; #Delta#eta; track p_T;",4,nVarsBins,nMinimumEdgeBins,nMaximumEdgeBins);
	hCorrLead->Sumw2();
	fOutputCorr->Add(hCorrLead);	
      }
    }
    
   
  } 
}


//____________________| Terminate 
void AliAnalysisTaskSEDplusCorrelations::Terminate(Option_t */*option*/) {
  
  fOutputCorr = dynamic_cast<TList*> (GetOutputData(2));
  if (!fOutputCorr) {
    printf("ERROR: fOutputCorr not available\n");
    return;
  }  
  
  fDplusCuts = dynamic_cast<AliRDHFCutsDplustoKpipi*>(GetOutputData(3));
  if(!fDplusCuts){
    printf("ERROR: fDplusCuts(Dplus Track Cuts) are not available\n");
    return;
  }
  
  fAssoCuts = dynamic_cast<AliHFAssociatedTrackCuts*>(GetOutputData(4));
  if(!fAssoCuts){
    printf("ERROR: fAssoCuts(Associated Track Cuts) are not available\n");
    return;
  }
  
  fCounter = dynamic_cast<AliNormalizationCounter*>(GetOutputData(5));
  if (!fCounter) {
    printf("ERROR: fCounter not available\n");
    return;
  }	  
  
  return;	  
}


//EOF Jitendra
