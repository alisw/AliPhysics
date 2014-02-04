/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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
//
//
//             Base class for DStar - Hadron Correlations Analysis
//
//-----------------------------------------------------------------------
//          
//
//						   Author S.Bjelogrlic
//                         Utrecht University 
//                      sandro.bjelogrlic@cern.ch
//
//-----------------------------------------------------------------------

/* $Id$ */

//#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TVector3.h>
#include <TChain.h>
#include "TROOT.h"

#include "AliAnalysisTaskDStarCorrelations.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliHFAssociatedTrackCuts.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODPidHF.h"
#include "AliVParticle.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliAODHandler.h"
#include "AliESDtrack.h"
#include "AliAODMCParticle.h"
#include "AliNormalizationCounter.h"
#include "AliReducedParticle.h"
#include "AliHFCorrelator.h"
#include "AliAODMCHeader.h"
#include "AliEventPoolManager.h"
#include "AliVertexingHFUtils.h"

using std::cout;
using std::endl;


ClassImp(AliAnalysisTaskDStarCorrelations)


//__________________________________________________________________________
AliAnalysisTaskDStarCorrelations::AliAnalysisTaskDStarCorrelations() :
AliAnalysisTaskSE(),
fhandler(0x0),
fmcArray(0x0),
fCounter(0x0),
fCorrelator(0x0),
fselect(0),
fmontecarlo(kFALSE),
fmixing(kFALSE),
fFullmode(kFALSE),
fSystem(pp),
fEfficiencyVariable(kNone),
fReco(kTRUE),
fUseEfficiencyCorrection(kFALSE),
fUseDmesonEfficiencyCorrection(kFALSE),
fUseCentrality(kFALSE),
fUseHadronicChannelAtKineLevel(kFALSE),
fPhiBins(32),
fEvents(0),
fDebugLevel(0),
fDisplacement(0),
fDim(0),
fNofPtBins(0),
fMaxEtaDStar(0.9),
fDMesonSigmas(0),

fD0Window(0),

fOutput(0x0),
fOutputMC(0x0),
fCuts(0),
fCuts2(0),
fUtils(0),
fTracklets(0),
fDeffMapvsPt(0),
fDeffMapvsPtvsMult(0),
fDeffMapvsPtvsEta(0)
{
    SetDim();
    // default constructor
}

//__________________________________________________________________________
AliAnalysisTaskDStarCorrelations::AliAnalysisTaskDStarCorrelations(const Char_t* name,AliRDHFCutsDStartoKpipi* cuts, AliHFAssociatedTrackCuts *AsscCuts,AliAnalysisTaskDStarCorrelations::CollSyst syst,Bool_t mode) :
AliAnalysisTaskSE(name),

fhandler(0x0),
fmcArray(0x0),
fCounter(0x0),
fCorrelator(0x0),
fselect(0),
fmontecarlo(kFALSE),
fmixing(kFALSE),
fFullmode(mode),
fSystem(syst),
fEfficiencyVariable(kNone),
fReco(kTRUE),
fUseEfficiencyCorrection(kFALSE),
fUseDmesonEfficiencyCorrection(kFALSE),
fUseCentrality(kFALSE),
fUseHadronicChannelAtKineLevel(kFALSE),
fPhiBins(32),
fEvents(0),
fDebugLevel(0),
fDisplacement(0),
fDim(0),
fNofPtBins(0),
fMaxEtaDStar(0.9),
fDMesonSigmas(0),
fD0Window(0),

fOutput(0x0),
fOutputMC(0x0),
fCuts(0),
fCuts2(AsscCuts),
fUtils(0),
fTracklets(0),
fDeffMapvsPt(0),
fDeffMapvsPtvsMult(0),
fDeffMapvsPtvsEta(0)
{
     Info("AliAnalysisTaskDStarCorrelations","Calling Constructor");
  
    SetDim();
    if(fSystem == AA)  fUseCentrality = kTRUE; else fUseCentrality = kFALSE;
  
    fCuts=cuts;
    fNofPtBins= fCuts->GetNPtBins();
    //cout << "Enlarging the DZero window " << endl;
    EnlargeDZeroMassWindow();
   // cout << "Done" << endl;
    
   
  DefineInput(0, TChain::Class());
  
  DefineOutput(1,TList::Class()); // histos from data and MC
  DefineOutput(2,TList::Class()); // histos from MC
  DefineOutput(3,AliRDHFCutsDStartoKpipi::Class()); // my D meson cuts
  DefineOutput(4,AliHFAssociatedTrackCuts::Class()); // my associated tracks cuts
  DefineOutput(5,AliNormalizationCounter::Class());   // normalization
}

//__________________________________________________________________________

AliAnalysisTaskDStarCorrelations::~AliAnalysisTaskDStarCorrelations() {
  //
	// destructor
	//
	
	Info("AliAnalysisTaskDStarCorrelations","Calling Destructor");  
	
	if(fhandler) {delete fhandler; fhandler = 0;}
	//if(fPoolMgr) {delete fPoolMgr; fPoolMgr = 0;}    
	if(fmcArray) {delete fmcArray; fmcArray = 0;}
	if(fCounter) {delete fCounter; fCounter = 0;}
	if(fCorrelator) {delete fCorrelator; fCorrelator = 0;}
	if(fOutput) {delete fOutput; fOutput = 0;}
	if(fOutputMC) {delete fOutputMC; fOutputMC = 0;}
	if(fCuts) {delete fCuts; fCuts = 0;}
	if(fCuts2) {delete fCuts2; fCuts2=0;}
    if(fDeffMapvsPt){delete fDeffMapvsPt; fDeffMapvsPt=0;}
    if(fDeffMapvsPtvsMult){delete fDeffMapvsPtvsMult; fDeffMapvsPtvsMult=0;}
    if(fDeffMapvsPtvsEta){delete fDeffMapvsPtvsEta; fDeffMapvsPtvsEta=0;}

}

//___________________________________________________________
void AliAnalysisTaskDStarCorrelations::Init(){
	//
	// Initialization
	//
	if(fDebugLevel > 1) printf("AliAnalysisTaskDStarCorrelations::Init() \n");
	
	AliRDHFCutsDStartoKpipi* copyfCuts=new AliRDHFCutsDStartoKpipi(*fCuts);
	
	
    
	// Post the D* cuts
	PostData(3,copyfCuts);
	
	// Post the hadron cuts
	PostData(4,fCuts2);
    
	return;
}


//_________________________________________________
void AliAnalysisTaskDStarCorrelations::UserCreateOutputObjects(){
	Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());
	
	//slot #1
	//OpenFile(0);
	fOutput = new TList();
	fOutput->SetOwner();
	
	fOutputMC = new TList();
	fOutputMC->SetOwner();
	
	// define histograms
	DefineHistoForAnalysis();
    DefineThNSparseForAnalysis();
    

	fCounter = new AliNormalizationCounter(Form("%s",GetOutputSlot(5)->GetContainer()->GetName()));
	fCounter->Init();
	
    Double_t Pi = TMath::Pi();
	fCorrelator = new AliHFCorrelator("Correlator",fCuts2,fUseCentrality); // fCuts2 is the hadron cut object, fSystem to switch between pp or PbPb
	fCorrelator->SetDeltaPhiInterval(  -0.5*Pi, 1.5*Pi); // set correct phi interval
	//fCorrelator->SetDeltaPhiInterval((-0.5)*Pi,(1.5)*Pi); // set correct phi interval
	fCorrelator->SetEventMixing(fmixing); //set kFALSE/kTRUE for mixing Off/On
	fCorrelator->SetAssociatedParticleType(fselect); // set 1/2/3 for hadron/kaons/kzeros
	fCorrelator->SetApplyDisplacementCut(fDisplacement); //set kFALSE/kTRUE for using the displacement cut
	fCorrelator->SetUseMC(fmontecarlo);
	fCorrelator->SetUseReco(fReco);
    //	fCorrelator->SetKinkRemoval(kTRUE);
	Bool_t pooldef = fCorrelator->DefineEventPool();
	
	if(!pooldef) AliInfo("Warning:: Event pool not defined properly");
    
    fUtils = new AliAnalysisUtils();
    
    
	
	PostData(1,fOutput); // set the outputs
	PostData(2,fOutputMC); // set the outputs
	PostData(5,fCounter); // set the outputs
}

//_________________________________________________
void AliAnalysisTaskDStarCorrelations::UserExec(Option_t *){
  
  
  if(fDebugLevel){
    
    if(fReco) std::cout << "USING RECONSTRUCTION" << std::endl;
    if(!fReco) std::cout << "USING MC TRUTH" << std::endl;
    std::cout << " " << std::endl;
    std::cout << "=================================================================================" << std::endl;
    if(!fmixing){
      if(fselect==1) std::cout << "TASK::Correlation with hadrons on SE "<< std::endl;
      if(fselect==2) std::cout << "TASK::Correlation with kaons on SE "<< std::endl;
      if(fselect==3) std::cout << "TASK::Correlation with kzeros on SE "<< std::endl;
    }
    if(fmixing){
      if(fselect==1) std::cout << "TASK::Correlation with hadrons on ME "<< std::endl;
      if(fselect==2) std::cout << "TASK::Correlation with kaons on ME "<< std::endl;
      if(fselect==3) std::cout << "TASK::Correlation with kzeros on ME "<< std::endl;
    }
    
  }// end if debug
  
  
  
  if (!fInputEvent) {
    Error("UserExec","NO EVENT FOUND!");
		return;
  }
  
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
  if(!aodEvent){
    AliError("AOD event not found!");
    return;
  }
  
    fTracklets = aodEvent->GetTracklets();
    
    fEvents++; // event counter
    ((TH1D*)fOutput->FindObject("NofEvents"))->Fill(0);
    
    fCounter->StoreEvent(aodEvent,fCuts,fmontecarlo);
    
    // load MC array
    fmcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if(fmontecarlo && !fmcArray){
      AliError("Array of MC particles not found");
      return;
    }
    
    
    
    
    // ********************************************** START EVENT SELECTION ****************************************************
    
    Bool_t isEvSel=fCuts->IsEventSelected(aodEvent);
    
    if(!isEvSel) {
      
      if(fCuts->IsEventRejectedDueToPileup()) ((TH1D*)fOutput->FindObject("NofEvents"))->Fill(2);
      if(fCuts->IsEventRejectedDueToCentrality()) ((TH1D*)fOutput->FindObject("NofEvents"))->Fill(3);
      if(fCuts->IsEventRejectedDueToNotRecoVertex()) ((TH1D*)fOutput->FindObject("NofEvents"))->Fill(4);
      if(fCuts->IsEventRejectedDueToVertexContributors()) ((TH1D*)fOutput->FindObject("NofEvents"))->Fill(5);
      if(fCuts->IsEventRejectedDueToZVertexOutsideFiducialRegion()) ((TH1D*)fOutput->FindObject("NofEvents"))->Fill(6);
      if(fCuts->IsEventRejectedDueToTrigger()) ((TH1D*)fOutput->FindObject("NofEvents"))->Fill(7);
      if(fCuts->IsEventRejectedDuePhysicsSelection()) ((TH1D*)fOutput->FindObject("NofEvents"))->Fill(8);
      
      return;
    }
    
    // added event selection for pA
    
    if(fSystem == pA){
      
      if(fUtils->IsFirstEventInChunk(aodEvent)) {
	AliInfo("Rejecting the event - first in the chunk");
	((TH1D*)fOutput->FindObject("NofEvents"))->Fill(9);
	return;
        }
      if(!fUtils->IsVertexSelected2013pA(aodEvent)) {
	((TH1D*)fOutput->FindObject("NofEvents"))->Fill(10);
	AliInfo("Rejecting the event - bad vertex");
	return;
      }
    }
    // ******************************** END event selections **************************************************
    
    AliCentrality *centralityObj = 0;
    Double_t MultipOrCent = -1;
    
    if(fUseCentrality){
      /* if(fSystem == AA ){	*/	centralityObj = aodEvent->GetHeader()->GetCentralityP();
      MultipOrCent = centralityObj->GetCentralityPercentileUnchecked("V0M");
      //AliInfo(Form("Centrality is %f", MultipOrCent));
    }
    
    if(!fUseCentrality) MultipOrCent = AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aodEvent,-1.,1.);
    
        
    fCorrelator->SetAODEvent(aodEvent); // set the event to be processed
    
    ((TH1D*)fOutput->FindObject("NofEvents"))->Fill(1);
    
    Bool_t correlatorON = fCorrelator->Initialize(); //define the pool for mixing and check if event is in pool settings
	if(!correlatorON) {
	  AliInfo("AliHFCorrelator didn't initialize the pool correctly or processed a bad event");
	  return;
	}
	
	if(fmontecarlo) fCorrelator->SetMCArray(fmcArray);
		
	// check the event type
	// load MC header
	if(fmontecarlo){
	  AliAODMCHeader *mcHeader = dynamic_cast<AliAODMCHeader*>(aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
	  if (fmontecarlo && !mcHeader) {
	    AliError("Could not find MC Header in AOD");
	    return;
	  }
	  
	  Bool_t isMCeventgood = kFALSE;
	  
	  
	  Int_t eventType = mcHeader->GetEventType();
	  Int_t NMCevents = fCuts2->GetNofMCEventType();
	  
	  for(Int_t k=0; k<NMCevents; k++){
	    Int_t * MCEventType = fCuts2->GetMCEventType();
	    
	    if(eventType == MCEventType[k]) isMCeventgood= kTRUE;
	    ((TH1D*)fOutputMC->FindObject("EventTypeMC"))->Fill(eventType);
	  }
	  
	  if(NMCevents && !isMCeventgood){
            if(fDebugLevel)	std::cout << "The MC event " << eventType << " not interesting for this analysis: skipping" << std::endl;
	    return;
	  }
	  
	} // end if montecarlo
	
	
	// checks on vertex and multiplicity of the event
	AliAODVertex *vtx = aodEvent->GetPrimaryVertex();
	Double_t zVtxPosition = vtx->GetZ(); // zvertex
	
	
	if(fFullmode) ((TH2F*)fOutput->FindObject("EventPropsCheck"))->Fill(MultipOrCent,zVtxPosition);
	
	
	
	// D* reconstruction
	TClonesArray *arrayDStartoD0pi=0;
	if(!aodEvent && AODEvent() && IsStandardAOD()) {
	  // In case there is an AOD handler writing a standard AOD, use the AOD
	  // event in memory rather than the input (ESD) event.
	  aodEvent = dynamic_cast<AliAODEvent*> (AODEvent());
	  // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
	  // have to taken from the AOD event hold by the AliAODExtension
	  AliAODHandler* aodHandler = (AliAODHandler*)
	    ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
	  if(aodHandler->GetExtensions()) {
	    AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
	    AliAODEvent *aodFromExt = ext->GetAOD();
	    arrayDStartoD0pi=(TClonesArray*)aodFromExt->GetList()->FindObject("Dstar");
	  }
	} else {
	  arrayDStartoD0pi=(TClonesArray*)aodEvent->GetList()->FindObject("Dstar");
	}
	
	if(!aodEvent->GetPrimaryVertex() || TMath::Abs(aodEvent->GetMagneticField())<0.001) return;
	
    
    // get the poolbin
    
    Int_t poolbin = fCuts2->GetPoolBin(MultipOrCent, zVtxPosition);
  
    
	// initialize variables you will need for the D*
	
	Double_t ptDStar;//
	Double_t phiDStar;//
	Double_t etaDStar;//
	Bool_t isInPeak, isInDZeroSideBand, isInDStarSideBand, isDStarMCtag;
	Double_t invMassDZero;
	Double_t deltainvMDStar;
    
	
	Double_t mPDGD0=1.8648;//TDatabasePDG::Instance()->GetParticle(421)->Mass();
	Double_t mPDGDstar=2.01022;//TDatabasePDG::Instance()->GetParticle(413)->Mass();
	
	
	//MC tagging for DStar
	//D* and D0 prongs needed to MatchToMC method
	Int_t pdgDgDStartoD0pi[2]={421,211};
	Int_t pdgDgD0toKpi[2]={321,211};
	
	Bool_t isDStarCand = kFALSE;
    Bool_t isDfromB = kFALSE;
	Bool_t isEventMixingFilledPeak = kFALSE;
	Bool_t isEventMixingFilledSB = kFALSE;
    Bool_t EventHasDStarCandidate = kFALSE;
    Bool_t EventHasDZeroSideBandCandidate = kFALSE;
    Bool_t EventHasDStarSideBandCandidate = kFALSE;
	//loop on D* candidates
	
	Int_t looponDCands = 0;
	if(fReco) looponDCands = arrayDStartoD0pi->GetEntriesFast();
	if(!fReco) looponDCands = fmcArray->GetEntriesFast();
	
	Int_t nOfDStarCandidates = 0;
	Int_t nOfSBCandidates = 0;
	
	Double_t DmesonEfficiency = 1.;
	Double_t DmesonWeight = 1.;
    Double_t efficiencyvariable = -999;
    
	
	
	for (Int_t iDStartoD0pi = 0; iDStartoD0pi<looponDCands; iDStartoD0pi++) {
	  isInPeak = kFALSE;
	  isInDStarSideBand = kFALSE;
	  isInDZeroSideBand = kFALSE;
	  isDStarMCtag = kFALSE;
	  isDfromB = kFALSE;
	  ptDStar = -123.4;
	  phiDStar = -999;
	  etaDStar = -56.;
	  invMassDZero = - 999;
	  deltainvMDStar = -998;
	  AliAODRecoCascadeHF* dstarD0pi;
	  AliAODRecoDecayHF2Prong* theD0particle;
	  AliAODMCParticle* DStarMC=0x0;
      Short_t daughtercharge = -2;
	  Int_t trackiddaugh0 = -1; // track id if it is reconstruction - label if it is montecarlo info
	  Int_t trackiddaugh1 = -1;
	  Int_t trackidsoftPi = -1;
	  
	  // start the if reconstructed candidates from here ************************************************
	  
	  if(fReco){//// if reconstruction is applied
	    dstarD0pi = (AliAODRecoCascadeHF*)arrayDStartoD0pi->At(iDStartoD0pi);
	    if(!dstarD0pi->GetSecondaryVtx()) continue;
	    theD0particle = (AliAODRecoDecayHF2Prong*)dstarD0pi->Get2Prong();
	    if (!theD0particle) continue;
            
	    
	    // track quality cuts
	    Int_t isTkSelected = fCuts->IsSelected(dstarD0pi,AliRDHFCuts::kTracks); // quality cuts on tracks
	    // region of interest + topological cuts + PID
	    Int_t isSelected=fCuts->IsSelected(dstarD0pi,AliRDHFCuts::kCandidate); //selected
            
            //apply track selections
	    if(!isTkSelected) continue;
	    if(!isSelected) continue;
        if(!fCuts->IsInFiducialAcceptance(dstarD0pi->Pt(),dstarD0pi->YDstar())) continue;
          
          
          ptDStar = dstarD0pi->Pt();
          phiDStar = dstarD0pi->Phi();
          etaDStar = dstarD0pi->Eta();
          if(TMath::Abs(etaDStar) > fMaxEtaDStar) continue;
          if(fEfficiencyVariable == kMult || fEfficiencyVariable == kCentr)  efficiencyvariable = MultipOrCent;
          if(fEfficiencyVariable == kEta) efficiencyvariable = etaDStar;
          if(fEfficiencyVariable == kRapidity) efficiencyvariable = dstarD0pi->YDstar();
          if(fEfficiencyVariable == kNone) efficiencyvariable = 0;
          
        // get the D meson efficiency
        DmesonEfficiency = fCuts2->GetTrigWeight(dstarD0pi->Pt(),efficiencyvariable);
            
         

	    if(fUseDmesonEfficiencyCorrection){
	      if(DmesonEfficiency>1.e-5) DmesonWeight = 1./DmesonEfficiency;
	      else {// THIS ELSE STATEMENT MUST BE REFINED: THE EFFICIENCY MAP HAS TO BE REPLACED WITH A WEIGHT MAP COOKED A PRIORI
		if(ptDStar>2.) DmesonWeight = 0.5; // at high pt a zero value in the efficiency can come only from stat fluctutations in MC -> 0.5 is an arbitrary asymptotic value
		else DmesonWeight = 1.e+5; // at low pt it can be that the efficiency is really low
	      }
	    }
            else DmesonWeight = 1.; 
            
            // continue;
          
	    Int_t mcLabelDStar = -999;
          if(fmontecarlo){
	      // find associated MC particle for D* ->D0toKpi
	      mcLabelDStar = dstarD0pi->MatchToMC(413,421,pdgDgDStartoD0pi,pdgDgD0toKpi,fmcArray/*,kFALSE*/);
	      if(mcLabelDStar>=0) isDStarMCtag = kTRUE;
              if(!isDStarMCtag) continue;
            AliAODMCParticle *MCDStar = (AliAODMCParticle*)fmcArray->At(mcLabelDStar);
            //check if DStar from B
            Int_t labelMother = MCDStar->GetMother();
            AliAODMCParticle * mother = dynamic_cast<AliAODMCParticle*>(fmcArray->At(labelMother));
            if(!mother) continue;
            Int_t motherPDG =TMath::Abs(mother->PdgCode());
            if((motherPDG>=500 && motherPDG <600) || (motherPDG>=5000 && motherPDG<6000 )) isDfromB = kTRUE;
	    }
                        
	    
	    phiDStar = fCorrelator->SetCorrectPhiRange(phiDStar);
            
            // set the phi of the D meson in the correct range
            
	    Int_t ptbin=fCuts->PtBin(dstarD0pi->Pt());
            
	    
	    Double_t dmDStarWindow = 0.0019/3;// 0.0019 = 3 sigma
	    //	Double_t mD0Window=0.074/3;
            
            Double_t mD0Window= fD0Window[ptbin]/3;
            //cout << "Check with new window " << fD0Window[ptbin]/3 << endl;
            
	    
            
	    invMassDZero = dstarD0pi->InvMassD0();
	    if(!fmixing && fFullmode) ((TH2F*)fOutput->FindObject("D0InvMass"))->Fill(ptDStar,invMassDZero);
            
	    deltainvMDStar = dstarD0pi->DeltaInvMass();
            
	    //good D0 candidates
	    if (TMath::Abs(invMassDZero-mPDGD0)<fDMesonSigmas[1]*mD0Window){
	      
	      if(!fmixing)	((TH2F*)fOutput->FindObject("DeltaInvMass"))->Fill(ptDStar,deltainvMDStar,DmesonWeight);
	      // good D*?
	      if(TMath::Abs(deltainvMDStar-(mPDGDstar-mPDGD0))<fDMesonSigmas[0]*dmDStarWindow){
		
		if(!fmixing)	((TH1F*)fOutput->FindObject("RecoPtDStar"))->Fill(ptDStar,DmesonWeight);
		if(!fmixing)	((TH2F*)fOutput->FindObject("PhiEtaTrigger"))->Fill(phiDStar,etaDStar);
		isInPeak = kTRUE;
        EventHasDStarCandidate = kTRUE;
		nOfDStarCandidates++;
	      } // end Good D*
              
                //  D* sideband?
	      if((deltainvMDStar-(mPDGDstar-mPDGD0)>fDMesonSigmas[2]*dmDStarWindow) && (deltainvMDStar-(mPDGDstar-mPDGD0)<fDMesonSigmas[3]*dmDStarWindow)){
		isInDStarSideBand = kTRUE;
              EventHasDStarSideBandCandidate = kTRUE;
	      } // end D* sideband
              
            }// end good D0 candidates
            
            //D0 sidebands
	    if (TMath::Abs(invMassDZero-mPDGD0)>fDMesonSigmas[2]*mD0Window && TMath::Abs(invMassDZero-mPDGD0)<fDMesonSigmas[3]*mD0Window ){
	      if(!fmixing)((TH2F*)fOutput->FindObject("bkgDeltaInvMass"))->Fill(ptDStar,deltainvMDStar,DmesonWeight);
	      if(!fmixing && fFullmode)((TH2F*)fOutput->FindObject("D0InvMassinSB"))->Fill(ptDStar,invMassDZero,DmesonWeight);
              
	      if(TMath::Abs(deltainvMDStar-(mPDGDstar-mPDGD0))<fDMesonSigmas[0] *dmDStarWindow){ // is in DStar peak region?
		if(!fmixing)	((TH1F*)fOutput->FindObject("RecoPtBkg"))->Fill(ptDStar,DmesonWeight);
		isInDZeroSideBand = kTRUE;
        EventHasDZeroSideBandCandidate = kTRUE;
		nOfSBCandidates++;
		if(!fmixing)	((TH2F*)fOutput->FindObject("PhiEtaSideBand"))->Fill(phiDStar,etaDStar);
	      }
              
	    }//end if sidebands
            
	   
            
            
	    if(!isInPeak && !isInDStarSideBand && !isInDZeroSideBand) continue; // skip if it is not side band or peak event - SAVE CPU TIME
	    
	    
            // check properties of the events containing the D*

     
            
            isDStarCand = kTRUE;
            
            // charge of the daughter od the
	    daughtercharge = ((AliAODTrack*)dstarD0pi->GetBachelor())->Charge();
            
	    
	    trackiddaugh0 = ((AliAODTrack*)theD0particle->GetDaughter(0))->GetID();
	    trackiddaugh1 = ((AliAODTrack*)theD0particle->GetDaughter(1))->GetID();
	    trackidsoftPi = ((AliAODTrack*)dstarD0pi->GetBachelor())->GetID();
	    
	    // end here the reco
            
	    
	  }// end of if for applied reconstruction to D*
	  
	  Int_t DStarLabel = -1;
	  
	  if(!fReco){ // use pure MC information
          
        // get the DStar Particle
	    DStarMC = dynamic_cast<AliAODMCParticle*>(fmcArray->At(iDStartoD0pi));
	    if (!DStarMC) {
	      AliWarning("Careful: DStar MC Particle not found in tree, skipping");
	      continue;
	    }
	    DStarLabel = DStarMC->GetLabel();
        if(DStarLabel>0)isDStarMCtag = kTRUE;
	    
	    Int_t PDG =TMath::Abs(DStarMC->PdgCode());
	    if(PDG !=413) continue; // skip if it is not a DStar
	    // check fiducial acceptance
        if(!fCuts->IsInFiducialAcceptance(DStarMC->Pt(),DStarMC->Y())) continue;
            
            //check if DStar from B
            Int_t labelMother = DStarMC->GetMother();
            AliAODMCParticle * mother = dynamic_cast<AliAODMCParticle*>(fmcArray->At(labelMother));
	         if(!mother) continue;
            Int_t motherPDG =TMath::Abs(mother->PdgCode());
            if((motherPDG>=500 && motherPDG <600) || (motherPDG>=5000 && motherPDG<6000 )) isDfromB = kTRUE;
          
          Bool_t isDZero = kFALSE;
          Bool_t isSoftPi = kFALSE;
          
          if(fUseHadronicChannelAtKineLevel){
          //check decay channel on MC ************************************************
                Int_t NDaugh = DStarMC->GetNDaughters();
                if(NDaugh != 2) continue; // skip decay channels w/0 2 prongs
                
                for(Int_t i=0; i<NDaugh;i++){ // loop on daughters
                        Int_t daugh_label = DStarMC->GetDaughter(i);
                        AliAODMCParticle* mcDaughter = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daugh_label));
                        if(!mcDaughter) continue;
                        Int_t daugh_pdg = TMath::Abs(mcDaughter->GetPdgCode());
                        if(fDebugLevel) std::cout << "Daughter " << i << " pdg code is " << daugh_pdg << std::endl;
	      
                        if(daugh_pdg == 421) {isDZero = kTRUE;
                            Int_t NDaughD0 = mcDaughter->GetNDaughters();
                            if(NDaughD0 != 2) continue; // skip decay channels w/0 2 prongs
                            trackiddaugh0 = mcDaughter->GetDaughter(0);
                            trackiddaugh1 = mcDaughter->GetDaughter(1);
                            Bool_t isKaon = kFALSE;
                            Bool_t isPion = kFALSE;
		
                            for(Int_t k=0;k<NDaughD0;k++){
                                Int_t labelD0daugh = mcDaughter->GetDaughter(k);
                                AliAODMCParticle* mcGrandDaughter = dynamic_cast<AliAODMCParticle*>(fmcArray->At(labelD0daugh));
                                if(!mcGrandDaughter) continue;
                                Int_t granddaugh_pdg = TMath::Abs(mcGrandDaughter->GetPdgCode());
                                if(granddaugh_pdg==321) isKaon = kTRUE;
                                if(granddaugh_pdg==211) isPion = kTRUE;
                            }
                            if(!isKaon || !isKaon) continue; // skip if not correct decay channel of D0
                        }
	      
                        if(daugh_pdg == 211) {
                            isSoftPi = kTRUE;
                            daughtercharge = mcDaughter->Charge();
                            trackidsoftPi = daugh_label;}
                    }
              if(!isDZero || !isSoftPi) continue; // skip if not correct decay channel
          } // end check decay channel
	    
	    ptDStar = DStarMC->Pt();
	    phiDStar = DStarMC->Phi();
	    etaDStar = DStarMC->Eta();
          
          if(TMath::Abs(etaDStar) > fMaxEtaDStar) continue;
	    
	  } // end use pure MC information
	  
    
        // getting the number of triggers in the MCtag D* case
        if(fmontecarlo && isDStarMCtag) ((TH1F*)fOutput->FindObject("MCtagPtDStar"))->Fill(ptDStar);
	  if(fmontecarlo && isDStarMCtag && !isDfromB) ((TH1D*)fOutputMC->FindObject("MCtagPtDStarfromCharm"))->Fill(ptDStar);
	  if(fmontecarlo && isDStarMCtag && isDfromB) ((TH1D*)fOutputMC->FindObject("MCtagPtDStarfromBeauty"))->Fill(ptDStar);
	  
	  
	  fCorrelator->SetTriggerParticleProperties(ptDStar,phiDStar,etaDStar); // pass to the object the necessary trigger part parameters
	  fCorrelator->SetTriggerParticleDaughterCharge(daughtercharge);
	  
	  
	  // ************************************************ CORRELATION ANALYSIS STARTS HERE
	  
	  
	  Bool_t execPool = fCorrelator->ProcessEventPool();
	  
	  
	  if(fmixing && !execPool) {
	    AliInfo("Mixed event analysis: pool is not ready");
            if(!isEventMixingFilledPeak && isInPeak)  {
	      ((TH1D*)fOutput->FindObject("CheckPoolReadiness"))->Fill(1);
	      isEventMixingFilledPeak = kTRUE;
            }
            if (!isEventMixingFilledSB && isInDZeroSideBand)  {
	      ((TH1D*)fOutput->FindObject("CheckPoolReadiness"))->Fill(3);
	      isEventMixingFilledSB=kTRUE;
            }
            
	    continue;
	  }
	  
	  // check event topology
	  if(fmixing&&execPool){
            // pool is ready - run checks on bins filling
            if(!isEventMixingFilledPeak && isInPeak)  {
	      ((TH1D*)fOutput->FindObject("CheckPoolReadiness"))->Fill(0);
	      if(fFullmode) EventMixingChecks(aodEvent);
	      isEventMixingFilledPeak = kTRUE;
            }
            
            if(!isEventMixingFilledSB && isInDZeroSideBand) {
	      ((TH1D*)fOutput->FindObject("CheckPoolReadiness"))->Fill(2);
	      isEventMixingFilledSB=kTRUE;
            }
	  }
	  
	  Int_t NofEventsinPool = 1;
	  if(fmixing) NofEventsinPool = fCorrelator->GetNofEventsInPool();
	  
	  
	  for (Int_t jMix =0; jMix < NofEventsinPool; jMix++){// loop on events in the pool; if it is SE analysis, stops at one
	    
	    Bool_t analyzetracks = fCorrelator->ProcessAssociatedTracks(jMix);
	    if(!analyzetracks) {
	      AliInfo("AliHFCorrelator::Cannot process the track array");
	      continue;
	    }
            
            //initialization of variables for correlations with leading particles
            Double_t DeltaPhiLeading = -999.;
          Double_t DeltaEtaLeading = -999.;
	
	    
	    Int_t NofTracks = fCorrelator->GetNofTracks();
            
	    
	    if(isInPeak && fFullmode) ((TH1D*)fOutput->FindObject("NofTracksInPeak"))->Fill(NofTracks);
	    if(isInDZeroSideBand && fFullmode) ((TH1D*)fOutput->FindObject("NofTracksInSB"))->Fill(NofTracks);
	    
	    
	    
            Double_t arraytofill[5];
            Double_t MCarraytofill[7];
            
            
            Double_t weight;
            
          for(Int_t iTrack = 0; iTrack<NofTracks; iTrack++){ // looping on track candidates
              Bool_t runcorrelation = fCorrelator->Correlate(iTrack);
              if(!runcorrelation) continue;
              
              Double_t DeltaPhi = fCorrelator->GetDeltaPhi();
              Double_t DeltaEta = fCorrelator->GetDeltaEta();
	      
              AliReducedParticle * hadron = fCorrelator->GetAssociatedParticle();
              if(!hadron) {/*cout << "No Hadron" << endl;*/ continue;}
              
              Double_t ptHad = hadron->Pt();
              Double_t phiHad = hadron->Phi();
              Double_t etaHad = hadron->Eta();
              Int_t label = hadron->GetLabel();
              Int_t trackid = hadron->GetID();
              Double_t efficiency = hadron->GetWeight();
              
              weight = 1;
              if(fUseEfficiencyCorrection && efficiency){
                  weight = DmesonWeight * (1./efficiency);
              }
        
              phiHad = fCorrelator->SetCorrectPhiRange(phiHad);
              
              
              if(fFullmode) ((TH2F*)fOutput->FindObject("WeightChecks"))->Fill(ptHad,efficiency);
              
              arraytofill[0] = DeltaPhi;
              arraytofill[1] = DeltaEta;
              arraytofill[2] = ptDStar;
              arraytofill[3] = ptHad;
              arraytofill[4] = poolbin;
                
              
              MCarraytofill[0] = DeltaPhi;
              MCarraytofill[1] = DeltaEta;
              MCarraytofill[2] = ptDStar;
              MCarraytofill[3] = ptHad;
              MCarraytofill[4] = poolbin;
	      
	     
	      if(fmontecarlo){
		if(label<0 && fFullmode) ((TH2D*)fOutputMC->FindObject("TrackLabels"))->Fill(0.,NofTracks);
		if(label>=0 && fFullmode) ((TH2D*)fOutputMC->FindObject("TrackLabels"))->Fill(1.,NofTracks);
		if(label<0) continue; // skip track with wrong label
	      }
              
               Bool_t isDdaughter = kFALSE;
            // skip the D daughters in the correlation
              if(!fmixing && fReco){
                  if(trackid == trackiddaugh0) continue;
                  if(trackid == trackiddaugh1) continue;
                  if(trackid == trackidsoftPi) continue;
              }
              
              if(!fmixing && !fReco){
                  AliAODMCParticle *part = (AliAODMCParticle*)fmcArray->At(label);
                  if(!part) continue;
                  if(IsDDaughter(DStarMC, part)) continue;
                  cout << "Skipping DStar  daugheter " << endl;
              }
              if(!fmixing && !fReco && fmontecarlo){  // skip D* Daughetrs if it is Pure MCDStar
                    Int_t hadronlabel = label;
                    for(Int_t k=0; k<4;k++){ // go back 4 generations and check the mothers
                        if(DStarLabel<0){ break;}
                        if(hadronlabel<0) { break;}
                        AliAODMCParticle* mcParticle = dynamic_cast<AliAODMCParticle*>(fmcArray->At(hadronlabel));
                        if(!mcParticle) {AliInfo("NO MC PARTICLE"); break;}
                        hadronlabel = mcParticle->GetMother();
                        if(hadronlabel == DStarLabel) isDdaughter = kTRUE;
                    }
		
                  if(isDdaughter && fDebugLevel){
                      std::cout << "It is the D* daughter with label " << label << std::endl;
                      std::cout << "Daughter 0 label = " << trackiddaugh0 << std::endl;
                      std::cout << "Daughter 1 label = " << trackiddaugh1 << std::endl;
                      std::cout << "Soft pi label = " << trackidsoftPi << std::endl;
                  }
                    
					if(isDdaughter) continue; // skip if track is from DStar
				}
                
				// ================ FILL CORRELATION HISTOGRAMS ===============================
				
                // monte carlo case (mc tagged D*)
	      if((fmontecarlo && isDStarMCtag) || (fmontecarlo && !fReco)){ // check correlations of MC tagged DStars in MonteCarlo
                
		Bool_t* PartSource = fCuts2->IsMCpartFromHF(label,fmcArray); // check source of associated particle (hadron/kaon/K0)
	                      
		MCarraytofill[5] = 0;
		if(PartSource[0]) MCarraytofill[5] = 1;
		if(PartSource[1]) MCarraytofill[5] = 2;
		if(PartSource[2]&&PartSource[0]) MCarraytofill[5] = 3;
                    if(PartSource[2]&&PartSource[1]) MCarraytofill[5] = 4;
                    if(PartSource[3]) MCarraytofill[5] = 5;
                    if(!isDfromB) MCarraytofill[6] = 0;
                    if(isDfromB) MCarraytofill[6] = 1;
		    if(!fReco && TMath::Abs(etaHad)>0.8) {
		      delete [] PartSource;
		      continue; // makes sure you study the correlation on MC  truth only if particles are in acceptance
		    }
                    ((THnSparseF*)fOutputMC->FindObject("MCDStarCorrelationsDStarHadron"))->Fill(MCarraytofill);
                    
                    delete[] PartSource;
				}
              
                // Good DStar canidates
				if(isInPeak)  {
                    
					if(!fReco && TMath::Abs(etaHad)>0.8) continue; // makes sure you study the correlation on MC  truth only if particles are in acceptance
                    if(fselect==1)  ((THnSparseF*)fOutput->FindObject("CorrelationsDStarHadron"))->Fill(arraytofill,weight);
                    if(fselect==2)  ((THnSparseF*)fOutput->FindObject("CorrelationsDStarKaon"))->Fill(arraytofill,weight);
                    if(fselect==3)  ((THnSparseF*)fOutput->FindObject("CorrelationsDStarKZero"))->Fill(arraytofill,weight);
					
				    ((TH3F*)fOutput->FindObject("PhiEtaPart"))->Fill(phiHad,etaHad,MultipOrCent);
					if(fFullmode)((TH1D*)fOutput->FindObject("TracksInPeakSpectra"))->Fill(ptHad);
                                    
				}
              
                // Sidebands from D0 candidate
				if(isInDZeroSideBand) {
					
					if(!fReco && TMath::Abs(etaHad)>0.8) continue; // makes sure you study the correlation on MC  truth only if particles are in acceptance
                    if(fselect==1)  ((THnSparseF*)fOutput->FindObject("DZeroBkgCorrelationsDStarHadron"))->Fill(arraytofill,weight);
                    if(fselect==2)  ((THnSparseF*)fOutput->FindObject("DZeroBkgCorrelationsDStarKaon"))->Fill(arraytofill,weight);
                    if(fselect==3)  ((THnSparseF*)fOutput->FindObject("DZeroBkgCorrelationsDStarKZero"))->Fill(arraytofill,weight);
                    
					if(fFullmode) ((TH1D*)fOutput->FindObject("TracksInSBSpectra"))->Fill(ptHad);
                    
                }
              
              // Sidebands from D* candidate
                if(isInDStarSideBand) {
					
					if(!fReco && TMath::Abs(etaHad)>0.8) continue; // makes sure you study the correlation on MC  truth only if particles are in acceptance
                    if(fselect==1 && fFullmode)  ((THnSparseF*)fOutput->FindObject("DStarBkgCorrelationsDStarHadron"))->Fill(arraytofill,weight);
                    if(fselect==2 && fFullmode)  ((THnSparseF*)fOutput->FindObject("DStarBkgCorrelationsDStarKaon"))->Fill(arraytofill,weight);
                    if(fselect==3 && fFullmode)  ((THnSparseF*)fOutput->FindObject("DStarBkgCorrelationsDStarKZero"))->Fill(arraytofill,weight);

				}
                
                
			} // end loop on track candidates
            
			
			
            // fill the leading particle histograms
            
            if(isInPeak && fFullmode) ((TH3D*)fOutput->FindObject("LeadingCand"))->Fill(DeltaPhiLeading,ptDStar,DeltaEtaLeading);
			if(isInDZeroSideBand && fFullmode) ((TH3D*)fOutput->FindObject("LeadingSB"))->Fill(DeltaPhiLeading,ptDStar,DeltaEtaLeading);
            
		} // end loop on events in the pool
        
	}// end loop on D* candidates
	
    
 
        // check events with D* or SB canidates
    if(fFullmode && EventHasDStarCandidate)  ((TH2F*)fOutput->FindObject("EventPropsCheckifDStar"))->Fill(MultipOrCent,zVtxPosition);
     if(fFullmode && EventHasDZeroSideBandCandidate)  ((TH2F*)fOutput->FindObject("EventPropsCheckifDZeroSB"))->Fill(MultipOrCent,zVtxPosition);
    
    if(fFullmode && EventHasDStarCandidate)  ((TH2F*)fOutput->FindObject("EventPropsCheckifDStarSB"))->Fill(MultipOrCent,zVtxPosition);
    
    
	if(fFullmode) ((TH2F*)fOutput->FindObject("DStarCandidates"))->Fill(nOfDStarCandidates,MultipOrCent);
    if(fFullmode) ((TH2F*)fOutput->FindObject("SBCandidates"))->Fill(nOfSBCandidates,MultipOrCent);
	
    // update event pool
    Bool_t updated = fCorrelator->PoolUpdate();
	
    //	if(updated) EventMixingChecks(aodEvent);
	if(!updated) AliInfo("Pool was not updated");
	
    
} //end the exec

//________________________________________ terminate ___________________________
void AliAnalysisTaskDStarCorrelations::Terminate(Option_t*)
{    
	// The Terminate() function is the last function to be called during
	// a query. It always runs on the client, it can be used to present
	// the results graphically or save the results to file.
	
	AliAnalysisTaskSE::Terminate();
	
	fOutput = dynamic_cast<TList*> (GetOutputData(1));
	if (!fOutput) {     
		printf("ERROR: fOutput not available\n");
		return;
	}

	return;
}
//_____________________________________________________
Bool_t AliAnalysisTaskDStarCorrelations::IsDDaughter(AliAODMCParticle* d, AliAODMCParticle* track) const {
    
    //Daughter removal in MCKine case
    Bool_t isDaughter = kFALSE;
    Int_t labelD0 = d->GetLabel();
    
    Int_t mother = track->GetMother();
    
    //Loop on the mothers to find the D0 label (it must be the trigger D0, not a generic D0!)
    while (mother > 0){
        AliAODMCParticle* mcMoth = dynamic_cast<AliAODMCParticle*>(fmcArray->At(mother)); //it's the mother of the track!
        if (mcMoth){
            if (mcMoth->GetLabel() == labelD0) isDaughter = kTRUE;
            mother = mcMoth->GetMother(); //goes back by one
        } else{
            AliError("Failed casting the mother particle!");
            break;
        }
    }
    
    return isDaughter;
}

//_____________________________________________________
void AliAnalysisTaskDStarCorrelations::DefineThNSparseForAnalysis(){
    
    //cout << "DEFINING THNSPARSES "<< endl;
    
    Double_t Pi = TMath::Pi();
	Int_t nbinscorr = fPhiBins;
	Double_t lowcorrbin = -0.5*Pi;
	Double_t upcorrbin = 1.5*Pi;
    // define the THnSparseF
    
    //sparse bins
    
    //1 delta_phi
    //2 delta_eta
    //3 D* pt
    //4 multiplicity
    //5 track pt
    //6 zVtx position
    
    Int_t nbinsPool = (fCuts2->GetNZvtxPoolBins())*(fCuts2->GetNCentPoolBins());
    
    
    Int_t nbinsSparse[5]=         {nbinscorr,   32,30,250,nbinsPool};
    Double_t binLowLimitSparse[5]={lowcorrbin,-1.6, 0,  0,-0.5};
    Double_t binUpLimitSparse[5]= {upcorrbin,  1.6,30,25,nbinsPool-0.5};
  
    Int_t MCnbinsSparse[7]=         {nbinscorr,   32,30,250,nbinsPool,10,2};    
    Double_t MCbinLowLimitSparse[7]={lowcorrbin,-1.6, 0,  0,-0.5,-0.5,-0.5};      //  
    Double_t MCbinUpLimitSparse[7]= {upcorrbin,  1.6,30,25,nbinsPool-0.5,9.5,1.5};
    
    TString sparsename = "CorrelationsDStar";
    if(fselect==1) sparsename += "Hadron";
	if(fselect==2) sparsename += "Kaon";
	if(fselect==3) sparsename += "KZero";
    
    TString D0Bkgsparsename = "DZeroBkg";
    D0Bkgsparsename += sparsename;
    
    TString DStarBkgsparsename = "DStarBkg";
    DStarBkgsparsename += sparsename;
    
    TString MCSparseName = "MCDStar";
    MCSparseName += sparsename;
    // signal correlations
    THnSparseF * Correlations = new THnSparseF(sparsename.Data(),"Correlations for signal",5,nbinsSparse,binLowLimitSparse,binUpLimitSparse);
    
    // bkg correlations from D0 sidebands
    THnSparseF * DZeroBkgCorrelations = new THnSparseF(D0Bkgsparsename.Data(),"Bkg Correlations estimated with D0 sidebands",5,nbinsSparse,binLowLimitSparse,binUpLimitSparse);
    
    // bkg correlations from D* sidebands
    THnSparseF * DStarBkgCorrelations = new THnSparseF(DStarBkgsparsename.Data(),"Bkg Correlations estimated with D* sidebands",5,nbinsSparse,binLowLimitSparse,binUpLimitSparse);
    
    
    THnSparseF * MCCorrelations = new THnSparseF(MCSparseName.Data(),"MC Correlations",7,MCnbinsSparse,MCbinLowLimitSparse,MCbinUpLimitSparse);
    
    MCCorrelations->GetAxis(5)->SetBinLabel(1," All ");
	MCCorrelations->GetAxis(5)->SetBinLabel(2," from hadron Heavy flavour");
	MCCorrelations->GetAxis(5)->SetBinLabel(3," from c->D");
	MCCorrelations->GetAxis(5)->SetBinLabel(4," from b->D");
	MCCorrelations->GetAxis(5)->SetBinLabel(5," from b->B");
	MCCorrelations->GetAxis(5)->SetBinLabel(6," from quark Heavy flavour");
	MCCorrelations->GetAxis(5)->SetBinLabel(7," from c");
	MCCorrelations->GetAxis(5)->SetBinLabel(8," from b");
    
    MCCorrelations->GetAxis(6)->SetBinLabel(1," if D* from c");
    MCCorrelations->GetAxis(6)->SetBinLabel(2," if D* from b");
    
    Correlations->Sumw2();
    DZeroBkgCorrelations->Sumw2();
    DStarBkgCorrelations->Sumw2();
    
    fOutput->Add(Correlations);
    fOutput->Add(DZeroBkgCorrelations);
    if(fFullmode) fOutput->Add(DStarBkgCorrelations);
    if(fmontecarlo) fOutputMC->Add(MCCorrelations);
    
    
}
//__________________________________________________________________________________________________
void AliAnalysisTaskDStarCorrelations::DefineHistoForAnalysis(){
	
	Double_t Pi = TMath::Pi();
	Int_t nbinscorr = fPhiBins;
	Double_t lowcorrbin = -0.5*Pi ; // shift the bin by half the width so that at 0 is it the bin center
	Double_t upcorrbin = 1.5*Pi ;
	
	// ========================= histograms for both Data and MonteCarlo
	
	
	TH1D * NofEvents = new TH1D("NofEvents","NofEvents",12,-0.5,11.5);
    NofEvents->GetXaxis()->SetBinLabel(1," All events");
	NofEvents->GetXaxis()->SetBinLabel(2," Selected events");
	NofEvents->GetXaxis()->SetBinLabel(3," Rejected - SPD Pileup");
	NofEvents->GetXaxis()->SetBinLabel(4," Rejected - Centrality");
	NofEvents->GetXaxis()->SetBinLabel(5," Rejected - No Reco Vtx");
	NofEvents->GetXaxis()->SetBinLabel(6," Rejected - Vtx Contr.");
	NofEvents->GetXaxis()->SetBinLabel(7," Rejected - Vtx outside fid.acc.");
	NofEvents->GetXaxis()->SetBinLabel(8," Rejected - Trigger");
    NofEvents->GetXaxis()->SetBinLabel(9," Rejected - Phys.Sel");
    NofEvents->GetXaxis()->SetBinLabel(10," Rejected - pA - 1st in chunk");
    NofEvents->GetXaxis()->SetBinLabel(11," Rejected - pA - bad vtx");
    fOutput->Add(NofEvents);
	
	
	
	
	TH2F *D0InvMass = new TH2F("D0InvMass","K#pi invariant mass distribution",300,0,30,1500,0.5,3.5);
	if(!fmixing && fFullmode) fOutput->Add(D0InvMass);
	
	TH2F *D0InvMassinSB = new TH2F("D0InvMassinSB","K#pi invariant mass distribution in sb",300,0,30,1500,0.5,3.5);
	if(!fmixing && fFullmode) fOutput->Add(D0InvMassinSB);
	
	//TH2F *DeltaInvMass = new TH2F("DeltaInvMass","K#pi#pi - K#pi invariant mass distribution",300,0,30,750,0.1,0.2);
	//if(!fmixing) fOutput->Add(DeltaInvMass);
    TH2F *DeltaInvMass = new TH2F("DeltaInvMass","K#pi#pi - K#pi invariant mass distribution; D* p_{T}; #DeltaInvMass",30,0,30,750,0.1,0.2);
	if(!fmixing) fOutput->Add(DeltaInvMass);
	
	TH2F *bkgDeltaInvMass = new TH2F("bkgDeltaInvMass","K#pi#pi - K#pi invariant mass distribution; SB p_{T}; #DeltaInvMass",30,0,30,750,0.1,0.2);
	if(!fmixing) fOutput->Add(bkgDeltaInvMass);
    
    DeltaInvMass->Sumw2();
    bkgDeltaInvMass->Sumw2();
	
	TH1F *RecoPtDStar = new TH1F("RecoPtDStar","RECO DStar pt distribution",50,0,50);
	if(!fmixing) fOutput->Add(RecoPtDStar);
	
	TH1F *RecoPtBkg = new TH1F("RecoPtBkg","RECO pt distribution side bands",50,0,50);
	if(!fmixing) fOutput->Add(RecoPtBkg);
    
    TH1D *MCtagPtDStarfromCharm = new TH1D("MCtagPtDStarfromCharm","RECO pt of MCtagged DStars from charm",50,0,50);
    if(fmontecarlo) fOutputMC->Add(MCtagPtDStarfromCharm);
    
    TH1D *MCtagPtDStarfromBeauty = new TH1D("MCtagPtDStarfromBeauty","RECO pt of MCtagged DStars from beauty",50,0,50);
    if(fmontecarlo) fOutputMC->Add(MCtagPtDStarfromBeauty);
	
	TH1F *MCtagPtDStar = new TH1F("MCtagPtDStar","RECO pt of MCtagged DStars side bands",50,0,50);
	if(!fmixing) fOutput->Add(MCtagPtDStar);
	
	TH2F *KZeroSpectra = new TH2F("KZeroSpectra","Spectra of K0s",500,0.3,0.8,250,0,25);
    if(fselect==3 && fFullmode) fOutput->Add(KZeroSpectra);
	
	TH2F *KZeroSpectraifHF = new TH2F("KZeroSpectraifHF","Spectra of K0s in association with a D*",500,0.3,0.8,250,0,25);
	if(fselect==3 && fFullmode) fOutput->Add(KZeroSpectraifHF);
	
	TH1D * NofTracksInPeak = new TH1D("NofTracksInPeak","N of associated tracks per D trigger; Nof tracks; Entries",500,-0.5,499.5);
	if(fFullmode) fOutput->Add(NofTracksInPeak);
	
	TH1D * NofTracksInSB = new TH1D("NofTracksInSB","N of associated tracks per SideBand trigger; Nof tracks; Entries",500,-0.5,499.5);
	if(fFullmode) fOutput->Add(NofTracksInSB);
	
	TH1D * TracksInPeakSpectra = new TH1D("TracksInPeakSpectra","Pt Spectra tracks with D trigger; p_{T} GeV/c; Entries",500,-0.5,49.5);
	if(fFullmode)fOutput->Add(TracksInPeakSpectra);
	
	TH1D * TracksInSBSpectra = new TH1D("TracksInSBSpectra","Pt Spectra tracks with SideBand trigger; p_{T} GeV/c; Entries",500,-0.5,49.5);
	if(fFullmode)fOutput->Add(TracksInSBSpectra);
	
	
	//TH2I * EventMixingCheck = new TH2I("EventMixingCheck","EventMixingCheck",5,-0.5,4.5,7,-0.5,6.5);
	//if(fmixing) fOutput->Add(EventMixingCheck);
    
    
    TH2F * EventPropsCheck = new TH2F("EventPropsCheck","Properties of the event; Multiplicity; ZVtx Position [cm]",1000,0,1000,40,-10,10);
	if(fFullmode)fOutput->Add(EventPropsCheck);
    
    TH2F * EventPropsCheckifDStar = new TH2F("EventPropsCheckifDStar","Properties of the event with D* Cand; Multiplicity; ZVtx Position [cm]",1000,0,1000,40,-10,10);
	if(fFullmode)fOutput->Add(EventPropsCheckifDStar);
    
    TH2F * EventPropsCheckifDZeroSB = new TH2F("EventPropsCheckifDZeroSB","Properties of the event with D* Cand; Multiplicity; ZVtx Position [cm]",1000,0,1000,40,-10,10);
	if(fFullmode)fOutput->Add(EventPropsCheckifDZeroSB);
    
    TH2F * EventPropsCheckifDStarSB = new TH2F("EventPropsCheckifDStarSB","Properties of the event with D* Cand; Multiplicity; ZVtx Position [cm]",1000,0,1000,40,-10,10);
	if(fFullmode)fOutput->Add(EventPropsCheckifDStarSB);
    
	
    TH2F * WeightChecks = new TH2F("WeightChecks","Checks on efficiency correction",300,0,30,100,0.005,1.005);
	if(fFullmode)fOutput->Add(WeightChecks);
    
	
	
	TH2F * PhiEtaTrigger = new TH2F("PhiEtaTrigger","#phi distribution of the trigger particle",nbinscorr,lowcorrbin,upcorrbin,18,-0.9,0.9);
	fOutput->Add(PhiEtaTrigger);
	
	TH2F * PhiEtaSideBand = new TH2F("PhiEtaSideBand","#phi distribution of the sideband particle",nbinscorr,lowcorrbin,upcorrbin,18,-0.9,0.9);
	fOutput->Add(PhiEtaSideBand);
	
	TH3F * PhiEtaPart = new TH3F("PhiEtaPart","#phi distribution of the associated particle; #phi; #eta; multiplicity",nbinscorr,lowcorrbin,upcorrbin,18,-0.9,0.9,100,0,1000);
	fOutput->Add(PhiEtaPart);
    
    TH2F * DStarCandidates = new TH2F("DStarCandidates","# of D* candidates per event vs multiplicity",6,-0.5,5.5,50,0,500);
	if(fFullmode)fOutput->Add(DStarCandidates);
    
    TH2F * SBCandidates = new TH2F("SBCandidates","# of SB candidates per event vs multiplicity",6,-0.5,5.5,50,0,500);
	if(fFullmode)fOutput->Add(SBCandidates);
	
	
	//correlations histograms
	TString histoname1 = "DPhiDStar";
	if(fselect==1) histoname1 += "Hadron";
	if(fselect==2) histoname1 += "Kaon";
	if(fselect==3) histoname1 += "KZero";
	
	/*
     TH3D * DPhiDStar = new TH3D(histoname1.Data(),histoname1.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-2,2);
     
     TH3D * DPhiDStarKZero1 = new TH3D("DPhiDStarKZero1","DPhiDStarKZero1",nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-2,2);
     
     //side band background histograms
     TString histoname2 = "bkg";
     histoname2 += histoname1;
     TH3D * bkgDPhiDStar = new TH3D(histoname2.Data(),histoname2.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-2,2);
     TH3D * bkgDPhiDStarKZero1 = new TH3D("bkgDPhiDStarKZero1","bkgDPhiDStarKZero1",nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-2,2);
     
     
     fOutput->Add(DPhiDStar);
     
     if(fselect==3){fOutput->Add(DPhiDStarKZero1);}
     
     fOutput->Add(bkgDPhiDStar);
     
     if(fselect==3){fOutput->Add(bkgDPhiDStarKZero1);}
     
     */
	// leading particle
	TH3D * leadingcand = new TH3D("LeadingCand","LeadingCand",nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-2,2);
	TH3D * leadingsidebands = new TH3D("LeadingSB","LeadingSB",nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-2,2);
	
	if(fFullmode)fOutput->Add(leadingcand);
	if(fFullmode)fOutput->Add(leadingsidebands);
	
	// ========================= histos for analysis on MC only
	
	TH1D * EventTypeMC = new TH1D("EventTypeMC","EventTypeMC",100,-0.5,99.5);
	if(fmontecarlo) fOutputMC->Add(EventTypeMC);
	
	TH1F * MCSources = new TH1F("MCSources","Origin of associated particles in MC", 10, -0.5, 9.5);
	MCSources->GetXaxis()->SetBinLabel(1," All ");
	MCSources->GetXaxis()->SetBinLabel(2," from hadron Heavy flavour");
	MCSources->GetXaxis()->SetBinLabel(3," from c->D");
	MCSources->GetXaxis()->SetBinLabel(4," from b->D");
	MCSources->GetXaxis()->SetBinLabel(5," from b->B");
	MCSources->GetXaxis()->SetBinLabel(6," from quark Heavy flavour");
	MCSources->GetXaxis()->SetBinLabel(7," from c");
	MCSources->GetXaxis()->SetBinLabel(8," from b");
	
	if(fmontecarlo) fOutputMC->Add(MCSources);
    
    // leading particle from mc source
    TH1F * LeadingMCSources = new TH1F("LeadingMCSources","Origin of associated leading particles in MC", 10, -0.5, 9.5);
	LeadingMCSources->GetXaxis()->SetBinLabel(1," All ");
	LeadingMCSources->GetXaxis()->SetBinLabel(2," from hadron Heavy flavour");
	LeadingMCSources->GetXaxis()->SetBinLabel(3," from c->D");
	LeadingMCSources->GetXaxis()->SetBinLabel(4," from b->D");
	LeadingMCSources->GetXaxis()->SetBinLabel(5," from b->B");
	LeadingMCSources->GetXaxis()->SetBinLabel(6," from quark Heavy flavour");
	LeadingMCSources->GetXaxis()->SetBinLabel(7," from c");
	LeadingMCSources->GetXaxis()->SetBinLabel(8," from b");
	
	if(fmontecarlo && fFullmode) fOutputMC->Add(LeadingMCSources);
	
    // all hadrons
	TString histoname3 = "MCTag";
	histoname3 += histoname1;
	TH3D * MCTagDPhiDStar = new TH3D(histoname3.Data(),histoname3.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-2,2);
	
	TString histoname44 = "CharmDOrigin";
	histoname44 += histoname1;
	histoname44 += "MC";
	
	TH3D * CharmDOriginDPhiDStar = new TH3D(histoname44.Data(),histoname44.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-2,2);
	
	
	TString histoname54 = "BeautyDOrigin";
	histoname54 += histoname1;
	histoname54 += "MC";
	TH3D * BeautyDOriginDPhiDStar = new TH3D(histoname54.Data(),histoname54.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-2,2);
	
	TString histoname55 = "BeautyBOrigin";
	histoname55 += histoname1;
	histoname55 += "MC";
	TH3D * BeautyBOriginDPhiDStar = new TH3D(histoname55.Data(),histoname55.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-2,2);
	
	TString histoname4 = "CharmQuarkOrigin";
	histoname4 += histoname1;
	histoname4 += "MC";
	TH3D * CharmQuarkOriginDPhiDStar = new TH3D(histoname4.Data(),histoname4.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-2,2);
	
	TString histoname5 = "BeautyQuarkOrigin";
	histoname5 += histoname1;
	histoname5 += "MC";
	TH3D * BeautyQuarkOriginDPhiDStar = new TH3D(histoname5.Data(),histoname5.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-2,2);
	
	TString histoname6 = "NonHFOrigin";
	histoname6 += histoname1;
	histoname6 += "MC";
	TH3D * NonHFOriginDPhiDStar = new TH3D(histoname6.Data(),histoname6.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-2,2);
    
    if(fmontecarlo && fFullmode){
        
        fOutputMC->Add(MCTagDPhiDStar);
        fOutputMC->Add(CharmDOriginDPhiDStar);
        fOutputMC->Add(BeautyDOriginDPhiDStar);
        fOutputMC->Add(BeautyBOriginDPhiDStar);
        fOutputMC->Add(CharmQuarkOriginDPhiDStar);
        fOutputMC->Add(BeautyQuarkOriginDPhiDStar);
		fOutputMC->Add(NonHFOriginDPhiDStar);
        
	}
    
    // ========================= histos for analysis on MC
    // all leading hadron
	TString Leadinghistoname3 = "LeadingMCTag";
	Leadinghistoname3 += histoname1;
	TH3D * LeadingMCTagDPhiDStar = new TH3D(Leadinghistoname3.Data(),Leadinghistoname3.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-2,2);
    
	TString Leadinghistoname44 = "LeadingCharmDOrigin";
	Leadinghistoname44 += histoname1;
	Leadinghistoname44 += "MC";
	
	TH3D * LeadingCharmDOriginDPhiDStar = new TH3D(Leadinghistoname44.Data(),Leadinghistoname44.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-2,2);
	
	
	TString Leadinghistoname54 = "LeadingBeautyDOrigin";
	Leadinghistoname54 += histoname1;
	Leadinghistoname54 += "MC";
	TH3D * LeadingBeautyDOriginDPhiDStar = new TH3D(Leadinghistoname54.Data(),Leadinghistoname54.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-2,2);
	
	TString Leadinghistoname55 = "LeadingBeautyBOrigin";
	Leadinghistoname55 += histoname1;
	Leadinghistoname55 += "MC";
	TH3D * LeadingBeautyBOriginDPhiDStar = new TH3D(Leadinghistoname55.Data(),Leadinghistoname55.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-2,2);
	
	TString Leadinghistoname4 = "LeadingCharmQuarkOrigin";
	Leadinghistoname4 += histoname1;
	Leadinghistoname4 += "MC";
	TH3D * LeadingCharmQuarkOriginDPhiDStar = new TH3D(Leadinghistoname4.Data(),Leadinghistoname4.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-2,2);
	
	TString Leadinghistoname5 = "LeadingBeautyQuarkOrigin";
	Leadinghistoname5 += histoname1;
	Leadinghistoname5 += "MC";
	TH3D * LeadingBeautyQuarkOriginDPhiDStar = new TH3D(Leadinghistoname5.Data(),Leadinghistoname5.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-2,2);
    
    
	
	
	if(fmontecarlo && fFullmode){
		
		fOutputMC->Add(LeadingMCTagDPhiDStar);
		fOutputMC->Add(LeadingCharmDOriginDPhiDStar);
		fOutputMC->Add(LeadingBeautyDOriginDPhiDStar);
		fOutputMC->Add(LeadingBeautyBOriginDPhiDStar);
		fOutputMC->Add(LeadingCharmQuarkOriginDPhiDStar);
		fOutputMC->Add(LeadingBeautyQuarkOriginDPhiDStar);
		
	}
	
	TH3F * MCPhiEtaPart = new TH3F("MCPhiEtaPart","#phi distribution of the associated particle",nbinscorr,lowcorrbin,upcorrbin,50,-2.5,2.5,6,-0.5,6.5);
	MCPhiEtaPart->GetZaxis()->SetBinLabel(1,"All particles");
	MCPhiEtaPart->GetZaxis()->SetBinLabel(2,"from c quark");
	MCPhiEtaPart->GetZaxis()->SetBinLabel(3,"from b quark");
	MCPhiEtaPart->GetZaxis()->SetBinLabel(4,"from D from c");
	MCPhiEtaPart->GetZaxis()->SetBinLabel(5,"from D from b");
	MCPhiEtaPart->GetZaxis()->SetBinLabel(6,"from B from b");
	if(fmontecarlo) fOutputMC->Add(MCPhiEtaPart);
	
	TH2D * TrackLabels = new TH2D("TrackLabels","NofEvents;track label; multiplicity",2,-0.5,1.5,500,-0.5,499.5);
	if(fmontecarlo && fFullmode) fOutputMC->Add(TrackLabels);
	
	// ============================= EVENT MIXING CHECKS ======================================
	
	Int_t MaxNofEvents = fCuts2->GetMaxNEventsInPool();
	Int_t MinNofTracks = fCuts2->GetMinNTracksInPool();
	Int_t NofCentBins = fCuts2->GetNCentPoolBins();
	Double_t * CentBins = fCuts2->GetCentPoolBins();
	Int_t NofZVrtxBins = fCuts2->GetNZvtxPoolBins();
	Double_t *ZVrtxBins = fCuts2->GetZvtxPoolBins();
    
    
    
	Int_t k =0;
	
	if(fSystem == AA) k = 100; // PbPb centrality
	if(fSystem == pp || fSystem == pA) k = NofCentBins; // pp multiplicity
	
	
	//Double_t minvalue = CentBins[0];
	//Double_t maxvalue = CentBins[NofCentBins+1];
	//Double_t Zminvalue = ZVrtxBins[0];
	//Double_t Zmaxvalue = ZVrtxBins[NofCentBins+1];
	
	Double_t minvalue, maxvalue;
    Double_t Zminvalue, Zmaxvalue;
    
    Zminvalue = -15.;
    Zmaxvalue = 15;
    if(fSystem == AA) {minvalue = 0; maxvalue = 100;} // PbPb
    if(fSystem == pp || fSystem == pA) {minvalue = 0; maxvalue = 500;} // multilpicity
    
	//Double_t Nevents[]={0,2*MaxNofEvents/10,4*MaxNofEvents/10,6*MaxNofEvents/10,8*MaxNofEvents/10,MaxNofEvents};
    // Double_t Nevents[]={0,2*MaxNofEvents/10,4*MaxNofEvents/10,6*MaxNofEvents/10,8*MaxNofEvents/10,MaxNofEvents};
    //	Double_t * events = Nevents;
    Double_t eventsv[] ={0,1000000};
    //Double_t * events = new Double_t[2];
   // events[0] = 0;
//	events[1] = 1000000;
    Double_t *events = eventsv;
    Int_t Nevents = 1000000;
    //  TH3D * EventsPerPoolBin = new TH3D("EventsPerPoolBin","Number of events in bin pool",NofCentBins,CentBins,NofZVrtxBins,ZVrtxBins,Nevents,events);
    
    TH3D * EventsPerPoolBin = new TH3D("EventsPerPoolBin","Number of events in bin pool",NofCentBins,minvalue,maxvalue,NofZVrtxBins,-15,15,Nevents,events[0],events[1]);
    
	EventsPerPoolBin->GetXaxis()->SetTitle("Centrality/multiplicity ");
	EventsPerPoolBin->GetYaxis()->SetTitle("Z vertex [cm]");
	EventsPerPoolBin->GetZaxis()->SetTitle("Number of events in pool bin");
	if(fmixing && fFullmode) fOutput->Add(EventsPerPoolBin);
	
	Int_t MaxNofTracks = (MaxNofEvents+1)*MinNofTracks;
	//Int_t Diff = MaxNofTracks-MinNofTracks;
	
    //Double_t Ntracks[]={MinNofTracks,MinNofTracks+Diff/5,MinNofTracks+2*Diff/5,MinNofTracks+3*Diff/5,MinNofTracks+4*Diff/5,MaxNofTracks};
    //	Double_t  * trackN = Ntracks;
	
	TH3D * NofTracksPerPoolBin = new TH3D("NofTracksPerPoolBin","Number of tracks in bin pool",NofCentBins,minvalue,maxvalue,NofZVrtxBins,-15,15,MaxNofTracks,0,MaxNofTracks);
	NofTracksPerPoolBin->GetXaxis()->SetTitle("Centrality/multiplicity ");
	NofTracksPerPoolBin->GetYaxis()->SetTitle("Z vertex [cm]");
	NofTracksPerPoolBin->GetZaxis()->SetTitle("Number of tracks per bin");
	
	if(fmixing && fFullmode) fOutput->Add(NofTracksPerPoolBin);
	
	TH2D * NofPoolBinCalls = new TH2D("NofPoolBinCalls","Calls per pool bin",NofCentBins,CentBins,NofZVrtxBins,ZVrtxBins);
	NofPoolBinCalls->GetXaxis()->SetTitle("Centrality/multiplicity ");
	NofPoolBinCalls->GetYaxis()->SetTitle("Z vertex [cm]");
	if(fmixing && fFullmode) fOutput->Add(NofPoolBinCalls);
	
    
	
	TH2D * EventProps = new TH2D("EventProps","Event properties",100,minvalue,maxvalue,100,Zminvalue,Zmaxvalue);
	EventProps->GetXaxis()->SetTitle("Centrality/multiplicity ");
	EventProps->GetYaxis()->SetTitle("Z vertex [cm]");
	if(fmixing && fFullmode) fOutput->Add(EventProps);
    
    TH1D * CheckPoolReadiness = new TH1D("CheckPoolReadiness","Pool readiness",5,-0.5,4.5);
    CheckPoolReadiness->GetXaxis()->SetBinLabel(1,"Have a D cand, pool is ready");
	CheckPoolReadiness->GetXaxis()->SetBinLabel(2,"Have a D cand, pool is not ready");
    CheckPoolReadiness->GetXaxis()->SetBinLabel(3,"Have a SB cand, pool is ready");
	CheckPoolReadiness->GetXaxis()->SetBinLabel(4,"Have a SB cand, pool is not ready");
	
    if(fmixing) fOutput->Add(CheckPoolReadiness);
    
	
}

//__________________________________________________________________________________________________
void AliAnalysisTaskDStarCorrelations::EnlargeDZeroMassWindow(){
    

  //Float_t* ptbins = fCuts->GetPtBinLimits();
  if(fD0Window) delete fD0Window;
  fD0Window = new Float_t[fNofPtBins];
    
  AliInfo("Enlarging the D0 mass windows from cut object\n"); 
  Int_t nvars = fCuts->GetNVars();

  if(nvars<1){
    AliWarning("EnlargeDZeroMassWindow: 0 variables in cut object... check!");
    return;
  }
  Float_t** rdcutsvalmine=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++){
    rdcutsvalmine[iv]=new Float_t[fNofPtBins];
  }
    
    
    for (Int_t k=0;k<nvars;k++){
      for (Int_t j=0;j<fNofPtBins;j++){
            
	// enlarge D0 window
	if(k==0)	{
	  fD0Window[j] =fCuts->GetCutValue(0,j);
	  rdcutsvalmine[k][j] = 5.* fCuts->GetCutValue(0,j);
	  cout << "the set window = " << fD0Window[j] << " for ptbin " << j << endl;
	}
	else rdcutsvalmine[k][j] =fCuts->GetCutValue(k,j);
			
	// set same windows
	//rdcutsvalmine[k][j] =oldCuts->GetCutValue(k,j);
      }
    }
    
    fCuts->SetCuts(nvars,fNofPtBins,rdcutsvalmine);
    
    AliInfo("\n New windows set\n");     
    fCuts->PrintAll();
    
    
    for(Int_t iv=0;iv<nvars;iv++){
      delete rdcutsvalmine[iv];
    }
    delete [] rdcutsvalmine;
    
}


//____________________________  Run checks on event mixing ___________________________________________________
void AliAnalysisTaskDStarCorrelations::EventMixingChecks(AliAODEvent* AOD){
	
    
	AliCentrality *centralityObj = 0;
	Int_t multiplicity = -1;
	Double_t MultipOrCent = -1;
	
	// get the pool for event mixing
	if(fSystem != AA){ // pp
		multiplicity = AOD->GetNTracks();
		MultipOrCent = multiplicity; // convert from Int_t to Double_t
	}
	if(fSystem == AA){ // PbPb
		
		centralityObj = AOD->GetHeader()->GetCentralityP();
		MultipOrCent = centralityObj->GetCentralityPercentileUnchecked("V0M");
		AliInfo(Form("Centrality is %f", MultipOrCent));
	}
	
	AliAODVertex *vtx = AOD->GetPrimaryVertex();
	Double_t zvertex = vtx->GetZ(); // zvertex
	
	
	
	
	AliEventPool * pool = fCorrelator->GetPool();
	
    
	
	
	((TH2D*)fOutput->FindObject("NofPoolBinCalls"))->Fill(MultipOrCent,zvertex); // number of calls of pool
	((TH2D*)fOutput->FindObject("EventProps"))->Fill(MultipOrCent,zvertex); // event properties
	
	((TH3D*)fOutput->FindObject("EventsPerPoolBin"))->Fill(MultipOrCent,zvertex,pool->GetCurrentNEvents()); // number of events in the pool
	((TH3D*)fOutput->FindObject("NofTracksPerPoolBin"))->Fill(MultipOrCent,zvertex,pool->NTracksInPool()); // number of calls of pool
}
	




