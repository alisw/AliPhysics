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

/* $Id: AliAnalysisTaskDStarCorrelations.cxx 65139 2013-11-25 14:47:45Z fprino $ */

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
fBkgMethod(kDZeroSB),
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
fDmesonOutput(0x0),
fTracksOutput(0x0),
fEMOutput(0x0),
fCorrelationOutput(0x0),
fOutputMC(0x0),

fCuts(0),
fAssocCuts(0),
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
fBkgMethod(kDZeroSB),
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
fDmesonOutput(0x0),
fTracksOutput(0x0),
fEMOutput(0x0),
fCorrelationOutput(0x0),
fOutputMC(0x0),

fCuts(0),
fAssocCuts(AsscCuts),
fUtils(0),
fTracklets(0),
fDeffMapvsPt(0),
fDeffMapvsPtvsMult(0),
fDeffMapvsPtvsEta(0)
{
     Info("AliAnalysisTaskDStarCorrelations","Calling Constructor");
  
    SetDim();
    if(fSystem == pp)  fUseCentrality = kFALSE; else fUseCentrality = kTRUE;
    
    if(fSystem == AA) fBkgMethod = kDStarSB; else fBkgMethod = kDZeroSB;
  
    fCuts=cuts;
    fNofPtBins= fCuts->GetNPtBins();
    //cout << "Enlarging the DZero window " << endl;
    EnlargeDZeroMassWindow();
   // cout << "Done" << endl;
    
   
  DefineInput(0, TChain::Class());
  
  DefineOutput(1,TList::Class()); // histos from data and MC
  DefineOutput(2,TList::Class()); // histos from MC
    DefineOutput(3,TList::Class()); // histos from data and MC
    DefineOutput(4,TList::Class()); // histos from MC
    DefineOutput(5,TList::Class()); // histos from data and MC
      DefineOutput(6,TList::Class()); // histos from data and MC
    DefineOutput(7,AliNormalizationCounter::Class());   // normalization

  DefineOutput(8,AliRDHFCutsDStartoKpipi::Class()); // my D meson cuts
  DefineOutput(9,AliHFAssociatedTrackCuts::Class()); // my associated tracks cuts
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
	if(fAssocCuts) {delete fAssocCuts; fAssocCuts=0;}
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
	PostData(8,copyfCuts);
	
	// Post the hadron cuts
	PostData(9,fAssocCuts);
    
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
    
    fDmesonOutput = new TList();
	fDmesonOutput->SetOwner();
    
    fTracksOutput = new TList();
	fTracksOutput->SetOwner();
    
    fEMOutput = new TList();
	fEMOutput->SetOwner();
    
    fCorrelationOutput = new TList();
	fCorrelationOutput->SetOwner();
	
	// define histograms
	DefineHistoForAnalysis();
    DefineThNSparseForAnalysis();
    

    

	fCounter = new AliNormalizationCounter(Form("%s",GetOutputSlot(6)->GetContainer()->GetName()));
	fCounter->Init();
	
    Double_t Pi = TMath::Pi();
    //AliHFCorrelator(const Char_t* name, AliHFAssociatedTrackCuts *cuts, Bool_t useCentrality, AliRDHFCuts *cutObject)
	fCorrelator = new AliHFCorrelator("Correlator",fAssocCuts,fUseCentrality,fCuts); // fAssocCuts is the hadron cut object, fSystem to switch between pp or PbPb
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
	PostData(2,fDmesonOutput); // set the outputs
    PostData(3,fTracksOutput); // set the outputs
    PostData(4,fEMOutput); // set the outputs
    PostData(5,fCorrelationOutput); // set the outputs
    PostData(6,fOutputMC); // set the outputs
	PostData(7,fCounter); // set the outputs
}

//________________________________________  user exec ____________
void AliAnalysisTaskDStarCorrelations::UserExec(Option_t *){
    // cout << "Task debug check 1 " << endl;
     // ********************************************** LOAD THE EVENT ****************************************************
    
    if (!fInputEvent) {
        Error("UserExec","NO EVENT FOUND!");
		return;
    }
    
    AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
    if(!aodEvent){
        AliError("AOD event not found!");
        return;
    }
    
    fEvents++; // event counter
    ((TH1D*)fOutput->FindObject("NofEvents"))->Fill(0);
    
    fCounter->StoreEvent(aodEvent,fCuts,fmontecarlo); // store event in AliNormalizationCounter
    
    // load MC array
    fmcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if(fmontecarlo && !fmcArray){
        AliError("Array of MC particles not found");
        return;
    }
    
    // ********************************************** START EVENT SELECTION ****************************************************
    
    Bool_t isEvSel=fCuts->IsEventSelected(aodEvent);
    
    if(!isEvSel) {
        
        if(fCuts->IsEventRejectedDueToPileup()) {((TH1D*)fOutput->FindObject("NofEvents"))->Fill(2);                        cout << "Reject PILEUP" << endl;}
        if(fCuts->IsEventRejectedDueToCentrality()) {((TH1D*)fOutput->FindObject("NofEvents"))->Fill(3);                    cout << "Reject CENTRALITY" << endl;}
        if(fCuts->IsEventRejectedDueToNotRecoVertex()) {((TH1D*)fOutput->FindObject("NofEvents"))->Fill(4);                 cout << "Reject NOT RECO VTX" << endl;}
        if(fCuts->IsEventRejectedDueToVertexContributors()) {((TH1D*)fOutput->FindObject("NofEvents"))->Fill(5);            cout << "Reject VTX CONTRIB" << endl;}
        if(fCuts->IsEventRejectedDueToZVertexOutsideFiducialRegion()) {((TH1D*)fOutput->FindObject("NofEvents"))->Fill(6);  cout << "Reject VTX no fid reg " << endl;}
        if(fCuts->IsEventRejectedDueToTrigger()) {((TH1D*)fOutput->FindObject("NofEvents"))->Fill(7);                       cout << "Reject TRIGGER" << endl;}
        if(fCuts->IsEventRejectedDuePhysicsSelection()) {((TH1D*)fOutput->FindObject("NofEvents"))->Fill(8);                cout << "Reject PHYS SEL" << endl;}
        
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

    fCorrelator->SetAODEvent(aodEvent); // set the event in the correlator class
    
    Bool_t correlatorON = fCorrelator->Initialize(); //define the pool for mixing and check if event within the pool definition
	if(!correlatorON) {
        AliInfo("AliHFCorrelator didn't initialize the pool correctly or processed a bad event");
        return; // if not the case, skips the event
	}
    
    ((TH1D*)fOutput->FindObject("NofEvents"))->Fill(1); // count events that are passing the event selection
    
    
    
    // cout << "Task debug check 2 " << endl;
     // ********************************************** CENTRALITY/MULTIPLICITY  ****************************************************


    Double_t MultipOrCent = fCorrelator->GetCentrality(); // returns centrality or multiplicity
    
    
    
    // ********************************************** MC SETUP  ****************************************************
    
    if(fmontecarlo) {
        
        fCorrelator->SetMCArray(fmcArray); // export MC array into correlatior class
        
        AliAODMCHeader *mcHeader = dynamic_cast<AliAODMCHeader*>(aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
        if (fmontecarlo && !mcHeader) {
            AliError("Could not find MC Header in AOD");
            return;
        }
        
      // select MC event type (PP, GS etc) - those are defined in the associated tracks cut object
        Bool_t isMCeventgood = kFALSE;
        
        
        Int_t eventType = mcHeader->GetEventType();
        Int_t NMCevents = fAssocCuts->GetNofMCEventType();
        
        for(Int_t k=0; k<NMCevents; k++){
            Int_t * MCEventType = fAssocCuts->GetMCEventType();
            
            if(eventType == MCEventType[k]) isMCeventgood= kTRUE;
            ((TH1D*)fOutputMC->FindObject("EventTypeMC"))->Fill(eventType);
        }
        
        if(NMCevents && !isMCeventgood){
            if(fDebugLevel)	std::cout << "The MC event " << eventType << " not interesting for this analysis: skipping" << std::endl;
            return;
        }
        
        
    }// end if fmontecarlo
    
    
     // ********************************************** EVENT PROPERTIES  ****************************************************
    // checks on vertex and multiplicity of the event
	AliAODVertex *vtx = aodEvent->GetPrimaryVertex();
	Double_t zVtxPosition = vtx->GetZ(); // zvertex
	
 
  
    if(!aodEvent->GetPrimaryVertex() || TMath::Abs(aodEvent->GetMagneticField())<0.001) return;
    
    if(!fmixing) ((TH2F*)fOutput->FindObject("EventPropsCheck"))->Fill(MultipOrCent,zVtxPosition);
    if(!fmixing) ((TH1D*)fOutput->FindObject("MultiplicityOnlyCheck"))->Fill(AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aodEvent,-1.,1.));
    
    
   // Int_t poolbin = fAssocCuts->GetPoolBin(MultipOrCent, zVtxPosition); // get the event pool bin - commented for the moment - to be checked
    
    
    
     // ********************************************** D * selection  ****************************************************
    
    
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
    
    
    
    // D* related variables
    
    Double_t ptDStar;   // D* pt
	Double_t phiDStar;  // D* phi
	Double_t etaDStar;  // D* eta
	Bool_t isInPeak, isInDZeroSideBand, isInDStarSideBand, isDStarMCtag; // flags for selection
	Double_t invMassDZero; // D0 inv mass
	Double_t deltainvMDStar; // D* - D0 invarian mass
    
    Double_t mPDGD0=1.8648;//TDatabasePDG::Instance()->GetParticle(421)->Mass();
	Double_t mPDGDstar=2.01022;//TDatabasePDG::Instance()->GetParticle(413)->Mass();
	
	
	//MC tagging for DStar
	//D* and D0 prongs needed to MatchToMC method
	Int_t pdgDgDStartoD0pi[2]={421,211};
	Int_t pdgDgD0toKpi[2]={321,211};
	
	//Bool_t isDStarCand = kFALSE;
    Bool_t isDfromB = kFALSE;
	Bool_t isEventMixingFilledPeak = kFALSE;
	Bool_t isEventMixingFilledSB = kFALSE;
    Bool_t EventHasDStarCandidate = kFALSE;
    Bool_t EventHasDZeroSideBandCandidate = kFALSE;
    Bool_t EventHasDStarSideBandCandidate = kFALSE;
    Bool_t isTrackForPeakFilled = kFALSE;
    Bool_t isTrackForSBFilled = kFALSE;
	//loop on D* candidates
	
	Int_t looponDCands = 0;
	if(fReco) looponDCands = arrayDStartoD0pi->GetEntriesFast(); // load N of D* candidates from AOD array
	if(!fReco) looponDCands = fmcArray->GetEntriesFast(); // load N of D* candidates from MC array
	
	Int_t nOfDStarCandidates = 0;  // D candidates counter
	Int_t nOfSBCandidates = 0;     // SB candidates counter
	
	Double_t DmesonEfficiency = 1.; // Efficiency of D meson
	Double_t DmesonWeight = 1.;    // Applyed weight
    Double_t efficiencyvariable = -999; // Variable to be used (defined by the AddTask)
    Double_t dmDStarWindow = 0.0019/3; // DStar window
    
    // cout << "Task debug check 3 " << endl;
    // loop over D meson candidates
    for (Int_t iDStartoD0pi = 0; iDStartoD0pi<looponDCands; iDStartoD0pi++) {
    
        // initialize variables
        isInPeak = kFALSE;
        isInDStarSideBand = kFALSE;
        isInDZeroSideBand = kFALSE;
        
        EventHasDStarCandidate = kFALSE;
        EventHasDZeroSideBandCandidate = kFALSE;
        EventHasDStarSideBandCandidate = kFALSE;
     
        
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
        Int_t ptbin = -1;
        
        // ============================================ using reconstruction on Data or MC
        if(fReco){
            // cout << "Task debug check 4 " << endl;
            // get the D objects
            dstarD0pi = (AliAODRecoCascadeHF*)arrayDStartoD0pi->At(iDStartoD0pi);
            if(!dstarD0pi->GetSecondaryVtx()) continue;
            theD0particle = (AliAODRecoDecayHF2Prong*)dstarD0pi->Get2Prong();
            if (!theD0particle) continue;
            
            
            // apply topological cuts
            
            // track quality cuts
            Int_t isTkSelected = fCuts->IsSelected(dstarD0pi,AliRDHFCuts::kTracks); // quality cuts on tracks
            // region of interest + topological cuts + PID
            Int_t isSelected=fCuts->IsSelected(dstarD0pi,AliRDHFCuts::kCandidate); //selected
            
            //apply topological cuts
            if(!isTkSelected) continue;
            if(!isSelected) continue;
            if(!fCuts->IsInFiducialAcceptance(dstarD0pi->Pt(),dstarD0pi->YDstar())) continue;
            
            // get D candidate variables
            ptDStar = dstarD0pi->Pt();
            phiDStar = dstarD0pi->Phi();
            etaDStar = dstarD0pi->Eta();
            if(TMath::Abs(etaDStar) > fMaxEtaDStar) continue;
            if(fEfficiencyVariable == kMult || fEfficiencyVariable == kCentr)  efficiencyvariable = MultipOrCent;
            if(fEfficiencyVariable == kEta) efficiencyvariable = etaDStar;
            if(fEfficiencyVariable == kRapidity) efficiencyvariable = dstarD0pi->YDstar();
            if(fEfficiencyVariable == kNone) efficiencyvariable = 0;
           
            
             if(TMath::Abs(etaDStar) > fMaxEtaDStar) continue; // skip candidates outside the defined eta range
            
            phiDStar = fCorrelator->SetCorrectPhiRange(phiDStar); // set the Phi of the D* in the range defined a priori (-0.5 Pi - 1.5 Pi)
            ptbin=fCuts->PtBin(dstarD0pi->Pt()); // get the pt bin of the D*
            
            cout << "DStar pt = " << ptDStar << endl;
             cout << "pt bin = " << ptbin << endl;
            if(ptbin<1) continue;
            
             Double_t mD0Window= fD0Window[ptbin]/3;
             invMassDZero = dstarD0pi->InvMassD0();
             deltainvMDStar = dstarD0pi->DeltaInvMass();
            
            
            // get the D meson efficiency
            DmesonEfficiency = fAssocCuts->GetTrigWeight(dstarD0pi->Pt(),efficiencyvariable);
            
            // check this!
            if(fUseDmesonEfficiencyCorrection){
                if(DmesonEfficiency>1.e-5) DmesonWeight = 1./DmesonEfficiency;
                else {// THIS ELSE STATEMENT MUST BE REFINED: THE EFFICIENCY MAP HAS TO BE REPLACED WITH A WEIGHT MAP COOKED A PRIORI
                    if(ptDStar>2.) DmesonWeight = 0.5; // at high pt a zero value in the efficiency can come only from stat fluctutations in MC -> 0.5 is an arbitrary asymptotic value
                    else DmesonWeight = 1.e+5; // at low pt it can be that the efficiency is really low
                }
            }
            else DmesonWeight = 1.;
            
         
            // using montecarlo on reconstruction
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
            
            
            // fill mass histograms
            
            // cout << "Task debug check 5 " << endl;
            // fill D0 invariant mass
            if(!fmixing) {
                cout << " histo name = " << Form("histDZeroMass_%d",ptbin) << endl;
                ((TH1F*)fDmesonOutput->FindObject(Form("histDZeroMass_%d",ptbin)))->Fill(invMassDZero);
                // cout << "Task debug check 5.1 " << endl;
                if(fUseDmesonEfficiencyCorrection) ((TH1F*)fDmesonOutput->FindObject(Form("histDZeroMassWeight_%d",ptbin)))->Fill(invMassDZero,DmesonWeight);
            } // end if !mixing
            
        
            
            
            
            // good D0 canidates
            if (TMath::Abs(invMassDZero-mPDGD0)<fDMesonSigmas[1]*mD0Window){
                // fill D* invariant mass
                if(!fmixing){ ((TH1F*)fDmesonOutput->FindObject(Form("histDStarMass_%d",ptbin)))->Fill(deltainvMDStar);
                   // fill D* invariant mass if weighting
                    if(fUseDmesonEfficiencyCorrection) ((TH1F*)fDmesonOutput->FindObject(Form("histDStarMassWeight_%d",ptbin)))->Fill(deltainvMDStar,DmesonWeight);} // end if no mixing
                    isInPeak = kTRUE;
                    // fill other good candidate distributions - pt, phi eta
                    if(TMath::Abs(deltainvMDStar-(mPDGDstar-mPDGD0))<fDMesonSigmas[0]*dmDStarWindow)	{
                        ((TH1F*)fDmesonOutput->FindObject("RecoPtDStar"))->Fill(ptDStar,DmesonWeight); // fill pt
                        if(!fmixing)	((TH2F*)fDmesonOutput->FindObject("PhiInclusiveDStar"))->Fill(phiDStar,ptDStar); // fill phi, eta
                        if(!fmixing)	((TH2F*)fDmesonOutput->FindObject("EtaInclusiveDStar"))->Fill(etaDStar,ptDStar); // fill phi, eta
                        nOfDStarCandidates++;
                        EventHasDStarCandidate=kTRUE;
                    } // end if in good D* mass window
                    
                    // count D* sideband candidates
                    if(fBkgMethod == kDStarSB ){
                        if ((deltainvMDStar-(mPDGDstar-mPDGD0))>fDMesonSigmas[2]*dmDStarWindow && (deltainvMDStar-(mPDGDstar-mPDGD0))<fDMesonSigmas[3]*dmDStarWindow ){
                            ((TH1F*)fDmesonOutput->FindObject("RecoPtBkg"))->Fill(ptDStar,DmesonWeight); // fill pt
                            if(!fmixing)	((TH2F*)fDmesonOutput->FindObject("PhiSidebandDStar"))->Fill(phiDStar,ptDStar); // fill phi, eta
                            if(!fmixing)	((TH2F*)fDmesonOutput->FindObject("EtaSidebandDStar"))->Fill(etaDStar,ptDStar); // fill phi, eta
                            nOfSBCandidates++;
                            EventHasDStarSideBandCandidate = kTRUE;
                            }
                            
                        } // end if using DStar sidebans
                
                
            }// end good D0
           
            // cout << "Task debug check 6 " << endl;
            // Sideband D0
              if (TMath::Abs(invMassDZero-mPDGD0)>fDMesonSigmas[2]*mD0Window && TMath::Abs(invMassDZero-mPDGD0)<fDMesonSigmas[3]*mD0Window ){
                isInDZeroSideBand = kTRUE;
                 if(!fmixing){ ((TH1F*)fDmesonOutput->FindObject(Form("histDStarFromSBMass_%d",ptbin)))->Fill(deltainvMDStar);
                     if(fUseDmesonEfficiencyCorrection) ((TH1F*)fDmesonOutput->FindObject(Form("histDStarFromSBMassWeight_%d",ptbin)))->Fill(deltainvMDStar,DmesonWeight);
                     
                     if(fBkgMethod == kDZeroSB){
                         if(TMath::Abs(deltainvMDStar-(mPDGDstar-mPDGD0))<fDMesonSigmas[0]*dmDStarWindow)	{
                         
                             ((TH1F*)fDmesonOutput->FindObject("RecoPtBkg"))->Fill(ptDStar,DmesonWeight); // fill pt
                             if(!fmixing)	((TH2F*)fDmesonOutput->FindObject("PhiSidebandDStar"))->Fill(phiDStar,ptDStar); // fill phi, eta
                             if(!fmixing)	((TH2F*)fDmesonOutput->FindObject("EtaSidebandDStar"))->Fill(etaDStar,ptDStar); // fill phi, eta
                             nOfSBCandidates++;
                             EventHasDZeroSideBandCandidate = kTRUE;
                         }
                     }
                 }
            }// end SideBand D0
            // cout << "Task debug check 7 " << endl;
            
            if(! isInPeak && !isInDZeroSideBand) continue; // skip candidates that are not interesting
            if(TMath::Abs(deltainvMDStar)<0.142 || TMath::Abs(deltainvMDStar)>0.151) continue; // skip all D* canidates outside the interesting mass ranges
           // cout << "Good D* candidate" << endl;

            // get charge of soft pion
            daughtercharge = ((AliAODTrack*)dstarD0pi->GetBachelor())->Charge();
            
            // get ID of daughters used for reconstruction
            trackiddaugh0 = ((AliAODTrack*)theD0particle->GetDaughter(0))->GetID();
            trackiddaugh1 = ((AliAODTrack*)theD0particle->GetDaughter(1))->GetID();
            trackidsoftPi = ((AliAODTrack*)dstarD0pi->GetBachelor())->GetID();
        }// end if reconstruction
        
        
        
          // ============================================ using MC kinematics only
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
            ptbin=fCuts->PtBin(DStarMC->Pt());
            
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
                
                for(Int_t i=0; i<NDaugh;i++){ // loop on DStar daughters
                    Int_t daugh_label = DStarMC->GetDaughter(i);
                    AliAODMCParticle* mcDaughter = dynamic_cast<AliAODMCParticle*>(fmcArray->At(daugh_label));
                    if(!mcDaughter) continue;
                    Int_t daugh_pdg = TMath::Abs(mcDaughter->GetPdgCode());
                    if(fDebugLevel) std::cout << "Daughter " << i << " pdg code is " << daugh_pdg << std::endl;
                    
                    if(daugh_pdg == 421) {
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
                        if(!isKaon || !isPion) break; // skip if not correct decay channel of D0
                        isDZero = kTRUE;
                    } // end check on Dzero
                    
                    if(daugh_pdg == 211) {
                        isSoftPi = kTRUE;
                        daughtercharge = mcDaughter->Charge();
                        trackidsoftPi = daugh_label;}
                } // end loop on D* daughters
                if(!isDZero || !isSoftPi) continue; // skip if not correct decay channel
            } // end check decay channel
            
            ptDStar = DStarMC->Pt();
            phiDStar = DStarMC->Phi();
            etaDStar = DStarMC->Eta();
            
            if(TMath::Abs(etaDStar) > fMaxEtaDStar) continue;
            
            
        }// end if pure MC information
        
        
        
        
        
        // getting the number of triggers in the MCtag D* case
        if(fmontecarlo && isDStarMCtag) ((TH1F*)fOutputMC->FindObject("MCtagPtDStar"))->Fill(ptDStar);
        if(fmontecarlo && isDStarMCtag && !isDfromB) ((TH1D*)fOutputMC->FindObject("MCtagPtDStarfromCharm"))->Fill(ptDStar);
        if(fmontecarlo && isDStarMCtag && isDfromB) ((TH1D*)fOutputMC->FindObject("MCtagPtDStarfromBeauty"))->Fill(ptDStar);
        
        
        fCorrelator->SetTriggerParticleProperties(ptDStar,phiDStar,etaDStar); // pass to the object the necessary trigger part parameters
        fCorrelator->SetTriggerParticleDaughterCharge(daughtercharge);
        
        
        
        
         // ************************************************ CORRELATION ANALYSIS STARTS HERE *****************************************************
        
        cout << "Correlating " << endl;
        
        Bool_t execPool = fCorrelator->ProcessEventPool(); // checks pool readiness for mixed events
        
        // if analysis is on mixed event and pool settings are not satisfied, fill relevant histograms and skip
        if(fmixing && !execPool) {
            AliInfo("Mixed event analysis: pool is not ready");
            if(!isEventMixingFilledPeak && isInPeak)  {
                ((TH1D*)fEMOutput->FindObject("CheckPoolReadiness"))->Fill(1);
                isEventMixingFilledPeak = kTRUE;
            }
            if (!isEventMixingFilledSB && isInDZeroSideBand)  {
                ((TH1D*)fEMOutput->FindObject("CheckPoolReadiness"))->Fill(3);
                isEventMixingFilledSB=kTRUE;
            }
            continue;
        } // end if pool not ok
        // if analysis is on mixed event and pool settings are  satisfied, fill relevant histograms and continue
        if(fmixing&&execPool){
            // pool is ready - run checks on bins filling
            if(!isEventMixingFilledPeak && isInPeak)  {
                ((TH1D*)fEMOutput->FindObject("CheckPoolReadiness"))->Fill(0);
                if(fFullmode) EventMixingChecks(aodEvent);
                isEventMixingFilledPeak = kTRUE;
            }
            
            if(!isEventMixingFilledSB && isInDZeroSideBand) {
                ((TH1D*)fEMOutput->FindObject("CheckPoolReadiness"))->Fill(2);
                isEventMixingFilledSB=kTRUE;
            }
        } // end if pool ok
        
        
        
        
        Int_t NofEventsinPool = 1;
        if(fmixing) NofEventsinPool = fCorrelator->GetNofEventsInPool();
        
        //************************************************** LOOP ON EVENTS IN EVENT POOL *****************************************************
        
        for (Int_t jMix =0; jMix < NofEventsinPool; jMix++){// loop on events in the pool; if it is SE analysis, stops at one
            
            Bool_t analyzetracks = fCorrelator->ProcessAssociatedTracks(jMix); // process the associated tracks
            if(!analyzetracks) {
                AliInfo("AliHFCorrelator::Cannot process the track array");
                continue;
            }
            
            Double_t arraytofill[4];
          // Double_t MCarraytofill[4]; // think about this
            Double_t weight;
            
              Int_t NofTracks = fCorrelator->GetNofTracks(); // number of tracks in event
            
             //************************************************** LOOP ON TRACKS *****************************************************
            
             for(Int_t iTrack = 0; iTrack<NofTracks; iTrack++){ // looping on track candidates
                 
                 Bool_t runcorrelation = fCorrelator->Correlate(iTrack); // calculate correlations
                 if(!runcorrelation) continue;
                 
                 Double_t DeltaPhi = fCorrelator->GetDeltaPhi();
                 Double_t DeltaEta = fCorrelator->GetDeltaEta();
                 
                 AliReducedParticle * hadron = fCorrelator->GetAssociatedParticle();
                 if(!hadron) {/*cout << "No Hadron" << endl;*/ continue;}
                 
                 Double_t ptHad = hadron->Pt();
                 Double_t phiHad = hadron->Phi();
                    phiHad = fCorrelator->SetCorrectPhiRange(phiHad); // set phi in correct range
                 Double_t etaHad = hadron->Eta();
                 Int_t label = hadron->GetLabel();
                 Int_t trackid = hadron->GetID();
                 Double_t efficiency = hadron->GetWeight();
                 
                 
                 
                  
                 if(!isTrackForPeakFilled && !fmixing && EventHasDStarCandidate){
                     
                 	((TH2F*)fTracksOutput->FindObject("PhiInclusiveTracks"))->Fill(phiHad,ptHad); // fill phi, eta
                 	((TH2F*)fTracksOutput->FindObject("EtaInclusiveTracks"))->Fill(etaHad,ptHad); // fill phi, eta
                     isTrackForPeakFilled =  kTRUE; // makes sure you do not fill twice in case of more candidates
                 }
                 
                 if(!isTrackForSBFilled && !fmixing && (fBkgMethod == kDZeroSB) && EventHasDZeroSideBandCandidate){
                    ((TH2F*)fTracksOutput->FindObject("PhiSidebandTracks"))->Fill(phiHad,ptHad); // fill phi, eta
                    ((TH2F*)fTracksOutput->FindObject("EtaSidebandTracks"))->Fill(etaHad,ptHad); // fill phi, eta
                     isTrackForSBFilled = kTRUE;
                 }
                 
                 if(!isTrackForSBFilled && !fmixing && (fBkgMethod == kDStarSB) && EventHasDStarSideBandCandidate){
                     ((TH2F*)fTracksOutput->FindObject("PhiSidebandTracks"))->Fill(phiHad,ptHad); // fill phi, eta
                     ((TH2F*)fTracksOutput->FindObject("EtaSidebandTracks"))->Fill(etaHad,ptHad); // fill phi, eta
                     isTrackForSBFilled = kTRUE;
                 }
                 
                 
                 weight = 1;
                 if(fUseEfficiencyCorrection && efficiency){
                     weight = DmesonWeight * (1./efficiency);
                 }
                 
                 
                 arraytofill[0] = DeltaPhi;
                 arraytofill[1] = deltainvMDStar;
                 arraytofill[2] = DeltaEta;
                 arraytofill[3] = ptHad;
                // arraytofill[4] = poolbin;
                 
                 
                 // skip the D daughters in the correlation
              //   Bool_t isDdaughter = kFALSE;
                 if(!fmixing && fReco){ // for reconstruction
                     if(trackid == trackiddaugh0) continue;
                     if(trackid == trackiddaugh1) continue;
                     if(trackid == trackidsoftPi) continue;
                 }
                 
                 if(!fmixing && !fReco){ // for kinematic MC
                     AliAODMCParticle *part = (AliAODMCParticle*)fmcArray->At(label);
                     if(!part) continue;
                     if(IsDDaughter(DStarMC, part)) continue;
                     cout << "Skipping DStar  daugheter " << endl;
                 }
                 
                 
                 
                 // ============================================= FILL CORRELATION THNSparses ===============================
                 
                 // filling signal
                 if(isInPeak){
                     cout << "Filling signal " << endl;
                   //  if(!fReco && TMath::Abs(etaHad)>0.8) continue;
                     //cout ("CorrelationsDStarHadron_%d",ptbin)
                     if(fselect==1)((THnSparseF*)fCorrelationOutput->FindObject(Form("CorrelationsDStarHadron_%d",ptbin)))->Fill(arraytofill,weight);
                     if(fselect==2)((THnSparseF*)fCorrelationOutput->FindObject(Form("CorrelationsDStarKaon_%d",ptbin)))->Fill(arraytofill,weight);
                     if(fselect==3)((THnSparseF*)fCorrelationOutput->FindObject(Form("CorrelationsDStarKZero_%d",ptbin)))->Fill(arraytofill,weight);
                 }
                 
                  // filling bkg
                 if(fBkgMethod == kDStarSB && isInPeak){ // bkg from DStar
                   //  if(!fReco && TMath::Abs(etaHad)>0.8) continue;
                      cout << "Filling bkg from D* sidebands " << endl;
                     if(fselect==1)((THnSparseF*)fCorrelationOutput->FindObject(Form("CorrelationsBkgFromDStarSBHadron_%d",ptbin)))->Fill(arraytofill,weight);
                     if(fselect==2)((THnSparseF*)fCorrelationOutput->FindObject(Form("CorrelationsBkgFromDStarSBKaon_%d",ptbin)))->Fill(arraytofill,weight);
                     if(fselect==3)((THnSparseF*)fCorrelationOutput->FindObject(Form("CorrelationsBkgFromDStarSBKZero_%d",ptbin)))->Fill(arraytofill,weight);
                     
                 } // end if bkg from DStar
                 
                 if(fBkgMethod == kDZeroSB && isInDZeroSideBand){ // bkg from DStar
                   //  if(!fReco && TMath::Abs(etaHad)>0.8) continue;
                     cout << "Filling bkg from Dzero sidebands " << endl;
                     if(fselect==1)((THnSparseF*)fCorrelationOutput->FindObject(Form("CorrelationsBkgFromDZeroSBHadron_%d",ptbin)))->Fill(arraytofill,weight);
                     if(fselect==2)((THnSparseF*)fCorrelationOutput->FindObject(Form("CorrelationsBkgFromDZeroSBKaon_%d",ptbin)))->Fill(arraytofill,weight);
                     if(fselect==3)((THnSparseF*)fCorrelationOutput->FindObject(Form("CorrelationsBkgFromDZeroSBKZero_%d",ptbin)))->Fill(arraytofill,weight);
                     
                 } // end if bkg from DZero
                 
                                                                        
                 
             } // end loop on associated tracks
            
        } // end loop on events in event pool
        
        
    } // end loop on D mesons
    
    
    
    
     Bool_t updated = fCorrelator->PoolUpdate(); // update event pool
    
    
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
    
    Double_t Pi = TMath::Pi();
	Int_t nbinscorr = fPhiBins;
	Double_t lowcorrbin = -0.5*Pi ;
    Double_t upcorrbin = 1.5*Pi ;
    
    
    // create ThNSparses
    
    Int_t nofPtBins = fCuts->GetNPtBins();// number of ptbins
    
    
    //sparse bins
    
    //1 delta_phi
    //2 invariant mass D *
    //3 delta eta
    //4 track pt
    
    
    //Int_t nbinsPool = (fAssocCuts->GetNZvtxPoolBins())*(fAssocCuts->GetNCentPoolBins());
    
    
    // for reconstruction on Data and MC
    Int_t nbinsSparse[4]=         {nbinscorr ,     32 ,  32, 10};
    Double_t binLowLimitSparse[4]={lowcorrbin,0.14314 ,-1.6,  0};
    Double_t binUpLimitSparse[4]= {upcorrbin ,0.14794 , 1.6,  5};
    
    Int_t nbinsSparseDStarSB[4]=         {nbinscorr ,     40 ,  32, 10};
    Double_t binLowLimitSparseDStarSB[4]={lowcorrbin,0.14788 ,-1.6,  0};
    Double_t binUpLimitSparseDStarSB[4]= {upcorrbin ,0.1504 , 1.6,  5};
    
    TString signalSparseName = "";
    TString bkgSparseName = "";
    
    THnSparseF * CorrelationsSignal = NULL;
    THnSparseF * CorrelationsBkg = NULL;
    
    
    Float_t * ptbinlims = fCuts->GetPtBinLimits();
    
  
    
    
    for(Int_t iBin =0; iBin < nofPtBins; iBin++){ // create a mass histogram for each ptbin
        
        
        
        if(ptbinlims[iBin]<fCuts->GetMinPtCandidate() || ptbinlims[iBin]>fCuts->GetMaxPtCandidate())continue;


        
        signalSparseName = "CorrelationsDStar";
        if(fselect==1) signalSparseName += "Hadron_";
        if(fselect==2) signalSparseName += "Kaon_";
        if(fselect==3) signalSparseName += "KZero_";
        
        bkgSparseName = "CorrelationsBkg";
        if(fBkgMethod == kDStarSB ) bkgSparseName+="FromDStarSB";
        if(fBkgMethod == kDZeroSB ) bkgSparseName+="FromDZeroSB";
        if(fselect==1) bkgSparseName += "Hadron_";
        if(fselect==2) bkgSparseName += "Kaon_";
        if(fselect==3) bkgSparseName += "KZero_";
        
        signalSparseName+=Form("%d",iBin);
        bkgSparseName+=Form("%d",iBin);
        cout << "ThNSparses name = " << signalSparseName << endl;
        
        // define thnsparses for signal candidates
        CorrelationsSignal = new THnSparseF(signalSparseName.Data(),"Correlations for signal; #Delta#Phi; invariant mass;  #Delta #eta;AssocTrk p_{T}",4,nbinsSparse,binLowLimitSparse,binUpLimitSparse);
          CorrelationsSignal->Sumw2();
        fCorrelationOutput->Add(CorrelationsSignal);
        
        // define thnsparses for bkg candidates from DStar
        if(fBkgMethod == kDStarSB ){
        CorrelationsBkg = new THnSparseF(bkgSparseName.Data(),"Correlations for bkg from DStar; #Delta#Phi; invariant mass;  #Delta #eta;AssocTrk p_{T}",4,nbinsSparseDStarSB,binLowLimitSparseDStarSB,binUpLimitSparseDStarSB);
        CorrelationsBkg->Sumw2();
        fCorrelationOutput->Add(CorrelationsBkg);
        }
        
        // define thnsparses for bkg candidates from DZero
        if(fBkgMethod == kDZeroSB ){
            CorrelationsBkg = new THnSparseF(bkgSparseName.Data(),"Correlations for bkg from DZero; #Delta#Phi; invariant mass;  #Delta #eta;AssocTrk p_{T}",4,nbinsSparse,binLowLimitSparse,binUpLimitSparse);
            CorrelationsBkg->Sumw2();
            fCorrelationOutput->Add(CorrelationsBkg);
        }
        
    } // end loop on bins
    
    
    
}

//__________________________________________________________________________________________________
void AliAnalysisTaskDStarCorrelations::DefineHistoForAnalysis(){
    
    Double_t Pi = TMath::Pi();
	Int_t nbinscorr = fPhiBins;
	Double_t lowcorrbin = -0.5*Pi ;
    Double_t upcorrbin = 1.5*Pi ;
    
    // event counter
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
    
    //event properties
    TH2F * EventPropsCheck = new TH2F("EventPropsCheck","Properties of the event; Multiplicity/Centrality; ZVtx Position [cm]",1000,0,1000,40,-10,10);
    fOutput->Add(EventPropsCheck);
    
    //event properties
    TH1D * SPDMultiplicty = new TH1D("MultiplicityOnlyCheck","Properties of the event; SPD Multiplicity",1000,0,1000);
    fOutput->Add(SPDMultiplicty);
    
    
    // ===================================================   D* histograms
    TH1F * D0mass = NULL;
    TH1F * DStarMass = NULL;
    TH1F * DStarFromSBMass = NULL;
    
    TH1F * D0massWeighted = NULL;
    TH1F * DStarMassWeighted = NULL;
    TH1F * DStarFromSBMassWeighted = NULL;
    
    TString nameDZeroMass = "", nameDStarMass = "", nameDStarFromSBMass = "";
    
    Int_t nofPtBins = fCuts->GetNPtBins();// number of ptbins
    Float_t * ptbinlims = fCuts->GetPtBinLimits();
    
    //GetMinPtCandidate()

    
    for(Int_t iBin =0; iBin < nofPtBins; iBin++){ // create a mass histogram for each ptbin
        
      
        
        if(ptbinlims[iBin]<fCuts->GetMinPtCandidate() || ptbinlims[iBin]>fCuts->GetMaxPtCandidate())continue;
        
        
        std::cout << ">>> Ptbin = " << iBin << " limit = " << ptbinlims[iBin] << std::endl;
        
        nameDZeroMass = "histDZeroMass_";
        nameDStarMass = "histDStarMass_";
        nameDStarFromSBMass = "histDStarFromSBMass_";
        
        nameDZeroMass+=Form("%d",iBin);
        nameDStarMass+=Form("%d",iBin);
        nameDStarFromSBMass+=Form("%d",iBin);
        
        cout << "D zero histogram: " << nameDZeroMass << endl;
        
        D0mass = new TH1F(nameDZeroMass.Data(), Form("D^{0} invarians mass in bin %d; M(K#pi) GeV/c^{2};",iBin),200,1.75,1.95);
        DStarMass = new TH1F(nameDStarMass.Data(), Form("Delta invarians mass for candidates in bin %d; M(K#pi) GeV/c^{2};",iBin),200,0.1,0.2);
        DStarFromSBMass = new TH1F(nameDStarFromSBMass.Data(), Form("Delta invarians mass for sideband in bin %d; M(K#pi) GeV/c^{2};",iBin),200,0.1,0.2);
        
        if(!fmixing){
        fDmesonOutput->Add(D0mass);
        fDmesonOutput->Add(DStarMass);
        fDmesonOutput->Add(DStarFromSBMass);
        }
        
        // if using D meson efficiency, define weighted histos
        if(fUseDmesonEfficiencyCorrection){
           
            nameDZeroMass = "histDZeroMassWeight_";
            nameDStarMass = "histDStarMassWeight_";
            nameDStarFromSBMass = "histDStarFromSBMassWeight_";
            
            nameDZeroMass+=Form("%d",iBin);
            nameDStarMass+=Form("%d",iBin);
            nameDStarFromSBMass+=Form("%d",iBin);
            
            D0massWeighted = new TH1F(nameDZeroMass.Data(), Form("D^{0} invarians mass in bin %d eff weight; M(K#pi) GeV/c^{2};",iBin),200,1.75,1.95);
            DStarMassWeighted = new TH1F(nameDStarMass.Data(), Form("Delta invarians mass for candidates in bin %d eff weight; M(K#pi) GeV/c^{2};",iBin),200,0.1,0.2);
            DStarFromSBMassWeighted = new TH1F(nameDStarFromSBMass.Data(), Form("Delta invarians mass for sideband in bin %d eff weight; M(K#pi) GeV/c^{2};",iBin),200,0.1,0.2);
            
            if(!fmixing){
            fDmesonOutput->Add(D0massWeighted);
            fDmesonOutput->Add(DStarMassWeighted);
            fDmesonOutput->Add(DStarFromSBMassWeighted);
            }
        }
    }// end loop on pt bins
    
    
    // pt distributions
    TH1F *RecoPtDStar = new TH1F("RecoPtDStar","RECO DStar pt distribution",60,0,60);
    fDmesonOutput->Add(RecoPtDStar);
	TH1F *RecoPtBkg = new TH1F("RecoPtBkg","RECO pt distribution side bands",60,0,60);
    fDmesonOutput->Add(RecoPtBkg);
    
    // phi distribution
    TH2F * PhiInclusiveDStar = new TH2F("PhiInclusiveDStar","Azimuthal distributions of Inclusive Dmesons; #phi; pT;Entries",nbinscorr,lowcorrbin,upcorrbin,30,0,30);
    TH2F * PhiSidebandDStar = new TH2F("PhiSidebandDStar","Azimuthal distributions of Sideband Dmesons; #phi; pT;Entries",nbinscorr,lowcorrbin,upcorrbin,30,0,30);
    
    // eta distribution
    TH2F * EtaInclusiveDStar = new TH2F("EtaInclusiveDStar","Azimuthal distributions of Inclusive Dmesons; #eta; pT;Entries",20,-1,1,30,0,30);
    TH2F * EtaSidebandDStar = new TH2F("EtaSidebandDStar","Azimuthal distributions of Sideband Dmesons; #eta; pT;Entries",20,-1,1,30,0,30);
    
    if(!fmixing) fDmesonOutput->Add(PhiInclusiveDStar);
    if(!fmixing) fDmesonOutput->Add(PhiSidebandDStar);
    if(!fmixing) fDmesonOutput->Add(EtaInclusiveDStar);
    if(!fmixing) fDmesonOutput->Add(EtaSidebandDStar);
    
    
    // single track related histograms
    // phi distribution
    TH2F * PhiInclusiveTracks = new TH2F("PhiInclusiveTracks","Azimuthal distributions tracks if Inclusive Dmesons; #phi; pT;Entries",nbinscorr,lowcorrbin,upcorrbin,100,0,10);
    TH2F * PhiSidebandTracks = new TH2F("PhiSidebandTracks","Azimuthal distributions tracks if Sideband Dmesons; #phi; pT;Entries",nbinscorr,lowcorrbin,upcorrbin,100,0,10);
    
    // eta distribution
    TH2F * EtaInclusiveTracks = new TH2F("EtaInclusiveTracks","Azimuthal distributions of tracks if Inclusive Dmesons; #eta; pT;Entries",20,-1,1,100,0,10);
    TH2F * EtaSidebandTracks = new TH2F("EtaSidebandTracks","Azimuthal distributions of tracks if Sideband Dmesons; #eta; pT;Entries",20,-1,1,100,0,10);
    
    if(!fmixing) fTracksOutput->Add(PhiInclusiveTracks);
    if(!fmixing) fTracksOutput->Add(PhiSidebandTracks);
    if(!fmixing) fTracksOutput->Add(EtaInclusiveTracks);
    if(!fmixing) fTracksOutput->Add(EtaSidebandTracks);

    
    // Montecarlo for D*
    TH1D *MCtagPtDStarfromCharm = new TH1D("MCtagPtDStarfromCharm","RECO pt of MCtagged DStars from charm",50,0,50);
    if(fmontecarlo) fOutputMC->Add(MCtagPtDStarfromCharm);
    
    TH1D *MCtagPtDStarfromBeauty = new TH1D("MCtagPtDStarfromBeauty","RECO pt of MCtagged DStars from beauty",50,0,50);
    if(fmontecarlo) fOutputMC->Add(MCtagPtDStarfromBeauty);
	
	TH1F *MCtagPtDStar = new TH1F("MCtagPtDStar","RECO pt of MCtagged DStars side bands",50,0,50);
	if(fmontecarlo) fOutputMC->Add(MCtagPtDStar);
    
    
    // event mixing histograms
    TH1D * CheckPoolReadiness = new TH1D("CheckPoolReadiness","Pool readiness",5,-0.5,4.5);
    CheckPoolReadiness->GetXaxis()->SetBinLabel(1,"Have a D cand, pool is ready");
	CheckPoolReadiness->GetXaxis()->SetBinLabel(2,"Have a D cand, pool is not ready");
    CheckPoolReadiness->GetXaxis()->SetBinLabel(3,"Have a SB cand, pool is ready");
	CheckPoolReadiness->GetXaxis()->SetBinLabel(4,"Have a SB cand, pool is not ready");
	
    if(fmixing) fEMOutput->Add(CheckPoolReadiness);
    
    
   /* Int_t NofCentBins = fAssocCuts->GetNCentPoolBins();
    Int_t NofZVrtxBins = fAssocCuts->GetNZvtxPoolBins();
    Int_t nPoolBins = NofCentBins*NofZVrtxBins;
	
    
    TH1D * PoolBinDistribution  = new TH1D("PoolBinDistribution","Pool Bin Checks; PoolBin; Entry",nPoolBins5,-0.5,nPoolBins-0.5);
    fEMOutput->Add(PoolBinDistribution);
    
    TH2D * EventDistributionPerPoolBin  = new TH2D("EventDistributionPerPoolBin","Pool Bin Checks; PoolBin; Entry",nPoolBins5,-0.5,nPoolBins-0.5);
    fEMOutput->Add(EventDistributionPerPoolBin);
    */
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
	// check this function
    
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
	




