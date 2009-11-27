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
//                  Base class for DStar in Jets Analysis
//
//  The D* (+ and -) is reconstructed inside jets. Two different cuts are 
//  implemented:
//
//  1) C.Ivan D* cuts adapted for correlation analysis 
//  2) Topological cut enforcing the correlation D0-softPion pt + relaxed 
//     CosThetaStar. This second should be better for correlations.
//
//  USAGE:
//
//  The analysis is performed separately for D*+ and D*-. A flag in the .C is
//  taking care to decide which analysis.
//
//  The cut number 2 can be activaded with a flag in the .C (not active in this version)
//  Cuts 1) are the default. 
//
//  The task requires reconstructed jets in the AODs
//
//-----------------------------------------------------------------------
//                         Author A.Grelli 
//              ERC-QGP Utrecht University - a.grelli@uu.nl
//-----------------------------------------------------------------------

#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TVector3.h>

#include "AliAnalysisTaskSEDStarJets.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliLog.h"
#include "AliAODVertex.h"
#include "AliAODJet.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliESDtrack.h"
#include "AliAODMCParticle.h"

ClassImp(AliAnalysisTaskSEDStarJets)

//__________________________________________________________________________
AliAnalysisTaskSEDStarJets::AliAnalysisTaskSEDStarJets() :
  AliAnalysisTaskSE(),
  fCountMC(0),
  fCountAcc(0),
  fCountReco(0),
  fCountRecoAcc(0),
  fCountRecoITSClusters(0),
  fCountRecoPPR(0),
  fCountDStar(0),
  fCountDStarMC(0),
  fEvents(0),
  fMinITSClusters(0),
  fComputeD0(kTRUE),
  ftopologicalCut(kFALSE), 
  fRequireNormalization(kTRUE),
  fLorentzTrack1(0,0,0,0),
  fLorentzTrack2(0,0,0,0),
  fLorentzTrack3(0,0,0,0),
  fLorentzTrack4(0,0,0,0),
  fOutput(0),
  fD0ptvsSoftPtSignal(0),    
  fD0ptvsSoftPt(0),          
  ftrigger(0),   
  fPtPion(0),        
  fInvMass(0),       
  fRECOPtDStar(0),    
  fDStar(0),          
  fDiff(0),           
  fDiffSideBand(0),  
  fDStarMass(0),    
  fPhi(0),       
  fPhiBkg(0),        
  fTrueDiff(0),       
  fResZ(0),        
  fResZBkg(0),       
  fcharmpt(0),        
  fdstarE(0),        
  fEjet(0),        
  fPhijet(0),        
  fEtaJet(0),         
  fdstarpt(0),        
  fMCPionPt(0)    
{
  //
  // Default ctor
  //
}
//___________________________________________________________________________
AliAnalysisTaskSEDStarJets::AliAnalysisTaskSEDStarJets(const Char_t* name) :
  AliAnalysisTaskSE(name),
  fCountMC(0),
  fCountAcc(0),
  fCountReco(0),
  fCountRecoAcc(0),
  fCountRecoITSClusters(0),
  fCountRecoPPR(0),
  fCountDStar(0),
  fCountDStarMC(0),
  fEvents(0),
  fMinITSClusters(0),
  fComputeD0(kTRUE),
  ftopologicalCut(kFALSE),
  fRequireNormalization(kTRUE),
  fLorentzTrack1(0,0,0,0),
  fLorentzTrack2(0,0,0,0),
  fLorentzTrack3(0,0,0,0),
  fLorentzTrack4(0,0,0,0),
  fOutput(0),
  fD0ptvsSoftPtSignal(0),    
  fD0ptvsSoftPt(0),          
  ftrigger(0),   
  fPtPion(0),        
  fInvMass(0),       
  fRECOPtDStar(0),    
  fDStar(0),          
  fDiff(0),           
  fDiffSideBand(0),  
  fDStarMass(0),    
  fPhi(0),       
  fPhiBkg(0),        
  fTrueDiff(0),       
  fResZ(0),        
  fResZBkg(0),       
  fcharmpt(0),        
  fdstarE(0),        
  fEjet(0),        
  fPhijet(0),        
  fEtaJet(0),         
  fdstarpt(0),        
  fMCPionPt(0)    
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskSEDStarJets","Calling Constructor");
 
  DefineOutput(1,TList::Class());
}

//___________________________________________________________________________
AliAnalysisTaskSEDStarJets& AliAnalysisTaskSEDStarJets::operator=(const AliAnalysisTaskSEDStarJets& c) 
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c) ;
  }
 
  return *this;
}

//___________________________________________________________________________
AliAnalysisTaskSEDStarJets::AliAnalysisTaskSEDStarJets(const AliAnalysisTaskSEDStarJets& c) :
  AliAnalysisTaskSE(c),
  fCountMC(c.fCountMC),
  fCountAcc(c.fCountAcc),
  fCountReco(c.fCountReco),
  fCountRecoAcc(c.fCountRecoAcc),
  fCountRecoITSClusters(c.fCountRecoITSClusters),
  fCountRecoPPR(c.fCountRecoPPR),
  fCountDStar(c.fCountDStar),
  fCountDStarMC(c.fCountDStarMC),
  fEvents(c.fEvents),
  fMinITSClusters(c.fMinITSClusters),
  fComputeD0(c.fComputeD0),
  ftopologicalCut(c.ftopologicalCut),
  fRequireNormalization(c.fRequireNormalization),
  fLorentzTrack1(c.fLorentzTrack1),
  fLorentzTrack2(c.fLorentzTrack2),
  fLorentzTrack3(c.fLorentzTrack3),
  fLorentzTrack4(c.fLorentzTrack4),
  fOutput(c.fOutput),
  fD0ptvsSoftPtSignal(c.fD0ptvsSoftPtSignal),    
  fD0ptvsSoftPt(c.fD0ptvsSoftPt),          
  ftrigger(c.ftrigger),   
  fPtPion(c.fPtPion),        
  fInvMass(c.fInvMass),       
  fRECOPtDStar(c.fRECOPtDStar),    
  fDStar(c.fDStar),          
  fDiff(c.fDiff),           
  fDiffSideBand(c.fDiffSideBand),  
  fDStarMass(c.fDStarMass),    
  fPhi(c.fPhi),       
  fPhiBkg(c.fPhiBkg),        
  fTrueDiff(c.fTrueDiff),       
  fResZ(c.fResZ),        
  fResZBkg(c.fResZBkg),       
  fcharmpt(c.fcharmpt),        
  fdstarE(c.fdstarE),        
  fEjet(c.fEjet),        
  fPhijet(c.fPhijet),        
  fEtaJet(c.fEtaJet),         
  fdstarpt(c.fdstarpt),        
  fMCPionPt(c.fMCPionPt)       

{
  //
  // Copy Constructor
  //
}

//___________________________________________________________________________
AliAnalysisTaskSEDStarJets::~AliAnalysisTaskSEDStarJets() {
  // destructor
 
  Info("~AliAnalysisTaskSEDStarJets","Calling Destructor");  
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  } 
}

//_________________________________________________
void AliAnalysisTaskSEDStarJets::UserExec(Option_t *)
{
  // user exec
  if (!fInputEvent) {
    Error("UserExec","NO EVENT FOUND!");
    return;
  }
  
  // Load the event
  fEvents++;
  AliInfo(Form("Event %d",fEvents));
  if (fEvents%10000 ==0) AliInfo(Form("Event %d",fEvents));
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
  
  TClonesArray *arrayVerticesHF=0;
  
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
      arrayVerticesHF=(TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
    }
  } else {
    arrayVerticesHF=(TClonesArray*)aodEvent->GetList()->FindObject("D0toKpi");
  }
  
  if (!arrayVerticesHF){
    AliInfo("Could not find array of HF vertices, skipping the event");
    return;
  }else AliDebug(2, Form("Found %d vertices",arrayVerticesHF->GetEntriesFast()));   

  // Simulate a jet triggered sample
  TClonesArray *arrayofJets = (TClonesArray*)aodEvent->GetJets();
  if(aodEvent->GetNJets()<=0) return;
  AliInfo("found a jet: processing D* in jet analysis");
 
  // Get Prymary vertex --- be careful for lhc09a5 it has larger uncertanties
  
  AliAODVertex* prVtx = aodEvent->GetPrimaryVertex();
  Double_t primaryPos[3];
  prVtx->GetXYZ(primaryPos);

  //loop on the MC event - some basic MC info on D*, D0 and soft pion
  TClonesArray* mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!mcArray) AliError("Could not find Monte-Carlo in AOD");

  // counters for efficiencies
  Int_t icountMC   = 0;
  Int_t icountAcc  = 0;
  Int_t icountReco = 0;
  Int_t icountRecoAcc = 0;
  Int_t icountRecoITSClusters = 0;
  Int_t icountRecoPPR = 0;
  Int_t fiDstar    = 0;
  Int_t fDStarD0   = 0;
  Int_t fDStarMC   = 0;
   
  for (Int_t iPart=0; iPart<mcArray->GetEntriesFast(); iPart++) { 
    AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(mcArray->At(iPart));
    if (!mcPart) {
      AliWarning("Particle not found in tree, skipping"); 
      continue;
    }   
    
    // charm 
    if(TMath::Abs(mcPart->GetPdgCode())==4){
      fcharmpt->Fill(mcPart->Pt());
    }
    
    // fill energy and pt for D* in acceptance with correct prongs 
    Bool_t isOk = DstarInMC(mcPart,mcArray);
    
    if(TMath::Abs(mcPart->GetPdgCode())== 413 && isOk){ // DStar in MC
      fDStarMC++;
      AliDebug(2, "Found a DStar in MC with correct prongs and in acceptance");
      fdstarE ->Fill(mcPart->E());
      fdstarpt->Fill(mcPart->Pt());
    }
    
    if (GetGeneratedValuesFromMCParticle(mcPart, mcArray)){ //D0
   
      // D0s in montecrlo decaing in Kpi
      icountMC++;
      
      // check the MC-Acceptance level cut	
      Int_t daughter0 = mcPart->GetDaughter(0);
      Int_t daughter1 = mcPart->GetDaughter(1);
      
      AliDebug(2, Form("daughter0 = %d and daughter1 = %d",daughter0,daughter1));
     
      if(daughter0!=-1 && daughter1!=-1){
    
	AliAODMCParticle* mcPartDaughter0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(daughter0));
	AliAODMCParticle* mcPartDaughter1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(daughter1));
      
	Double_t eta0 = mcPartDaughter0->Eta();
	Double_t eta1 = mcPartDaughter1->Eta();
	
	Double_t pt0 = mcPartDaughter0->Pt();
	Double_t pt1 = mcPartDaughter1->Pt();
	
	AliDebug(2, Form("Daughter 0: eta = %f, pt = %f", eta0, pt0));
	AliDebug(2, Form("Daughter 1: eta = %f, pt = %f", eta1, pt1));
	
        //daughter D0 in acceptance
	Bool_t daught0inAcceptance = (TMath::Abs(eta0) <= 0.9 && pt0 >= 0.1); 
	Bool_t daught1inAcceptance = (TMath::Abs(eta1) <= 0.9 && pt1 >= 0.1); 
	
	if (daught0inAcceptance && daught1inAcceptance) {
	  AliDebug(2, "D0 Daughter particles in acceptance");
	  
	  icountAcc++;
	  
	  Int_t motherMCD0 = mcPart->GetMother();
	  if(motherMCD0==-1) continue;
	  AliAODMCParticle* mcMothD0 = (AliAODMCParticle*)mcArray->At(motherMCD0); 
	  if(mcMothD0->GetPdgCode()==413 || mcMothD0->GetPdgCode()== -413 ) fDStarD0++;
	}
      }
    }
    else {
      AliDebug(3,"Problems in the task");
      continue;
    }
  }    
  
  AliDebug(2, Form("Found %i MC particles that are D0 in Kpi!!",icountMC));
  AliDebug(2, Form("Found %i MC particles that are D0 in Kpi and satisfy Acc cuts!!",icountAcc));
  
  // Now perform the D* in jet reconstruction
  
  // fill statistic
  fCountMC += icountMC;
  fCountAcc += icountAcc;
  fCountDStar += fDStarD0;
  fCountDStarMC +=fDStarMC;
  
  Int_t efficiencyCeck = 0;
  Int_t efficiency = 0;
  
  // Normalization factor
  if(fRequireNormalization){       
    ftrigger->Fill(1);
  }

  Int_t nJets = 0; // for reco D0 countings

  for (Int_t iJets = 0; iJets<arrayofJets->GetEntriesFast(); iJets++) { // main loop on jets
    
    AliAODJet* jet = (AliAODJet*)arrayofJets->At(iJets);

    //jets variables
    Double_t ejet   = jet->E();
    Double_t phiJet = jet->Phi();
    Double_t etaJet = jet->Eta();

    // fill energy, eta and phi of the jet
    fEjet   ->Fill(ejet);
    fPhijet ->Fill(phiJet);
    fEtaJet ->Fill(etaJet);
    
    nJets++;

    // loop over the tracks to search for candidates soft pions
    for (Int_t i=0; i<aodEvent->GetNTracks(); i++){ 
      
      AliAODTrack* aodTrack = aodEvent->GetTrack(i);
      if(efficiencyCeck == 1)   efficiency = -999; 

      Double_t vD0[4] = {0.,0.,0.,0.};   
      Double_t vD0bar[4] ={0.,0.,0.,0.};

      Int_t pCharge = aodTrack->Charge();

      // few selections on soft pion
      Bool_t tPrimCand =  aodTrack->IsPrimaryCandidate(); // is it primary? 

      if(aodTrack->Pt()>= 5 || aodTrack->Pt()< 0.08 ) continue; //cut on soft pion pt VERY large 
                                                                  //~ D*s of pt >80GeV with a soft pion of 5GeV! 
      if(TMath::Abs(aodTrack->Eta())>0.9) continue;

      // if D*+ analysis tha D0 and pi+         
      if(fComputeD0 && pCharge!= 1 ) continue; 
      // if D*- analysis tha D0bar and pi-
      if(!fComputeD0 && pCharge!= -1 ) continue; 
      
      if(tPrimCand && arrayVerticesHF->GetEntriesFast()>0){ // isPion and is Primary, no PID for now
	
        // label to the candidate soft pion
        Int_t pLabel = aodTrack->GetLabel();
        
	// prepare the TLorentz vector for the pion	
	Float_t pionMass = TDatabasePDG::Instance()->GetParticle(211)->Mass(); 
	fLorentzTrack3.SetPxPyPzE(aodTrack->Px(),aodTrack->Py(),aodTrack->Pz(),aodTrack->E(pionMass)); 
	
        // search for the D0
	for (Int_t iVertex = 0; iVertex<arrayVerticesHF->GetEntriesFast(); iVertex++) { 
	  
	  Double_t invM      = 0;          
	  Double_t invMDStar = 0; 
          Double_t dPhi = 0;

	  efficiencyCeck = 1; 
	  
	  AliAODRecoDecayHF2Prong* vtx = (AliAODRecoDecayHF2Prong*)arrayVerticesHF->At(iVertex);
	  
	  Double_t pt = vtx->Pt();

	  Bool_t acceptanceProng0 = (TMath::Abs(vtx->EtaProng(0))<= 0.9 && vtx->PtProng(0) > 0.1);
	  Bool_t acceptanceProng1 = (TMath::Abs(vtx->EtaProng(1))<= 0.9 && vtx->PtProng(1) > 0.1);
	  
          Int_t pdgDgD0toKpi[2]={321,211};

          Int_t mcLabel = -1; 
          if(nJets == 1) mcLabel = vtx->MatchToMC(421, mcArray,2,pdgDgD0toKpi) ;   //MC D0
         
	  if (acceptanceProng0 && acceptanceProng1) {
              
	    AliDebug(2,"D0 reco daughters in acceptance");
	    if(mcLabel!=-1) icountRecoAcc++;   
	    
            AliAODTrack *track0 = (AliAODTrack*)vtx->GetDaughter(0);
	    AliAODTrack *track1 = (AliAODTrack*)vtx->GetDaughter(1);

            // check for ITS refit (already required at the ESD filter level )
	     Bool_t kRefitITS = kTRUE;
	    
	    if((!(track0->GetStatus()&AliESDtrack::kITSrefit)|| (!(track1->GetStatus()&AliESDtrack::kITSrefit)))) {
	      kRefitITS = kFALSE;
	     }
	    
	    Int_t ncls0=0,ncls1=0;
	    for(Int_t l=0;l<6;l++) {
	      if(TESTBIT(track0->GetITSClusterMap(),l)) ncls0++;
	      if(TESTBIT(track1->GetITSClusterMap(),l)) ncls1++;
	     }
	    AliDebug(2, Form("n clusters = %d", ncls0));
	    if (ncls0 >= fMinITSClusters && ncls1 >= fMinITSClusters) {
	      
	      if(mcLabel!=-1) icountRecoITSClusters++; 
	      
	      // D0 cutting varibles
	      Double_t cosThetaStar = 9999.;
	      Double_t pTpi = 0.;
	      Double_t pTK = 0.;
	      
	      // D0, D0bar
	      if (fComputeD0){
		cosThetaStar = vtx->CosThetaStarD0();
	 	    
		pTpi = vtx->PtProng(0);
		pTK =  vtx->PtProng(1);    
	      }else{
		cosThetaStar = vtx->CosThetaStarD0bar();
	          
		pTpi = vtx->PtProng(1);
		pTK =  vtx->PtProng(0);    
	      }
  
	      AliDebug(2,"D0 reco daughters in acceptance");
	    
	      Double_t dca =   vtx->GetDCA()*1E4;  	    
	      Double_t d01 =   vtx->Getd0Prong(0)*1E4;
	      Double_t d00 =   vtx->Getd0Prong(1)*1E4;
	      Double_t d0d0 =  vtx->Prodd0d0()*1E8; 

	      TVector3 mom(vtx->Px(),vtx->Py(),vtx->Pz());
      	      TVector3 flight((vtx->Xv())- primaryPos[0],(vtx->Yv())- primaryPos[1],(vtx->Zv())- primaryPos[2]); 
	      Double_t pta = mom.Angle(flight); 
	      Double_t cosPointingAngle = TMath::Cos(pta);
	      	      
	      // D* cuts for correlation analysis 	    
	      Double_t cuts[6] = {9999999., 1.1, 0., 9999999., 9999999., 0.};
	      
	      if (pt <= 1){ // first bin not optimized
		cuts[0] = 400;
		cuts[1] = 0.8;
		cuts[2] = 0.21;
		cuts[3] = 500;
		cuts[4] = -1000;
		cuts[5] = 0.8;  
	      }
	      else if (pt > 1 && pt <= 2){
		cuts[0] = 300; 
		cuts[1] = 0.7; 
		cuts[2] = 0.8; 
		cuts[3] = 210; 
		cuts[4] = -2000;
		cuts[5] = 0.8;  
	      }
	      else if (pt > 2 && pt <= 3){
		cuts[0] = 400;
		cuts[1] = 0.8; 
		cuts[2] = 0.8;
		cuts[3] = 420; // looser for correlations
		cuts[4] = -1000;
		cuts[5] = 0.8;   
	      }
	      else if (pt > 3 && pt <= 5){
		cuts[0] = 200;  
		cuts[1] = 0.8; 
		cuts[2] = 1.2;  
		cuts[3] = 560; //looser for correlations 
		cuts[4] = -1000;
		cuts[5] = 0.8;  
	      }
	      else if (pt > 5){
		cuts[0] = 800;
		cuts[1] = 1.0;
		cuts[2] = 1.2; 
		cuts[3] = 700;  
		cuts[4] = -1000;
		cuts[5] = 0.8;  
	      }
	      // apply D0 cuts
	      
	      if (dca < cuts[0] 
		  && TMath::Abs(cosThetaStar) < cuts[1]  
		  && pTpi > cuts[2] 
		  && pTK > cuts[2]  
		  && TMath::Abs(d01) < cuts[3] 
		  && TMath::Abs(d00) < cuts[3]  
		  && d0d0 < cuts[4] 
		  && cosPointingAngle > cuts[5]
		){
		
		Int_t v0quality = -1;
		
		if(fComputeD0){ // D0 from default
		  
		  if(vtx->ChargeProng(0)==-1){
		    fLorentzTrack1.SetPxPyPzE( vtx->PxProng(0),vtx->PyProng(0), vtx->PzProng(0),vtx->EProng(0,321) );
		    fLorentzTrack2.SetPxPyPzE( vtx->PxProng(1),vtx->PyProng(1), vtx->PzProng(1),vtx->EProng(1,211) );
		  }else{		  
		    fLorentzTrack1.SetPxPyPzE( vtx->PxProng(0),vtx->PyProng(0), vtx->PzProng(0),vtx->EProng(0,211) );
		    fLorentzTrack2.SetPxPyPzE( vtx->PxProng(1),vtx->PyProng(1), vtx->PzProng(1),vtx->EProng(1,321) );
		  }
		  
		  v0quality = mcLabel; // is a true D0? 
		  
		  vD0[0] =  (fLorentzTrack1+fLorentzTrack2).Px();
		  vD0[1] =  (fLorentzTrack1+fLorentzTrack2).Py();
		  vD0[2] =  (fLorentzTrack1+fLorentzTrack2).Pz();		
		  vD0[3] =  (fLorentzTrack1+fLorentzTrack2).E();
	      
		  fLorentzTrack4.SetPxPyPzE(vD0[0],vD0[1],vD0[2],vD0[3]);
		  
		}else{ // D0bar analysis
		  
		  if(vtx->ChargeProng(0)==-1){		    
		    fLorentzTrack1.SetPxPyPzE( vtx->PxProng(0),vtx->PyProng(0), vtx->PzProng(0),vtx->EProng(0,211) );
		    fLorentzTrack2.SetPxPyPzE( vtx->PxProng(1),vtx->PyProng(1), vtx->PzProng(1),vtx->EProng(1,321) );
		  }else{
		    fLorentzTrack1.SetPxPyPzE( vtx->PxProng(0),vtx->PyProng(0), vtx->PzProng(0),vtx->EProng(0,321) );
		    fLorentzTrack2.SetPxPyPzE( vtx->PxProng(1),vtx->PyProng(1), vtx->PzProng(1),vtx->EProng(1,211) );
		  }
		  
		  v0quality = mcLabel; // only abs allowed
		  
		  vD0bar[0] = (fLorentzTrack1+fLorentzTrack2).Px();
		  vD0bar[1] = (fLorentzTrack1+fLorentzTrack2).Py();
		  vD0bar[2] = (fLorentzTrack1+fLorentzTrack2).Pz();
		  vD0bar[3] = (fLorentzTrack1+fLorentzTrack2).E();
		
		  fLorentzTrack4.SetPxPyPzE(vD0bar[0],vD0bar[1],vD0bar[2],vD0bar[3]);		  
		}
		
		// D0-D0bar related variables
		
		invM = GetInvariantMass(fLorentzTrack1,fLorentzTrack2);
		fInvMass->Fill(invM); 
		
		if(v0quality>=0){
		  icountRecoPPR++;             
		}
		
		//DStar invariant mass
		invMDStar = GetInvariantMassDStar(fLorentzTrack3,fLorentzTrack4);
		
		//conversion from phi TLorentzVerctor to phi aliroot
		Double_t kconvert = 0;
		Double_t phiDStar = (fLorentzTrack3 + fLorentzTrack4).Phi();	      
		kconvert = phiDStar;
		if(phiDStar<0) kconvert = phiDStar + 2*(TMath::Pi());
		phiDStar = kconvert;
		              
		// dphi between jet and D* 
		dPhi = phiJet - phiDStar;
		
		//Just for plotting pourposal
		if(dPhi<=-1.58) dPhi = TMath::Abs(dPhi);
		if(dPhi>4.72){	
		  dPhi = dPhi-2*(TMath::Pi());
		}
		
                //Alternative cutting procedure 
                //the cut on cosThetaStar is to reduce near collinear
                //background from jet fragmentation

                Bool_t tCut = EvaluateCutOnPiD0pt(vtx,aodTrack);

		if(ftopologicalCut && tCut) continue;
		if(ftopologicalCut && TMath::Abs(cosThetaStar)>cuts[1]) continue;

		if(invM >= 1.829 && invM <= 1.901){ // ~4 sigma cut on D0 mass
		  
		  if(v0quality !=-1 && pLabel!=-1) {

		    AliAODMCParticle* mcPion = (AliAODMCParticle*)mcArray->At(pLabel);
                    Int_t motherMCPion = mcPion->GetMother();
		    
                    if(motherMCPion!=-1){ //mother of soft pion cand
		      AliAODMCParticle* mcMother = (AliAODMCParticle*)mcArray->At(motherMCPion); 
		      if(TMath::Abs(mcMother->GetPdgCode()) == 413){
			
			fMCPionPt->Fill(mcPion->Pt()); 
			
			fDStarMass->Fill(invMDStar); 
			fTrueDiff->Fill(invMDStar-invM);
			fiDstar++;
		      }
		    }
		  }
		  
		  fDStar->Fill(invMDStar);
		  fDiff->Fill(invMDStar-invM); // M(Kpipi)-M(Kpi)
		  
		  // now the dphi signal and the fragmentation function 
		  if((invMDStar-invM)<=0.150 && (invMDStar-invM)>=0.140) { 
		    
		    //fill candidates D* and soft pion reco pt
		    fPtPion->Fill(aodTrack->Pt());		  
		    fRECOPtDStar->Fill((fLorentzTrack3 + fLorentzTrack4).Pt());		    
		    fPhi ->Fill(dPhi);
		    
         	    if(dPhi>=-0.5 && dPhi<=0.5){  // evaluate in the near side		      
		      Double_t dStarMom = (fLorentzTrack3 + fLorentzTrack4).P();
		      Double_t zFrag = (TMath::Cos(dPhi)*dStarMom)/ejet;                               
		      fResZ->Fill(TMath::Abs(zFrag));		      
		    }		    
		  }
		}
		
		// evaluate side band background
		SideBandBackground(invM, invMDStar, ejet, dPhi);
		
		invM      = 0;      
		invMDStar = 0;          
		
	      }    // end PPR cuts
	    }    // end ITS cluster
	  }    // end acceptance
	}    // D0 cand
      }   // softpion
    }  // tracks
  } // jets

  AliDebug(2, Form("Found %i Reco particles that are D0!!",icountReco));
   
  fCountReco+= fiDstar;
  fCountRecoAcc+= icountRecoAcc;
  fCountRecoITSClusters+= icountRecoITSClusters;
  fCountRecoPPR+= icountRecoPPR;

  PostData(1,fOutput);
}

//________________________________________ terminate ___________________________
void AliAnalysisTaskSEDStarJets::Terminate(Option_t*)
{    
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  
  Info("Terminate","");
  AliAnalysisTaskSE::Terminate();
  
  AliInfo(Form("Found %i MC particles that are D0->kpi, in %d events",fCountMC,fEvents));
  AliInfo(Form("Found %i of that MC D0->kpi are in acceptance, in %d events",fCountAcc,fEvents));
  AliInfo(Form("Found %i MC particles that are D*->D0pi, in %d events",fCountDStarMC,fEvents));
  AliInfo(Form("Found %i of that MC D*->D0pi(D0->kpi) are in acceptance, in %d events",fCountDStar,fEvents));

  AliInfo(Form("Found %i RECO particles after cuts that are DStar, in %d events",fCountReco,fEvents));
  AliInfo(Form("Among the above, found %i reco D0 that are decaying in K+pi and are in the requested acceptance, in %d events",fCountRecoAcc,fEvents));
  AliInfo(Form("Among the above, found %i reco D0 that are decaying in K+pi and have at least %d clusters in ITS, in %d events",fCountRecoITSClusters,fMinITSClusters,fEvents));
  AliInfo(Form("Among the above, found %i reco D0 that are decaying in K+pi and satisfy PPR cuts, in %d events",fCountRecoPPR,fEvents));

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }
  
  fMCPionPt     = dynamic_cast<TH1F*>(fOutput->FindObject("fMCPionPt"));
  fcharmpt      = dynamic_cast<TH1F*>(fOutput->FindObject("fcharmpt"));
  fdstarE       = dynamic_cast<TH1F*>(fOutput->FindObject("fdstarE"));
  fdstarpt      = dynamic_cast<TH1F*>(fOutput->FindObject("fdstarpt"));
  fDStarMass    = dynamic_cast<TH1F*>(fOutput->FindObject("fDStarMass"));
  fTrueDiff     = dynamic_cast<TH1F*>(fOutput->FindObject("fTrueDiff"));
  fInvMass      = dynamic_cast<TH1F*>(fOutput->FindObject("fInvMass"));
  fPtPion       = dynamic_cast<TH1F*>(fOutput->FindObject("fPtPion "));
  fDStar        = dynamic_cast<TH1F*>(fOutput->FindObject("fDStar"));
  fDiff         = dynamic_cast<TH1F*>(fOutput->FindObject("fDiff"));
  fDiffSideBand = dynamic_cast<TH1F*>(fOutput->FindObject("fDiffSideBand"));
  ftrigger      = dynamic_cast<TH1F*>(fOutput->FindObject("ftrigger"));
  fRECOPtDStar  = dynamic_cast<TH1F*>(fOutput->FindObject("fRECOPtDStar"));
  fEjet         = dynamic_cast<TH1F*>(fOutput->FindObject("fEjet"));
  fPhijet       = dynamic_cast<TH1F*>(fOutput->FindObject("fPhijet"));
  fEtaJet       = dynamic_cast<TH1F*>(fOutput->FindObject("fEtaJet"));
  fPhi          = dynamic_cast<TH1F*>(fOutput->FindObject("fPhi"));
  fResZ         = dynamic_cast<TH1F*>(fOutput->FindObject("fResZ"));
  fResZBkg      = dynamic_cast<TH1F*>(fOutput->FindObject("fResZBkg"));
  fPhiBkg       = dynamic_cast<TH1F*>(fOutput->FindObject("fPhiBkg"));
  fD0ptvsSoftPtSignal = dynamic_cast<TH2F*>(fOutput->FindObject("fD0ptvsSoftPtSignal"));
  fD0ptvsSoftPt       = dynamic_cast<TH2F*>(fOutput->FindObject("fD0ptvsSoftPt"));
 
}
//___________________________________________________________________________

void AliAnalysisTaskSEDStarJets::UserCreateOutputObjects() { 
 // output
  
  Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());
  
  //slot #1  
  OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner();
  // define histograms
  DefineHistoFroAnalysis();

  return;
}
//_______________________________________D0 in MC ______________________________________

Bool_t AliAnalysisTaskSEDStarJets::GetGeneratedValuesFromMCParticle(AliAODMCParticle* const mcPart, TClonesArray* const mcArray) const {
  // D0 in MC
 
  Bool_t isOk = kFALSE;

  // is a D0?  
  if(TMath::Abs(mcPart->GetPdgCode())!=421) return isOk;

  // check whether the D0 decays in pi+K
  
  // getting the daughters
  Int_t daughter0 = mcPart->GetDaughter(0);
  Int_t daughter1 = mcPart->GetDaughter(1);
  
  AliDebug(2, Form("daughter0 = %d and daughter1 = %d",daughter0,daughter1));
  if (daughter0 == 0 || daughter1 == 0) {
    AliDebug(2, "Error! the D0 MC doesn't have correct daughters!!");
    return isOk;  
  }
  if (TMath::Abs(daughter1 - daughter0)!= 1) {
    AliDebug(2, "The D0 MC doesn't come from a 2-prong decay, skipping!!");
    return isOk;  
  }

  AliAODMCParticle* mcPartDaughter0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(daughter0));
  AliAODMCParticle* mcPartDaughter1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(daughter1)); 
  
  if (!mcPartDaughter0 || !mcPartDaughter1) {
    AliWarning("At least one Daughter Particle not found in tree, skipping"); 
    return isOk;  
  }
  
  // check for the correct daughters
  if((TMath::Abs(mcPartDaughter0->GetPdgCode())!=211) && (TMath::Abs(mcPartDaughter0->GetPdgCode())!=321)) return isOk;
  if((TMath::Abs(mcPartDaughter1->GetPdgCode())!=211) && (TMath::Abs(mcPartDaughter1->GetPdgCode())!=321)) return isOk;
  
  isOk = kTRUE;
  
  return isOk;

}

//_______________________________D0 invariant mass________________________________
 
Double_t AliAnalysisTaskSEDStarJets::GetInvariantMass(TLorentzVector LorentzTrack1, TLorentzVector LorentzTrack2){
  
  return TMath::Abs((LorentzTrack1+LorentzTrack2).M());   // invariant mass of two tracks   
}
//______________________________D* invariant mass_________________________________

Double_t AliAnalysisTaskSEDStarJets::GetInvariantMassDStar(TLorentzVector LorentzTrack4, TLorentzVector LorentzTrack3){
  
  return TMath::Abs((LorentzTrack4+LorentzTrack3).M());   // invariant mass of two tracks   
}                  
//_________________________________D* in MC _______________________________________

Bool_t  AliAnalysisTaskSEDStarJets::DstarInMC(AliAODMCParticle* const mcPart, TClonesArray* mcArray){
  // D* in MC
  Bool_t dStarKpi = kFALSE;

  // is a D*?
  if (mcPart->GetPdgCode() != 413 && mcPart->GetPdgCode() != -413 ) {
    AliDebug(2, "warning! This is not a Dstar!!");
    return dStarKpi;  
  }
  
  // getting the daughters
  Int_t daughter0 = mcPart->GetDaughter(0);
  Int_t daughter1 = mcPart->GetDaughter(1);
  
  AliDebug(2, Form("daughter0 = %d and daughter1 = %d",daughter0,daughter1));
  
  // check if the daughter are correct. Should be everytime the case!
  if (daughter0 == 0 || daughter1 == 0) {
    AliDebug(2, "Error! the D* MC doesn't have correct daughters!!");
    return dStarKpi;  
  }
  if (TMath::Abs(daughter0 - daughter1) != 1) {
    AliDebug(2, "The D* MC doesn't come from a 2-prong decay, skipping!!");
    return dStarKpi;  
  }
 
  AliAODMCParticle* mcPartDaughter0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(daughter0));
  AliAODMCParticle* mcPartDaughter1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(daughter1));

  if (!mcPartDaughter0 || !mcPartDaughter1) {
    AliWarning("At least one Daughter Particle not found in tree, skipping"); 
    return dStarKpi;  
   }

  if((TMath::Abs(mcPartDaughter0->GetPdgCode())!=211) && (TMath::Abs(mcPartDaughter0->GetPdgCode())!=421)) return dStarKpi;
  if((TMath::Abs(mcPartDaughter1->GetPdgCode())!=211) && (TMath::Abs(mcPartDaughter1->GetPdgCode())!=421)) return dStarKpi;
  
  // are the daughters in acceptance? 
  if((TMath::Abs(mcPartDaughter0->Eta())<=0.9) && (TMath::Abs(mcPartDaughter1->Eta())<=0.9)){
    AliDebug(2, "The D* MC is in acceptance");
  } 

  dStarKpi = kTRUE;
  return dStarKpi;

}

//___________________________________ hiostograms _______________________________________

Bool_t  AliAnalysisTaskSEDStarJets::DefineHistoFroAnalysis(){
  
  // Invariant mass related histograms
  fInvMass = new TH1F("invMass","Kpi invariant mass distribution",1500,.5,3.5);
  fInvMass->SetStats(kTRUE);
  fInvMass->GetXaxis()->SetTitle("GeV/c");
  fInvMass->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fInvMass);
  
  fDStar = new TH1F("invMassDStar","DStar invariant mass after D0 cuts ",600,1.8,2.4);
  fDStar->SetStats(kTRUE);
  fDStar->GetXaxis()->SetTitle("GeV/c");
  fDStar->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fDStar);

  fDiff = new TH1F("Diff","M(kpipi)-M(kpi)",750,0.1,0.2);
  fDiff->SetStats(kTRUE);
  fDiff->GetXaxis()->SetTitle("M(kpipi)-M(kpi) GeV/c^2");
  fDiff->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fDiff);
  
  fDiffSideBand = new TH1F("DiffSide","M(kpipi)-M(kpi) Side Band Background",750,0.1,0.2);
  fDiffSideBand->SetStats(kTRUE);
  fDiffSideBand->GetXaxis()->SetTitle("M(kpipi)-M(Kpi) GeV/c^2");
  fDiffSideBand->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fDiffSideBand); 
 
  fDStarMass = new TH1F("RECODStar2","RECO DStar invariant mass distribution",750,1.5,2.5);
  fOutput->Add(fDStarMass);

  fTrueDiff  = new TH1F("dstar","True Reco diff",750,0,0.2);
  fOutput->Add(fTrueDiff);

  // trigger normalization
  ftrigger = new TH1F("Normalization","Normalization factor for correlations",1,0,10);
  ftrigger->SetStats(kTRUE);
  fOutput->Add(ftrigger);

  //correlation fistograms
  fPhi = new TH1F("phi","Delta phi between Jet axis and DStar ",25,-1.57,4.72);
  fPhi->SetStats(kTRUE);
  fPhi->GetXaxis()->SetTitle("#Delta #phi (rad)");
  fPhi->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fPhi);

  fPhiBkg = new TH1F("phiBkg","Delta phi between Jet axis and DStar background ",25,-1.57,4.72);
  fPhiBkg->SetStats(kTRUE);
  fPhiBkg->GetXaxis()->SetTitle("#Delta #phi (rad)");
  fPhiBkg->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fPhiBkg);

  fRECOPtDStar = new TH1F("RECODStar1","RECO DStar pt distribution",600,0,15);
  fRECOPtDStar->SetStats(kTRUE);
  fRECOPtDStar->SetLineColor(2);
  fRECOPtDStar->GetXaxis()->SetTitle("GeV/c");
  fRECOPtDStar->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fRECOPtDStar);

  fPtPion = new TH1F("pionpt","Primary pions candidates pt ",500,0,5);
  fPtPion->SetStats(kTRUE);
  fPtPion->GetXaxis()->SetTitle("GeV/c");
  fPtPion->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fPtPion);
  
  fcharmpt = new TH1F("charmpt","MC Charm pt distribution",    10000,0,5000);
  fdstarE  = new TH1F("dstare", "MC D* energy distribution",   10000,0,5000);
  fdstarpt = new TH1F("dstarpt","MC D* pt distribution",       10000,0,5000);
  fOutput->Add(fcharmpt);
  fOutput->Add(fdstarE);
  fOutput->Add(fdstarpt);
  
  // jet related fistograms
  fEjet      = new TH1F("ejet",  "UA1 algorithm jet energy distribution",1000,0,500);
  fPhijet    = new TH1F("Phijet","UA1 algorithm jet #phi distribution",  200,-7,7);
  fEtaJet    = new TH1F("Etajet","UA1 algorithm jet #eta distribution",  200,-7,7); 
  fOutput->Add(fEjet);
  fOutput->Add(fPhijet);
  fOutput->Add(fEtaJet);

  fMCPionPt = new TH1F("pionptMC2","Primary pions pt from MC ",500,0,5);
  fMCPionPt->SetStats(kTRUE);
  fMCPionPt->GetXaxis()->SetTitle("GeV/c");
  fMCPionPt->GetYaxis()->SetTitle("Entries");
  fOutput->Add(fMCPionPt);
  
  fResZ      = new TH1F("FragFunc","Fragmentation function ",50,0,1);
  fResZBkg   = new TH1F("FragFuncBkg","Fragmentation function background",50,0,1);  
  fOutput->Add(fResZ);
  fOutput->Add(fResZBkg);

  fD0ptvsSoftPt       = new TH2F("D0piRec","Candidates (background + signal)",100,0,3,100,0,15);
  fD0ptvsSoftPtSignal = new TH2F("D0PiSignal","Signal",100,0,3,100,0,15);
  fOutput->Add(fD0ptvsSoftPt);
  fOutput->Add(fD0ptvsSoftPtSignal);

  return kTRUE;
  
}

//______________________________Phase space cut alternative to PPR _______________________________________

Bool_t AliAnalysisTaskSEDStarJets::EvaluateCutOnPiD0pt(AliAODRecoDecayHF2Prong* const vtx,  AliAODTrack* const aodTrack) {

  // The soft pion pt and DO pt are strongly correlated. It can be shown that ~ 95% of the signal in constrained
  // into a narrow band defined by  10 < D0pt/softPt < 18. This cut can be used with a relaxed CosThetaStar cut
  // to reconstruct the D*. 
  
  Double_t softPt   = 0;
  Double_t d0ptReco = 0;
  
  softPt = aodTrack->Pt();
  d0ptReco = vtx->Pt();
  fD0ptvsSoftPt->Fill(softPt,d0ptReco);
  
  if(softPt>0){
    Double_t ratio = d0ptReco/softPt;
    if( ratio <=10. || ratio>=18. ) return kFALSE;
  }
 
  fD0ptvsSoftPtSignal->Fill(softPt,d0ptReco); 
  return kTRUE;
}

//______________________________ side band background for D*___________________________________

void AliAnalysisTaskSEDStarJets::SideBandBackground(Double_t invM, Double_t invMDStar, Double_t ejet, Double_t dPhi){

  //  D* side band background method. Two side bands, in M(Kpi) are taken at ~6 sigmas 
  // (expected detector resolution) on the left and right frm the D0 mass. Each band
  //  has a width of ~5 sigmas. Two band needed  for opening angle considerations   
  
  if((invM>=1.763 && invM<=1.811) || (invM>=1.919 && invM<=1.963)){
    
    fDiffSideBand->Fill(invMDStar-invM); // M(Kpipi)-M(Kpi) side band background
    
    if ((invMDStar-invM)<=0.150 && (invMDStar-invM)>=0.140) {                                                  
      fPhiBkg->Fill(dPhi);
      
      if(dPhi>=-0.5 && dPhi<=0.5){  // evaluate in the near side
	
	Double_t dStarMomBkg = (fLorentzTrack3 + fLorentzTrack4).P();
	Double_t zFragBkg = (TMath::Cos(dPhi)*dStarMomBkg)/ejet;                               
	fResZBkg->Fill(TMath::Abs(zFragBkg));
	
      }
    }
  }
}
