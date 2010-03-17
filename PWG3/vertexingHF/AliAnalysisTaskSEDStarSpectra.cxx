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
//                  Base class for DStar Analysis
//
//  Side Band and like sign background are implemented in the macro
//
//  Utrech Optimized cuts used and PID is on request (flag in che .C) 
//  The D* spectra study is done in pt bins:
//
//  [0,1] [1,2] [2,3] [3,5] [5,8] [8,14]
//
//-----------------------------------------------------------------------
//
//                         Author A.Grelli 
//              ERC-QGP Utrecht University - a.grelli@uu.nl
//           Contributors: C.Ivan  - Utrecht University
//
//-----------------------------------------------------------------------

#include <TSystem.h>
#include <TParticle.h>
#include <TH1I.h>
#include "TROOT.h"

#include "AliPID.h"
#include "AliTPCPIDResponse.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliLog.h"
#include "AliAODVertex.h"
#include "AliAODJet.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAnalysisVertexingHF.h"
#include "AliESDtrack.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSEDStarSpectra.h"

ClassImp(AliAnalysisTaskSEDStarSpectra)

//__________________________________________________________________________
AliAnalysisTaskSEDStarSpectra::AliAnalysisTaskSEDStarSpectra():
  AliAnalysisTaskSE(),
  fEvents(0),
  fVHF(0),  
  fMinITSClusters(0),
  fMinITSClustersSoft(0),
  fUseMCInfo(kTRUE), 
  fOutput(0),
  fNSigma(3),
  fPID(kTRUE),
  fAODTrack(0),
  fMCDStarPt(0),    
  fCEvents(0),     
  fDStarMass(0),  
  fTrueDiff(0),   
  fTrueDiff2(0),  
  fInvMass(0),    
  fInvMass1(0),   
  fInvMass2(0),  
  fInvMass3(0),    
  fInvMass4(0),   
  fInvMass5(0),   
  fPtDStar(0),     
  fDStar(0),   
  fDiff(0),   
  fDiff1(0),    
  fDiff2(0),    
  fDiff3(0),  
  fDiff4(0),  
  fDiff5(0), 
  fDiffSideBand(0),
  fDiffSideBand1(0),
  fDiffSideBand2(0),
  fDiffSideBand3(0),
  fDiffSideBand4(0),
  fDiffSideBand5(0),
  fDiffWrongSign(0),
  fDiffWrongSign1(0),
  fDiffWrongSign2(0),
  fDiffWrongSign3(0),
  fDiffWrongSign4(0),
  fDiffWrongSign5(0)

{
  //
  // Default ctor
  //
}
//___________________________________________________________________________
AliAnalysisTaskSEDStarSpectra::AliAnalysisTaskSEDStarSpectra(const Char_t* name) :
  AliAnalysisTaskSE(name),
  fEvents(0),
  fVHF(0),  
  fMinITSClusters(0),
  fMinITSClustersSoft(0),
  fUseMCInfo(kTRUE),
  fOutput(0),
  fNSigma(3),
  fPID(kTRUE),
  fAODTrack(0),
  fMCDStarPt(0),    
  fCEvents(0),     
  fDStarMass(0),  
  fTrueDiff(0),   
  fTrueDiff2(0),  
  fInvMass(0),    
  fInvMass1(0),   
  fInvMass2(0),  
  fInvMass3(0),    
  fInvMass4(0),   
  fInvMass5(0),   
  fPtDStar(0),    
  fDStar(0),   
  fDiff(0),   
  fDiff1(0),    
  fDiff2(0),    
  fDiff3(0),  
  fDiff4(0),  
  fDiff5(0), 
  fDiffSideBand(0),
  fDiffSideBand1(0),
  fDiffSideBand2(0),
  fDiffSideBand3(0),
  fDiffSideBand4(0),
  fDiffSideBand5(0),
  fDiffWrongSign(0),
  fDiffWrongSign1(0),
  fDiffWrongSign2(0),
  fDiffWrongSign3(0),
  fDiffWrongSign4(0),
  fDiffWrongSign5(0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskSEDStarSpectra","Calling Constructor");
  
  DefineOutput(1,TList::Class());
}

//___________________________________________________________________________
AliAnalysisTaskSEDStarSpectra& AliAnalysisTaskSEDStarSpectra::operator=(const AliAnalysisTaskSEDStarSpectra& c) 
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
AliAnalysisTaskSEDStarSpectra::AliAnalysisTaskSEDStarSpectra(const AliAnalysisTaskSEDStarSpectra& c) :
  AliAnalysisTaskSE(c),
  fEvents(c.fEvents),
  fVHF(c.fVHF),  
  fMinITSClusters(c.fMinITSClusters),
  fMinITSClustersSoft(c.fMinITSClustersSoft),
  fUseMCInfo(c.fUseMCInfo),
  fOutput(c.fOutput),
  fNSigma(c.fNSigma),
  fPID(c.fPID),
  fAODTrack(c.fAODTrack),
  fMCDStarPt(c.fMCDStarPt),    
  fCEvents(c.fCEvents),     
  fDStarMass(c.fDStarMass),  
  fTrueDiff(c.fTrueDiff),   
  fTrueDiff2(c.fTrueDiff2),  
  fInvMass(c.fInvMass),    
  fInvMass1(c.fInvMass1),   
  fInvMass2(c.fInvMass2),  
  fInvMass3(c.fInvMass3),    
  fInvMass4(c.fInvMass4),   
  fInvMass5(c.fInvMass5),   
  fPtDStar(c.fPtDStar),    
  fDStar(c.fDStar),   
  fDiff(c.fDiff),   
  fDiff1(c.fDiff1),    
  fDiff2(c.fDiff2),    
  fDiff3(c.fDiff3),  
  fDiff4(c.fDiff4),  
  fDiff5(c.fDiff5), 
  fDiffSideBand(c.fDiffSideBand),
  fDiffSideBand1(c.fDiffSideBand1),
  fDiffSideBand2(c.fDiffSideBand2),
  fDiffSideBand3(c.fDiffSideBand3),
  fDiffSideBand4(c.fDiffSideBand4),
  fDiffSideBand5(c.fDiffSideBand5),
  fDiffWrongSign(c.fDiffWrongSign),
  fDiffWrongSign1(c.fDiffWrongSign1),
  fDiffWrongSign2(c.fDiffWrongSign2),
  fDiffWrongSign3(c.fDiffWrongSign3),
  fDiffWrongSign4(c.fDiffWrongSign4),
  fDiffWrongSign5(c.fDiffWrongSign5)
{
  //
  // Copy Constructor
  //
}

//___________________________________________________________________________
AliAnalysisTaskSEDStarSpectra::~AliAnalysisTaskSEDStarSpectra() {
  //
  // destructor
  //
  Info("~AliAnalysisTaskSEDStarSpectra","Calling Destructor");
  
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
  if (fVHF) {
    delete fVHF;
    fVHF = 0;
  } 
}
//_________________________________________________
void AliAnalysisTaskSEDStarSpectra::Init(){
  //
  // Initialization
  //

  if(fDebug > 1) printf("AnalysisTaskSEDStarSpectra::Init() \n");

  gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/ConfigVertexingHF.C");
  fVHF = (AliAnalysisVertexingHF*)gROOT->ProcessLine("ConfigVertexingHF()");
  //fVHF->PrintStatus();

  return;
}

//_________________________________________________
void AliAnalysisTaskSEDStarSpectra::UserExec(Option_t *)
{
  // user exec
  if (!fInputEvent) {
    Error("UserExec","NO EVENT FOUND!");
    return;
  }
  
  fCEvents->Fill(1);
  // Load the event
  fEvents++;
  AliInfo(Form("Event %d",fEvents));
  if (fEvents%10000 ==0) AliInfo(Form("Event %d",fEvents));
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
  TClonesArray *arrayDStartoD0pi=0;
  Init();
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
 
 // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aodEvent->GetPrimaryVertex();

  // counters for efficiencies
  Int_t icountReco = 0;
  
  //D* and D0 prongs needed to MatchToMC method
  Int_t pdgDgDStartoD0pi[2]={421,211};
  Int_t pdgDgD0toKpi[2]={321,211};
  
  if (!arrayDStartoD0pi){
    AliInfo("Could not find array of HF vertices, skipping the event");
    return;
  }else AliDebug(2, Form("Found %d vertices",arrayDStartoD0pi->GetEntriesFast())); 
    
  // loop over the tracks to search for candidates soft pion
  
  for (Int_t iDStartoD0pi = 0; iDStartoD0pi<arrayDStartoD0pi->GetEntriesFast(); iDStartoD0pi++) {
    
    // D* candidates
    AliAODRecoCascadeHF* dstarD0pi = (AliAODRecoCascadeHF*)arrayDStartoD0pi->At(iDStartoD0pi);
    
    // D0 from the reco cascade
    AliAODRecoDecayHF2Prong* theD0particle = (AliAODRecoDecayHF2Prong*)dstarD0pi->Get2Prong();
    Bool_t unsetvtx=kFALSE;
    
    Double_t finvM =0;
    Double_t finvMDStar = 0;
    // needed for pointing angle
    if(!theD0particle->GetOwnPrimaryVtx()) {
      theD0particle->SetOwnPrimaryVtx(vtx1);
      unsetvtx=kTRUE;
    }    
    Bool_t isDStar = 0;

    // mc analysis 
    if(fUseMCInfo){
    //MC array need for maching
      TClonesArray* mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!mcArray) AliError("Could not find Monte-Carlo in AOD");
      // find associated MC particle for D* ->D0toKpi
      Int_t mcLabel = dstarD0pi->MatchToMC(413,421,pdgDgDStartoD0pi,pdgDgD0toKpi,mcArray);
      if(mcLabel>=0) isDStar = 1;
    }

    // soft pion
    AliAODTrack *track2 = (AliAODTrack*)dstarD0pi->GetBachelor(); 

    //D0tokpi
    AliAODTrack *track0 = (AliAODTrack*)theD0particle->GetDaughter(0);
    AliAODTrack *track1 = (AliAODTrack*)theD0particle->GetDaughter(1);
    
    Double_t pt = dstarD0pi->Pt();
    
    //kaon PID on low pt D*
    if(pt<3 && fPID){
      if (dstarD0pi->Charge()>0){
	if(!SelectPID(track1,fNSigma)) continue;
      }else{
	if(!SelectPID(track0,fNSigma)) continue;
      }
    }

    // reft in ITS for soft pion
    //if((!(track2->GetStatus()&AliESDtrack::kITSrefit))) continue;
       
    // cut in acceptance for the soft pion and for the D0 daughters      
    Bool_t acceptanceProng0 = (TMath::Abs(theD0particle->EtaProng(0))<= 0.9 && theD0particle->PtProng(0) >= 0.1);
    Bool_t acceptanceProng1 = (TMath::Abs(theD0particle->EtaProng(1))<= 0.9 && theD0particle->PtProng(1) >= 0.1);
    // soft pion acceptance ... is it fine 0.9?????
    Bool_t acceptanceProng2 = (TMath::Abs(track2->Eta())<= 1.0 && track2->Pt() >= 0.05);
    
    if (acceptanceProng0 && acceptanceProng1 && acceptanceProng2) {
      AliDebug(2,"D* reco daughters in acceptance");

      // cut on the min n. of clusters in ITS for the D0 and soft pion
      Int_t ncls0=0,ncls1=0,ncls2=0;
      for(Int_t l=0;l<6;l++) {
	if(TESTBIT(track0->GetITSClusterMap(),l)) ncls0++;
	if(TESTBIT(track1->GetITSClusterMap(),l)) ncls1++;
	if(TESTBIT(track2->GetITSClusterMap(),l)) ncls2++;
      }       
      // see AddTask for soft pion and D0 prongs ITS clusters request
      if (ncls0 >= fMinITSClusters && ncls1 >= fMinITSClusters && ncls2>=fMinITSClustersSoft) { 
 
        // tag the D* decay do not double count background
        Int_t decayTag = track2->Charge(); 
        // D0 pt needed for the cuts       
	Double_t ptD0 = theD0particle->Pt();
       
	Int_t okD0 = 0;
	Int_t okD0bar = 0;
        Int_t okD0WrongSign =0;
	Int_t okD0barWrongSign =0;

	// optimized D* cuts from UU
	Bool_t tCutDStar = SetUtrechtSelections(ptD0,fVHF);  
	Bool_t tCutOk = kFALSE;
	if(tCutDStar) tCutOk = theD0particle->SelectD0(fVHF->GetD0toKpiCuts(),okD0,okD0bar);

        // flags for wrong sign
	okD0WrongSign=okD0;
	okD0barWrongSign=okD0bar;

	// search for D*
	if(tCutOk){	
	  // correct decay 
          if(decayTag>0 ? (okD0bar = 0) : (okD0 = 0));

	  finvM = dstarD0pi->InvMassD0();

	  if(okD0 == 1) fInvMass->Fill(finvM); 
	  if(okD0bar == 1) fInvMass->Fill(finvM); 

	  //DStar invariant mass
	  finvMDStar = dstarD0pi->InvMassDstarKpipi();
	  
	  if(finvM >= 1.829 && finvM <= 1.901){ // ~3 sigma cut on D0 mass
	   
	    if(isDStar == 1) { 
	      fDStarMass->Fill(finvMDStar); 
	      fTrueDiff ->Fill(dstarD0pi->DeltaInvMass());
	      fTrueDiff2->Fill(pt,dstarD0pi->DeltaInvMass());
	      fMCDStarPt->Fill(pt);  
	    }

            // D* candidates - pt integrated
	    if(okD0==1){
	      fDStar->Fill(finvMDStar);
	      fDiff->Fill(dstarD0pi->DeltaInvMass()); // M(Kpipi)-M(Kpi) 
	    }else if(okD0bar==1){
	      fDStar->Fill(finvMDStar);
	      fDiff->Fill(dstarD0pi->DeltaInvMass()); // M(Kpipi)-M(Kpi) 
	    }

	    //D* candidates - pt bins
	    if(pt>1 && pt<=2){ // [1-2]
	      if(okD0==1){
		fDiff1->Fill(dstarD0pi->DeltaInvMass());
		fInvMass1->Fill(finvM);
	      }else if(okD0bar==1){
		fDiff1->Fill(dstarD0pi->DeltaInvMass());
		fInvMass1->Fill(finvM);
	      }
	    }
	    if( pt>2 && pt<=3){  // [2-3]
	      if(okD0==1){
		fDiff2->Fill(dstarD0pi->DeltaInvMass());
		fInvMass2->Fill(finvM);
	      }else if(okD0bar==1){
		fDiff2->Fill(dstarD0pi->DeltaInvMass());
		fInvMass2->Fill(finvM);
	      }
	    }
	    if(pt>3  && pt<=5){ // [3-5]
	      if(okD0==1){
		fDiff3->Fill(dstarD0pi->DeltaInvMass());
		fInvMass3->Fill(finvM);
	      }else if(okD0bar==1){
		fDiff3->Fill(dstarD0pi->DeltaInvMass());
		fInvMass3->Fill(finvM);
	      }
	    }
	    if( pt>5 && pt<8){ // [5-8]
	      if(okD0==1){
		fDiff4->Fill(dstarD0pi->DeltaInvMass());
		fInvMass4->Fill(finvM);
	      }else if(okD0bar==1){
		fDiff4->Fill(dstarD0pi->DeltaInvMass());
		fInvMass4->Fill(finvM);	
	      }
	    }	    
	    if(pt>=8){ // [>8]
	      if(okD0==1){
		fDiff5->Fill(dstarD0pi->DeltaInvMass());
		fInvMass5->Fill(finvM);
	      }else if(okD0bar==1){
		fDiff5->Fill(dstarD0pi->DeltaInvMass());
		fInvMass5->Fill(finvM);
	      }
	    }
	    
	    // D* pt distribution
	    if((dstarD0pi->DeltaInvMass())>=0.143920 && (dstarD0pi->DeltaInvMass())<=0.14702){
	      if(okD0==1) fPtDStar->Fill(pt);
	      if(okD0bar==1) fPtDStar->Fill(pt);
	    }
 
	  }	
	  // evaluate side band background
	  SideBandBackground(finvM,finvMDStar,pt,okD0,okD0bar);
	  
	  finvM      = 0;      
	  finvMDStar = 0;          
	  
	} // end cuts
	
        // wrong sign background
	if(tCutOk){	
	  // correct decay 
          if(decayTag>0 ? ( okD0WrongSign = 0) : ( okD0barWrongSign = 0));

          //wrong D0 inv mass
          if(okD0WrongSign==1){
	    finvM = theD0particle->InvMassD0();
	  }else if(okD0WrongSign==0){
	    finvM = theD0particle->InvMassD0bar();
	  }
	  // wrong D* inv mass	  
	  Double_t px[3],py[3],pz[3];
	  UInt_t pdg[3]={321,211,211};
	  pdg[0] = (decayTag>0 ? 321 : 211); // positive daughter of D0
	  px[0] = theD0particle->PxProng(0);
	  py[0] = theD0particle->PyProng(0);
	  pz[0] = theD0particle->PzProng(0);
	  pdg[1] = (decayTag>0 ? 211 : 321); // negative daughter of D0
	  px[1] = theD0particle->PxProng(1);
	  py[1] = theD0particle->PyProng(1);
	  pz[1] = theD0particle->PzProng(1);
	  pdg[2] = 211; // soft pion
	  px[2] = track2->Px();
	  py[2] = track2->Py();
	  pz[2] = track2->Pz();

	  Short_t dummycharge=0;
	  Double_t dummyd0[3]={0,0,0};
	  AliAODRecoDecay *rd = new AliAODRecoDecay(0x0,3,dummycharge,px,py,pz,dummyd0);
	  
	  finvMDStar = rd->InvMass(3,pdg);
	  
	  delete rd; rd=NULL;
	  
	  if(finvM >= 1.829 && finvM <= 1.901){ // ~3 sigma cut on D0 mass
	    WrongSignForDStar(finvM,finvMDStar,pt,okD0WrongSign,okD0barWrongSign);
	  }
	}// end of wrong sign
      }
    }
  }
  
  AliDebug(2, Form("Found %i Reco particles that are D*!!",icountReco));
  
  PostData(1,fOutput);
}
//________________________________________ terminate ___________________________
void AliAnalysisTaskSEDStarSpectra::Terminate(Option_t*)
{    
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  
  Info("Terminate","");
  AliAnalysisTaskSE::Terminate();
  
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }
  
  fMCDStarPt      = dynamic_cast<TH1F*>(fOutput->FindObject("fMCDStarPt"));
  fCEvents        = dynamic_cast<TH1F*>(fOutput->FindObject("fCEvents"));
  fDStarMass      = dynamic_cast<TH1F*>(fOutput->FindObject("fDStarMass"));
  fTrueDiff       = dynamic_cast<TH1F*>(fOutput->FindObject("fTrueDiff"));
  fTrueDiff2      = dynamic_cast<TH2F*>(fOutput->FindObject("fTrueDiff2"));
  fInvMass        = dynamic_cast<TH1F*>(fOutput->FindObject("fInvMass"));
  fInvMass1       = dynamic_cast<TH1F*>(fOutput->FindObject("fInvMass1"));
  fInvMass2       = dynamic_cast<TH1F*>(fOutput->FindObject("fInvMass2"));
  fInvMass3       = dynamic_cast<TH1F*>(fOutput->FindObject("fInvMass3"));
  fInvMass4       = dynamic_cast<TH1F*>(fOutput->FindObject("fInvMass4"));
  fInvMass5       = dynamic_cast<TH1F*>(fOutput->FindObject("fInvMass5"));
  fPtDStar        = dynamic_cast<TH1F*>(fOutput->FindObject("fPtDStar "));
  fDStar          = dynamic_cast<TH1F*>(fOutput->FindObject("fDStar"));
  fDiff           = dynamic_cast<TH1F*>(fOutput->FindObject("fDiff"));
  fDiff1          = dynamic_cast<TH1F*>(fOutput->FindObject("fDiff1"));
  fDiff2          = dynamic_cast<TH1F*>(fOutput->FindObject("fDiff2"));
  fDiff3          = dynamic_cast<TH1F*>(fOutput->FindObject("fDiff3"));
  fDiff4          = dynamic_cast<TH1F*>(fOutput->FindObject("fDiff4"));
  fDiff5          = dynamic_cast<TH1F*>(fOutput->FindObject("fDiff5"));
  fDiffSideBand   = dynamic_cast<TH1F*>(fOutput->FindObject("fDiffSideBand"));
  fDiffSideBand1  = dynamic_cast<TH1F*>(fOutput->FindObject("fDiffSideBand1"));
  fDiffSideBand2  = dynamic_cast<TH1F*>(fOutput->FindObject("fDiffSideBand2"));
  fDiffSideBand3  = dynamic_cast<TH1F*>(fOutput->FindObject("fDiffSideBand3"));
  fDiffSideBand4  = dynamic_cast<TH1F*>(fOutput->FindObject("fDiffSideBand4"));
  fDiffSideBand5  = dynamic_cast<TH1F*>(fOutput->FindObject("fDiffSideBand5"));
  fDiffWrongSign  = dynamic_cast<TH1F*>(fOutput->FindObject("fDiffWrongSign"));
  fDiffWrongSign1 = dynamic_cast<TH1F*>(fOutput->FindObject("fDiffWrongSign1"));
  fDiffWrongSign2 = dynamic_cast<TH1F*>(fOutput->FindObject("fDiffWrongSign2"));
  fDiffWrongSign3 = dynamic_cast<TH1F*>(fOutput->FindObject("fDiffWrongSign3"));
  fDiffWrongSign4 = dynamic_cast<TH1F*>(fOutput->FindObject("fDiffWrongSign4"));
  fDiffWrongSign5 = dynamic_cast<TH1F*>(fOutput->FindObject("fDiffWrongSign5"));

}
//___________________________________________________________________________
void AliAnalysisTaskSEDStarSpectra::UserCreateOutputObjects() { 
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

//___________________________________ hiostograms _______________________________________
Bool_t  AliAnalysisTaskSEDStarSpectra::DefineHistoFroAnalysis(){

  // Invariant mass related histograms
  fInvMass = new TH1F("invMass","K#pi invariant mass distribution",1500,.5,3.5);
  fInvMass->SetStats(kTRUE);
  fInvMass->GetXaxis()->SetTitle("M(K#pi) GeV/c");
  fInvMass->GetYaxis()->SetTitle("Entries");
  
  fOutput->Add(fInvMass);
  
  fInvMass1 = new TH1F("invMass1","K#pi invariant mass distribution",3500,1.5,2.2);
  fInvMass1->SetStats(kTRUE);
  fInvMass1->GetXaxis()->SetTitle("M(K#pi) GeV/c");
  fInvMass1->GetYaxis()->SetTitle("Entries");
  
  fOutput->Add(fInvMass1);
  
  fInvMass2 = new TH1F("invMass2","K#pi invariant mass distribution",350,1.5,2.2);
  fInvMass2->SetStats(kTRUE);
  fInvMass2->GetXaxis()->SetTitle("M(K#pi) GeV/c");
  fInvMass2->GetYaxis()->SetTitle("Entries");
  
  fOutput->Add(fInvMass2);
  
  fInvMass3 = new TH1F("invMass3","K#pi invariant mass distribution",350,1.5,2.2);
  fInvMass3->SetStats(kTRUE);
  fInvMass3->GetXaxis()->SetTitle("M(K#pi) GeV/c");
  fInvMass3->GetYaxis()->SetTitle("Entries");
  
  fOutput->Add(fInvMass3);
  
  fInvMass4 = new TH1F("invMass4","K#pi invariant mass distribution",350,1.5,2.2);
  fInvMass4->SetStats(kTRUE);
  fInvMass4->GetXaxis()->SetTitle("M(K#pi) GeV/c");
  fInvMass4->GetYaxis()->SetTitle("Entries");

  fOutput->Add(fInvMass4);
  
  fInvMass5 = new TH1F("invMass5","K#pi invariant mass distribution",350,1.5,2.2);
  fInvMass5->SetStats(kTRUE);
  fInvMass5->GetXaxis()->SetTitle("M(K#pi) GeV/c");
  fInvMass5->GetYaxis()->SetTitle("Entries");

  fOutput->Add(fInvMass5);

  fDStar = new TH1F("invMassDStar","DStar invariant mass after D0 cuts ",600,1.8,2.4);
  fDStar->SetStats(kTRUE);
  fDStar->GetXaxis()->SetTitle("GeV/c");
  fDStar->GetYaxis()->SetTitle("Entries");

  fOutput->Add(fDStar);

  fCEvents = new TH1F("fCEvents","conter",10,0,10);
  fCEvents->SetStats(kTRUE);
  fCEvents->GetXaxis()->SetTitle("1");
  fCEvents->GetYaxis()->SetTitle("counts");

  fOutput->Add(fCEvents);

  fDiff = new TH1F("Diff","M(K#pi#pi)-M(K#pi)",750,0.1,0.2);
  fDiff->SetStats(kTRUE);
  fDiff->GetXaxis()->SetTitle("M(K#pi#pi)-M(K#pi) GeV/c^2");
  fDiff->GetYaxis()->SetTitle("Entries");

  fOutput->Add(fDiff);

  fDiff1 = new TH1F("Diff1","M(K#pi#pi)-M(K#pi)",750,0.1,0.2);
  fDiff1->SetStats(kTRUE);
  fDiff1->GetXaxis()->SetTitle("M(K#pi#pi)-M(K#pi) GeV/c^2");
  fDiff1->GetYaxis()->SetTitle("Entries");

  fOutput->Add(fDiff1);

  fDiff2 = new TH1F("Diff2","M(K#pi#pi)-M(K#pi)",750,0.1,0.2);
  fDiff2->SetStats(kTRUE);
  fDiff2->GetXaxis()->SetTitle("M(K#pi#pi)-M(K#pi) GeV/c^2");
  fDiff2->GetYaxis()->SetTitle("Entries");

  fOutput->Add(fDiff2);

  fDiff3 = new TH1F("Diff3","M(K#pi#pi)-M(K#pi)",750,0.1,0.2);
  fDiff3->SetStats(kTRUE);
  fDiff3->GetXaxis()->SetTitle("M(K#pi#pi)-M(K#pi) GeV/c^2");
  fDiff3->GetYaxis()->SetTitle("Entries");

  fOutput->Add(fDiff3);

  fDiff4 = new TH1F("Diff4","M(K#pi#pi)-M(K#pi)",750,0.1,0.2);
  fDiff4->SetStats(kTRUE);
  fDiff4->GetXaxis()->SetTitle("M(K#pi#pi)-M(K#pi) GeV/c^2");
  fDiff4->GetYaxis()->SetTitle("Entries");

  fOutput->Add(fDiff4);
  
  fDiff5 = new TH1F("Diff5","M(K#pi#pi)-M(K#pi)",750,0.1,0.2);
  fDiff5->SetStats(kTRUE);
  fDiff5->GetXaxis()->SetTitle("M(K#pi#pi)-M(K#pi) GeV/c^2");
  fDiff5->GetYaxis()->SetTitle("Entries");

  fOutput->Add(fDiff5);

  fDiffSideBand = new TH1F("DiffSide","M(kpipi)-M(kpi) Side Band Background",750,0.1,0.2);
  fDiffSideBand->SetStats(kTRUE);
  fDiffSideBand->GetXaxis()->SetTitle("M(kpipi)-M(Kpi) GeV/c^2");
  fDiffSideBand->GetYaxis()->SetTitle("Entries");

  fOutput->Add(fDiffSideBand); 

  fDiffSideBand1 = new TH1F("DiffSide1","M(kpipi)-M(kpi) Side Band Background",750,0.1,0.2);
  fDiffSideBand1->SetStats(kTRUE);
  fDiffSideBand1->GetXaxis()->SetTitle("M(kpipi)-M(Kpi) GeV/c^2");
  fDiffSideBand1->GetYaxis()->SetTitle("Entries");

  fOutput->Add(fDiffSideBand1); 

  fDiffSideBand2 = new TH1F("DiffSide2","M(kpipi)-M(kpi) Side Band Background",750,0.1,0.2);
  fDiffSideBand2->SetStats(kTRUE);
  fDiffSideBand2->GetXaxis()->SetTitle("M(kpipi)-M(Kpi) GeV/c^2");
  fDiffSideBand2->GetYaxis()->SetTitle("Entries");

  fOutput->Add(fDiffSideBand2); 

  fDiffSideBand3 = new TH1F("DiffSide3","M(kpipi)-M(kpi) Side Band Background",750,0.1,0.2);
  fDiffSideBand3->SetStats(kTRUE);
  fDiffSideBand3->GetXaxis()->SetTitle("M(kpipi)-M(Kpi) GeV/c^2");
  fDiffSideBand3->GetYaxis()->SetTitle("Entries");

  fOutput->Add(fDiffSideBand3); 

  fDiffSideBand4 = new TH1F("DiffSide4","M(kpipi)-M(kpi) Side Band Background",750,0.1,0.2);
  fDiffSideBand4->SetStats(kTRUE);
  fDiffSideBand4->GetXaxis()->SetTitle("M(kpipi)-M(Kpi) GeV/c^2");
  fDiffSideBand4->GetYaxis()->SetTitle("Entries");

  fOutput->Add(fDiffSideBand4); 

  fDiffSideBand5 = new TH1F("DiffSide5","M(kpipi)-M(kpi) Side Band Background",750,0.1,0.2);
  fDiffSideBand5->SetStats(kTRUE);
  fDiffSideBand5->GetXaxis()->SetTitle("M(kpipi)-M(Kpi) GeV/c^2");
  fDiffSideBand5->GetYaxis()->SetTitle("Entries");

  fOutput->Add(fDiffSideBand5); 

  fDiffWrongSign = new TH1F("DiffWrong","M(kpipi)-M(kpi) Side Band Background",750,0.1,0.2);
  fDiffWrongSign->SetStats(kTRUE);
  fDiffWrongSign->GetXaxis()->SetTitle("M(kpipi)-M(Kpi) GeV/c^2");
  fDiffWrongSign->GetYaxis()->SetTitle("Entries");

  fOutput->Add(fDiffWrongSign); 

  fDiffWrongSign1 = new TH1F("DiffWrong1","M(kpipi)-M(kpi) Side Band Background",750,0.1,0.2);
  fDiffWrongSign1->SetStats(kTRUE);
  fDiffWrongSign1->GetXaxis()->SetTitle("M(kpipi)-M(Kpi) GeV/c^2");
  fDiffWrongSign1->GetYaxis()->SetTitle("Entries");

  fOutput->Add(fDiffWrongSign1); 

  fDiffWrongSign2 = new TH1F("DiffWrong2","M(kpipi)-M(kpi) Side Band Background",750,0.1,0.2);
  fDiffWrongSign2->SetStats(kTRUE);
  fDiffWrongSign2->GetXaxis()->SetTitle("M(kpipi)-M(Kpi) GeV/c^2");
  fDiffWrongSign2->GetYaxis()->SetTitle("Entries");

  fOutput->Add(fDiffWrongSign2); 

  fDiffWrongSign3 = new TH1F("DiffWrong3","M(kpipi)-M(kpi) Side Band Background",750,0.1,0.2);
  fDiffWrongSign3->SetStats(kTRUE);
  fDiffWrongSign3->GetXaxis()->SetTitle("M(kpipi)-M(Kpi) GeV/c^2");
  fDiffWrongSign3->GetYaxis()->SetTitle("Entries");

  fOutput->Add(fDiffWrongSign3); 

  fDiffWrongSign4 = new TH1F("DiffWrong4","M(kpipi)-M(kpi) Side Band Background",750,0.1,0.2);
  fDiffWrongSign4->SetStats(kTRUE);
  fDiffWrongSign4->GetXaxis()->SetTitle("M(kpipi)-M(Kpi) GeV/c^2");
  fDiffWrongSign4->GetYaxis()->SetTitle("Entries");

  fOutput->Add(fDiffWrongSign4); 

  fDiffWrongSign5 = new TH1F("DiffWrong5","M(kpipi)-M(kpi) Side Band Background",750,0.1,0.2);
  fDiffWrongSign5->SetStats(kTRUE);
  fDiffWrongSign5->GetXaxis()->SetTitle("M(kpipi)-M(Kpi) GeV/c^2");
  fDiffWrongSign5->GetYaxis()->SetTitle("Entries");

  fOutput->Add(fDiffWrongSign5); 
 
  fDStarMass = new TH1F("RECODStar2","RECO DStar invariant mass distribution",750,1.5,2.5);

  fOutput->Add(fDStarMass);
 
  fTrueDiff  = new TH1F("dstar","True Reco diff",750,0,0.2);

  fOutput->Add(fTrueDiff);
  fTrueDiff2 = new TH2F("DiffDstar_pt","True Reco diff vs pt",100,0,15,900,0,0.3);

  fOutput->Add(fTrueDiff2);

  fPtDStar = new TH1F("DStarpt","Reconstructed D* candidates momentum ",500,0,25);
  fPtDStar->SetStats(kTRUE);
  fPtDStar->GetXaxis()->SetTitle("GeV/c");
  fPtDStar->GetYaxis()->SetTitle("Entries");

  fOutput->Add(fPtDStar);

  fMCDStarPt = new TH1F("DStarMC2","Reconstructed D* momentum MC tagged ",500,0,25);
  fMCDStarPt->SetStats(kTRUE);
  fMCDStarPt->GetXaxis()->SetTitle("GeV/c");
  fMCDStarPt->GetYaxis()->SetTitle("Entries");

  fOutput->Add(fMCDStarPt);

  return kTRUE;
  
}
//______________________________ side band background for D*___________________________________
void AliAnalysisTaskSEDStarSpectra::SideBandBackground(Double_t finvM, Double_t finvMDStar,  Double_t pt, Int_t okD0, Int_t okD0bar ){

  //  D* side band background method. Two side bands, in M(Kpi) are taken at ~6 sigmas 
  // (expected detector resolution) on the left and right frm the D0 mass. Each band
  //  has a width of ~5 sigmas. Two band needed  for opening angle considerations   
  
  if((finvM>=1.740 && finvM<=1.829) || (finvM>=1.901 && finvM<=1.990)){
    
    if(okD0==1)       fDiffSideBand->Fill(finvMDStar-finvM); // M(Kpipi)-M(Kpi) side band background
    if(okD0bar==1)    fDiffSideBand->Fill(finvMDStar-finvM); // M(Kpipi)-M(Kpi) side band background

    // pt bins
    if( pt <= 1){
      if(okD0==1)     fDiffSideBand1->Fill(finvMDStar-finvM);
      if(okD0bar==1)  fDiffSideBand1->Fill(finvMDStar-finvM);
    }
    if( pt > 1 && pt <=3){
      if(okD0==1 )    fDiffSideBand2->Fill(finvMDStar-finvM);
      if(okD0bar==1 ) fDiffSideBand2->Fill(finvMDStar-finvM);
    }
    if( pt >3  && pt <=5){
      if(okD0==1)     fDiffSideBand3->Fill(finvMDStar-finvM);
      if(okD0bar==1)  fDiffSideBand3->Fill(finvMDStar-finvM);
    }
    if( pt > 5  && pt <8){
      if(okD0==1)     fDiffSideBand4->Fill(finvMDStar-finvM);
      if(okD0bar==1)  fDiffSideBand4->Fill(finvMDStar-finvM);
    }  
    if( pt >= 8){
      if(okD0==1)     fDiffSideBand5->Fill(finvMDStar-finvM);
      if(okD0bar==1)  fDiffSideBand5->Fill(finvMDStar-finvM);
    }
    
  }
}

//________________________________________________________________________________________________________________
void AliAnalysisTaskSEDStarSpectra::WrongSignForDStar(Double_t finvM, Double_t finvMDStar,  Double_t pt, Int_t okD0, Int_t okD0bar){
  //
  // assign the wrong charge to the soft pion to create background
  //
  if(okD0==1 )       fDiffWrongSign->Fill(finvMDStar-finvM); // M(Kpipi)-M(Kpi) side band background
  if(okD0bar==1 )     fDiffWrongSign->Fill(finvMDStar-finvM); // M(Kpipi)-M(Kpi) side band background
  
  // pt bins
  if( pt <= 1){
    if(okD0==1)      fDiffWrongSign1->Fill(finvMDStar-finvM);
    if(okD0bar==1)   fDiffWrongSign1->Fill(finvMDStar-finvM);
  }
  if( pt > 1 && pt <=3){
    if(okD0==1)      fDiffWrongSign2->Fill(finvMDStar-finvM);
    if(okD0bar==1)   fDiffWrongSign2->Fill(finvMDStar-finvM);
  }
  if( pt >3  && pt <=5){
    if(okD0==1)      fDiffWrongSign3->Fill(finvMDStar-finvM);
    if(okD0bar==1)   fDiffWrongSign3->Fill(finvMDStar-finvM);
  }
  if( pt > 5  && pt <8){
    if(okD0==1)      fDiffWrongSign4->Fill(finvMDStar-finvM);
    if(okD0bar==1)   fDiffWrongSign4->Fill(finvMDStar-finvM);
  }  
  if( pt >= 8){
    if(okD0==1)      fDiffWrongSign5->Fill(finvMDStar-finvM);
    if(okD0bar==1)   fDiffWrongSign5->Fill(finvMDStar-finvM);
  }
  
}
//_____________________________________________________________________________________________
Bool_t AliAnalysisTaskSEDStarSpectra::SetUtrechtSelections(Double_t ptD0, AliAnalysisVertexingHF *fVHF){

  //cuts[0] = inv. mass half width [GeV]
  //cuts[1] = dca [cm]
  //cuts[2] = cosThetaStar
  //cuts[3] = pTK [GeV/c]
  //cuts[4] = pTPi [GeV/c]
  //cuts[5] = d0K [cm]   upper limit!
  //cuts[6] = d0Pi [cm]  upper limit!
  //cuts[7] = d0d0 [cm^2]
  //cuts[8] = cosThetaPoint

  if(ptD0<1.){
    fVHF->SetD0toKpiCuts(0.450,0.04,0.8,0.21,0.21,0.021,0.021,-0.0002,0.9);
  }
  else if(ptD0 >=1. && ptD0 <=2.){
    fVHF->SetD0toKpiCuts(0.450,0.02,0.7,0.8,0.8,0.021,0.021,-0.0002,0.9);
  }  
  else if(ptD0 >2. && ptD0 <=3.){
    fVHF->SetD0toKpiCuts(0.450,0.04,0.8,0.8,0.8,0.035,0.042,-0.000085,0.9);
  } 
  else if(ptD0 >3. && ptD0 <=5.){
    fVHF->SetD0toKpiCuts(0.450,0.016,0.8,1.2,1.2,0.042,0.056,-0.000085,0.9);
  } 
  else if(ptD0 >5.){
    fVHF->SetD0toKpiCuts(0.450,0.08,1.0,1.2,1.2,0.07,0.07,0.0001,0.9);
  }  
  return kTRUE;
}
//_____________________________ pid _______________________________________-
Bool_t AliAnalysisTaskSEDStarSpectra::SelectPID(AliAODTrack *track, Double_t nsig){//pid - K = 3
  //TPC
  Bool_t isKaon=kTRUE;
  if ((track->GetStatus()&AliESDtrack::kTPCpid )==0) return isKaon;
  AliAODPid *pid = track->GetDetPid();
  static AliTPCPIDResponse theTPCpid;
  Double_t nsigma = theTPCpid.GetNumberOfSigmas(track->P(),pid->GetTPCsignal(),track->GetTPCClusterMap().CountBits(), (AliPID::kKaon));
  if (TMath::Abs(nsigma)>nsig) isKaon=kFALSE;
  //ITS
  //
  //
  //

  return isKaon;
}

