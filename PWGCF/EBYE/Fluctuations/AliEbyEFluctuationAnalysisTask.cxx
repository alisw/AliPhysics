#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TParticle.h"
#include "TObjArray.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"

#include "AliESDtrack.h"
#include "AliESDpid.h"

#include "AliEbyEFluctuationAnalysisTask.h"

// Event by event charge fluctuation analysis
// Authors: Satyajit Jena and Panos Cristakoglou
//          (PID by R+)
// 

ClassImp(AliEbyEFluctuationAnalysisTask)

//________________________________________________________________________
AliEbyEFluctuationAnalysisTask::AliEbyEFluctuationAnalysisTask(const char *name) 
  : AliAnalysisTaskSE(name), fESD(0), fOutputList(0), 
    fHistEventStats(0), fHistCentrality(0),
    fHistNMultMC(0), fHistNPlusNMinusMC(0), 
    fESDtrackCuts(0),
    fAnalysisType("ESD"), fAnalysisMode("TPC"),
    fCentralityEstimator("V0M"), fCentralityBins20(kFALSE),
  fVxMax(3.), fVyMax(3.), fVzMax(10.),
  fESDpid(NULL) {
  // Constructor
  
  fOutputList = new TList();
  fOutputList->SetOwner();

  for(Int_t iBin = 1; iBin <= nCentralityBins; iBin++) {
    fHistNMult[iBin-1] = NULL;
    fHistNPlusNMinus[iBin-1] = NULL;
  }

  // Define input and output slots here
  // Input slot #0 works with a TChain
  //  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(2, TList::Class());
}

//________________________________________________________________________
void AliEbyEFluctuationAnalysisTask::UserCreateOutputObjects() {
  // Create histograms
  // Called once

  printf("in UserCreateOutputObjects\n");

  //Event stats.
  TString gCutName[4] = {"Total","Offline trigger",
                         "Vertex","Analyzed"};
  fHistEventStats = new TH1F("fHistEventStats",
                             "Event statistics;;N_{events}",
                             4,0.5,4.5);
  for(Int_t i = 1; i <= 4; i++)
    fHistEventStats->GetXaxis()->SetBinLabel(i,gCutName[i-1].Data());
  fOutputList->Add(fHistEventStats);

  //ESD analysis
  if(fAnalysisType == "ESD") {
    fHistCentrality = new TH1F("fHistCentrality",";Centrality bin;Events",
			       nCentralityBins,0.5,nCentralityBins+0.5);
    fOutputList->Add(fHistCentrality);
    
    TString histName;
    for(Int_t iBin = 1; iBin <= nCentralityBins; iBin++) {
      histName = "fHistNMult"; histName += "Centrality"; histName += iBin; 
      fHistNMult[iBin-1] = new TH1F(histName.Data(), 
				    ";N_{mult.}",
				    500,0,3000);
      fOutputList->Add(fHistNMult[iBin-1]);
    }
    for(Int_t iBin = 1; iBin <= nCentralityBins; iBin++) {
      histName = "fHistNPlusNMinus"; histName += "Centrality"; histName += iBin; 
      fHistNPlusNMinus[iBin-1] = new TH2F(histName.Data(), 
					  ";N_{+};N_{-}",
					  2000,0.5,2000.5,2000,0.5,2000.5);  
      fOutputList->Add(fHistNPlusNMinus[iBin-1]);
    }
  }//ESD analysis
  else if(fAnalysisType == "MC") {
    TString histName = "fHistNMultMC";
    fHistNMultMC = new TH1F(histName.Data(), 
			    ";N_{mult.}",
			    600,0,6000);
    fOutputList->Add(fHistNMultMC);
    
    histName = "fHistNPlusNMinusMC";
    fHistNPlusNMinusMC = new TH2F(histName.Data(), 
				  ";N_{+};N_{-}",
				  3000,0.5,3000.5,3000,0.5,3000.5);  
    fOutputList->Add(fHistNPlusNMinusMC);
  }//MC analysis

  PostData(2, fOutputList);
}

//________________________________________________________________________
void AliEbyEFluctuationAnalysisTask::UserExec(Option_t *) {
  // Main loop
  // Called for each event
  
  Int_t nPlus = 0, nMinus = 0;

  // Post output data.
  //ESD analysis
  if(fAnalysisType == "ESD") {
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!fESD) {
      printf("ERROR: fESD not available\n");
      return;
    }
  
    fHistEventStats->Fill(1); //all events
    
    //Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
    //if(isSelected) {
      fHistEventStats->Fill(2); //triggered + centrality
      Printf("Event accepted");

      //Centrality stuff
      AliCentrality *centrality = fESD->GetCentrality();
      Int_t nCentrality = 0;
      if(fCentralityBins20) 
	nCentrality = centrality->GetCentralityClass5(fCentralityEstimator.Data());
      else
	nCentrality = centrality->GetCentralityClass10(fCentralityEstimator.Data());

      if((nCentrality < 0)||(nCentrality > 19)) return;
      
      if(fAnalysisMode == "TPC") {
	const AliESDVertex *vertex = fESD->GetPrimaryVertexTPC();
	if(vertex) {
	  if(vertex->GetNContributors() > 0) {
	    if(vertex->GetZRes() != 0) {
	      fHistEventStats->Fill(3); //events with a proper vertex
	      if(TMath::Abs(vertex->GetXv()) < fVxMax) {
		if(TMath::Abs(vertex->GetYv()) < fVyMax) {
		  if(TMath::Abs(vertex->GetZv()) < fVzMax) {
		    fHistEventStats->Fill(4); //analyzed events
		    
		    //Printf("Centrality percentile: %lf - Centrality: %d - Total tracks: %d",
		    //centrality->GetCentralityPercentile(fCentralityEstimator.Data()),
		    //nCentrality,fESD->GetNumberOfTracks());
		    		    
		    Int_t nAcceptedTracks = 0;
		    TObjArray *gTrackArray = 0;
		    if(fESDtrackCuts)
                      gTrackArray = fESDtrackCuts->GetAcceptedTracks(fESD,kTRUE);
                    if(gTrackArray) {
                      nAcceptedTracks = gTrackArray->GetEntries();

		      AliESDtrack* track = 0;
		      // Track loop to fill a pT spectrum
		      for (Int_t iTracks = 0; iTracks < nAcceptedTracks; iTracks++) {
			track = dynamic_cast<AliESDtrack *>(gTrackArray->At(iTracks));
			if (!track) {
			  printf("ERROR: Could not receive track %d\n", iTracks);
			  continue;
			}
			
			Short_t gCharge = track->Charge();
			if(gCharge > 0) nPlus += 1;
			if(gCharge < 0) nMinus += 1;

			
		      }//track loop
		    }//TObjArray valid object
		    //if((nCentrality >= 1)&&(nCentrality <= 20)) {
		    
		    fHistCentrality->Fill(nCentrality);
		    fHistNPlusNMinus[nCentrality-1]->Fill(nPlus,nMinus);
		    fHistNMult[nCentrality-1]->Fill(nPlus+nMinus);
		    //}
		  }//Vz cut
		}//Vy cut
	      }//Vx cut
	    }//Vz resolution
	  }//number of contributors
	}//valid vertex
      }//TPC analysis mode
      //}//physics selection
  }//ESD analysis level
  
  //MC analysis
  if(fAnalysisType == "MC") {
    AliMCEvent* mcEvent = MCEvent();
    if (!mcEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }
    AliStack* stack = mcEvent->Stack();
    if (!stack) {
      Printf("ERROR: Could not retrieve MC stack");
      return;
    }
    
    fHistEventStats->Fill(1); 
    fHistEventStats->Fill(2); 
    fHistEventStats->Fill(3); 
    fHistEventStats->Fill(4); //analyzed events
    for (Int_t iTracks = 0; iTracks < mcEvent->GetNumberOfTracks(); iTracks++) {
      AliVParticle* track = mcEvent->GetTrack(iTracks);
      if (!track) {
	Printf("ERROR: Could not receive track %d (mc loop)", iTracks);
	continue;
      }
      
      if(!stack->IsPhysicalPrimary(iTracks)) continue;
      if((track->Pt() < 0.3) && (track->Pt() > 1.5)) continue;
      if(TMath::Abs(track->Eta()) > 0.8) continue;
      
      Short_t gCharge = track->Charge();
      if(gCharge > 0) nPlus += 1;
      if(gCharge < 0) nMinus += 1;
    }//particle loop
    fHistNPlusNMinusMC->Fill(nPlus,nMinus);
    fHistNMultMC->Fill(nPlus+nMinus);
    
  }//MC analysis level


  PostData(2, fOutputList);
}      

//________________________________________________________________________
void AliEbyEFluctuationAnalysisTask::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}

//___________________________________________________________

Bool_t
AliEbyEFluctuationAnalysisTask::HasTPCPID(AliESDtrack *track) const
{
  /*
   * has TPC PID
   */

  /* check PID signal */
  if (track->GetTPCsignal() <= 0. ||
      track->GetTPCsignalN() == 0 ||
      !track->GetInnerParam()) return kFALSE;
  return kTRUE;
}
  
//___________________________________________________________

Bool_t
AliEbyEFluctuationAnalysisTask::HasTOFPID(AliESDtrack *track) const
{
  /*
   * has TOF PID
   */

  /* check TOF matched track */
  if (!(track->GetStatus() & AliESDtrack::kTOFout)||
      !(track->GetStatus() & AliESDtrack::kTIME)) return kFALSE;
  /* check integrated length */
  if (track->GetIntegratedLength() < 350.) return kFALSE;
  return kTRUE;
}

//___________________________________________________________

Double_t
AliEbyEFluctuationAnalysisTask::MakeTPCPID(AliESDtrack *track, Double_t *nSigma) const
{
  /*
   * make TPC PID
   * returns measured dEdx if PID available, otherwise -1.
   * fills nSigma array with TPC nsigmas for e, mu, pi, K, p
   */
  
  /* check TPC PID */
  if (!HasTPCPID(track)) return -1.;

  /* get TPC info */
  Double_t ptpc = track->GetInnerParam() ? track->GetInnerParam()->P() : 0.;
  Double_t dEdx = track->GetTPCsignal();
  Double_t dEdxN = track->GetTPCsignalN();

  /* loop over particles */
  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
    Double_t bethe = fESDpid->GetTPCResponse().GetExpectedSignal(ptpc, (AliPID::EParticleType)ipart);
    Double_t diff = dEdx - bethe;
    Double_t sigma = fESDpid->GetTPCResponse().GetExpectedSigma(ptpc, dEdxN, (AliPID::EParticleType)ipart);
    nSigma[ipart] = diff / sigma;
  }

  return dEdx;
}

//___________________________________________________________

Double_t
AliEbyEFluctuationAnalysisTask::MakeTOFPID(AliESDtrack *track, Double_t *nSigma) const
{
  /*
   * make TOF PID
   * returns measured beta if PID available, otherwise -1.
   * fills nSigma array with TOF nsigmas for e, mu, pi, K, p
   */
  
  /* check TOF PID */
  if (!HasTOFPID(track)) return -1.;

  /* get TOF info */
  Double_t p = track->P();
  Double_t time = track->GetTOFsignal() - fESDpid->GetTOFResponse().GetStartTime(p);
  Double_t length = track->GetIntegratedLength();
  Double_t timei[5];
  track->GetIntegratedTimes(timei);

  /* loop over particles */
  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
    Double_t timez = time - timei[ipart];
    Double_t sigma = fESDpid->GetTOFResponse().GetExpectedSigma(p, timei[ipart], AliPID::ParticleMass(ipart));
    nSigma[ipart] = timez / sigma;
  }

  return length / (time * 2.99792457999999984e-02);
}

//___________________________________________________________

void
AliEbyEFluctuationAnalysisTask::MakePID(AliESDtrack *track, Bool_t *pidFlag) const
{
  /*
   * make PID
   * fills PID QA plots
   * fills pidFlag array with PID flags for e, mu, pi, K, p
   */

  /* cut definitions
     (better put them as static variables so they can be changed from outside) */
  Double_t fgTPCPIDmomcut[AliPID::kSPECIES] = {0., 0., 0.5, 0.5, 0.7};
  Double_t fgTPCPIDsigmacut[AliPID::kSPECIES] = {0., 0., 2., 2., 2.};
  Double_t fgTPCPIDcompcut[AliPID::kSPECIES] = {0., 0., 3., 3., 3.};
  Double_t fgTOFPIDmomcut[AliPID::kSPECIES] = {0., 0., 1.5, 1.5, 2.0};
  Double_t fgTOFPIDsigmacut[AliPID::kSPECIES] = {0., 0., 2., 2., 2.};

  /* make pid and check if available */
  Double_t p = track->P();
  Double_t pt = track->Pt();
  Double_t ptpc = track->GetInnerParam() ? track->GetInnerParam()->P() : 0.;
  Double_t nsigmaTPC[AliPID::kSPECIES];
  Double_t nsigmaTOF[AliPID::kSPECIES];
  Double_t dEdx = MakeTPCPID(track, nsigmaTPC);
  Double_t beta = MakeTOFPID(track, nsigmaTPC);
  Bool_t hasTPCPID = dEdx > 0.;
  Bool_t hasTOFPID = beta > 0.;
  
  /* fill qa histos */
  if (hasTPCPID)
    fHistoTPCdEdx->Fill(ptpc, dEdx);
  if (hasTOFPID)
    fHistoTOFbeta->Fill(p, beta);
  
  /* loop over species */
  for (Int_t ipart = 0; ipart < AliPID::kSPECIES; ipart++) {
    
    /* reset PID flag */
    pidFlag[ipart] = kFALSE;
    
    /* fill qa histos */
    if (hasTPCPID)
      fHistoNSigmaTPC[ipart]->Fill(pt, nsigmaTPC[ipart]);
    if (hasTOFPID)
      fHistoNSigmaTOF[ipart]->Fill(pt, nsigmaTOF[ipart]);
    
    /* combined PID cuts */
    if (hasTPCPID && hasTOFPID) {
      if (pt < fgTOFPIDmomcut[ipart] &&
	  TMath::Abs(nsigmaTOF[ipart]) < fgTOFPIDsigmacut[ipart] &&
	  TMath::Abs(nsigmaTPC[ipart]) < fgTPCPIDcompcut[ipart]) 
	pidFlag[ipart] = kTRUE;
      
    }
    /* TPC-only PID cuts */
    else if (hasTPCPID && !hasTOFPID) {
      if (pt < fgTPCPIDmomcut[ipart] &&
	  TMath::Abs(nsigmaTPC[ipart]) < fgTPCPIDsigmacut[ipart]) 
	pidFlag[ipart] = kTRUE;
    }
    /* TOF-only PID cuts */
    else if (!hasTPCPID && hasTOFPID) {
      if (pt < fgTOFPIDmomcut[ipart] &&
	  TMath::Abs(nsigmaTOF[ipart]) < fgTOFPIDsigmacut[ipart]) 
	pidFlag[ipart] = kTRUE;
    }
    
  } /* end of loop over species */

}

//___________________________________________________________

void
AliEbyEFluctuationAnalysisTask::InitPID(AliESDEvent *event)
{
  /*
   * init PID
   */

  /* create ESDpid object if not there yet */
  if (!fESDpid) {
    /* instance object */
    Bool_t mcFlag = kFALSE; /*** WARNING: check whether is MC ***/
    fESDpid = new AliESDpid(mcFlag);
    /* set OADB path */
    fESDpid->SetOADBPath("$ALICE_ROOT/OADB");
  }

  /* init ESD PID */
  Int_t recoPass = 2; /*** WARNING: need to set the recoPass somehow better ***/
  fESDpid->InitialiseEvent(event, recoPass); /* warning: this apparently sets TOF time        
					      * resolution to some hardcoded value,     
					      * therefore we have to set correct   
					      * resolution value after this call */

  /* set TOF resolution */
  Double_t tofReso = 85.; /* ps */ /*** WARNING: need to set tofReso somehow better ***/
  fESDpid->GetTOFResponse().SetTimeResolution(tofReso);
  
}



