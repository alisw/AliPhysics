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

//-----------------------------------------------------------------
//         AliAnalysisTaskSpectraAOD class
//-----------------------------------------------------------------

#include "TChain.h"
#include "TTree.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliVParticle.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskSpectraAOD.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisDataContainer.h"
#include "AliSpectraAODHistoManager.h"
#include "AliSpectraAODTrackCuts.h"
#include "AliSpectraAODEventCuts.h"
#include "AliCentrality.h"
#include "TProof.h"
#include "AliPID.h"
#include "AliVEvent.h"
#include "AliPIDResponse.h"
#include "AliStack.h"
#include "AliSpectraAODPID.h"
#include <iostream>




using namespace AliSpectraNameSpace;
using namespace std;

ClassImp(AliAnalysisTaskSpectraAOD) // EX1 This stuff tells root to implement the streamer, inspector methods etc (we discussed about it today)

//________________________________________________________________________
AliAnalysisTaskSpectraAOD::AliAnalysisTaskSpectraAOD(const char *name) : AliAnalysisTaskSE(name), fAOD(0), fHistMan(0), fTrackCuts(0), fEventCuts(0),  fPID(0), fIsMC(0)
{
  // Default constructor

  DefineInput(0, TChain::Class());
  DefineOutput(1, AliSpectraAODHistoManager::Class());
  DefineOutput(2, AliSpectraAODEventCuts::Class());
  DefineOutput(3, AliSpectraAODTrackCuts::Class());
  DefineOutput(4, AliSpectraAODPID::Class());

}
//________________________________________________________________________
//________________________________________________________________________
void AliAnalysisTaskSpectraAOD::UserCreateOutputObjects()
{
  // create output objects
  fHistMan = new AliSpectraAODHistoManager("SpectraHistos");

  if (!fTrackCuts) AliFatal("Track Cuts should be set in the steering macro");
  if (!fEventCuts) AliFatal("Event Cuts should be set in the steering macro");
  if (!fPID)       AliFatal("PID object should be set in the steering macro");

  PostData(1, fHistMan  );
  PostData(2, fEventCuts);
  PostData(3, fTrackCuts);
  PostData(4, fPID      );

}
//________________________________________________________________________
void AliAnalysisTaskSpectraAOD::UserExec(Option_t *)
{
  // main event loop
  fAOD = dynamic_cast<AliAODEvent*>(fInputEvent);
  if (strcmp(fAOD->ClassName(), "AliAODEvent"))
    {
      AliFatal("Not processing AODs");
    }
    
  //check on centrality distribution
  fHistMan->GetPtHistogram("CentCheck")->Fill(fAOD->GetCentrality()->GetCentralityPercentile("V0M"),fAOD->GetHeader()->GetCentralityP()->GetCentralityPercentileUnchecked("V0M"));
  
  if(!fEventCuts->IsSelected(fAOD,fTrackCuts))return;//event selection
  
  //AliCentrality fAliCentral*;
  //   if ((fAOD->GetCentrality())->GetCentralityPercentile("V0M") > 5.) return;
  
  // First do MC to fill up the MC particle array, such that we can use it later
  TClonesArray *arrayMC = 0;
  if (fIsMC)
    {
      arrayMC = (TClonesArray*) fAOD->GetList()->FindObject(AliAODMCParticle::StdBranchName());
      if (!arrayMC) {
	AliFatal("Error: MC particles branch not found!\n");
      }
      Int_t nMC = arrayMC->GetEntries();
      for (Int_t iMC = 0; iMC < nMC; iMC++)
	{
	  AliAODMCParticle *partMC = (AliAODMCParticle*) arrayMC->At(iMC);
	  if(!partMC->Charge()) continue;//Skip neutrals
	  if(TMath::Abs(partMC->Eta()) > fTrackCuts->GetEta()){//charged hadron are filled inside the eta acceptance
	    fHistMan->GetPtHistogram(kHistPtGen)->Fill(partMC->Pt(),partMC->IsPhysicalPrimary());
	  }
	  //rapidity cut
	  if(TMath::Abs(partMC->Y())   > fTrackCuts->GetY()  ) continue;	    
	  // check for true PID + and fill P_t histos 
	  Int_t charge = partMC->Charge() > 0 ? kChPos : kChNeg ;
	  Int_t id = fPID->GetParticleSpecie(partMC);
	  if(id != kSpUndefined) {
	    fHistMan->GetHistogram2D(kHistPtGenTruePrimary,id,charge)->Fill(partMC->Pt(),partMC->IsPhysicalPrimary());
	  }
	}
    }
  
  //main loop on tracks
  for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++) {
    AliAODTrack* track = fAOD->GetTrack(iTracks);
    if (!fTrackCuts->IsSelected(track)) continue;
    
    fPID->FillQAHistos(fHistMan, track, fTrackCuts);
    
    // //cut on q vectors track-by-track FIXME
    // if ((qPos>fTrackCuts->GetQvecMin() && qPos<fTrackCuts->GetQvecMax() && track->Eta()<0) || (qNeg>fTrackCuts->GetQvecMin() && qNeg<fTrackCuts->GetQvecMax() && track->Eta()>=0)){
    
    //calculate DCA for AOD track
    Double_t d[2], covd[3];
    AliAODTrack* track_clone=(AliAODTrack*)track->Clone("track_clone"); // need to clone because PropagateToDCA updates the track parameters
    Bool_t isDCA = track_clone->PropagateToDCA(fAOD->GetPrimaryVertex(),fAOD->GetMagneticField(),9999.,d,covd);
    delete track_clone;
    if(!isDCA)d[0]=-999;
    
    fHistMan->GetPtHistogram(kHistPtRec)->Fill(track->Pt(),d[0]);  // PT histo
    
    // get identity and charge
    Int_t idRec  = fPID->GetParticleSpecie(track, fTrackCuts);
    
    Int_t charge = track->Charge() > 0 ? kChPos : kChNeg;
    
    // Fill histograms, only if inside y and nsigma acceptance
    if(idRec != kSpUndefined && fTrackCuts->CheckYCut ((AODParticleSpecies_t)idRec))fHistMan->GetHistogram2D(kHistPtRecSigma,idRec,charge)->Fill(track->Pt(),d[0]);
    //can't put a continue because we still have to fill allcharged primaries, done later
    
    /* MC Part */
    if (arrayMC) {
      AliAODMCParticle *partMC = (AliAODMCParticle*) arrayMC->At(TMath::Abs(track->GetLabel()));
      if (!partMC) { 
	AliError("Cannot get MC particle");
	continue; 
      }
      // Check if it is primary, secondary from material or secondary from weak decay
      Bool_t isPrimary           = partMC->IsPhysicalPrimary();
      Bool_t isSecondaryMaterial = kFALSE; 
      Bool_t isSecondaryWeak     = kFALSE; 
      if(!isPrimary) {
	Int_t mfl=-999,codemoth=-999;
	Int_t indexMoth=partMC->GetMother(); // FIXME ignore fakes? TO BE CHECKED, on ESD is GetFirstMother()
	if(indexMoth>=0){//is not fake
	  AliAODMCParticle* moth = (AliAODMCParticle*) arrayMC->At(indexMoth);
	  codemoth = TMath::Abs(moth->GetPdgCode());
	  mfl = Int_t (codemoth/ TMath::Power(10, Int_t(TMath::Log10(codemoth))));
	}
	if(mfl==3) isSecondaryWeak     = kTRUE; // add if(partMC->GetStatus() & kPDecay)? FIXME
	else       isSecondaryMaterial = kTRUE;
      }
      
      if (isPrimary)fHistMan->GetPtHistogram(kHistPtRecPrimary)->Fill(track->Pt(),d[0]);  // PT histo of primaries
      
      //nsigma cut (reconstructed nsigma)
      if(idRec == kSpUndefined) continue;
      
      // rapidity cut (reconstructed pt and identity)
      if(!fTrackCuts->CheckYCut ((AODParticleSpecies_t)idRec)) continue;
      
      // Get true ID
      Int_t idGen     = fPID->GetParticleSpecie(partMC);
      //if(TMath::Abs(partMC->Y())   > fTrackCuts->GetY()  ) continue;	    // FIXME: do we need a rapidity cut on the generated?
      // Fill histograms for primaries
      
      if (idRec == idGen) fHistMan->GetHistogram2D(kHistPtRecTrue,  idGen, charge)->Fill(track->Pt(),d[0]); 
      
      if (isPrimary) {
	fHistMan                    ->GetHistogram2D(kHistPtRecSigmaPrimary, idRec, charge)->Fill(track->Pt(),d[0]); 
	if(idGen != kSpUndefined) {
	  fHistMan                    ->GetHistogram2D(kHistPtRecPrimary,      idGen, charge)->Fill(track->Pt(),d[0]);
	  if (idRec == idGen) fHistMan->GetHistogram2D(kHistPtRecTruePrimary,  idGen, charge)->Fill(track->Pt(),d[0]); 
	}
      }
      //25th Apr - Muons are added to Pions -- FIXME
      if ( partMC->PdgCode() == 13 && idRec == kSpPion) { 
	fHistMan->GetPtHistogram(kHistPtRecTrueMuonPlus)->Fill(track->Pt(),d[0]); 
	if(isPrimary)
	  fHistMan->GetPtHistogram(kHistPtRecTruePrimaryMuonPlus)->Fill(track->Pt(),d[0]); 
      }
      if ( partMC->PdgCode() == -13 && idRec == kSpPion) { 
	fHistMan->GetPtHistogram(kHistPtRecTrueMuonMinus)->Fill(track->Pt(),d[0]); 
	if (isPrimary) {
	  fHistMan->GetPtHistogram(kHistPtRecTruePrimaryMuonMinus)->Fill(track->Pt(),d[0]); 
	}
      }
      
      ///..... END FIXME
      
      // Fill secondaries
      if(isSecondaryWeak    )  fHistMan->GetHistogram2D(kHistPtRecSigmaSecondaryWeakDecay, idRec, charge)->Fill(track->Pt(),d[0]);
      if(isSecondaryMaterial)  fHistMan->GetHistogram2D(kHistPtRecSigmaSecondaryMaterial , idRec, charge)->Fill(track->Pt(),d[0]);
      
    }//end if(arrayMC)
  } // end loop on tracks
  
  PostData(1, fHistMan  );
  PostData(2, fEventCuts);
  PostData(3, fTrackCuts);
  PostData(4, fPID      );
}

//_________________________________________________________________
void   AliAnalysisTaskSpectraAOD::Terminate(Option_t *)
{
  // Terminate
}
