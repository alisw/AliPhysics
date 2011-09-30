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
#include <iostream>

using namespace AliSpectraNameSpace;
using namespace std;

ClassImp(AliAnalysisTaskSpectraAOD) // EX1 This stuff tells root to implement the streamer, inspector methods etc (we discussed about it today)

//________________________________________________________________________
AliAnalysisTaskSpectraAOD::AliAnalysisTaskSpectraAOD(const char *name) : AliAnalysisTaskSE(name), fAOD(0), fHistMan(0), fTrackCuts(0), fEventCuts(0), fIsMC(0), fPIDResponse(0), fNSigmaPID(0), fYCut(0)
{
   // Default constructor
   fNSigmaPID = 3.;
   fYCut = .5;
   DefineInput(0, TChain::Class());
   DefineOutput(1, AliSpectraAODHistoManager::Class());
   DefineOutput(2, AliSpectraAODEventCuts::Class());
   DefineOutput(3, AliSpectraAODTrackCuts::Class());

}
//________________________________________________________________________
Bool_t AliAnalysisTaskSpectraAOD::CheckYCut(AODParticleSpecies_t species, AliAODTrack* track) const
{
    // check if the rapidity is within the set range
    // note: masses are hardcoded for now. we could look them up in the pdg database, but that would mean accecing it 100k+ times per run ...
    Double_t y;
    if (species == kProton) { y = track->Y(9.38271999999999995e-01); }
    if ( species == kKaon ) { y = track->Y(4.93676999999999977e-01); }
    if ( species == kPion)  { y = track->Y(1.39570000000000000e-01); }
    if (TMath::Abs(y) > fYCut || y < -998.) return kFALSE;
    return kTRUE;
}
//____________________________________________________________________________
Bool_t AliAnalysisTaskSpectraAOD::CheckYCut(AliAODMCParticle* particle) const
{
    // check if the rapidity is within the set range
    Double_t y = particle->Y();
    if (TMath::Abs(y) > fYCut || y < -998.) return kFALSE;
    return kTRUE;
}
//________________________________________________________________________
void AliAnalysisTaskSpectraAOD::UserCreateOutputObjects()
{
  // create output objects
    fHistMan = new AliSpectraAODHistoManager("SpectraHistos");

   if (!fTrackCuts) AliFatal("Track Cuts should be set in the steering macro");
   if (!fEventCuts) AliFatal("Event Cuts should be set in the steering macro");

   AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
   AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
   fPIDResponse = inputHandler->GetPIDResponse();

   PostData(1, fHistMan);
   PostData(2, fEventCuts);
   PostData(3, fTrackCuts);

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

   if(!fPIDResponse) {
     AliError("Cannot get pid response");
     return;
   }

   if (!fEventCuts->IsSelected(fAOD)) return;

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
         fHistMan->GetPtHistogram(kHistPtGen)->Fill(partMC->Pt());
         
	 if(TMath::Abs(partMC->Eta()) > fTrackCuts->GetEta()) continue;
         // check for true PID + and fill P_t histos 
         if (partMC->IsPhysicalPrimary() && CheckYCut(partMC) ) {// only primary vertices and y cut satisfied
            if ( partMC->PdgCode() == 2212) { fHistMan->GetHistogram(kHistPtGenTruePrimaryProtonPlus)->Fill(partMC->Pt()); } 
            if ( partMC->PdgCode() == -2212) { fHistMan->GetHistogram(kHistPtGenTruePrimaryProtonMinus)->Fill(partMC->Pt()); } 
            if ( partMC->PdgCode() == 321)  { fHistMan->GetHistogram(kHistPtGenTruePrimaryKaonPlus)->Fill(partMC->Pt()); } 
            if ( partMC->PdgCode() == -321) { fHistMan->GetHistogram(kHistPtGenTruePrimaryKaonMinus)->Fill(partMC->Pt()); } 
            if ( partMC->PdgCode() == 211)  { fHistMan->GetHistogram(kHistPtGenTruePrimaryPionPlus)->Fill(partMC->Pt()); } 
            if ( partMC->PdgCode() == -211) { fHistMan->GetHistogram(kHistPtGenTruePrimaryPionMinus)->Fill(partMC->Pt()); } 
         }
      }
   }
   
   for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++)
   {

      AliAODTrack* track = fAOD->GetTrack(iTracks);
      if (!fTrackCuts->IsSelected(track)) continue;

      fHistMan->GetHistogram(kHistPtRec)->Fill(track->Pt());  // PT histo
      fHistMan->GetPIDHistogram(kHistPIDTPC)->Fill(track->P(), track->GetTPCsignal()); // PID histo
      fHistMan->GetPIDHistogram(kHistPIDTPCPt)->Fill(track->Pt(), track->GetTPCsignal());
      
      // Sigma identify particles (for both MC and regular data):
      // note: as of now, very crude identification. how many sigmas are allowed? what to do with high pt?
      AliVParticle *inEvHMain = dynamic_cast<AliVParticle *>(track);
      Double_t nsigmaTPCkKaon = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kKaon)); 
      Double_t nsigmaTPCkProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kProton));
      Double_t nsigmaTPCkPion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kPion)); 
      if( ( nsigmaTPCkKaon < nsigmaTPCkPion ) && ( nsigmaTPCkKaon < nsigmaTPCkProton )) { 
	
          if ((nsigmaTPCkKaon > fNSigmaPID) || (!CheckYCut(kKaon, track) ) ) continue;
          if ( track->Charge() > 0 ) { fHistMan->GetHistogram(kHistPtRecSigmaKaonPlus)->Fill(track->Pt()); } 
          else { fHistMan->GetHistogram(kHistPtRecSigmaKaonMinus)->Fill(track->Pt()); } 
      }
      if( ( nsigmaTPCkProton < nsigmaTPCkKaon ) && ( nsigmaTPCkProton < nsigmaTPCkPion ) ) {
          if ( nsigmaTPCkProton > fNSigmaPID || (!CheckYCut(kProton, track) ) )  continue;
          if ( track->Charge() > 0 ) { fHistMan->GetHistogram(kHistPtRecSigmaProtonPlus)->Fill(track->Pt()); }
          else { fHistMan->GetHistogram(kHistPtRecSigmaProtonMinus)->Fill(track->Pt()); }
      }
      if( (nsigmaTPCkPion < nsigmaTPCkProton ) && ( nsigmaTPCkPion < nsigmaTPCkKaon ) ) {
          if (nsigmaTPCkPion > fNSigmaPID || (!CheckYCut(kPion, track) ) ) continue;
          if ( track->Charge() > 0 )  { fHistMan->GetHistogram(kHistPtRecSigmaPionPlus)->Fill(track->Pt()); }
          else  { fHistMan->GetHistogram(kHistPtRecSigmaPionMinus)->Fill(track->Pt()); }
      }
      /* MC Part */
      if (arrayMC) {
         AliAODMCParticle *partMC = (AliAODMCParticle*) arrayMC->At(TMath::Abs(track->GetLabel()));
         if (!partMC) { 
             AliError("Cannot get MC particle");
             continue; }
         if (CheckYCut(partMC)) {
         // primaries, true pid
         if ( partMC->PdgCode() == 2212) { fHistMan->GetHistogram(kHistPtRecTrueProtonPlus)->Fill(partMC->Pt()); 
            if (partMC->IsPhysicalPrimary()) {fHistMan->GetHistogram(kHistPtRecTruePrimaryProtonPlus)->Fill(partMC->Pt()); }}
         if ( partMC->PdgCode() == -2212) { fHistMan->GetHistogram(kHistPtRecTrueProtonMinus)->Fill(partMC->Pt()); 
            if (partMC->IsPhysicalPrimary()) {fHistMan->GetHistogram(kHistPtRecTruePrimaryProtonMinus)->Fill(partMC->Pt()); }}
         if ( partMC->PdgCode() == 321) { fHistMan->GetHistogram(kHistPtRecTrueKaonPlus)->Fill(partMC->Pt()); 
            if (partMC->IsPhysicalPrimary()) {fHistMan->GetHistogram(kHistPtRecTruePrimaryKaonPlus)->Fill(partMC->Pt()); }}
         if ( partMC->PdgCode() == -321) { fHistMan->GetHistogram(kHistPtRecTrueKaonMinus)->Fill(partMC->Pt()); 
            if (partMC->IsPhysicalPrimary()) {fHistMan->GetHistogram(kHistPtRecTruePrimaryKaonMinus)->Fill(partMC->Pt()); }}
         if ( partMC->PdgCode() == 211) { fHistMan->GetHistogram(kHistPtRecTruePionPlus)->Fill(partMC->Pt()); 
            if (partMC->IsPhysicalPrimary()) {fHistMan->GetHistogram(kHistPtRecTruePrimaryPionPlus)->Fill(partMC->Pt()); }}
         if ( partMC->PdgCode() == -211) { fHistMan->GetHistogram(kHistPtRecTruePionMinus)->Fill(partMC->Pt());  
            if (partMC->IsPhysicalPrimary()) {fHistMan->GetHistogram(kHistPtRecTruePrimaryPionMinus)->Fill(partMC->Pt()); }}
         }

         // primaries, sigma pid
         if (partMC->IsPhysicalPrimary()) { 
            if( ( nsigmaTPCkKaon < nsigmaTPCkPion ) && ( nsigmaTPCkKaon < nsigmaTPCkProton ) ) {
               if ( (nsigmaTPCkKaon > fNSigmaPID ) || (!CheckYCut(kKaon, track) ) ) continue; 
               if ( track->Charge() > 0 ) { fHistMan->GetHistogram(kHistPtRecSigmaPrimaryKaonPlus)->Fill(track->Pt()); } 
               else { fHistMan->GetHistogram(kHistPtRecSigmaPrimaryKaonMinus)->Fill(track->Pt()); } 
            }
            if( ( nsigmaTPCkProton < nsigmaTPCkKaon ) && ( nsigmaTPCkProton < nsigmaTPCkPion ) ) {
               if ( (nsigmaTPCkProton > fNSigmaPID ) || (!CheckYCut(kProton, track) ) ) continue;
               if ( track->Charge() > 0 ) { fHistMan->GetHistogram(kHistPtRecSigmaPrimaryProtonPlus)->Fill(track->Pt()); }
               else { fHistMan->GetHistogram(kHistPtRecSigmaPrimaryProtonMinus)->Fill(track->Pt()); }
            }
            if( (nsigmaTPCkPion < nsigmaTPCkProton ) && ( nsigmaTPCkPion < nsigmaTPCkKaon ) ) {
               if ( ( nsigmaTPCkPion > fNSigmaPID )  || (!CheckYCut(kPion, track) ) ) continue;
               if ( track->Charge() > 0 )  { fHistMan->GetHistogram(kHistPtRecSigmaPrimaryPionPlus)->Fill(track->Pt()); }
               else  { fHistMan->GetHistogram(kHistPtRecSigmaPrimaryPionMinus)->Fill(track->Pt()); }
            }
         }
         //secondaries
         else {
            if( ( nsigmaTPCkKaon < nsigmaTPCkPion ) && ( nsigmaTPCkKaon < nsigmaTPCkProton ) ) { 
               if ( (nsigmaTPCkKaon > fNSigmaPID )  || (!CheckYCut(kKaon, track) ) ) continue;
               if ( track->Charge() > 0 ) { fHistMan->GetHistogram(kHistPtRecSigmaSecondaryKaonPlus)->Fill(track->Pt()); } 
               else { fHistMan->GetHistogram(kHistPtRecSigmaSecondaryKaonMinus)->Fill(track->Pt()); } 
            }
            if( ( nsigmaTPCkProton < nsigmaTPCkKaon ) && ( nsigmaTPCkProton < nsigmaTPCkPion ) ) {
               if ( (nsigmaTPCkProton > fNSigmaPID )  || (!CheckYCut(kProton, track) ) ) continue;
               if ( track->Charge() > 0 ) { fHistMan->GetHistogram(kHistPtRecSigmaSecondaryProtonPlus)->Fill(track->Pt()); }
               else { fHistMan->GetHistogram(kHistPtRecSigmaSecondaryProtonMinus)->Fill(track->Pt()); }
            }
            if( (nsigmaTPCkPion < nsigmaTPCkProton ) && ( nsigmaTPCkPion < nsigmaTPCkKaon ) ) {
               if ( ( nsigmaTPCkPion > fNSigmaPID )  || (!CheckYCut(kPion, track) ) ) continue;
               if ( track->Charge() > 0 )  { fHistMan->GetHistogram(kHistPtRecSigmaSecondaryPionPlus)->Fill(track->Pt()); }
               else  { fHistMan->GetHistogram(kHistPtRecSigmaSecondaryPionMinus)->Fill(track->Pt()); }
            }
         }
      }
   }

   PostData(1, fHistMan);
   PostData(2, fEventCuts);
   PostData(3, fTrackCuts);
}
  
//_________________________________________________________________
void   AliAnalysisTaskSpectraAOD::Terminate(Option_t *)
{
   // Terminate
}
