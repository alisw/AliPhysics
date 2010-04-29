/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id:$ */

#include <TROOT.h>
#include <TSystem.h>
#include <TChain.h>
#include <TFile.h>
#include <TList.h>
#include <TH1I.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TTree.h>
#include <TBranch.h>
#include <TRandom.h>

#include "AliAnalysisTaskUE.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODHandler.h"
#include "AliStack.h"
#include "AliAODJet.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliKFVertex.h"

#include "AliGenPythiaEventHeader.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliInputEventHandler.h"
#include "AliStack.h"
#include "AliLog.h"

//__________________________________________________________________
// Analysis class for Underlying Event studies
//
// Look for correlations on the tranverse regions to 
// the leading charged jet
//
// This class needs as input AOD with track and Jets.
// The output is a list of histograms
//
// AOD can be either connected to the InputEventHandler  
// for a chain of AOD files 
// or 
// to the OutputEventHandler
// for a chain of ESD files, so this case class should be 
// in the train after the Jet finder
//
//    Arian.Abrahantes.Quintana@cern.ch 
//    Ernesto.Lopez.Torres@cern.ch
// 


ClassImp( AliAnalysisTaskUE)

////////////////////////////////////////////////////////////////////////


//____________________________________________________________________
AliAnalysisTaskUE:: AliAnalysisTaskUE(const char* name):
AliAnalysisTask(name,""),
fTrigger(0),
fDebug(0),
fDeltaAOD(kFALSE),
fDeltaAODBranch(""),
fAODBranch("jets"),
fArrayJets(0x0),           
fAOD(0x0),            
fAODjets(0x0),
fListOfHistos(0x0),  
fBinsPtInHist(30),     
fMinJetPtInHist(0.),
fMaxJetPtInHist(300.), 
fIsNorm2Area(kTRUE),
fUseMCParticleBranch(kFALSE),
fConstrainDistance(kTRUE),
fMinDistance(0.2),
fSimulateChJetPt(kFALSE),
fUseAliStack(kTRUE),
fnTracksVertex(3),  // QA tracks pointing to principal vertex (= 3 default) 
fZVertex(5.),
fAnaType(1),         
fRegionType(1),
fConeRadius(0.7),
fConePosition(1),
fAreaReg(1.5393), // Pi*0.7*0.7
fUseChPartJet(kFALSE),
fUseChPartMaxPt(kFALSE),
fUseChargeHadrons(kFALSE),
fUseSingleCharge(kFALSE),
fUsePositiveCharge(kTRUE),
fOrdering(1),
fFilterBit(0xFF),
fJetsOnFly(kFALSE),
fChJetPtMin(5.0),
fJet1EtaCut(0.2),
fJet2DeltaPhiCut(2.616),    // 150 degrees
fJet2RatioPtCut(0.8),
fJet3PtCut(15.),
fTrackPtCut(0.),
fTrackEtaCut(0.9),
fAvgTrials(1),
fhNJets(0x0),
fhEleadingPt(0x0),
fhMinRegPtDist(0x0),
fhRegionMultMin(0x0),
fhMinRegAvePt(0x0), 
fhMinRegSumPt(0x0),            
fhMinRegMaxPtPart(0x0),
fhMinRegSumPtvsMult(0x0),
fhdNdEtaPhiDist(0x0),        
fhFullRegPartPtDistVsEt(0x0), 
fhTransRegPartPtDistVsEt(0x0),
fhRegionSumPtMaxVsEt(0x0),
fhRegionMultMax(0x0),         
fhRegionMultMaxVsEt(0x0),     
fhRegionSumPtMinVsEt(0x0), //fhRegionMultMin(0x0),         
fhRegionMultMinVsEt(0x0),
fhRegionAveSumPtVsEt(0x0),
fhRegionDiffSumPtVsEt(0x0),
fhRegionAvePartPtMaxVsEt(0x0),
fhRegionAvePartPtMinVsEt(0x0),
fhRegionMaxPartPtMaxVsEt(0x0),
//fhRegForwardSumPtVsEt(0x0),
//fhRegForwardMultVsEt(0x0),
//fhRegBackwardSumPtVsEt(0x0),
//fhRegBackwardMultVsEt(0x0),
fhRegForwardMult(0x0),
fhRegForwardSumPtvsMult(0x0),
fhRegBackwardMult(0x0),
fhRegBackwardSumPtvsMult(0x0),
fhRegForwardPartPtDistVsEt(0x0),
fhRegBackwardPartPtDistVsEt(0x0),
fhRegTransMult(0x0),
fhRegTransSumPtVsMult(0x0),
fhMinRegSumPtJetPtBin(0x0),
fhMaxRegSumPtJetPtBin(0x0),
fhVertexMult(0x0),
fh1Xsec(0x0),
fh1Trials(0x0),
fSettingsTree(0x0)//,   fhValidRegion(0x0)
{
  // Default constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class());

}

//______________________________________________________________
Bool_t AliAnalysisTaskUE::Notify()
{
  //
  // Implemented Notify() to read the cross sections
  // and number of trials from pyxsec.root
  // Copy from AliAnalysisTaskJFSystematics
  fAvgTrials = 1;
  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  Float_t xsection = 0;
  Float_t ftrials  = 1;
  if(tree){
    TFile *curfile = tree->GetCurrentFile();
    if (!curfile) {
      Error("Notify","No current file");
      return kFALSE;
    }
    if(!fh1Xsec||!fh1Trials){
      Printf("%s%d No Histogram fh1Xsec",(char*)__FILE__,__LINE__);
      return kFALSE;
    }
    AliAnalysisHelperJetTasks::PythiaInfoFromFile(curfile->GetName(),xsection,ftrials); 
    fh1Xsec->Fill("<#sigma>",xsection);

   // construct average trials
   Float_t nEntries = (Float_t)tree->GetTree()->GetEntries();
   if(ftrials>=nEntries && nEntries>0.)fAvgTrials = ftrials/nEntries;
  }
  
  return kTRUE;
}

//____________________________________________________________________
void AliAnalysisTaskUE::ConnectInputData(Option_t* /*option*/)
{
  // Connect the input data  
  
  // We need AOD with tracks and jets.
  // Since AOD can be either connected to the InputEventHandler (input chain fron AOD files)
  // or to the OutputEventHandler ( AOD is create by a previus task in the train )
  // we need to check where it is and get the pointer to AODEvent in the right way
  
  // Delta AODs are also implemented
  
 
  if (fDebug > 1) AliInfo("ConnectInputData() ");
  
  TObject* handler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  
  if( handler && handler->InheritsFrom("AliAODInputHandler") ) { // input AOD
    fAOD = ((AliAODInputHandler*)handler)->GetEvent();
    if (fDebug > 1) AliInfo(" ==== Tracks from AliAODInputHandler");
    // Case when jets are reconstructed on the fly from AOD tracks
    // (the Jet Finder is using the AliJetAODReader) of InputEventHandler
    // and put in the OutputEventHandler AOD. Useful whe you want to reconstruct jets with
    // different parameters to default ones stored in the AOD or to use a different algorithm
    if( fJetsOnFly ) {
      handler = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
      if( handler && handler->InheritsFrom("AliAODHandler") ) {
        fAODjets = ((AliAODHandler*)handler)->GetAOD();
        if (fDebug > 1) AliInfo(" ==== Jets from AliAODHandler (on the fly)");
      }
    } else {
      fAODjets = fAOD;
      if (fDebug > 1) AliInfo(" ==== Jets from AliAODInputHandler");
    }
  } else {  //output AOD
    handler = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
    if( handler && handler->InheritsFrom("AliAODHandler") ) {
      fAOD = ((AliAODHandler*)handler)->GetAOD();
      fAODjets = fAOD;
      if (fDebug > 1) AliInfo(" ==== Tracks and Jets from AliAODHandler");
    } else {
      AliFatal("I can't get any AOD Event Handler");
      return;
    }
  }
}

//____________________________________________________________________
void  AliAnalysisTaskUE::CreateOutputObjects()
{
  // Create the output container
  //
  if (fDebug > 1) AliInfo("CreateOutPutData()");
  //
  //  Histograms

  // OpenFile(0);
  CreateHistos();
  fListOfHistos->SetOwner(kTRUE);  
  PostData(0,fListOfHistos);

}

//____________________________________________________________________
void  AliAnalysisTaskUE::Exec(Option_t */*option*/)
{
  //Trigger selection ************************************************
  /*AliInputEventHandler* inputHandler = (AliInputEventHandler*)
         ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
  if (inputHandler->IsEventSelected()) {
    if (fDebug > 1) AliInfo(" Trigger Selection: event ACCEPTED ... ");
  } else {
    if (fDebug > 1) AliInfo(" Trigger Selection: event REJECTED ... ");
    return;
  }
   */                                     
  //Event selection (vertex) *****************************************

  Int_t nVertex = fAOD->GetNumberOfVertices();
  if( nVertex > 0 ) { // Only one vertex (reject pileup??)
     AliAODVertex* vertex = (AliAODVertex*)fAOD->GetPrimaryVertex();
     TString tpcvertex("TPCVertex");
     if (vertex->GetName() == tpcvertex){
	if (fDebug > 1) AliInfo(Form("Primary vertex selection: %s event REJECTED ...",vertex->GetName()));
	return;
        }
     Int_t nTracksPrim = vertex->GetNContributors();
     Double_t zVertex = vertex->GetZ();
     if (fDebug > 1) AliInfo(Form(" Vertex in = %f with %d particles by  %s data ...",zVertex,nTracksPrim,vertex->GetName()));
     // Select a quality vertex by number of tracks?
     if( nTracksPrim < fnTracksVertex || TMath::Abs(zVertex) > fZVertex ) {
        if (fDebug > 1) AliInfo(" Primary-vertex Selection: event REJECTED ...");
        return;
     }
     if (fDebug > 1) AliInfo(" Primary-vertex Selection: event ACCEPTED...");
  } else {
     if (fDebug > 1) AliInfo(" Primary-vertex Selection: event REJECTED ...");
     return;
  } 
  /*AliKFVertex primVtx(*(fAOD->GetPrimaryVertex()));
  Int_t nTracksPrim=primVtx.GetNContributors();
  if (fDebug > 1) AliInfo(Form(" Primary-vertex Selection: %d",nTracksPrim));
  if(!nTracksPrim){
    if (fDebug > 1) AliInfo(" Primary-vertex Selection: event REJECTED ...");
    return;
  }
  if (fDebug > 1) AliInfo(" Primary-vertex Selection: event ACCEPTED ...");
  */
  // Execute analysis for current event
  //
  if ( fDebug > 3 ) AliInfo( " Processing event..." );
  // fetch the pythia header info and get the trials
  AliMCEventHandler* mcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  Float_t nTrials = 1;
  if (mcHandler) {
       AliMCEvent* mcEvent = mcHandler->MCEvent();
       if (mcEvent) {
         AliGenPythiaEventHeader*  pythiaGenHeader = AliAnalysisHelperJetTasks::GetPythiaEventHeader(mcEvent);
         if(pythiaGenHeader){
           nTrials = pythiaGenHeader->Trials();
         }
      }
    }
  fh1Trials->Fill("#sum{ntrials}",fAvgTrials);
  AnalyseUE();

  // Post the data
  PostData(0, fListOfHistos);
}

//____________________________________________________________________
void  AliAnalysisTaskUE::AnalyseUE()
{
   Double_t const kMyTolerance = 0.0000001; 
  //
  // Look for correlations on the tranverse regions to 
  // the leading charged jet
  // 
  
  // ------------------------------------------------
  // Find Leading Jets 1,2,3 
  // (could be skipped if Jets are sort by Pt...)
  Double_t maxPtJet1 = 0.; 
  Int_t    index1 = -1;
  Double_t maxPtJet2 = 0.;   // jet 2 need for back to back inclusive
  Int_t    index2 = -1;
  Double_t maxPtJet3 = 0.;   // jet 3 need for back to back exclusive
  Int_t    index3 = -1;
  TVector3 jetVect[3];
  Int_t nJets = 0;
  
  
  if( !fUseChPartJet && fAnaType != 4 ) {
    
    // Use AOD Jets
    if(fDeltaAOD){
      if (fDebug > 1) AliInfo(" ==== Jets From  Delta-AODs !");
      if (fDebug > 1) AliInfo(Form(" ====  Reading Branch: %s  ",fDeltaAODBranch.Data()));
 	   fArrayJets = (TClonesArray*)fAODjets->GetList()->FindObject(fDeltaAODBranch.Data());
    } else {
      if (fDebug > 1) AliInfo(" ==== Read Standard-AODs  !");
      if (fDebug > 1) AliInfo(Form(" ====  Reading Branch: %s  ",fAODBranch.Data()));
      fArrayJets = (TClonesArray*)fAODjets->GetList()->FindObject(fAODBranch.Data());
    }
    if (!fArrayJets){
       AliFatal(" No jet-array! ");
       return;
    }
    nJets=fArrayJets->GetEntries();
    //printf("AOD %d jets \n", nJets);

    for( Int_t i=0; i<nJets; ++i ) {
      AliAODJet* jet = (AliAODJet*)fArrayJets->At(i);
      Double_t jetPt = jet->Pt();//*1.666; // FIXME Jet Pt Correction ?????!!!
 
      if( jetPt > maxPtJet1 ) {
	     maxPtJet3 = maxPtJet2; index3 = index2;
	     maxPtJet2 = maxPtJet1; index2 = index1;
	     maxPtJet1 = jetPt; index1 = i;
      } else if( jetPt > maxPtJet2 ) {
	     maxPtJet3 = maxPtJet2; index3 = index2;
	     maxPtJet2 = jetPt; index2 = i;
      } else if( jetPt > maxPtJet3 ) {
	     maxPtJet3 = jetPt; index3 = i;
      }
    }

    if( index1 != -1 ) {
      AliAODJet *jet =(AliAODJet*) fArrayJets->At(index1);
      if(jet)jetVect[0].SetXYZ(jet->Px(), jet->Py(), jet->Pz());
    }
    if( index2 != -1 ) {
      AliAODJet* jet= (AliAODJet*) fArrayJets->At(index2);
      if(jet)jetVect[1].SetXYZ(jet->Px(), jet->Py(), jet->Pz());
    }
    if( index3 != -1 ) {
       AliAODJet* jet = (AliAODJet*) fArrayJets->At(index3);
      if(jet)jetVect[2].SetXYZ(jet->Px(), jet->Py(), jet->Pz());
    }
    
  }

  if( fUseChPartJet )  {

    // Use "Charged Particle Jets"
    TObjArray* jets = FindChargedParticleJets();
    if( jets ) {
      nJets = jets->GetEntriesFast();
      if( nJets > 0 ) {
	     index1 = 0;
	     AliAODJet* jet = (AliAODJet*)jets->At(0);
	     maxPtJet1 = jet->Pt();
	     jetVect[0].SetXYZ(jet->Px(), jet->Py(), jet->Pz());
      }
      if( nJets > 1 ) {
	     index2 = 1;
	     AliAODJet* jet = (AliAODJet*)jets->At(1);
        maxPtJet2 = jet->Pt();
	     jetVect[1].SetXYZ(jet->Px(), jet->Py(), jet->Pz());
      }
      if( nJets > 2 ) {
	     index3 = 2;
	     AliAODJet* jet = (AliAODJet*)jets->At(2);
        maxPtJet3 = jet->Pt();
	     jetVect[2].SetXYZ(jet->Px(), jet->Py(), jet->Pz());
      }
      
      jets->Delete();
      delete jets;
      if (fDebug > 1) AliInfo(" ==== Jets Created in OUR class  !: CDF algorithm ");
    }
  }


  if( fAnaType == 4 )  {

    // Use "Max Pt Charge Particle"
    TObjArray* tracks = SortChargedParticles();
    if( tracks ) {
      nJets = tracks->GetEntriesFast();
      if( nJets > 0 ) {
        index1 = 0;
        AliAODTrack* jet = (AliAODTrack*)tracks->At(0);
        maxPtJet1 = jet->Pt();
        jetVect[0].SetXYZ(jet->Px(), jet->Py(), jet->Pz());
      }
      tracks->Clear();
      delete tracks; 
    }
  }


  fhNJets->Fill(nJets);
  
  if( fDebug > 1 ) {
    if( index1 < 0 ) {
      AliInfo("\n   Skipping Event, not jet found...");
      return;
    } else
      AliInfo(Form("\n   Pt Leading Jet = %6.1f eta=%5.3f ",  maxPtJet1, jetVect[0].Eta() ));
  }
  
  
  // ----------------------------------------------
  // Cut events by jets topology
  // fAnaType:
  //     1 = inclusive,
  //         - Jet1 |eta| < fJet1EtaCut
  //     2 = back to back inclusive
  //         - fulfill case 1
  //         - |Jet1.Phi - Jet2.Phi| > fJet2DeltaPhiCut
  //         - Jet2.Pt/Jet1Pt > fJet2RatioPtCut
  //     3 = back to back exclusive
  //         - fulfill case 2
  //         - Jet3.Pt < fJet3PtCut
  
  if( index1 < 0 || TMath::Abs(jetVect[0].Eta()) > fJet1EtaCut) {
    if( fDebug > 1 ) AliInfo("\n   Skipping Event...Jet1 |eta| > fJet1EtaCut");
    return;
  }
  // back to back inclusive
  if( fAnaType > 1 && fAnaType < 4 && index2 == -1 ) {
    if( fDebug > 1 ) AliInfo("\n   Skipping Event... no second Jet found");
    return;
  }
  if( fAnaType > 1 && fAnaType < 4 && index2 > -1 ) {
    if( TMath::Abs(jetVect[0].DeltaPhi(jetVect[1])) < fJet2DeltaPhiCut ||
       maxPtJet2/maxPtJet1 < fJet2RatioPtCut ) {
      if( fDebug > 1 ) AliInfo("\n   Skipping Event... |Jet1.Phi - Jet2.Phi| < fJet2DeltaPhiCut");
      return;
    }
  }
  // back to back exclusive
  if( fAnaType > 2 && fAnaType < 4 && index3 > -1 ) {
    if( maxPtJet3 > fJet3PtCut ) {
      if( fDebug > 1 ) AliInfo("\n   Skipping Event... Jet3.Pt > fJet3PtCut ");
      return;
    }
  }
  
  //fhEleadingPt->Fill( maxPtJet1 );

  // ----------------------------------------------
  // Find max and min regions
  Double_t sumPtRegionPosit = 0.;
  Double_t sumPtRegionNegat = 0.;
  Double_t sumPtRegionForward = 0;
  Double_t sumPtRegionBackward = 0;
  Double_t maxPartPtRegion  = 0.;
  Int_t    nTrackRegionPosit = 0;
  Int_t    nTrackRegionNegat = 0;
  Int_t    nTrackRegionForward = 0;
  Int_t    nTrackRegionBackward = 0;
  static Double_t const  kPI     = TMath::Pi();
  static Double_t const  kTWOPI  = 2.*kPI;
  static Double_t const  k270rad = 270.*kPI/180.;
  static Double_t const  k120rad = 120.*kPI/180.;
  
  //Area for Normalization Purpose at Display histos
  // Forward and backward
  Double_t normArea = 1.;
  // Transverse
  if (fIsNorm2Area) {
    SetRegionArea(jetVect);
    normArea =  2.*fTrackEtaCut*k120rad ;
  } else fAreaReg = 1.;
  
  if (!fUseMCParticleBranch){
    AliAODVertex* vertex = (AliAODVertex*)fAOD->GetPrimaryVertex();
    fhVertexMult->Fill( vertex->GetNContributors() );
    
    fhEleadingPt->Fill( maxPtJet1 );
    Int_t nTracks = fAOD->GetNTracks();
    if (fDebug > 1) AliInfo(Form(" ==== AOD tracks = %d \n ",nTracks));
    
    for (Int_t ipart=0; ipart<nTracks; ++ipart) {
      
      AliAODTrack* part = fAOD->GetTrack( ipart );
      if (fDebug > 1) AliInfo(Form(" ==== AOD track = %d pt = %f charge = %d \n ",ipart,part->Pt(),part->Charge()));
      if ( !part->TestFilterBit(fFilterBit) ) continue; // track cut selection
      if (!part->IsPrimaryCandidate()) continue; // reject whatever is not linked to collision point
      // PID Selection: Reject everything but hadrons
      Bool_t isHadron = part->GetMostProbablePID()==AliAODTrack::kPion || 
                        part->GetMostProbablePID()==AliAODTrack::kKaon || 
                        part->GetMostProbablePID()==AliAODTrack::kProton;
      if ( fUseChargeHadrons && !isHadron ) continue;
      
      if ( !part->Charge() ) continue; //Only charged
      if ( fUseSingleCharge ) { // Charge selection
        if ( fUsePositiveCharge && part->Charge() < 0.) continue; // keep Positives
        if ( !fUsePositiveCharge && part->Charge() > 0.) continue; // keep Negatives
      }
      
      if ( part->Pt() < fTrackPtCut ) continue;
      if( TMath::Abs(part->Eta()) > fTrackEtaCut ) continue;
      
      TVector3 partVect(part->Px(), part->Py(), part->Pz());
      Bool_t isFlagPart = kTRUE;
      Double_t deltaPhi = jetVect[0].DeltaPhi(partVect)+k270rad;
      if( deltaPhi > kTWOPI )  deltaPhi-= kTWOPI;
      if (fAnaType != 4 ) fhdNdEtaPhiDist->Fill( deltaPhi, maxPtJet1 );
      else if (TMath::Abs(deltaPhi-k270rad) >= kMyTolerance && TMath::Abs(jetVect[0].Eta()-partVect.Eta()) >= kMyTolerance){
         fhdNdEtaPhiDist->Fill( deltaPhi, maxPtJet1 );
         isFlagPart = kFALSE;
      }
      
      fhFullRegPartPtDistVsEt->Fill( part->Pt(), maxPtJet1 );
      
      Int_t region = IsTrackInsideRegion( jetVect, &partVect );  
      
      if (region == 1) {
        if( maxPartPtRegion < part->Pt() ) maxPartPtRegion = part->Pt();
        sumPtRegionPosit += part->Pt();
        nTrackRegionPosit++;
        fhTransRegPartPtDistVsEt->Fill( part->Pt(), maxPtJet1 );
      }
      if (region == -1) {
        if( maxPartPtRegion < part->Pt() ) maxPartPtRegion = part->Pt();
        sumPtRegionNegat += part->Pt();
        nTrackRegionNegat++;
        fhTransRegPartPtDistVsEt->Fill( part->Pt(), maxPtJet1 );
      }
      if (region == 2){ //forward
        sumPtRegionForward += part->Pt();
        nTrackRegionForward++;
        fhRegForwardPartPtDistVsEt->Fill( part->Pt(), maxPtJet1 );
      }
      if (region == -2){ //backward
        sumPtRegionBackward += part->Pt();
        nTrackRegionBackward++;
        fhRegBackwardPartPtDistVsEt->Fill( part->Pt(), maxPtJet1 );
      }
    }//end loop AOD tracks
  }
  else {
    
    // this is the part we only use when we have MC information
    // More than a test for values of it also resumes the reconstruction efficiency of jets
    // As commented bellow if available for the data, we try to pair reconstructed jets with simulated ones
    // afterwards we kept angular variables of MC jet to perform UE analysis over MC particles
    // TODO: Handle Multiple jet environment. 06/2009 just suited for inclusive jet condition ( fAnaType = 1 ) 
    
      AliMCEventHandler* mcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
      if (!mcHandler) {
        Printf("ERROR: Could not retrieve MC event handler");
        return;
      }
      
      AliMCEvent* mcEvent = mcHandler->MCEvent();
      if (!mcEvent) {
        Printf("ERROR: Could not retrieve MC event");
        return;
      }
    AliGenPythiaEventHeader*  pythiaGenHeader = AliAnalysisHelperJetTasks::GetPythiaEventHeader(mcEvent);
    if(!pythiaGenHeader){
      return;
    }
    
    //Get Jets from MC header
    Int_t nPythiaGenJets = pythiaGenHeader->NTriggerJets();
    AliAODJet pythiaGenJets[4];
    TVector3 jetVectnew[4];
    Int_t iCount = 0;
    for(int ip = 0;ip < nPythiaGenJets;++ip){
      if (iCount>3) break;
      Float_t p[4];
      pythiaGenHeader->TriggerJet(ip,p);
      TVector3 tempVect(p[0],p[1],p[2]);
      if ( TMath::Abs(tempVect.Eta())>fJet1EtaCut ) continue;
      pythiaGenJets[iCount].SetPxPyPzE(p[0],p[1],p[2],p[3]);
      jetVectnew[iCount].SetXYZ(pythiaGenJets[iCount].Px(), pythiaGenJets[iCount].Py(), pythiaGenJets[iCount].Pz());
      iCount++;
    }
    
    if (!iCount) return;// no jet in eta acceptance
    
    //Search the index of the nearest MC jet to the leading jet reconstructed from the input data
    Int_t index = 0;
    if (fConstrainDistance){
      Float_t deltaR = 0.;
      Float_t dRTemp = 0.;
      for (Int_t i=0; i<iCount; i++){
         if (!i) {
         dRTemp = jetVectnew[i].DeltaR(jetVect[0]);
         index = i;
         }
         deltaR = jetVectnew[i].DeltaR(jetVect[0]);
         if (deltaR < dRTemp){
         index = i;
         dRTemp = deltaR;
         }
      }
   
      if (jetVectnew[index].DeltaR(jetVect[0]) > fMinDistance) return;
    }
    //Let's add some taste to jet and simulate pt of charged alone 
    //eta and phi are kept as original
    //Play a Normal Distribution
    Float_t random = 1.;  
    if (fSimulateChJetPt){
      while(1){
        random = gRandom->Gaus(0.6,0.25);
        if (random > 0. && random < 1. && 
            (random * jetVectnew[index].Pt()>6.)) break;
      }
    }
    
    //Set new Pt & Fill histogram accordingly
    maxPtJet1 = random * jetVectnew[index].Pt();
    
    
    fhEleadingPt->Fill( maxPtJet1 );

    if (fUseAliStack){//Try Stack Information to perform UE analysis
    
      AliStack* mcStack = mcEvent->Stack();//Load Stack
      Int_t nTracksMC = mcStack->GetNtrack();
      for (Int_t iTracks = 0; iTracks < nTracksMC; iTracks++) {
        //Cuts
        if(!(mcStack->IsPhysicalPrimary(iTracks))) continue;
        
        TParticle* mctrk = mcStack->Particle(iTracks);
        
        Double_t charge = mctrk->GetPDG()->Charge();
        if (charge == 0) continue;
        
        if ( fUseSingleCharge ) { // Charge selection
          if ( fUsePositiveCharge && charge < 0.) continue; // keep Positives
          if ( !fUsePositiveCharge && charge > 0.) continue; // keep Negatives
        }
        
        //Kinematics cuts on particle
        if ((mctrk->Pt() < fTrackPtCut) || (TMath::Abs(mctrk->Eta()) > fTrackEtaCut )) continue;
        
        Bool_t isHadron = TMath::Abs(mctrk->GetPdgCode())==211 ||
                          TMath::Abs(mctrk->GetPdgCode())==2212 ||
                          TMath::Abs(mctrk->GetPdgCode())==321;
        
        if ( fUseChargeHadrons && !isHadron ) continue;
        
        TVector3 partVect(mctrk->Px(), mctrk->Py(), mctrk->Pz());

        Double_t deltaPhi = jetVectnew[index].DeltaPhi(partVect)+k270rad;
        if( deltaPhi > 2.*TMath::Pi() )  deltaPhi-= 2.*TMath::Pi();
        fhdNdEtaPhiDist->Fill( deltaPhi, maxPtJet1 );
        
        fhFullRegPartPtDistVsEt->Fill( mctrk->Pt(), maxPtJet1 );
        
        //We are not interested on stack organization but don't loose track of info
        TVector3 tempVector =  jetVectnew[0];
        jetVectnew[0] = jetVectnew[index];
        jetVectnew[index] = tempVector;
        
        Int_t region = IsTrackInsideRegion( jetVectnew, &partVect );  
        
        if (region == 1) {
          if( maxPartPtRegion < mctrk->Pt() ) maxPartPtRegion = mctrk->Pt();
          sumPtRegionPosit += mctrk->Pt();
          nTrackRegionPosit++;
          fhTransRegPartPtDistVsEt->Fill( mctrk->Pt(), maxPtJet1 );
        }
        if (region == -1) {
          if( maxPartPtRegion < mctrk->Pt() ) maxPartPtRegion = mctrk->Pt();
          sumPtRegionNegat += mctrk->Pt();
          nTrackRegionNegat++;
          fhTransRegPartPtDistVsEt->Fill( mctrk->Pt(), maxPtJet1 );
        }
        if (region == 2){ //forward
          sumPtRegionForward += mctrk->Pt();
          nTrackRegionForward++;
          fhRegForwardPartPtDistVsEt->Fill( mctrk->Pt(), maxPtJet1 );
        }
        if (region == -2){ //backward
          sumPtRegionBackward += mctrk->Pt();
          nTrackRegionBackward++;
          fhRegBackwardPartPtDistVsEt->Fill( mctrk->Pt(), maxPtJet1 );
        }
        
      }  // end loop stack Particles
      
    }else{//Try mc Particle

      TClonesArray* farray = (TClonesArray*)fAOD->FindListObject("mcparticles");
       
      Int_t ntrks = farray->GetEntries();
      if (fDebug>1) AliInfo(Form("In UE MC analysis tracks %d \n",ntrks));
      for(Int_t i =0 ; i < ntrks; i++){   
        AliAODMCParticle* mctrk = (AliAODMCParticle*)farray->At(i);
        //Cuts
        if (!(mctrk->IsPhysicalPrimary())) continue;
        //if (!(mctrk->IsPrimary())) continue;
        
        if (mctrk->Charge() == 0 || mctrk->Charge()==-99) continue;
        
        if (mctrk->Pt() < fTrackPtCut ) continue;
        if( TMath::Abs(mctrk->Eta()) > fTrackEtaCut ) continue;
        
        Bool_t isHadron = TMath::Abs(mctrk->GetPdgCode())==211 ||
                          TMath::Abs(mctrk->GetPdgCode())==2212 ||
                          TMath::Abs(mctrk->GetPdgCode())==321;
        
        if ( fUseChargeHadrons && !isHadron ) continue;
        
        TVector3 partVect(mctrk->Px(), mctrk->Py(), mctrk->Pz());

        Double_t deltaPhi = jetVectnew[index].DeltaPhi(partVect)+k270rad;
        if( deltaPhi > 2.*TMath::Pi() )  deltaPhi-= 2.*TMath::Pi();
        fhdNdEtaPhiDist->Fill( deltaPhi, maxPtJet1 );

        fhFullRegPartPtDistVsEt->Fill( mctrk->Pt(), maxPtJet1 );
        
        //We are not interested on stack organization but don't loose track of info
        TVector3 tempVector =  jetVectnew[0];
        jetVectnew[0] = jetVectnew[index];
        jetVectnew[index] = tempVector;
        
        Int_t region = IsTrackInsideRegion( jetVectnew, &partVect );  
        
        if (region == 1) { //right
          if( maxPartPtRegion < mctrk->Pt() ) maxPartPtRegion = mctrk->Pt();
          sumPtRegionPosit += mctrk->Pt();
          nTrackRegionPosit++;
          fhTransRegPartPtDistVsEt->Fill( mctrk->Pt(), maxPtJet1 );
        }
        if (region == -1) { //left
          if( maxPartPtRegion < mctrk->Pt() ) maxPartPtRegion = mctrk->Pt();
          sumPtRegionNegat += mctrk->Pt();
          nTrackRegionNegat++;
          fhTransRegPartPtDistVsEt->Fill( mctrk->Pt(), maxPtJet1 );
        }
        if (region == 2){ //forward
          sumPtRegionForward += mctrk->Pt();
          nTrackRegionForward++;
          fhRegForwardPartPtDistVsEt->Fill( mctrk->Pt(), maxPtJet1 );
        }
        if (region == -2){ //backward
          sumPtRegionBackward += mctrk->Pt();
          nTrackRegionBackward++;
          fhRegBackwardPartPtDistVsEt->Fill( mctrk->Pt(), maxPtJet1 );
        }
        
      }//end loop AliAODMCParticle tracks
    }
  }
  
  Double_t avePosRegion = (nTrackRegionPosit) ? sumPtRegionPosit/nTrackRegionPosit : 0.;
  Double_t aveNegRegion = (nTrackRegionNegat) ? sumPtRegionNegat/nTrackRegionNegat : 0.;
  if( avePosRegion > aveNegRegion ) {
     FillAvePartPtRegion( maxPtJet1, avePosRegion/fAreaReg, aveNegRegion/fAreaReg );
  } else {
     FillAvePartPtRegion( maxPtJet1, aveNegRegion/fAreaReg, avePosRegion/fAreaReg );
  }
  
  //How quantities will be sorted before Fill Min and Max Histogram
  //  1=Plots will be CDF-like
  //  2=Plots will be Marchesini-like
  //  3=Minimum zone is selected as the one having lowest pt per track 
  if( fOrdering == 1 ) {
    if( sumPtRegionPosit > sumPtRegionNegat ) {
      FillSumPtRegion( maxPtJet1, sumPtRegionPosit/fAreaReg, sumPtRegionNegat/fAreaReg );
    } else {
      FillSumPtRegion( maxPtJet1, sumPtRegionNegat/fAreaReg, sumPtRegionPosit/fAreaReg );
    }
    if (nTrackRegionPosit > nTrackRegionNegat ) {
      FillMultRegion( maxPtJet1, nTrackRegionPosit/fAreaReg, nTrackRegionNegat/fAreaReg, sumPtRegionNegat/fAreaReg );
    } else {
      FillMultRegion( maxPtJet1, nTrackRegionNegat/fAreaReg, nTrackRegionPosit/fAreaReg, sumPtRegionPosit/fAreaReg );
    }
  } else if( fOrdering == 2 ) {
    if (sumPtRegionPosit > sumPtRegionNegat) {
      FillSumPtRegion( maxPtJet1, sumPtRegionPosit/fAreaReg, sumPtRegionNegat/fAreaReg );
      FillMultRegion( maxPtJet1, nTrackRegionPosit/fAreaReg, nTrackRegionNegat/fAreaReg, sumPtRegionNegat/fAreaReg );
    } else {
      FillSumPtRegion( maxPtJet1, sumPtRegionNegat/fAreaReg, sumPtRegionPosit/fAreaReg );
      FillMultRegion( maxPtJet1, nTrackRegionNegat/fAreaReg, nTrackRegionPosit/fAreaReg, sumPtRegionPosit/fAreaReg );
    }
  } else if( fOrdering == 3 ){
     if (avePosRegion > aveNegRegion) {
        FillSumPtRegion( maxPtJet1, sumPtRegionPosit/fAreaReg, sumPtRegionNegat/fAreaReg );
        FillMultRegion( maxPtJet1, nTrackRegionPosit/fAreaReg, nTrackRegionNegat/fAreaReg, sumPtRegionNegat/fAreaReg );
     }else{
        FillSumPtRegion( maxPtJet1, sumPtRegionNegat/fAreaReg, sumPtRegionPosit/fAreaReg );
        FillMultRegion( maxPtJet1, nTrackRegionNegat/fAreaReg, nTrackRegionPosit/fAreaReg, sumPtRegionPosit/fAreaReg );
     }
  }
  fhRegionMaxPartPtMaxVsEt->Fill(maxPtJet1, maxPartPtRegion );
  
  // Compute pedestal like magnitudes
  fhRegionDiffSumPtVsEt->Fill(maxPtJet1, TMath::Abs(sumPtRegionPosit-sumPtRegionNegat)/(2.0*fAreaReg));
  fhRegionAveSumPtVsEt->Fill(maxPtJet1, (sumPtRegionPosit+sumPtRegionNegat)/(2.0*fAreaReg));
  // Transverse as a whole
  fhRegTransMult->Fill( maxPtJet1, nTrackRegionPosit + nTrackRegionNegat, (nTrackRegionPosit + nTrackRegionNegat)/(2.0*fAreaReg));
  fhRegTransSumPtVsMult->Fill(maxPtJet1, nTrackRegionPosit + nTrackRegionNegat , (sumPtRegionNegat + sumPtRegionPosit)/(2.0*fAreaReg) );
  
  // Fill Histograms for Forward and away side w.r.t. leading jet direction
  // Pt dependence
  //fhRegForwardSumPtVsEt->Fill( maxPtJet1, sumPtRegionForward/normArea );
  //fhRegForwardMultVsEt->Fill( maxPtJet1, nTrackRegionForward/normArea );
  //fhRegBackwardSumPtVsEt->Fill( maxPtJet1, sumPtRegionBackward/normArea );
  //fhRegBackwardMultVsEt->Fill( maxPtJet1, nTrackRegionBackward/normArea );
  // Multiplicity dependence
  fhRegForwardMult->Fill(maxPtJet1, nTrackRegionForward, nTrackRegionForward/normArea);
  fhRegForwardSumPtvsMult->Fill(maxPtJet1, nTrackRegionForward,sumPtRegionForward/normArea);
  fhRegBackwardMult->Fill(maxPtJet1, nTrackRegionBackward, nTrackRegionBackward/normArea );
  fhRegBackwardSumPtvsMult->Fill(maxPtJet1, nTrackRegionBackward,sumPtRegionBackward/normArea);
  
}

//____________________________________________________________________
void AliAnalysisTaskUE::FillSumPtRegion( Double_t leadingE, Double_t ptMax, Double_t ptMin  )
{
  // Fill sumPt of control regions
  
  // Max cone
  fhRegionSumPtMaxVsEt->Fill( leadingE, ptMax );
  // Min cone
  fhRegionSumPtMinVsEt->Fill( leadingE, ptMin );
  // MAke distributions for UE comparison with MB data
  fhMinRegSumPt->Fill(ptMin);
  fhMinRegSumPtJetPtBin->Fill(leadingE, ptMin);
  fhMaxRegSumPtJetPtBin->Fill(leadingE, ptMax);
}

//____________________________________________________________________
void AliAnalysisTaskUE::FillAvePartPtRegion( Double_t leadingE, Double_t ptMax, Double_t ptMin  )
{
  // Fill average particle Pt of control regions
  
  // Max cone
  fhRegionAvePartPtMaxVsEt->Fill( leadingE, ptMax );
  // Min cone
  fhRegionAvePartPtMinVsEt->Fill( leadingE, ptMin );
  // MAke distributions for UE comparison with MB data
  fhMinRegAvePt->Fill(ptMin);
}

//____________________________________________________________________
void AliAnalysisTaskUE::FillMultRegion( Double_t leadingE, Double_t nTrackPtmax, Double_t nTrackPtmin, Double_t ptMin  )
{
  // Fill Nch multiplicity of control regions
  
  // Max cone
  fhRegionMultMaxVsEt->Fill( leadingE, nTrackPtmax );
  fhRegionMultMax->Fill( nTrackPtmax );
  // Min cone
  fhRegionMultMinVsEt->Fill( leadingE, nTrackPtmin );
  fhRegionMultMin->Fill( nTrackPtmin );
  // MAke distributions for UE comparison with MB data
  fhMinRegSumPtvsMult->Fill(nTrackPtmin,ptMin);
}

//____________________________________________________________________
Int_t AliAnalysisTaskUE::IsTrackInsideRegion(TVector3 *jetVect, TVector3 *partVect) 
{  
  // return de region in delta phi
  // -1 negative delta phi 
  //  1 positive delta phi
  //  0 outside region
  static const Double_t k60rad  = 60.*TMath::Pi()/180.;
  static const Double_t k120rad = 120.*TMath::Pi()/180.;
  
  Int_t region = 0;
  if( fRegionType == 1 ) {
    if( TMath::Abs(partVect->Eta()) > fTrackEtaCut ) return 0;
    // transverse regions
    if (jetVect[0].DeltaPhi(*partVect) < -k60rad && jetVect[0].DeltaPhi(*partVect) > -k120rad ) region = -1; //left
    if (jetVect[0].DeltaPhi(*partVect) > k60rad && jetVect[0].DeltaPhi(*partVect) < k120rad ) region = 1;    //right
    if (TMath::Abs(jetVect[0].DeltaPhi(*partVect)) < k60rad ) region = 2;    //forward
    if (TMath::Abs(jetVect[0].DeltaPhi(*partVect)) > k120rad ) region = -2; //backward
    
  } else if( fRegionType == 2 ) {
    // Cone regions
    Double_t deltaR = 0.;
    
    TVector3 positVect,negatVect;
    if (fConePosition==1){
      positVect.SetMagThetaPhi(1, 2.*atan(exp(-jetVect[0].Eta())), jetVect[0].Phi()+TMath::PiOver2());
      negatVect.SetMagThetaPhi(1, 2.*atan(exp(-jetVect[0].Eta())), jetVect[0].Phi()-TMath::PiOver2());
    }else if (fConePosition==2){
       if(fAnaType<2) AliFatal("Prevent error in Analysis type there might be only 1 jet. To avoid overflow better Correct UE config");
       positVect.SetMagThetaPhi(1, 2.*atan(exp(-(jetVect[0].Eta()+jetVect[1].Eta())/2.)), jetVect[0].Phi()+TMath::PiOver2());
       negatVect.SetMagThetaPhi(1, 2.*atan(exp(-(jetVect[0].Eta()+jetVect[1].Eta())/2.)), jetVect[0].Phi()-TMath::PiOver2());
    }else if (fConePosition==3){
       if(fAnaType<2) AliFatal("Prevent error in Analysis type there might be only 1 jet. To avoid overflow better Correct UE config");
       Double_t weightEta = jetVect[0].Eta() * jetVect[0].Pt()/(jetVect[0].Pt() + jetVect[1].Pt()) + 
                            jetVect[1].Eta() * jetVect[1].Pt()/(jetVect[0].Pt() + jetVect[1].Pt());
       //Double_t weightEta = jetVect[0].Eta() * jetVect[0].Mag()/(jetVect[0].Mag() + jetVect[1].Mag()) + 
       //                     jetVect[1].Eta() * jetVect[1].Mag()/(jetVect[0].Mag() + jetVect[1].Mag());
       positVect.SetMagThetaPhi(1, 2.*atan(exp(-weightEta)), jetVect[0].Phi()+TMath::PiOver2());
       negatVect.SetMagThetaPhi(1, 2.*atan(exp(-weightEta)), jetVect[0].Phi()-TMath::PiOver2());
    }
    if (TMath::Abs(positVect.DeltaPhi(*partVect)) < fConeRadius ) { 
      region = 1;  
      deltaR = positVect.DrEtaPhi(*partVect);
    } else if (TMath::Abs(negatVect.DeltaPhi(*partVect)) < fConeRadius) { 
      region = -1;  
      deltaR = negatVect.DrEtaPhi(*partVect);
    }
    
    if (deltaR > fConeRadius) region = 0;
    
  } else 
    AliError("Unknow region type");
  
  // For debug (to be removed)
  if( fDebug > 5 && region != 0 ) AliInfo(Form("Delta Phi = %3.2f region = %d \n", jetVect[0].DeltaPhi(*partVect), region)); //fhValidRegion->Fill( partVect->Eta()-jetVect[0].Eta(), jetVect[0].DeltaPhi(*partVect) );
  
  return region;
}

//____________________________________________________________________
TObjArray*  AliAnalysisTaskUE::SortChargedParticles()
{
  //  return an array with all charged particles ordered according to their pT .
  Int_t nTracks = fAOD->GetNTracks();
  if( !nTracks ) return 0;
  TObjArray* tracks = new TObjArray(nTracks);

  for (Int_t ipart=0; ipart<nTracks; ++ipart) {
    AliAODTrack* part = fAOD->GetTrack( ipart );
    if( !part->TestFilterBit( fFilterBit ) ) continue; // track cut selection
    if( !part->Charge() ) continue;
    if( part->Pt() < fTrackPtCut ) continue;
    tracks->AddLast( part );
  }
  QSortTracks( *tracks, 0, tracks->GetEntriesFast() );

  nTracks = tracks->GetEntriesFast();
  if( !nTracks ) return 0;

  return tracks;
}

//____________________________________________________________________
TObjArray*  AliAnalysisTaskUE::FindChargedParticleJets()
{
  // Return a TObjArray of "charged particle jets"
  //
  // Charged particle jet definition from reference:
  //
  // "Charged jet evolution and the underlying event
  //  in proton-antiproton collisions at 1.8 TeV"
  //  PHYSICAL REVIEW D 65 092002, CDF Collaboration
  //
  // We defined "jets" as circular regions in eta-phi space with
  // radius defined by R = sqrt( (eta-eta0)^2 +(phi-phi0)^2 ).
  // Our jet algorithm is as follows:
  //   1- Order all charged particles according to their pT .
  //   2- Start with the highest pT particle and include in the jet all
  //      particles within the radius R=0.7 considering each particle
  //      in the order of decreasing pT and recalculating the centroid
  //      of the jet after each new particle is added to the jet .
  //   3- Go to the next highest pT particle not already included in
  //      a jet and add to the jet all particles not already included in
  //      a jet within R=0.7.
  //   4- Continue until all particles are in a jet.
  // We defined the transverse momentum of the jet to be
  // the scalar pT sum of all the particles within the jet, where pT
  // is measured with respect to the beam axis
  
  //  1 - Order all charged particles according to their pT .
  Int_t nTracks = fAOD->GetNTracks();
  if( !nTracks ) return 0;
  TObjArray tracks(nTracks);
  
  for (Int_t ipart=0; ipart<nTracks; ++ipart) {
    AliAODTrack* part = fAOD->GetTrack( ipart );
    if( !part->TestFilterBit(fFilterBit) ) continue; // track cut selection
    if( !part->Charge() ) continue;
    if( part->Pt() < fTrackPtCut ) continue;
    tracks.AddLast(part);
  }
  QSortTracks( tracks, 0, tracks.GetEntriesFast() );
  
  nTracks = tracks.GetEntriesFast();
  if( !nTracks ) return 0;

  TObjArray *jets = new TObjArray(nTracks);
  TIter itrack(&tracks);
  while( nTracks ) {
    // 2- Start with the highest pT particle ...
    Float_t px,py,pz,pt; 
    AliAODTrack* track = (AliAODTrack*)itrack.Next();
    if( !track ) continue;
    px = track->Px();
    py = track->Py();
    pz = track->Pz();
    pt = track->Pt(); // Use the energy member to store Pt
    jets->AddLast( new TLorentzVector(px, py, pz, pt) );
    tracks.Remove( track );
    TLorentzVector* jet = (TLorentzVector*)jets->Last();
    jet->SetPtEtaPhiE( 1., jet->Eta(), jet->Phi(), pt );
    // 3- Go to the next highest pT particle not already included...
    AliAODTrack* track1;
    Double_t fPt = jet->E();
    while ( (track1  = (AliAODTrack*)(itrack.Next())) ) {
      Double_t tphi = track1->Phi(); // here Phi is from 0 <-> 2Pi
      if (tphi > TMath::Pi()) tphi -= 2. * TMath::Pi(); // convert to  -Pi <-> Pi
      Double_t dphi = TVector2::Phi_mpi_pi(jet->Phi()-tphi);
      Double_t r = TMath::Sqrt( (jet->Eta()-track1->Eta())*(jet->Eta()-track1->Eta()) +
                               dphi*dphi );
      if( r < fConeRadius ) {
        fPt   = jet->E()+track1->Pt();  // Scalar sum of Pt
        // recalculating the centroid
        Double_t eta = jet->Eta()*jet->E()/fPt + track1->Eta()*track1->Pt()/fPt;
        Double_t phi = jet->Phi()*jet->E()/fPt + tphi*track1->Pt()/fPt;
        jet->SetPtEtaPhiE( 1., eta, phi, fPt );
        tracks.Remove( track1 );
      }
    }
    
    tracks.Compress();
    nTracks = tracks.GetEntries();
    //   4- Continue until all particles are in a jet.
    itrack.Reset();
  } // end while nTracks
  
  // Convert to AODjets....
  Int_t njets = jets->GetEntriesFast();
  TObjArray* aodjets = new TObjArray(njets);
  aodjets->SetOwner(kTRUE);
  for(Int_t ijet=0; ijet<njets; ++ijet) {
    TLorentzVector* jet = (TLorentzVector*)jets->At(ijet);
    if (jet->E() < fChJetPtMin) continue;
    Float_t px, py,pz,en; // convert to 4-vector
    px = jet->E() * TMath::Cos(jet->Phi());  // Pt * cos(phi)
    py = jet->E() * TMath::Sin(jet->Phi());  // Pt * sin(phi)
    pz = jet->E() / TMath::Tan(2.0 * TMath::ATan(TMath::Exp(-jet->Eta())));
    en = TMath::Sqrt(px * px + py * py + pz * pz);

    aodjets->AddLast( new AliAODJet(px, py, pz, en) );
  }
  jets->Delete();
  delete jets;
  
  // Order jets according to their pT .
  QSortTracks( *aodjets, 0, aodjets->GetEntriesFast() );
  
  // debug
  if (fDebug>3) AliInfo(Form(" %d Charged jets found\n",njets));
  
  return aodjets;
}

//____________________________________________________________________
void  AliAnalysisTaskUE::QSortTracks(TObjArray &a, Int_t first, Int_t last)
{
  // Sort array of TObjArray of tracks by Pt using a quicksort algorithm.
  
  static TObject *tmp;
  static int i;           // "static" to save stack space
  int j;
  
  while (last - first > 1) {
    i = first;
    j = last;
    for (;;) {
      while (++i < last && ((AliVParticle*)a[i])->Pt() > ((AliVParticle*)a[first])->Pt() )
        ;
      while (--j > first && ((AliVParticle*)a[j])->Pt() < ((AliVParticle*)a[first])->Pt() )
        ;
      if (i >= j)
        break;
      
      tmp  = a[i];
      a[i] = a[j];
      a[j] = tmp;
    }
    if (j == first) {
      ++first;
      continue;
    }
    tmp = a[first];
    a[first] = a[j];
    a[j] = tmp;
    if (j - first < last - (j + 1)) {
      QSortTracks(a, first, j);
      first = j + 1;   // QSortTracks(j + 1, last);
    } else {
      QSortTracks(a, j + 1, last);
      last = j;        // QSortTracks(first, j);
    }
  }
}

//____________________________________________________________________
void AliAnalysisTaskUE::SetRegionArea(TVector3 *jetVect)
{
  Double_t fAreaCorrFactor=0.;
  Double_t deltaEta = 0.;
  if (fRegionType==1) fAreaReg = 2.*fTrackEtaCut*60.*TMath::Pi()/180.;
  else if (fRegionType==2){ 
    deltaEta = 0.9-TMath::Abs(jetVect[0].Eta());
    if (deltaEta>fConeRadius) fAreaReg = TMath::Pi()*fConeRadius*fConeRadius;
    else{
      fAreaCorrFactor = fConeRadius*fConeRadius*TMath::ACos(deltaEta/fConeRadius) -
      (deltaEta*fConeRadius)*TMath::Sqrt( 1. - deltaEta*deltaEta/(fConeRadius*fConeRadius));
      fAreaReg=TMath::Pi()*fConeRadius*fConeRadius-fAreaCorrFactor;
    }
  }else AliWarning("Unknown Rgion Type");
  if (fDebug>10) AliInfo(Form("\n dEta=%5.3f Angle =%5.3f Region Area = %5.3f Corr Factor=%5.4f \n",deltaEta,TMath::ACos(deltaEta/fConeRadius),fAreaReg,fAreaCorrFactor));
}

//____________________________________________________________________
void  AliAnalysisTaskUE::CreateHistos()
{
  fListOfHistos = new TList();
  
  
  fhNJets = new TH1F("fhNJets", "n Jet",  20, 0, 20);
  fhNJets->SetXTitle("# of jets");
  fhNJets->Sumw2();
  fListOfHistos->Add( fhNJets );                 // At(0) 
  
  fhEleadingPt  = new TH1F("hEleadingPt",   "Leading Jet P_{T}",  fBinsPtInHist, fMinJetPtInHist,   fMaxJetPtInHist);
  fhEleadingPt->SetXTitle("P_{T} (GeV/c)");
  fhEleadingPt->SetYTitle("dN/dP_{T}");
  fhEleadingPt->Sumw2();
  fListOfHistos->Add( fhEleadingPt );            // At(1)
  
  fhMinRegPtDist = new TH1F("hMinRegPtDist",   "P_{T} distribution in Min zone",  50,0.,20.);
  fhMinRegPtDist->SetXTitle("P_{T} (GeV/c)");
  fhMinRegPtDist->SetYTitle("dN/dP_{T}");
  fhMinRegPtDist->Sumw2();
  fListOfHistos->Add( fhMinRegPtDist );          // At(2)
  
  fhRegionMultMin = new TH1F("hRegionMultMin",      "N_{ch}^{90, min}",  21, -0.5,   20.5);
  fhRegionMultMin->SetXTitle("N_{ch tracks}");
  fhRegionMultMin->Sumw2();
  fListOfHistos->Add( fhRegionMultMin );         // At(3)            
  
  fhMinRegAvePt = new TH1F("hMinRegAvePt", "#LTp_{T}#GT",  50, 0.,   20.);
  fhMinRegAvePt->SetXTitle("P_{T} (GeV/c)");
  fhMinRegAvePt->Sumw2();
  fListOfHistos->Add( fhMinRegAvePt );           // At(4)
  
  fhMinRegSumPt = new TH1F("hMinRegSumPt", "#Sigma p_{T} ",  50, 0.,   20.);
  fhMinRegSumPt->SetYTitle("Ed^{3}N_{tracks}/dp^{3} (c^{3}/GeV^{2})");  
  fhMinRegSumPt->SetXTitle("#Sigma p_{T} (GeV/c)");
  fhMinRegSumPt->Sumw2();
  fListOfHistos->Add( fhMinRegSumPt );           // At(5)
  
  fhMinRegMaxPtPart = new TH1F("hMinRegMaxPtPart", "max(p_{T})|_{event} ",  50, 0.,   20.);
  fhMinRegMaxPtPart->SetYTitle("Ed^{3}N_{tracks}/dp^{3} (c^{3}/GeV^{2})");  
  fhMinRegMaxPtPart->SetXTitle("p_{T} (GeV/c)");
  fhMinRegMaxPtPart->Sumw2();
  fListOfHistos->Add( fhMinRegMaxPtPart );       // At(6)
  
  fhMinRegSumPtvsMult = new TH1F("hMinRegSumPtvsMult", "#Sigma p_{T} vs. Multiplicity ",  21, -0.5,   20.5);
  fhMinRegSumPtvsMult->SetYTitle("#Sigma p_{T} (GeV/c)");  
  fhMinRegSumPtvsMult->SetXTitle("N_{charge}");
  fhMinRegSumPtvsMult->Sumw2();
  fListOfHistos->Add( fhMinRegSumPtvsMult );     // At(7);
  
  fhdNdEtaPhiDist  = new TH2F("hdNdEtaPhiDist", Form("Charge particle density |#eta|<%3.1f vs #Delta#phi", fTrackEtaCut),  
                               62, 0.,   2.*TMath::Pi(), fBinsPtInHist, fMinJetPtInHist, fMaxJetPtInHist);
  fhdNdEtaPhiDist->SetXTitle("#Delta#phi");
  fhdNdEtaPhiDist->SetYTitle("Leading Jet P_{T}");
  fhdNdEtaPhiDist->Sumw2();
  fListOfHistos->Add( fhdNdEtaPhiDist );        // At(8)
  
  // Can be use to get part pt distribution for differente Jet Pt bins
  fhFullRegPartPtDistVsEt = new TH2F("hFullRegPartPtDistVsEt", Form( "dN/dP_{T} |#eta|<%3.1f vs Leading Jet P_{T}", fTrackEtaCut),
                                     200,0.,100., fBinsPtInHist, fMinJetPtInHist, fMaxJetPtInHist);
  fhFullRegPartPtDistVsEt->SetYTitle("Leading Jet P_{T}");
  fhFullRegPartPtDistVsEt->SetXTitle("p_{T}");
  fhFullRegPartPtDistVsEt->Sumw2();
  fListOfHistos->Add( fhFullRegPartPtDistVsEt );  // At(9) 
  
  // Can be use to get part pt distribution for differente Jet Pt bins
  fhTransRegPartPtDistVsEt = new TH2F("hTransRegPartPtDistVsEt", Form( "dN/dP_{T} in tranvese regions |#eta|<%3.1f vs Leading Jet P_{T}", fTrackEtaCut),  
                                     200,0.,100., fBinsPtInHist, fMinJetPtInHist,   fMaxJetPtInHist);
  fhTransRegPartPtDistVsEt->SetYTitle("Leading Jet P_{T}");
  fhTransRegPartPtDistVsEt->SetXTitle("p_{T}");
  fhTransRegPartPtDistVsEt->Sumw2();
  fListOfHistos->Add( fhTransRegPartPtDistVsEt );  // At(10)
  
  
  fhRegionSumPtMaxVsEt = new TH1F("hRegionSumPtMaxVsEt",  "P_{T}^{90, max} vs Leading Jet P_{T}",  fBinsPtInHist, fMinJetPtInHist,   fMaxJetPtInHist);
  fhRegionSumPtMaxVsEt->SetXTitle("P_{T} (GeV/c)");
  fhRegionSumPtMaxVsEt->Sumw2();
  fListOfHistos->Add( fhRegionSumPtMaxVsEt );      // At(11)
  
  fhRegionSumPtMinVsEt = new TH1F("hRegionSumPtMinVsEt",   "P_{T}^{90, min} vs Leading Jet P_{T}",  fBinsPtInHist, fMinJetPtInHist,   fMaxJetPtInHist);
  fhRegionSumPtMinVsEt->SetXTitle("P_{T} (GeV/c)");
  fhRegionSumPtMinVsEt->Sumw2();
  fListOfHistos->Add( fhRegionSumPtMinVsEt );      // At(12)
  
  fhRegionMultMax = new TH1I("hRegionMultMax",      "N_{ch}^{90, max}",  21, -0.5,   20.5);
  fhRegionMultMax->SetXTitle("N_{ch tracks}");
  fhRegionMultMax->Sumw2();
  fListOfHistos->Add( fhRegionMultMax );           // At(13)
  
  fhRegionMultMaxVsEt = new TH1F("hRegionMultMaxVsEt",  "N_{ch}^{90, max} vs Leading Jet P_{T}",  fBinsPtInHist, fMinJetPtInHist,   fMaxJetPtInHist);
  fhRegionMultMaxVsEt->SetXTitle("E (GeV hRegionAveSumPtVsEt/c)");
  fhRegionMultMaxVsEt->Sumw2();
  fListOfHistos->Add( fhRegionMultMaxVsEt );       // At(14)
  
  fhRegionMultMinVsEt = new TH1F("hRegionMultMinVsEt",  "N_{ch}^{90, min} vs Leading Jet P_{T}",  fBinsPtInHist, fMinJetPtInHist,   fMaxJetPtInHist);
  fhRegionMultMinVsEt->SetXTitle("E (GeV/c)");
  fhRegionMultMinVsEt->Sumw2();
  fListOfHistos->Add( fhRegionMultMinVsEt );      // At(15)
  
  fhRegionAveSumPtVsEt = new TH1F("hRegionAveSumPtVsEt", "(P_{T}^{90, max} + P_{T}^{90, min})/2 vs Leading Jet P_{T}",  fBinsPtInHist, fMinJetPtInHist,   fMaxJetPtInHist);
  fhRegionAveSumPtVsEt->SetXTitle("P_{T} (GeV/c)");
  fhRegionAveSumPtVsEt->Sumw2();
  fListOfHistos->Add( fhRegionAveSumPtVsEt );     // At(16)
  
  fhRegionDiffSumPtVsEt= new TH1F("hRegionPtDiffVsEt", "(P_{T}^{90, max} - P_{T}^{90, min}) vs Leading Jet P_{T}",  fBinsPtInHist, fMinJetPtInHist,   fMaxJetPtInHist);
  fhRegionDiffSumPtVsEt->SetXTitle("P_{T} (GeV/c)");
  fhRegionDiffSumPtVsEt->Sumw2();
  fListOfHistos->Add( fhRegionDiffSumPtVsEt );    // At(17)
  
  fhRegionAvePartPtMaxVsEt = new TH1F("hRegionAvePartPtMaxVsEt", "#LTp_{T}#GT^{90, max} vs Leading Jet P_{T}",  fBinsPtInHist, fMinJetPtInHist,   fMaxJetPtInHist);
  fhRegionAvePartPtMaxVsEt->SetXTitle("P_{T} (GeV/c)");
  fhRegionAvePartPtMaxVsEt->Sumw2();
  fListOfHistos->Add( fhRegionAvePartPtMaxVsEt );  // At(18)
  
  fhRegionAvePartPtMinVsEt = new TH1F("hRegionAvePartPtMinVsEt", "#LTp_{T}#GT^{90, min} vs Leading Jet P_{T}",  fBinsPtInHist, fMinJetPtInHist,   fMaxJetPtInHist);
  fhRegionAvePartPtMinVsEt->SetXTitle("P_{T} (GeV/c)");
  fhRegionAvePartPtMinVsEt->Sumw2();
  fListOfHistos->Add( fhRegionAvePartPtMinVsEt );   // At(19)
  
  fhRegionMaxPartPtMaxVsEt = new TH1F("hRegionMaxPartPtMaxVsEt", "max(p_{T})^{90} vs Leading Jet P_{T}",  fBinsPtInHist, fMinJetPtInHist,   fMaxJetPtInHist);
  fhRegionMaxPartPtMaxVsEt->SetXTitle("P_{T} (GeV/c)");
  fhRegionMaxPartPtMaxVsEt->Sumw2();
  fListOfHistos->Add( fhRegionMaxPartPtMaxVsEt );    // At(20)
  /*
  fhRegForwardSumPtVsEt = new TH1F("hRegForwardSumPtVsEt", "Forward #sum{p_{T}} vs Leading Jet P_{T}",  fBinsPtInHist, fMinJetPtInHist,   fMaxJetPtInHist);
  fhRegForwardSumPtVsEt->SetXTitle("P_{T} (GeV/c)");
  fhRegForwardSumPtVsEt->Sumw2();
  fListOfHistos->Add( fhRegForwardSumPtVsEt );    // At(21)
  
  fhRegForwardMultVsEt = new TH1F("hRegForwardMultVsEt", "Forward #sum{N_{ch}} vs Leading Jet P_{T}",  fBinsPtInHist, fMinJetPtInHist,   fMaxJetPtInHist);
  fhRegForwardMultVsEt->SetXTitle("P_{T} (GeV/c)");
  fhRegForwardMultVsEt->Sumw2();
  fListOfHistos->Add( fhRegForwardMultVsEt );    // At(22)
  
  fhRegBackwardSumPtVsEt = new TH1F("hRegBackwardSumPtVsEt", "Backward #sum{p_{T}} vs Leading Jet P_{T}",  fBinsPtInHist, fMinJetPtInHist,   fMaxJetPtInHist);
  fhRegBackwardSumPtVsEt->SetXTitle("P_{T} (GeV/c)");
  fhRegBackwardSumPtVsEt->Sumw2();
  fListOfHistos->Add( fhRegBackwardSumPtVsEt );    // At(23)
  
  fhRegBackwardMultVsEt = new TH1F("hRegBackwardMultVsEt", "Backward #sum{N_{ch}} vs Leading Jet P_{T}",  fBinsPtInHist, fMinJetPtInHist,   fMaxJetPtInHist);
  fhRegBackwardMultVsEt->SetXTitle("P_{T} (GeV/c)");
  fhRegBackwardMultVsEt->Sumw2();
  fListOfHistos->Add( fhRegBackwardMultVsEt );    // At(24)
                                  */
  fhRegForwardMult = new TH2F("hRegForwardMult",      "N_{ch}^{forward}",  
                              fBinsPtInHist, fMinJetPtInHist, fMaxJetPtInHist, 21, -0.5,   20.5);
  fhRegForwardMult->SetXTitle("N_{ch tracks}");
  fhRegForwardMult->Sumw2();
  fListOfHistos->Add( fhRegForwardMult );           // At(25)
  
  fhRegForwardSumPtvsMult = new TH2F("hRegForwardSumPtvsMult", "Forward #Sigma p_{T} vs. Multiplicity ",  
                                     fBinsPtInHist, fMinJetPtInHist, fMaxJetPtInHist, 21, -0.5,   20.5);
  fhRegForwardSumPtvsMult->SetYTitle("#Sigma p_{T} (GeV/c)");  
  fhRegForwardSumPtvsMult->SetXTitle("N_{charge}");
  fhRegForwardSumPtvsMult->Sumw2();
  fListOfHistos->Add( fhRegForwardSumPtvsMult );     // At(26);
  
  fhRegBackwardMult = new TH2F("hRegBackwardMult",      "N_{ch}^{backward}",  
                               fBinsPtInHist, fMinJetPtInHist, fMaxJetPtInHist, 21, -0.5,   20.5);
  fhRegBackwardMult->SetXTitle("N_{ch tracks}");
  fhRegBackwardMult->Sumw2();
  fListOfHistos->Add( fhRegBackwardMult );           // At(27)
  
  fhRegBackwardSumPtvsMult = new TH2F("hRegBackwardSumPtvsMult", "Backward #Sigma p_{T} vs. Multiplicity ",  
                                      fBinsPtInHist, fMinJetPtInHist, fMaxJetPtInHist, 21, -0.5,   20.5);
  fhRegBackwardSumPtvsMult->SetYTitle("#Sigma p_{T} (GeV/c)");  
  fhRegBackwardSumPtvsMult->SetXTitle("N_{charge}");
  fhRegBackwardSumPtvsMult->Sumw2();
  fListOfHistos->Add( fhRegBackwardSumPtvsMult );     // At(28);
  
  fhRegForwardPartPtDistVsEt = new TH2F("hRegForwardPartPtDistVsEt", Form( "dN/dP_{T} |#eta|<%3.1f vs Leading Jet P_{T}", fTrackEtaCut),
                                       200,0.,100., fBinsPtInHist, fMinJetPtInHist, fMaxJetPtInHist);
  fhRegForwardPartPtDistVsEt->SetYTitle("Leading Jet P_{T}");
  fhRegForwardPartPtDistVsEt->SetXTitle("p_{T}");
  fhRegForwardPartPtDistVsEt->Sumw2();
  fListOfHistos->Add( fhRegForwardPartPtDistVsEt );  // At(29) 
  
  fhRegBackwardPartPtDistVsEt = new TH2F("hRegBackwardPartPtDistVsEt", Form( "dN/dP_{T} |#eta|<%3.1f vs Leading Jet P_{T}", fTrackEtaCut),
                                         200,0.,100., fBinsPtInHist, fMinJetPtInHist, fMaxJetPtInHist);
  fhRegBackwardPartPtDistVsEt->SetYTitle("Leading Jet P_{T}");
  fhRegBackwardPartPtDistVsEt->SetXTitle("p_{T}");
  fhRegBackwardPartPtDistVsEt->Sumw2();
  fListOfHistos->Add( fhRegBackwardPartPtDistVsEt );  // At(30) 
  
  fhRegTransMult  = new TH2F("hRegTransMult",      "N_{ch}^{transv}",  
                             fBinsPtInHist, fMinJetPtInHist, fMaxJetPtInHist, 21, -0.5,   20.5);
  fhRegTransMult->SetXTitle("N_{ch tracks}");
  fhRegTransMult->Sumw2();
  fListOfHistos->Add( fhRegTransMult );           // At(31)
  
  fhRegTransSumPtVsMult = new TH2F("hRegTransSumPtVsMult", "Transverse #Sigma p_{T} vs. Multiplicity ",  
                                   fBinsPtInHist, fMinJetPtInHist, fMaxJetPtInHist, 21, -0.5,   20.5);
  fhRegTransSumPtVsMult->SetYTitle("#Sigma p_{T} (GeV/c)");  
  fhRegTransSumPtVsMult->SetXTitle("N_{charge}");
  fhRegTransSumPtVsMult->Sumw2();
  fListOfHistos->Add( fhRegTransSumPtVsMult );     // At(32);
  
  fhMinRegSumPtJetPtBin = new TH2F("hMinRegSumPtJetPtBin",      "Transverse Min Reg #Sigma p_{T} per jet Pt Bin",  
                                   fBinsPtInHist, fMinJetPtInHist, fMaxJetPtInHist, 50, 0.,   20.);
  fhMinRegSumPtJetPtBin->SetXTitle("Leading Jet P_{T}");
  fhMinRegSumPtJetPtBin->Sumw2();
  fListOfHistos->Add( fhMinRegSumPtJetPtBin );           // At(33)
  
  fhMaxRegSumPtJetPtBin = new TH2F("hMaxRegSumPtJetPtBin",      "Transverse Max Reg #Sigma p_{T} per jet Pt Bin",  
                                   fBinsPtInHist, fMinJetPtInHist, fMaxJetPtInHist, 50, 0.,   20.);
  fhMaxRegSumPtJetPtBin->SetXTitle("Leading Jet P_{T}");
  fhMaxRegSumPtJetPtBin->Sumw2();
  fListOfHistos->Add( fhMaxRegSumPtJetPtBin );           // At(34)
  
  fhVertexMult = new TH1F("hVertexMult",      "Multiplicity in Main Vertex", 81, -0.5 , 80.5);
  fhVertexMult->SetXTitle("Main Vertex Multiplicity");
  fhVertexMult->Sumw2();
  fListOfHistos->Add( fhVertexMult ); //At(35)
  
  fh1Xsec = new TProfile("fh1Xsec","xsec from pyxsec.root",1,0,1); 
  fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");
  fh1Xsec->Sumw2();
  fListOfHistos->Add( fh1Xsec );            //At(36)
  
  fh1Trials = new TH1F("fh1Trials","trials from pyxsec.root",1,0,1);
  fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
  fh1Trials->Sumw2();
  fListOfHistos->Add( fh1Trials ); //At(37)
  
  fSettingsTree   = new TTree("UEAnalysisSettings","Analysis Settings in UE estimation");
  fSettingsTree->Branch("fFilterBit", &fFilterBit,"FilterBit/I");
  fSettingsTree->Branch("fTrigger", &fTrigger,"TriggerFlag/I");
  fSettingsTree->Branch("fConeRadius", &fConeRadius,"Rad/D");
  fSettingsTree->Branch("fJet1EtaCut", &fJet1EtaCut, "LeadJetEtaCut/D");
  fSettingsTree->Branch("fJet2DeltaPhiCut", &fJet2DeltaPhiCut, "DeltaPhi/D");
  fSettingsTree->Branch("fJet2RatioPtCut", &fJet2RatioPtCut, "Jet2Ratio/D");
  fSettingsTree->Branch("fJet3PtCut", &fJet3PtCut, "Jet3PtCut/D");
  fSettingsTree->Branch("fTrackPtCut", &fTrackPtCut, "TrackPtCut/D");
  fSettingsTree->Branch("fTrackEtaCut", &fTrackEtaCut, "TrackEtaCut/D");
  fSettingsTree->Branch("fAnaType", &fAnaType, "Ana/I");        
  fSettingsTree->Branch("fRegionType", &fRegionType,"Reg/I");
  fSettingsTree->Branch("fOrdering", &fOrdering,"OrderMeth/I");
  fSettingsTree->Branch("fUseChPartJet", &fUseChPartJet,"UseChPart/O");
  fSettingsTree->Branch("fUseChargeHadrons", &fUseChargeHadrons,"UseChHadrons/O");
  fSettingsTree->Branch("fUseSingleCharge", &fUseSingleCharge,"UseSingleCh/O");
  fSettingsTree->Branch("fUsePositiveCharge", &fUsePositiveCharge,"UsePositiveCh/O");
  fSettingsTree->Fill();

  
  fListOfHistos->Add( fSettingsTree );    // At(38)
  
  /*   
   // For debug region selection
   fhValidRegion = new TH2F("hValidRegion", "dN/d#eta/d#phi",      
   20, -2.,2., 62, -TMath::Pi(),   TMath::Pi());
   fhValidRegion->SetYTitle("#Delta#phi");
   fhValidRegion->Sumw2();
   fListOfHistos->Add( fhValidRegion );  // At(15)
   */
}



//____________________________________________________________________
void  AliAnalysisTaskUE::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  
  //Save Analysis Settings
  //  WriteSettings();
  
  // Normalize histos to region area TODO: 
  // Normalization done at Analysis, taking into account 
  // area variations on per-event basis (cone case)
   
   //HIGH WARNING!!!!!: DO NOT SCALE ANY OF THE ORIGINAL HISTOGRAMS
   //MAKE A COPY, DRAW IT, And later sacale that copy. CAF Issue!!!!!
  
  // Draw some Test plot to the screen
  if (!gROOT->IsBatch()) {
    
    // Update pointers reading them from the output slot
    fListOfHistos = dynamic_cast<TList*> (GetOutputData(0));
    if( !fListOfHistos ) {
      AliError("Histogram List is not available");
      return;
    }
    fhNJets              = (TH1F*)fListOfHistos->At(0);
    fhEleadingPt         = (TH1F*)fListOfHistos->At(1);
    fhdNdEtaPhiDist      = (TH2F*)fListOfHistos->At(8);
    fhRegionSumPtMaxVsEt = (TH1F*)fListOfHistos->At(11);
    fhRegionSumPtMinVsEt = (TH1F*)fListOfHistos->At(12);
    fhRegionMultMaxVsEt  = (TH1F*)fListOfHistos->At(14);
    fhRegionMultMinVsEt  = (TH1F*)fListOfHistos->At(15);
    fhRegionAveSumPtVsEt = (TH1F*)fListOfHistos->At(16);
    //fhRegForwardSumPtVsEt = (TH1F*)fListOfHistos->At(21);
    //fhRegBackwardSumPtVsEt = (TH1F*)fListOfHistos->At(23);
    
    //fhValidRegion  = (TH2F*)fListOfHistos->At(21);
    
    // Canvas
    TCanvas* c1 = new TCanvas("c1",Form("sumPt dist (%s)", GetTitle()),60,60,1100,700);
    c1->Divide(2,2);
    c1->cd(1);
    TH1F *h1r = new TH1F("hRegionEtvsSumPtMax" , "", fBinsPtInHist,  fMinJetPtInHist, fMaxJetPtInHist);
    h1r->Divide(fhRegionSumPtMaxVsEt,fhEleadingPt,1,1);
    //h1r->Scale( areafactor );
    h1r->SetMarkerStyle(20);
    h1r->SetXTitle("P_{T} of Leading Jet (GeV/c)");
    h1r->SetYTitle("P_{T}^{90, max}");
    h1r->DrawCopy("p");
    
    c1->cd(2);
    
    TH1F *h2r = new TH1F("hRegionEtvsSumPtMin" , "", fBinsPtInHist,  fMinJetPtInHist, fMaxJetPtInHist);
    h2r->Divide(fhRegionSumPtMinVsEt,fhEleadingPt,1,1);
    //h2r->Scale( areafactor );
    h2r->SetMarkerStyle(20);
    h2r->SetXTitle("P_{T} of Leading Jet (GeV/c)");
    h2r->SetYTitle("P_{T}^{90, min}");
    h2r->DrawCopy("p");
    
    c1->cd(3);
    TH1F *h4r = new TH1F("hRegionEtvsDiffPt" , "", fBinsPtInHist,  fMinJetPtInHist, fMaxJetPtInHist);
    //TH1F *h41r = new TH1F("hRegForwvsDiffPt" , "", fBinsPtInHist,  fMinJetPtInHist, fMaxJetPtInHist);
    //TH1F *h42r = new TH1F("hRegBackvsDiffPt" , "", fBinsPtInHist,  fMinJetPtInHist, fMaxJetPtInHist);
    //h41r->Divide(fhRegForwardSumPtVsEt,fhEleadingPt,1,1);
    //h42r->Divide(fhRegBackwardSumPtVsEt,fhEleadingPt,1,1);
    h4r->Divide(fhRegionAveSumPtVsEt,fhEleadingPt,1,1);
    //h4r->Scale(2.); // make average
    //h4r->Scale( areafactor );
    h4r->SetYTitle("#DeltaP_{T}^{90}");
    h4r->SetXTitle("P_{T} of Leading Jet (GeV/c)");
    h4r->SetMarkerStyle(20);
    h4r->DrawCopy("p");
    /*h41r->SetXTitle("P_{T} of Leading Jet (GeV/c)");
    h41r->SetMarkerStyle(22);
    h41r->DrawCopy("p same");
    h42r->SetXTitle("P_{T} of Leading Jet (GeV/c)");
    h42r->SetMarkerStyle(23);
    h42r->SetMarkerColor(kRed);
    h42r->DrawCopy("p same");
    */
    c1->cd(4);
    TH1F *h5r = new TH1F("hRegionMultMaxVsEtleading",   "",  fBinsPtInHist, fMinJetPtInHist,   fMaxJetPtInHist);
    TH1F *h6r = new TH1F("hRegionMultMinVsEtleading",   "",  fBinsPtInHist, fMinJetPtInHist,   fMaxJetPtInHist);
    
    h5r->Divide(fhRegionMultMaxVsEt,fhEleadingPt,1,1);
    h6r->Divide(fhRegionMultMinVsEt,fhEleadingPt,1,1);
    //h5r->Scale( areafactor );
    //h6r->Scale( areafactor );
    h5r->SetYTitle("N_{Tracks}^{90}");
    h5r->SetXTitle("P_{T} of Leading Jet (GeV/c)");
    h5r->SetMarkerStyle(20);
    h5r->DrawCopy("p");
    h6r->SetMarkerStyle(21);
    h6r->SetMarkerColor(2);
    h6r->SetYTitle("N_{Tracks}^{90}");
    h6r->SetXTitle("P_{T} of Leading Jet (GeV/c)");
    h6r->DrawCopy("p same");
    c1->Update();
    
    //Get Normalization
    fh1Xsec                = (TProfile*)fListOfHistos->At(21);
    fh1Trials              = (TH1F*)fListOfHistos->At(22);
    
    Double_t xsec = fh1Xsec->GetBinContent(1);
    Double_t ntrials = fh1Trials->GetBinContent(1);
    Double_t normFactor = xsec/ntrials;
    if(fDebug > 1)Printf("xSec %f nTrials %f Norm %f \n",xsec,ntrials,normFactor);
    
    
    
    TCanvas* c2 = new TCanvas("c2","Jet Pt dist",160,160,1200,800);
    TH1 * copy = 0;
    c2->Divide(2,2);
    c2->cd(1);
    fhEleadingPt->SetMarkerStyle(20);
    fhEleadingPt->SetMarkerColor(2);
    //if( normFactor > 0.) fhEleadingPt->Scale(normFactor);
    //fhEleadingPt->Draw("p");
    copy = fhEleadingPt->DrawCopy("p");
    if( normFactor > 0.) copy->Scale(normFactor);
    gPad->SetLogy();
    
    c2->cd(2);
    Int_t xbin1 = fhdNdEtaPhiDist->GetYaxis()->FindFixBin(fMinJetPtInHist);
    Int_t xbin2 = fhdNdEtaPhiDist->GetYaxis()->FindFixBin(fMaxJetPtInHist);
    TH1D* dNdEtaPhiDistAllJets = fhdNdEtaPhiDist->ProjectionX("dNdEtaPhiDistAllJets",xbin1,xbin2);
    dNdEtaPhiDistAllJets->SetMarkerStyle(20);
    dNdEtaPhiDistAllJets->SetMarkerColor(2);
    dNdEtaPhiDistAllJets->DrawCopy("p");
    gPad->SetLogy();
    
    c2->cd(3);      
    fhNJets->DrawCopy();
    
    //c2->cd(4);      
    //fhValidRegion->DrawCopy("p");
    
    //fhTransRegPartPtDist  = (TH1F*)fListOfHistos->At(2);
    fhRegionMultMin          = (TH1F*)fListOfHistos->At(3);
    fhMinRegAvePt     = (TH1F*)fListOfHistos->At(4);
    fhMinRegSumPt     = (TH1F*)fListOfHistos->At(5);   
    //fhMinRegMaxPtPart     = (TH1F*)fListOfHistos->At(6);
    fhMinRegSumPtvsMult   = (TH1F*)fListOfHistos->At(7);
    
    
    // Canvas
    TCanvas* c3 = new TCanvas("c3"," pT dist",160,160,1200,800);
    c3->Divide(2,2);
    c3->cd(1);
    /*fhTransRegPartPtDist->SetMarkerStyle(20);
     fhTransRegPartPtDist->SetMarkerColor(2); 
     fhTransRegPartPtDist->Scale(areafactor/fhTransRegPartPtDist->GetEntries());
     fhTransRegPartPtDist->DrawCopy("p");
     gPad->SetLogy();
     */
    c3->cd(2); 
    fhMinRegSumPt->SetMarkerStyle(20);
    fhMinRegSumPt->SetMarkerColor(2);  
    //fhMinRegSumPt->Scale(areafactor);
    fhMinRegSumPt->DrawCopy("p");
    gPad->SetLogy();
    
    c3->cd(3);
    fhMinRegAvePt->SetMarkerStyle(20);
    fhMinRegAvePt->SetMarkerColor(2);  
    //fhMinRegAvePt->Scale(areafactor);
    fhMinRegAvePt->DrawCopy("p");
    gPad->SetLogy();
    
    c3->cd(4);
    
    TH1F *h7r = new TH1F("hRegionMultMinVsMult",   "",  21, -0.5,   20.5);
    
    h7r->Divide(fhMinRegSumPtvsMult,fhRegionMultMin,1,1);
    
    h7r->SetMarkerStyle(20);
    h7r->SetMarkerColor(2);   
    h7r->DrawCopy("p");
    
    c3->Update();
    
    
    /*    c2->cd(3);      
     fhValidRegion->SetMarkerStyle(20);
     fhValidRegion->SetMarkerColor(2);
     fhValidRegion->DrawCopy("p");
     */    
    c2->Update();
  } else {
    AliInfo(" Batch mode, not histograms will be shown...");
  }
  
  if( fDebug > 1 ) 
    AliInfo("End analysis");
  
}

void  AliAnalysisTaskUE::WriteSettings()
{ 
  if (fDebug>5){
    AliInfo(" All Analysis Settings in Saved Tree");
    fSettingsTree->Scan();
  }
}
