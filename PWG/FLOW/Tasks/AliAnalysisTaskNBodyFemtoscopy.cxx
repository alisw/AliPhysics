/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/************************************** 
* Femtoscopy task for N bodies * 
*Laura Serksnyte
*laura.serksnyte@cern.ch
**************************************/ 
  
#include "Riostream.h"
#include "AliAnalysisTaskNBodyFemtoscopy.h"
#include "AliLog.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "TFile.h"
#include "AliMultSelection.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskNBodyFemtoscopy)

//================================================================================================================

AliAnalysisTaskNBodyFemtoscopy::AliAnalysisTaskNBodyFemtoscopy(const char *name, Bool_t useParticleWeights): 
 AliAnalysisTaskSE(name), 
 fHistList(NULL),
 // Control histograms:
 fControlHistogramsList(NULL),
 fPtHist(NULL),
 fCentralityHist(NULL),
 fNumberOfTracksHist(NULL),
 fCentralityHistVzCut(NULL),
 fPtHistEtaCut(NULL),
 fPtHistEtaCutPTCut(NULL),
 fPtHistEtaCutPTCutPhiCut(NULL),
 fNbinsPt(1000),
 fMinBinPt(0.),
 fMaxBinPt(10.),
 fCentralityEstimator("V0M"),
 fNbinsCentrality(100),
 fMinCentrality(0.),
 fMaxCentrality(100.),
 // Final results:
 fFinalResultsList(NULL),
 fRejectEventsNoPrimaryVertex(kTRUE),
 fCutOnVertexZ(kTRUE),
  // Global track cuts:
 fApplyCommonTrackCuts(kTRUE),
 fFilterBit(128)
 
 
 {
  // Constructor.
 
  AliDebug(2,"AliAnalysisTaskNBodyFemtoscopy::AliAnalysisTaskNBodyFemtoscopy(const char *name, Bool_t useParticleWeights)");

  // Base list:
  fHistList = new TList();
  fHistList->SetName("FemtoNBody");
  fHistList->SetOwner(kTRUE);

  // Initialize all arrays:
  this->InitializeArrays();

  // Define input and output slots here
  // Input slot #0 works with an AliFlowEventSimple
  //DefineInput(0, AliFlowEventSimple::Class());  
  // Input slot #1 is needed for the weights input file:
  //if(useParticleWeights)
  //{
  // DefineInput(1, TList::Class());   
  //}  
  // Output slot #0 is reserved              
  // Output slot #1 writes into a TList container

  DefineOutput(1, TList::Class());  

  if(useParticleWeights)
  {
   // not needed for the time being
  }

} // AliAnalysisTaskNBodyFemtoscopy::AliAnalysisTaskNBodyFemtoscopy(const char *name, Bool_t useParticleWeights): 

//================================================================================================================

AliAnalysisTaskNBodyFemtoscopy::AliAnalysisTaskNBodyFemtoscopy(): 
 AliAnalysisTaskSE(),
 fHistList(NULL),
 // Control histograms:
 fControlHistogramsList(NULL),
 fPtHist(NULL),
 fCentralityHist(NULL),
 fNumberOfTracksHist(NULL),
 fCentralityHistVzCut(NULL),
 fPtHistEtaCut(NULL),
 fPtHistEtaCutPTCut(NULL),
 fPtHistEtaCutPTCutPhiCut(NULL),
 fNbinsPt(1000),
 fMinBinPt(0.),
 fMaxBinPt(10.),
 fCentralityEstimator("V0M"),
 fNbinsCentrality(100),
 fMinCentrality(0.),
 fMaxCentrality(100.),
 // Final results:
 fFinalResultsList(NULL),
 fRejectEventsNoPrimaryVertex(kTRUE),
 fCutOnVertexZ(kTRUE),
  // Global track cuts:
 fApplyCommonTrackCuts(kTRUE),
 fFilterBit(128)
{
  // Dummy constructor.
 
  AliDebug(2,"AliAnalysisTaskNBodyFemtoscopy::AliAnalysisTaskNBodyFemtoscopy()");
  this->InitializeArrays();

} // AliAnalysisTaskNBodyFemtoscopy::AliAnalysisTaskNBodyFemtoscopy():

//================================================================================================================

AliAnalysisTaskNBodyFemtoscopy::~AliAnalysisTaskNBodyFemtoscopy()
{
 // Destructor.

 if(fHistList) delete fHistList;
  
} // AliAnalysisTaskNBodyFemtoscopy::~AliAnalysisTaskNBodyFemtoscopy()

//================================================================================================================

void AliAnalysisTaskNBodyFemtoscopy::UserCreateOutputObjects() 
{
 // Called at every worker node to initialize.

 // a) Trick to avoid name clashes, part 1;
 // b) Book and nest all lists;
 // c) Book all objects;
 // *) Trick to avoid name clashes, part 2.
  
 // a) Trick to avoid name clashes, part 1:
 Bool_t oldHistAddStatus = TH1::AddDirectoryStatus(); 
 TH1::AddDirectory(kFALSE);

 // b) Book and nest all lists:
 this->BookAndNestAllLists();

 // c) Book all objects:
 this->BookControlHistograms();
 this->BookFinalResultsHistograms();

 // *) Trick to avoid name clashes, part 2:
 TH1::AddDirectory(oldHistAddStatus);

 PostData(1,fHistList);

} // void AliAnalysisTaskNBodyFemtoscopy::UserCreateOutputObjects() 

//================================================================================================================

void AliAnalysisTaskNBodyFemtoscopy::UserExec(Option_t *) 
{
 // Main loop (called for each event).
 // a) Get pointer to AOD event, chech multiplicity, apply event cuts;
 // b) Start analysis over AODs;
 // c) Reset event-by-event objects;
 // d) PostData.

 // a) Get pointer to AOD event:
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(InputEvent()); // from TaskSE
 if(!aAOD){return;}
 // Check Multiplicity
 AliMultSelection *ams = (AliMultSelection*)aAOD->FindListObject("MultSelection");
 if(!ams){return;}
 if(ams->GetMultiplicityPercentile(Form("%s", fCentralityEstimator.Data())) >= fMinCentrality && ams->GetMultiplicityPercentile(Form("%s", fCentralityEstimator.Data())) < fMaxCentrality)
 {
  fCentralityHist->Fill(ams->GetMultiplicityPercentile(Form("%s", fCentralityEstimator.Data())));
 }
 else
 {
  return; 
 }
 //Apply other event cuts
 if(!this->CommonEventCuts(aAOD, ams->GetMultiplicityPercentile(Form("%s", fCentralityEstimator.Data())))){
  return;
 }

 // b) Start analysis over AODs:
 Int_t nTracks = aAOD->GetNumberOfTracks(); // number of all tracks in current event 
 fNumberOfTracksHist->Fill(nTracks);
 for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks
 {
  AliAODTrack *aTrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack)); // getting a pointer to a track
  if(!aTrack){ continue;} // protection against NULL pointers

  Double_t pt = aTrack->Pt(); // Pt
  // Fill pt control histogram before cuts
  fPtHist->Fill(pt); // filling pt distribution
  // Apply analysis cuts, fill histograms inside the cut function
  if(!this->CommonTrackCuts(aTrack)){continue;}

 } // for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks
  
 // c) Reset event-by-event objects:
 // ...

 // d) PostData:
 PostData(1,fHistList);

} // void AliAnalysisTaskNBodyFemtoscopy::UserExec(Option_t *)




//================================================================================================================

Bool_t AliAnalysisTaskNBodyFemtoscopy::CommonEventCuts(AliVEvent *ave, Float_t centrality)
{
 // Apply event cuts.

 // a) Check which event type
 AliMCEvent *aMC = dynamic_cast<AliMCEvent*>(ave);
 AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(ave);
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(ave);

 if(aMC)
 {
  // TBI
 }
 else if(aESD)
 {
  // TBI
 }
 else if(aAOD)
 {
  // a) Cuts on AliAODEvent:
  if(fRejectEventsNoPrimaryVertex && !aAOD->GetPrimaryVertex()) return kFALSE;

  // b) Cuts on AliAODVertex:
  AliAODVertex *avtx = (AliAODVertex*)aAOD->GetPrimaryVertex();
  if(fCutOnVertexZ)
  {
   if(avtx->GetZ() < fVertexZ[0]) return kFALSE;
   if(avtx->GetZ() > fVertexZ[1]) return kFALSE;
  }

 } // else if(aAOD)
 fCentralityHistVzCut->Fill(centrality);

 return kTRUE;

} // Bool_t AliAnalysisTaskNBodyFemtoscopy::CommonEventCuts(AliVEvent *ave)



//================================================================================================================
Bool_t AliAnalysisTaskNBodyFemtoscopy::CommonTrackCuts(AliAODTrack *atrack)
{
 // Apply track cuts.
 if(fApplyCommonTrackCuts) {

  TString sMethodName = "Bool_t AliAnalysisTaskNBodyFemtoscopy::CommonTrackCuts(AliAODTrack *atrack)";
  if(!atrack){Fatal(sMethodName.Data(),"!atrack");}

  if(!atrack->TestFilterBit(128)) return kFALSE; 

  Double_t pt = atrack->Pt(); // Pt
  if(atrack->Eta()<fEtaRange[0]) return kFALSE;
  if(atrack->Eta()>=fEtaRange[1]) return kFALSE;
  fPtHistEtaCut->Fill(pt); // filling pt distribution after eta cuts
  if(pt<fPtRange[0]) return kFALSE;
  if(pt>=fPtRange[1]) return kFALSE;
  fPtHistEtaCutPTCut->Fill(pt); // filling pt distribution after eta and pt cuts
  if(atrack->Phi()<fPhiRange[0]) return kFALSE;
  if(atrack->Phi()>=fPhiRange[1]) return kFALSE;
  fPtHistEtaCutPTCutPhiCut->Fill(pt); // filling pt distribution after eta, pt and phi cuts
  
 }
 return kTRUE;

} // Bool_t AliAnalysisTaskNBodyFemtoscopy::CommonTrackCuts(AliAODTrack *atrack)
//================================================================================================================

void AliAnalysisTaskNBodyFemtoscopy::Terminate(Option_t *)
{
 // Accessing the merged output list.

 fHistList = (TList*)GetOutputData(1);
 if(!fHistList){exit(1);}

 // Do some calculation in offline mode here:

 // ... your code for offline calculations ...

 // Update the output file with new results:
 TFile *f = new TFile("AnalysisResults.root","RECREATE");
 fHistList->Write(fHistList->GetName(),TObject::kSingleKey);

 delete f;

} // end of void AliAnalysisTaskNBodyFemtoscopy::Terminate(Option_t *)

//================================================================================================================

void AliAnalysisTaskNBodyFemtoscopy::InitializeArrays()
{
 // Initialize all data members which are arrays in this method.
  fPtRange[0]=0.2;
  fPtRange[1]=5.0;
  fEtaRange[0]=-0.8;
  fEtaRange[1]=0.8;
  fVertexZ[0]=-10.0;
  fVertexZ[1]=10.0;
  fPhiRange[0]=0.0;
  fPhiRange[1]=TMath::TwoPi();


} // void AliAnalysisTaskNBodyFemtoscopy::InitializeArrays()

//=======================================================================================================================

void AliAnalysisTaskNBodyFemtoscopy::BookAndNestAllLists()
{
 // Book and nest all lists nested in the base list fHistList.

 // a) Book and nest lists for control histograms;
 // b) Book and nest lists for final results.

 TString sMethodName = "void AliAnalysisTaskNBodyFemtoscopy::BookAndNestAllLists()";
 if(!fHistList){Fatal(sMethodName.Data(),"fHistList is NULL");}

 // a) Book and nest lists for control histograms:
 fControlHistogramsList = new TList();
 fControlHistogramsList->SetName("ControlHistograms");
 fControlHistogramsList->SetOwner(kTRUE);
 fHistList->Add(fControlHistogramsList);

 // b) Book and nest lists for final results:
 fFinalResultsList = new TList();
 fFinalResultsList->SetName("FinalResults");
 fFinalResultsList->SetOwner(kTRUE);
 fHistList->Add(fFinalResultsList);

} // void AliAnalysisTaskNBodyFemtoscopy::BookAndNestAllLists()

//=======================================================================================================================

void AliAnalysisTaskNBodyFemtoscopy::BookControlHistograms()
{
 // Book all control histograms.

 // a) Book histogram to hold pt spectra;
 // b) Book histogram for centrality distribution
 // c) Book histogram for track number distribution
 // d) Book histogram for centrality distribution after event cuts
 // e) Book histogram for pt spectra after track cuts

 // a) Book histogram to hold pt spectra before cuts:
 fPtHist = new TH1F("fPtHist","atrack->Pt()",fNbinsPt,fMinBinPt,fMaxBinPt);
 fPtHist->SetFillColor(kBlue-10);
 fPtHist->GetXaxis()->SetTitle("p_{t}");
 fControlHistogramsList->Add(fPtHist);
 
 // b) Book histogram for centrality distribution before cuts:
 fCentralityHist = new TH1F("fCentralityHist","ams->GetMultiplicityPercentile()",fNbinsCentrality,fMinCentrality,fMaxCentrality);
 fCentralityHist->SetFillColor(kRed-10);
 fCentralityHist->GetXaxis()->SetTitle("Centrality percentile");
 fControlHistogramsList->Add(fCentralityHist);

 // c) Book histogram for track number distribution
 fNumberOfTracksHist = new TH1F("fNumberOfTracksHist","aAOD->GetNumberOfTracks()",100,0,20000);
 fNumberOfTracksHist->SetFillColor(kRed-10);
 fNumberOfTracksHist->GetXaxis()->SetTitle("Number of tracks");
 fControlHistogramsList->Add(fNumberOfTracksHist);

 // c) Book histogram for centrality distribution after event cuts
 fCentralityHistVzCut = new TH1F("fCentralityHistVzCut","ams->GetMultiplicityPercentile() after vz cut",fNbinsCentrality,fMinCentrality,fMaxCentrality);
 fCentralityHistVzCut->SetFillColor(kRed-10);
 fCentralityHistVzCut->GetXaxis()->SetTitle("Centrality percentile");
 fControlHistogramsList->Add(fCentralityHistVzCut);

 // d) Book histogram for pt spectra after track cuts
 fPtHistEtaCut = new TH1F("fPtHistEtaCut","atrack->Pt() after Eta cuts",fNbinsPt,fMinBinPt,fMaxBinPt);
 fPtHistEtaCut->SetFillColor(kBlue-10);
 fPtHistEtaCut->GetXaxis()->SetTitle("p_{t}");
 fControlHistogramsList->Add(fPtHistEtaCut);

 fPtHistEtaCutPTCut = new TH1F("fPtHistEtaCutPTCut","atrack->Pt() after Eta and pT cuts",fNbinsPt,fMinBinPt,fMaxBinPt);
 fPtHistEtaCutPTCut->SetFillColor(kBlue-10);
 fPtHistEtaCutPTCut->GetXaxis()->SetTitle("p_{t}");
 fControlHistogramsList->Add(fPtHistEtaCutPTCut);

 fPtHistEtaCutPTCutPhiCut = new TH1F("fPtHistEtaCutPTCutPhiCut","atrack->Pt() after Eta, pT and Phi cuts",fNbinsPt,fMinBinPt,fMaxBinPt);
 fPtHistEtaCutPTCutPhiCut->SetFillColor(kBlue-10);
 fPtHistEtaCutPTCutPhiCut->GetXaxis()->SetTitle("p_{t}");
 fControlHistogramsList->Add(fPtHistEtaCutPTCutPhiCut);


 
} // void AliAnalysisTaskNBodyFemtoscopy::BookControlHistograms()

//=======================================================================================================================

void AliAnalysisTaskNBodyFemtoscopy::BookFinalResultsHistograms()
{
 // Book all histograms to hold the final results.

} // void AliAnalysisTaskNBodyFemtoscopy::BookFinalResultsHistograms()

//=======================================================================================================================

