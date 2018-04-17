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

/********************************************** 
* template class for student projects         *
* Marcel Lesch (marcel.lesch@cern.ch)         *
**********************************************/ 
  
#include "Riostream.h"
#include "AliAnalysisTaskStudentsML.h"
#include "AliLog.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"
#include "TRandom.h" 
#include "TComplex.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1F.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskStudentsML)

//================================================================================================================

AliAnalysisTaskStudentsML::AliAnalysisTaskStudentsML(const char *name, Bool_t useParticleWeights): 
 AliAnalysisTaskSE(name), 
 fHistList(NULL),
 // Control histograms:
 fControlHistogramsList(NULL),
 fPtHist(NULL),
 fNbins(1000),
 fMinBin(0.),
 fMaxBin(10.),
 fPhiHist(NULL),
 fEtaHist(NULL),
 fMultiHisto(NULL),
 fMaxCorrelator(8),
 // Final results:
 fFinalResultsList(NULL)
 {
  // Constructor.
 
  AliDebug(2,"AliAnalysisTaskStudentsML::AliAnalysisTaskStudentsML(const char *name, Bool_t useParticleWeights)");

  // Base list:
  fHistList = new TList();
  fHistList->SetName("outputStudentAnalysis");
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

} // AliAnalysisTaskStudentsML::AliAnalysisTaskStudentsML(const char *name, Bool_t useParticleWeights): 

//================================================================================================================

AliAnalysisTaskStudentsML::AliAnalysisTaskStudentsML(): 
 AliAnalysisTaskSE(),
 fHistList(NULL),
 // Control histograms:
 fControlHistogramsList(NULL),
 fPtHist(NULL),
 fNbins(1000),
 fMinBin(0.),
 fMaxBin(10.),
 fPhiHist(NULL),
 fEtaHist(NULL),
 fMultiHisto(NULL),
 fMaxCorrelator(8),
 // Final results:
 fFinalResultsList(NULL)
{
  // Dummy constructor.
 
  AliDebug(2,"AliAnalysisTaskStudentsML::AliAnalysisTaskStudentsML()");

} // AliAnalysisTaskStudentsML::AliAnalysisTaskStudentsML():

//================================================================================================================

AliAnalysisTaskStudentsML::~AliAnalysisTaskStudentsML()
{
 // Destructor.

 if(fHistList) delete fHistList;
  
} // AliAnalysisTaskStudentsML::~AliAnalysisTaskStudentsML()

//================================================================================================================

void AliAnalysisTaskStudentsML::UserCreateOutputObjects() 
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

} // void AliAnalysisTaskStudentsML::UserCreateOutputObjects() 

//================================================================================================================

void AliAnalysisTaskStudentsML::UserExec(Option_t *) 
{
 // Main loop (called for each event).

 // a) Get pointer to AOD event:
 // b) Start analysis over AODs;
 // c) Reset event-by-event objects;
 // d) PostData.

 // a) Get pointer to AOD event:
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(InputEvent()); // from TaskSE
 if(!aAOD){return;}





 // b) Start analysis over AODs:
 Int_t nTracks = aAOD->GetNumberOfTracks(); // number of all tracks in current event 
 fMultiHisto->Fill(nTracks);

 for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks
 {
  AliAODTrack *aTrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack)); // getting a pointer to a track
  if(!aTrack){continue;} // protection against NULL pointers
  if(!aTrack->TestFilterBit(128)){continue;} // filter bit 128 denotes TPC-only tracks, use only them for the analysis

  // example variables for each track:
  Double_t px = aTrack->Px(); // x-component of momenta
  Double_t py = aTrack->Py(); // y-component of momenta
  Double_t pz = aTrack->Pz(); // z-component of momenta
  Double_t e = aTrack->E();  // energy
  Double_t phi = aTrack->Phi(); // azimuthal angle
  Double_t eta = aTrack->Eta(); // pseudorapidity
  Double_t charge = aTrack->Charge(); // charge
  Double_t pt = aTrack->Pt(); // Pt (transverse momentum)


  // apply some cuts: e.g. take for the analysis only particles in -0.8 < eta < 0.8, and 0.2 < pT < 5.0
  // ... implementation of particle cuts ...
  if( (-0.8 < eta) && (eta < 0.8) && (0.2 < pt) && (pt < 5.0) ) 
  {
   // Fill some control histograms with the particles which passed cuts:
   fPtHist->Fill(pt);
   fPhiHist->Fill(phi); 
   fEtaHist->Fill(eta);
   

   // Do some analysis only with the particles which passed the cuts
   // ... your analysis code ... 

	//fill arrays for angles (weights and bin)

  } // if( (-0.8 < eta) && (eta < 0.8) && (0.2 < pt) && (pt < 5.0) )  
 } // for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks
   
 // c) Reset event-by-event objects:
 // ...


 // d) PostData:
 PostData(1,fHistList);

} // void AliAnalysisTaskStudentsML::UserExec(Option_t *)

//================================================================================================================

void AliAnalysisTaskStudentsML::Terminate(Option_t *)
{
 // Accessing the merged output list. 

 fHistList = (TList*)GetOutputData(1);
 if(!fHistList){exit(1);}

 // Do some calculation in offline mode here:
 // ...

 TFile *f = new TFile("AnalysisResults.root","RECREATE");
 fHistList->Write(fHistList->GetName(),TObject::kSingleKey);

 delete f;

} // end of void AliAnalysisTaskStudentsML::Terminate(Option_t *)

//================================================================================================================

void AliAnalysisTaskStudentsML::InitializeArrays()
{
 // Initialize all data members which are arrays in this method.

 // ...

} // void AliAnalysisTaskStudentsML::InitializeArrays()

//=======================================================================================================================

void AliAnalysisTaskStudentsML::BookAndNestAllLists()
{
 // Book and nest all lists nested in the base list fHistList.

 // a) Book and nest lists for control histograms;
 // b) Book and nest lists for final results.

 TString sMethodName = "void AliAnalysisTaskStudentsML::BookAndNestAllLists()";
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

} // void AliAnalysisTaskStudentsML::BookAndNestAllLists()

//=======================================================================================================================

void AliAnalysisTaskStudentsML::BookControlHistograms()
{
 // Book all control histograms.

 // a) Book histogram to hold pt distribution;
 // b) Book histogram to hold phi distribution;
 // c) Book histogram to hold eta distribution;
 // d) Book histogam to hold multiplicty distribution;
 // e) ...

 
 

 // a) Book histogram to hold pt spectra:
 fPtHist = new TH1F("fPtHist","atrack->Pt()",fNbins,fMinBin,fMaxBin);
 fPtHist->SetStats(kFALSE);
 fPtHist->SetFillColor(kBlue-10);
 fPtHist->GetXaxis()->SetTitle("p_{t}");
 fControlHistogramsList->Add(fPtHist);

 // b) Book histogram to hold phi distribution:
 fPhiHist = new TH1F("fPhiHist","Phi Distribution",1000,0.,6.3);
 fPhiHist->GetXaxis()->SetTitle("Phi");
 fPhiHist->SetLineColor(4);
 fControlHistogramsList->Add(fPhiHist);

 // c) Book histogram to hold eta distribution:
 fEtaHist = new TH1F("fEtaHist","Eta Distribution",1000,-1.,1.);
 fEtaHist->GetXaxis()->SetTitle("Eta");
 fEtaHist->SetLineColor(4);
 fControlHistogramsList->Add(fEtaHist);

 // d) Book histogam to hold multiplicty distribution:
 fMultiHisto = new TH1F("fMultiHisto","Multiplicity",300,0.,3.); //histogram for multiplicity
 fMultiHisto->GetXaxis()->SetTitle("Multiplicity M");
 fControlHistogramsList->Add(fMultiHisto);

 // e) ...

} // void AliAnalysisTaskStudentsML::BookControlHistograms()

//=======================================================================================================================

void AliAnalysisTaskStudentsML::BookFinalResultsHistograms()
{
 // Book all histograms to hold the final results.

 // fFinalResultsList->Add(...);

} // void AliAnalysisTaskStudentsML::BookFinalResultsHistograms()

//=======================================================================================================================

