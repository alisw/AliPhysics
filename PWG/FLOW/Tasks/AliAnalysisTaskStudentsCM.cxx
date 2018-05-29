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
* template class for student projects * 
**************************************/ 
  
#include "Riostream.h"
#include "AliAnalysisTaskStudentsCM.h"
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

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskStudentsCM)

//================================================================================================================

AliAnalysisTaskStudentsCM::AliAnalysisTaskStudentsCM(const char *name, Bool_t useParticleWeights):		// {Definition of the constructor of the class}
 AliAnalysisTaskSE(name), 
 fHistList(NULL),		// {These are pointers. There is no address to point to when they are initialise by default, so NULL}
 // Control histograms:
 fControlHistogramsList(NULL),
 fTestHistogramsList(NULL),	// {TEST}
 fPtHistNoCut(NULL),
 fPtHist(NULL),
 fNbins(1000),
 fMinBin(0.),
 fMaxBin(20.),
 fPhiHist(NULL),
 fEtaHist(NULL),
 fEnergyHist(NULL),	// TEST
 fEnergyHistNoCut(NULL),
 fPtPhiHist(NULL),	// TEST
 fMassSquareHist(NULL), // TEST
 // Final results:
 fFinalResultsList(NULL),
 fEtaMassSquareHist(NULL)
 {
  // Constructor. {The data members must initialise in the EXACT SAME order than they are declared in the header file ! If a data member is added, it must be put at the same place in the list in the .h and the .cxx}
 
  AliDebug(2,"AliAnalysisTaskStudentsCM::AliAnalysisTaskStudentsCM(const char *name, Bool_t useParticleWeights)");

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

} // AliAnalysisTaskStudentsCM::AliAnalysisTaskStudentsCM(const char *name, Bool_t useParticleWeights): 

//================================================================================================================

AliAnalysisTaskStudentsCM::AliAnalysisTaskStudentsCM():		// {Dummy constructor. Simply paste the initialised data members again here.}
 AliAnalysisTaskSE(),
 fHistList(NULL),
 // Control histograms:
 fControlHistogramsList(NULL),
 fTestHistogramsList(NULL),		// {TEST}
 fPtHistNoCut(NULL),
 fPtHist(NULL),
 fNbins(1000),
 fMinBin(0.),
 fMaxBin(20.),
 fPhiHist(NULL),
 fEtaHist(NULL),
 fEnergyHist(NULL),	// {TEST}
 fEnergyHistNoCut(NULL),
 fPtPhiHist(NULL),	// TEST
 fMassSquareHist(NULL),	// TEST
 // Final results:
 fFinalResultsList(NULL),
 fEtaMassSquareHist(NULL)
{
  // Dummy constructor.
 
  AliDebug(2,"AliAnalysisTaskStudentsCM::AliAnalysisTaskStudentsCM()");

} // AliAnalysisTaskStudentsCM::AliAnalysisTaskStudentsCM():

//================================================================================================================

AliAnalysisTaskStudentsCM::~AliAnalysisTaskStudentsCM()
{
 // Destructor. {Delete the main TList deletes everything in it automatically}

 if(fHistList) delete fHistList;
  
} // AliAnalysisTaskStudentsCM::~AliAnalysisTaskStudentsCM()

//================================================================================================================

void AliAnalysisTaskStudentsCM::UserCreateOutputObjects() 
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
 this->BookTestHistograms();		// {TEST}
 this->BookFinalResultsHistograms();

 // *) Trick to avoid name clashes, part 2:
 TH1::AddDirectory(oldHistAddStatus);

 PostData(1,fHistList);

} // void AliAnalysisTaskStudentsCM::UserCreateOutputObjects() 

//================================================================================================================

void AliAnalysisTaskStudentsCM::UserExec(Option_t *) 
{
 // Main loop (called for each event). {Put the tasks that are applied to each event there.}

 // a) Get pointer to AOD event:
 // b) Start analysis over AODs;
 // c) Reset event-by-event objects;
 // d) PostData.

 // a) Get pointer to AOD event:
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(InputEvent()); // from TaskSE
 if(!aAOD){return;}

 // b) Start analysis over AODs:
 Int_t nTracks = aAOD->GetNumberOfTracks(); // number of all tracks in current event 
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
  Double_t mSquare = 0.;	// mass^2

  // Fill the control histogram containing the transverse momentum before the application of any cut
  fPtHistNoCut->Fill(pt);
  fEnergyHistNoCut->Fill(e);
 
  // apply some cuts: e.g. take for the analysis only particles in -0.8 < eta < 0.8, and 0.2 < pT < 5.0
  // ... implementation of particle cuts ...
  if( (-0.8 < eta) && (eta < 0.8) && (0.2 < pt) && (pt < 5.0) )
  {
   // Fill some control histograms with the particles which passed cuts:
   fPtHist->Fill(pt);
   fPhiHist->Fill(phi); 
   //fEtaHist->Fill(eta);	// moved to final results as a test

   // Fill the personal histograms into the test subfolder of the output
   fEnergyHist->Fill(e);		// {TEST}
   fPtPhiHist->Fill(pt,phi);	// TEST


   // Do some analysis only with the particles which passed the cuts
   // ... your analysis code ...
   // TEST
   mSquare = e*e - (px*px + py*py + pz*pz);

   // Fill some final histograms
   fEtaHist->Fill(eta);
   fMassSquareHist->Fill(mSquare);
   fEtaMassSquareHist->Fill(eta,mSquare); 

  } // if( (-0.8 < eta) && (eta < 0.8) && (0.2 < pt) && (pt < 5.0) )  
 } // for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks
   
 // c) Reset event-by-event objects:
 // ...

 // d) PostData:
 PostData(1,fHistList);

} // void AliAnalysisTaskStudentsCM::UserExec(Option_t *)

//================================================================================================================

void AliAnalysisTaskStudentsCM::Terminate(Option_t *)
{
 // Accessing the merged output list. 

 fHistList = (TList*)GetOutputData(1);
 if(!fHistList){exit(1);}

 // Do some calculation in offline mode here: {i.e. once the run over the events has been done}
 // ...


 TFile *f = new TFile("AnalysisResults.root","RECREATE");		// {Name of the output file}
 fHistList->Write(fHistList->GetName(),TObject::kSingleKey);	// {Writing of the main TList in the output file}

 delete f;

} // end of void AliAnalysisTaskStudentsCM::Terminate(Option_t *)

//================================================================================================================

void AliAnalysisTaskStudentsCM::InitializeArrays()
{
 // Initialize all data members which are arrays in this method.

 // ...

} // void AliAnalysisTaskStudentsCM::InitializeArrays()

//=======================================================================================================================

void AliAnalysisTaskStudentsCM::BookAndNestAllLists()
{
 // Book and nest all lists nested in the base list fHistList.

 // a) Book and nest lists for control histograms;
 // b) Book and nest lists for final results.

 TString sMethodName = "void AliAnalysisTaskStudentsCM::BookAndNestAllLists()";
 if(!fHistList){Fatal(sMethodName.Data(),"fHistList is NULL");}

 // a) Book and nest lists for control histograms:
 fControlHistogramsList = new TList();
 fControlHistogramsList->SetName("ControlHistograms");
 fControlHistogramsList->SetOwner(kTRUE);
 fHistList->Add(fControlHistogramsList);

 // {TEST} Book and nest lists for personal test histograms:
 fTestHistogramsList = new TList();
 fTestHistogramsList->SetName("TestHistograms");
 fTestHistogramsList->SetOwner(kTRUE);
 fHistList->Add(fTestHistogramsList);

 // b) Book and nest lists for final results:
 fFinalResultsList = new TList();
 fFinalResultsList->SetName("FinalResults");
 fFinalResultsList->SetOwner(kTRUE);
 fHistList->Add(fFinalResultsList);

} // void AliAnalysisTaskStudentsCM::BookAndNestAllLists()

//=======================================================================================================================

void AliAnalysisTaskStudentsCM::BookControlHistograms()
{
 // Book all control histograms.

 // a) Book histogram to hold pt distribution before the application of cuts;
 // b) Book histogram to hold pt distribution;
 // c) Book histogram to hold phi distribution;
 // d) Book histogram to hold eta distribution
 // e) ...


 // a) Book histogram to hold pt distribution before the application of cuts;
 fPtHistNoCut = new TH1F("fPtHistNoCut", "atrack->Pt()", fNbins, fMinBin, fMaxBin);
 fPtHistNoCut->SetStats(kTRUE);
 fPtHistNoCut->SetFillColor(kRed);
 fPtHistNoCut->GetXaxis()->SetTitle("p_{t}");
 fControlHistogramsList->Add(fPtHistNoCut);

 // b) Book histogram to hold pt spectra:
 fPtHist = new TH1F("fPtHist","atrack->Pt()",fNbins,fMinBin,fMaxBin);
 fPtHist->SetStats(kTRUE);
 fPtHist->SetFillColor(kBlue-10);
 fPtHist->GetXaxis()->SetTitle("p_{t}");
 fControlHistogramsList->Add(fPtHist);

 // c) Book histogram to hold phi distribution:
 fPhiHist = new TH1F("fPhiHist","Phi Distribution",1000,0.,6.3);
 fPhiHist->GetXaxis()->SetTitle("Phi");
 fPhiHist->SetLineColor(4);
 fControlHistogramsList->Add(fPhiHist);

 // c) Book histogram to hold eta distribution:
 /*fEtaHist = new TH1F("fEtaHist","Eta Distribution",1000,-1.,1.);
 fEtaHist->GetXaxis()->SetTitle("Eta");
 fEtaHist->SetLineColor(4);
 fControlHistogramsList->Add(fEtaHist);*/
 
} // void AliAnalysisTaskStudentsCM::BookControlHistograms()

//=======================================================================================================================

void AliAnalysisTaskStudentsCM::BookTestHistograms()
{
 // Book all personal histograms.

 // a) Book histogram to hold energy distribution
 // b) Book histogram to hold energy distribution before the application of cuts
 // c) Book 2D histogram to hold (pt, phi) distribution

 // a) Book histogram to hold energy distribution {TEST}
 fEnergyHist = new TH1F("fEnergyHist", "Energy Distribution",fNbins,fMinBin,fMaxBin);
 fEnergyHist->SetStats(kTRUE);
 fEnergyHist->SetLineColor(3);
 fEnergyHist->GetXaxis()->SetTitle("E");
 fTestHistogramsList->Add(fEnergyHist);

 // b) Book histogram to hold energy distribution before the application of cuts
 fEnergyHistNoCut = new TH1F("fEnergyHistNoCut", "Energy distribution before cuts", fNbins, fMinBin, fMaxBin);
 fEnergyHistNoCut->SetStats(kTRUE);
 fEnergyHistNoCut->GetXaxis()->SetTitle("E");
 fTestHistogramsList->Add(fEnergyHistNoCut);

 // c) Book 2D histogram to hold (pt, phi) distribution TEST
 fPtPhiHist = new TH2F("fPtPhiHist", "(Pt, Phi) Distribution", fNbins,fMinBin,fMaxBin,1000,0.,6.3);
 fPtPhiHist->SetStats(kFALSE);
 fPtPhiHist->GetXaxis()->SetTitle("p_{t}");
 fPtPhiHist->GetYaxis()->SetTitle("Phi");
 fTestHistogramsList->Add(fPtPhiHist);

} // void AliAnalysisTaskStudentsCM::BookTestHistograms()

// ======================================================================================================================

void AliAnalysisTaskStudentsCM::BookFinalResultsHistograms()
{
 // Book all histograms to hold the final results.

 // a) Book histogram to hold eta (previously into control histogram) TEST
 // b) Book histogram to hold mSquare TEST
 // c) Book histogram to hold mSquare as a function of eta TEST

 // a) Book histogram to hold eta (previously into control histogram) TEST
 fEtaHist = new TH1F("fEtaHist","Eta Distribution",1000,-1.,1.);
 fEtaHist->GetXaxis()->SetTitle("Eta");
 fEtaHist->SetLineColor(4);
 fFinalResultsList->Add(fEtaHist);

 // b) Book histogram to hold mSquare TEST
 fMassSquareHist = new TH1F("fMassSquareHist", "Mass squared distribution",1000,0.,15.);
 fMassSquareHist->GetXaxis()->SetTitle("m^{2}");
 fMassSquareHist->SetLineColor(1);
 fFinalResultsList->Add(fMassSquareHist);

 // c) Book histogram to hold mSquare as a function of eta TEST
 fEtaMassSquareHist = new TH2F("fEtaMassSquareHist", "(eta, m^2) distribution",1000,-1.,1.,1000,0.,5.);
 fEtaMassSquareHist->GetXaxis()->SetTitle("Eta");
 fEtaMassSquareHist->GetYaxis()->SetTitle("m^2");
 fFinalResultsList->Add(fEtaMassSquareHist);

 // fFinalResultsList->Add(...);

} // void AliAnalysisTaskStudentsCM::BookFinalResultsHistograms()

//=======================================================================================================================

