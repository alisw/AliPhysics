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
#include "AliAnalysisTaskStudentsMW.h"
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
#include <TArrayF.h>
#include <vector>
#include "TMath.h"
#include "TF1.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskStudentsMW)

//================================================================================================================

AliAnalysisTaskStudentsMW::AliAnalysisTaskStudentsMW(const char *name, Bool_t useParticleWeights): 
 AliAnalysisTaskSE(name), 
 fHistList(NULL),
 // Control histograms:
 fControlHistogramsList(NULL),
 fPtHist(NULL),
 fNbins(1000),
 fMinBin(0.),
 fMaxBin(10.),
 fCentralityHist(NULL),
 fMinCentrality(0.),
 fMaxCentrality(100.),
 fPhiHistBeforeTrackSeletion(NULL),
 fEtaHistBeforeTrackSeletion(NULL),
 fTotalMultBeforeTrackSeletion(NULL),
 fMultiHistoBeforeTrackSeletion(NULL),
 fPhiHistAfterTrackSeletion(NULL),
 fEtaHistAfterTrackSeletion(NULL),
 fTotalMultAfterTrackSeletion(NULL),
 fMultiHistoAfterTrackSeletion(NULL),
 fMultiHistoBeforeMultCut(NULL),
 //SelectionCuts
 fMainFilter(0),

 //Global
 bCutOnVertexZ(kFALSE), 
 fMinVertexZ(-10.),
 fMaxVertexZ(10.),
 fVertexZBefore(NULL),
 fVertexZAfter(NULL),
   //Physics
 bCutOnEta(kTRUE),
 bCutOnPt(kTRUE),
 fMinEtaCut(-0.8),
 fMaxEtaCut(0.8),
 fMinPtCut(0.2),
 fMaxPtCut(5.0),
 //Variables for the correlation
 fMaxCorrelator(8),
 bUseWeights(kFALSE),
 fNumber(6),  //number of correlation first correlator
 fNumberSecond(6), //number of correlation second correlator
 fMinNumberPart(8),
 fh1(0), fh2(0), fh3(0), fh4(0), fh5(0), fh6(0), fh7(0), fh8(0),  //harmonics
 fa1(0), fa2(0), fa3(0), fa4(0), fa5(0), fa6(0), fa7(0), fa8(0),  //second set of harmonics
 fParticles(0),
 fAngles(NULL),
 fWeights(NULL),
 fBin(NULL),
 fCentralityres(NULL),
 fCentralitySecond(NULL),
 fCentralitySecondSquare(NULL),
 fCov(NULL),
 fCounterHistogram(NULL),
 // Final results:
 fFinalResultsList(NULL)
 {
  // Constructor.
 
  AliDebug(2,"AliAnalysisTaskStudentsMW::AliAnalysisTaskStudentsMW(const char *name, Bool_t useParticleWeights)");

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

  //num =1;

  DefineOutput(1, TList::Class());  

  if(useParticleWeights)
  {
   // not needed for the time being
  }

} // AliAnalysisTaskStudentsMW::AliAnalysisTaskStudentsMW(const char *name, Bool_t useParticleWeights): 

//================================================================================================================

AliAnalysisTaskStudentsMW::AliAnalysisTaskStudentsMW(): 
 AliAnalysisTaskSE(),
 fHistList(NULL),
  // Control histograms:
 fControlHistogramsList(NULL),
 fPtHist(NULL),
 fNbins(1000),
 fMinBin(0.),
 fMaxBin(10.),
 fCentralityHist(NULL),
 fMinCentrality(0.),
 fMaxCentrality(100.),
 fPhiHistBeforeTrackSeletion(NULL),
 fEtaHistBeforeTrackSeletion(NULL),
 fTotalMultBeforeTrackSeletion(NULL),
 fMultiHistoBeforeTrackSeletion(NULL),
 fPhiHistAfterTrackSeletion(NULL),
 fEtaHistAfterTrackSeletion(NULL),
 fTotalMultAfterTrackSeletion(NULL),
 fMultiHistoAfterTrackSeletion(NULL),
 fMultiHistoBeforeMultCut(NULL),
 //SelectionCuts
 fMainFilter(0),

 //Global
 bCutOnVertexZ(kFALSE), 
 fMinVertexZ(-10.),
 fMaxVertexZ(10.),
 fVertexZBefore(NULL),
 fVertexZAfter(NULL),
   //Physics
 bCutOnEta(kTRUE),
 bCutOnPt(kTRUE),
 fMinEtaCut(-0.8),
 fMaxEtaCut(0.8),
 fMinPtCut(0.2),
 fMaxPtCut(5.0),
 //Variables for the correlation
 fMaxCorrelator(8),
 bUseWeights(kFALSE),
 fNumber(6),  //number of correlation first correlator
 fNumberSecond(6), //number of correlation second correlator
 fMinNumberPart(8),
 fh1(0), fh2(0), fh3(0), fh4(0), fh5(0), fh6(0), fh7(0), fh8(0),  //harmonics
 fa1(0), fa2(0), fa3(0), fa4(0), fa5(0), fa6(0), fa7(0), fa8(0),  //second set of harmonics
 fParticles(0),
 fAngles(NULL),
 fWeights(NULL),
 fBin(NULL),
 fCentralityres(NULL),
 fCentralitySecond(NULL),
 fCentralitySecondSquare(NULL),
 fCov(NULL),
 fCounterHistogram(NULL),
fFinalResultsList(NULL)
{
  // Dummy constructor.
 
  AliDebug(2,"AliAnalysisTaskStudentsMW::AliAnalysisTaskStudentsMW()");

} // AliAnalysisTaskStudentsMW::AliAnalysisTaskStudentsMW():

//================================================================================================================

AliAnalysisTaskStudentsMW::~AliAnalysisTaskStudentsMW()
{
 // Destructor.

 if(fHistList) delete fHistList;
  
} // AliAnalysisTaskStudentsMW::~AliAnalysisTaskStudentsMW()

//================================================================================================================

void AliAnalysisTaskStudentsMW::UserCreateOutputObjects() 
{
 // Called at every worker node to initialize.

 // a) Trick to avoid name clashes, part 1;
 // b) Book and nest all lists;
 // c) Book all objects;
 // *) Trick to avoid name clashes, part 2.
 // d) 

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

} // void AliAnalysisTaskStudentsMW::UserCreateOutputObjects() 

//================================================================================================================

void AliAnalysisTaskStudentsMW::UserExec(Option_t *) 
{
 // Main loop (called for each event).

 // a) Get pointer to AOD event:
 // b) Start analysis over AODs;
 // c) Reset event-by-event objects;
 // d) PostData.


 // a) Get pointer to AOD event:
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(InputEvent()); // from TaskSE
 if(!aAOD){return;}

 AliMultSelection *ams = (AliMultSelection*)aAOD->FindListObject("MultSelection");
 if(!ams){return;}
 if(ams->GetMultiplicityPercentile("V0M") >= fMinCentrality && ams->GetMultiplicityPercentile("V0M") < fMaxCentrality)
 {
  fCentralityHist->Fill(ams->GetMultiplicityPercentile("V0M"));
 }
 else
 {
  return; // this event do not belong to the centrality class specified for this particular analysis
 }

 AliAODVertex *avtx = (AliAODVertex*)aAOD->GetPrimaryVertex();

 fVertexZBefore->Fill(avtx->GetZ());

 if(avtx->GetZ() < fMinVertexZ) return;
 if(avtx->GetZ() > fMaxVertexZ) return;
 
 fVertexZAfter->Fill(avtx->GetZ());

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

 //b.0) Start analysis over AODs:
 Int_t nTracks = aAOD->GetNumberOfTracks(); // number of all tracks in current event 
 fAngles = new TArrayD(nTracks); //new Array
 fParticles=0; //number of particles after selections

 if(nTracks>0){fMultiHistoBeforeTrackSeletion->Fill(nTracks);} //multiplicity distribution before track selection
 for(Int_t u=0;u<nTracks;u++){fTotalMultBeforeTrackSeletion->Fill(0.5);} //total number of particles in whole centrality class before track sel.


  //b.1) Loop over the tracks in the event with PhysicsSelection(Eta Cut, Pt Cut)
 for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks
 {
    AliAODTrack *aTrack = dynamic_cast<AliAODTrack*>(aAOD->GetTrack(iTrack)); // getting a pointer to a track
    if(!aTrack){continue;} // protection against NULL pointers
    if(!aTrack->TestFilterBit(fMainFilter)){continue;} //Check if in Filter
    
     Double_t phi = aTrack->Phi(); // azimuthal angle
     Double_t eta = aTrack->Eta(); // pseudorapidity
     Double_t pt = aTrack->Pt(); // Pt (transverse momentum)

      // Fill some control histograms with the particles before track selection:
     fPhiHistBeforeTrackSeletion->Fill(phi); 
     fEtaHistBeforeTrackSeletion->Fill(eta);

      if(!TrackSelection(aTrack)){continue;} //Track did not pass physics selection 
	
      // Fill some control histograms with the particles after track selection:
     fPtHist->Fill(pt);
     fPhiHistAfterTrackSeletion->Fill(phi); 
     fEtaHistAfterTrackSeletion->Fill(eta);

     //Fill angle array  
     fAngles->AddAt(phi,fParticles);
     fParticles += 1;


 } // for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks
 
 fAngles->Set(fParticles); 


 if(fParticles>0){fMultiHistoAfterTrackSeletion->Fill(fParticles);} //multiplicity distribution after track selection
 for(Int_t u=0;u<fParticles;u++){fTotalMultAfterTrackSeletion->Fill(0.5);} //total number of particles in whole centrality class after track sel.

 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

 //b.2) Multi-Particle Correlation;

 if(fParticles>=fMinNumberPart) //do the correlation only if there are more than 8 particles in the event
 { 

    // Calculate Q-vectors for available angles and weights;
    this->CalculateQvectors();

    Double_t FirstCorrelation=0.;
    Double_t Weight_FirstCorrelation=0.;
    Double_t SecondCorrelation=0.;
    Double_t Weight_SecondCorrelation=0.;
    
    //~~~~~~~~~~~~~~~~~

    this->Correlation(fNumber,fh1, fh2, fh3, fh4, fh5, fh6, fh7, fh8);  //do the correlation for the first set

    FirstCorrelation=fRecursion[0][fNumber-2]->GetBinContent(1);
    Weight_FirstCorrelation=fRecursion[0][fNumber-2]->GetBinContent(2);

    fRecursion[0][fNumber-2]->Reset(); //Reset
    fRecursion[1][fNumber-2]->Reset(); //Reset

    //~~~~~~~~~~~~~~~~~

    this->Correlation(fNumberSecond,fa1, fa2, fa3, fa4, fa5, fa6, fa7, fa8);  //do the correlation for the second set

    SecondCorrelation=fRecursion[0][fNumberSecond-2]->GetBinContent(1);
    Weight_SecondCorrelation=fRecursion[0][fNumberSecond-2]->GetBinContent(2);
    
    fRecursion[0][fNumberSecond-2]->Reset(); //Reset
    fRecursion[1][fNumberSecond-2]->Reset(); //Reset

    //~~~~~~~~~~~~~~~~~
    
    fCentralityres->Fill(0.5,FirstCorrelation,Weight_FirstCorrelation); //safe output first set of harmonics
    fCentralitySecond->Fill(0.5,SecondCorrelation,Weight_SecondCorrelation); //safe output second set of harmonics

   fCentralitySecondSquare->Fill(0.5,SecondCorrelation*SecondCorrelation,Weight_SecondCorrelation*Weight_SecondCorrelation);
   fCov->Fill(0.5,FirstCorrelation*SecondCorrelation,Weight_FirstCorrelation*Weight_SecondCorrelation);

  } //if(fParticles>=fMinNumberPart)

 
 // c) Reset event-by-event objects:
 fParticles=0;
 delete fAngles;

 

 // d) PostData:
 PostData(1,fHistList);

} // void AliAnalysisTaskStudentsMW::UserExec(Option_t *)

//================================================================================================================

void AliAnalysisTaskStudentsMW::Terminate(Option_t *)
{

 // Accessing the merged output list. 

 fHistList = (TList*)GetOutputData(1);
 if(!fHistList){exit(1);}

 // Do some calculation in offline mode here:
 // ...

 TFile *f = new TFile("AnalysisResults.root","RECREATE");
 fHistList->Write(fHistList->GetName(),TObject::kSingleKey);

 delete f;

} // end of void AliAnalysisTaskStudentsMW::Terminate(Option_t *)

//================================================================================================================

void AliAnalysisTaskStudentsMW::InitializeArrays()
{
 // Initialize all data members which are arrays in this method.

   for(Int_t cs=0;cs<2;cs++) 
   {
     for(Int_t c=0;c<fMaxCorrelator;c++)
     {
   
      fRecursion[cs][c] = NULL; //! [cs]: real (0) or imaginary part (1) ....
   
     }  
    }  //for(Int_t cs=0;cs<2;cs++)

   for(Int_t js=0;js<49;js++) 
   {
     for(Int_t j=0;j<9;j++)
     {
   
      fQvector[js][j] = TComplex(0.,0.); //! 
   
     } 
   } 

} // void AliAnalysisTaskStudentsMW::InitializeArrays()

//=======================================================================================================================

void AliAnalysisTaskStudentsMW::BookAndNestAllLists()
{
 // Book and nest all lists nested in the base list fHistList.

 // a) Book and nest lists for control histograms;
 // b) Book and nest lists for final results.

 TString sMethodName = "void AliAnalysisTaskStudentsMW::BookAndNestAllLists()";
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

} // void AliAnalysisTaskStudentsMW::BookAndNestAllLists()

//=======================================================================================================================

void AliAnalysisTaskStudentsMW::BookControlHistograms()
{
 // Book all control histograms.

 // a) Book histogram to hold pt spectra;
 // b) ...

 // a) Book histogram to hold pt spectra and centrality:
 fPtHist = new TH1F("fPtHist","atrack->Pt()",fNbins,fMinBin,fMaxBin);
 fPtHist->SetStats(kFALSE);
 fPtHist->SetFillColor(kBlue-10);
 fPtHist->GetXaxis()->SetTitle("p_{t}");
 fControlHistogramsList->Add(fPtHist);

 fCentralityHist = new TH1F("fCentralityHist","ams->GetMultiplicityPercentile()",100,0.,100.);
 fControlHistogramsList->Add(fCentralityHist);

 // b) Book histogram to hold phi distribution before track selection:
 fPhiHistBeforeTrackSeletion = new TH1F("fPhiHistBeforeTrackSeletion","Phi Distribution",1000,0.,6.3);
 fPhiHistBeforeTrackSeletion->GetXaxis()->SetTitle("Phi");
 fPhiHistBeforeTrackSeletion->SetLineColor(4);
 fControlHistogramsList->Add(fPhiHistBeforeTrackSeletion);

 // c) Book histogram to hold eta distribution before track selection:
 fEtaHistBeforeTrackSeletion = new TH1F("fEtaHistBeforeTrackSeletion","Eta Distribution",1000,-1.,1.);
 fEtaHistBeforeTrackSeletion->GetXaxis()->SetTitle("Eta");
 fEtaHistBeforeTrackSeletion->SetLineColor(4);
 fControlHistogramsList->Add(fEtaHistBeforeTrackSeletion);

 // d) Book Mult. Histo before before track selection
 fTotalMultBeforeTrackSeletion = new TH1F("fTotalMultBeforeTrackSeletion","Mult. Counts per Class before brute cut",1,0.,1.);
 fTotalMultBeforeTrackSeletion->GetYaxis()->SetTitle("Counts");
 fTotalMultBeforeTrackSeletion->SetLineColor(4);
 fControlHistogramsList->Add(fTotalMultBeforeTrackSeletion);
 
 // e) Book histogam to hold multiplicty distribution before track selection:
 fMultiHistoBeforeTrackSeletion = new TH1F("fMultiHistoBeforeTrackSeletion","Multiplicity",5000,0.,5000.); 
 fMultiHistoBeforeTrackSeletion->GetXaxis()->SetTitle("Multiplicity M");
 fControlHistogramsList->Add(fMultiHistoBeforeTrackSeletion);

 // f) Book histogram to hold phi distribution after track selection:
 fPhiHistAfterTrackSeletion = new TH1F("fPhiHistAfterTrackSeletion","Phi Distribution",1000,0.,6.3);
 fPhiHistAfterTrackSeletion->GetXaxis()->SetTitle("Phi");
 fPhiHistAfterTrackSeletion->SetLineColor(4);
 fControlHistogramsList->Add(fPhiHistAfterTrackSeletion);

 // g) Book histogram to hold eta distribution after track selection:
  fEtaHistAfterTrackSeletion = new TH1F("fEtaHistAfterTrackSeletion","Eta Distribution",1000,-1.,1.);
 fEtaHistAfterTrackSeletion->GetXaxis()->SetTitle("Eta");
 fEtaHistAfterTrackSeletion->SetLineColor(4);
 fControlHistogramsList->Add(fEtaHistAfterTrackSeletion);

 // h) Book Mult. Histo before after track selection
  fTotalMultAfterTrackSeletion = new TH1F("fTotalMultAfterTrackSeletion","Mult. Counts per Class before brute cut",1,0.,1.);
 fTotalMultAfterTrackSeletion->GetYaxis()->SetTitle("Counts");
 fTotalMultAfterTrackSeletion->SetLineColor(4);
 fControlHistogramsList->Add(fTotalMultAfterTrackSeletion);
 
 // i) Book histogam to hold multiplicty distribution after track selection:
 fMultiHistoAfterTrackSeletion = new TH1F("fMultiHistoAfterTrackSeletion","Multiplicity",5000,0.,5000.);
 fMultiHistoAfterTrackSeletion->GetXaxis()->SetTitle("Multiplicity M");
 fControlHistogramsList->Add(fMultiHistoAfterTrackSeletion);

  // o) Book histogam for Vertex Z after Cut
 fVertexZBefore = new TH1F("fVertexZBefore","fVertexZBefore",1000,-20.,20.); 
 fVertexZBefore->GetXaxis()->SetTitle("");
 fControlHistogramsList->Add(fVertexZBefore);

  // p) Book histogam for Vertex Z after Cut
 fVertexZAfter = new TH1F("fVertexZAfter","fVertexZAfter",1000,-20.,20.); 
 fVertexZAfter->GetXaxis()->SetTitle("");
 fControlHistogramsList->Add(fVertexZAfter);

 // q) Book histogram to debug
 fCounterHistogram = new TH1F("fCounterHistogram","Histogram for some checks",3,0.,3.); 
 fControlHistogramsList->Add(fCounterHistogram);

} // void AliAnalysisTaskStudentsMW::BookControlHistograms()

//=======================================================================================================================

void AliAnalysisTaskStudentsMW::BookFinalResultsHistograms()
{
 // Book all histograms to hold the final results.
  
 fCentralityres = new TProfile("fCentrality","Result Analysis First Set Correlators",1,0.,1.); //centrality dependet output
 fCentralityres->GetXaxis()->SetTitle("");
 fCentralityres->GetYaxis()->SetTitle("flow");
 fCentralityres->Sumw2();
 fFinalResultsList->Add(fCentralityres);

 fCentralitySecond = new TProfile("fCentralitySecond","Result Analysis Second Set Correlators",1,0.,1.); //centrality dependet output
 fCentralitySecond->GetXaxis()->SetTitle("");
 fCentralitySecond->GetYaxis()->SetTitle("flow");
 fCentralitySecond->Sumw2(); 
 fFinalResultsList->Add(fCentralitySecond);

 fCentralitySecondSquare = new TProfile("fCentralitySecondSquare","Result Analysis Second Set Correlators Squared",1,0.,1.); //centrality dependet output
 fCentralitySecondSquare->GetXaxis()->SetTitle("");
 fCentralitySecondSquare->Sumw2();  
 fFinalResultsList->Add(fCentralitySecondSquare);

 fCov = new TProfile("fCov","Result Analysis Covariance Term",1,0.,1.); //centrality dependet output
 fCov->GetXaxis()->SetTitle("");
 fCov->Sumw2();  
 fFinalResultsList->Add(fCov);

 Cosmetics();
 
} // void AliAnalysisTaskStudentsMW::BookFinalResultsHistograms()


//=======================================================================================================================

void AliAnalysisTaskStudentsMW::Cosmetics()
{
 // Book everything here.
  
 for(Int_t cs=0;cs<2;cs++) 
 {
  for(Int_t c=0;c<fMaxCorrelator;c++)
  {
   
   fRecursion[cs][c] = new TProfile("","",2,0.,2.); 
   //fRecursion[cs][c]->Sumw2();
 
   //NOTE for fRecursion: 1.) [cs] will say if its the real (0) or imaginary part (1) 
   // 2.) [c] gives gives the kind of correlation. [n] is the (n+2)-particle correlation
   //3.) first bin holds value of correlation (single event). Second bin holds the weight! 
   
  } // end of for(Int_t c=0;c<fMaxCorrelator;c++) 
 } // end of for(Int_t cs=0;cs<2;cs++) 

} // void Cosmetics()

//=======================================================================================================================

 Bool_t AliAnalysisTaskStudentsMW::TrackSelection(AliAODTrack *aTrack)
 {

        // example variables for each track:
 	/*Double_t px = aTrack->Px(); // x-component of momenta
 	Double_t py = aTrack->Py(); // y-component of momenta
 	Double_t pz = aTrack->Pz(); // z-component of momenta
 	Double_t e = aTrack->E();  // energy
        Double_t charge = aTrack->Charge(); // charge
 	Double_t phi = aTrack->Phi(); // azimuthal angle*/
 	Double_t eta = aTrack->Eta(); // pseudorapidity
 	Double_t pt = aTrack->Pt(); // Pt (transverse momentum)

	if(bCutOnEta) 
	{
	  if(eta<fMinEtaCut) return kFALSE;
	  if(eta>fMaxEtaCut) return kFALSE;
	}

	if(bCutOnPt) 
	{
	  if(pt<fMinPtCut) return kFALSE;
	  if(pt>fMaxPtCut) return kFALSE;
	}

    return kTRUE;

 }// end AliAnalysisTaskStudentsMW::PhysicsSelection(AliAODTrack *aTrack)

//=======================================================================================================================

void AliAnalysisTaskStudentsMW::CalculateQvectors()
{
 // Calculate Q-vectors.

 // a) Make sure all Q-vectors are initially zero;
 // b) Calculate Q-vectors for available angles and weights. 

 // a) Make sure all Q-vectors are initially zero:
 for(Int_t h=0;h<49;h++)
 {
  for(Int_t p=0;p<9;p++)
  {
   fQvector[h][p] = TComplex(0.,0.);
  } //  for(Int_t p=0;p<kMaxPower;p++)
 } // for(Int_t h=0;h<kMaxHarmonic;h++)

 // b) Calculate Q-vectors for available angles and weights: 
 Double_t dPhi2 = 0.; // particle angle
 Double_t wPhi = 1.; // particle weight
 Double_t wPhiToPowerP = 1.; // particle weight raised to power p
 for(Int_t i=0;i<fParticles;i++) // loop over particles
 {
  dPhi2 = fAngles->GetAt(i);
  if(bUseWeights){wPhi = fWeights->GetAt(i);}
  for(Int_t h=0;h<49;h++)
  {
   for(Int_t p=0;p<9;p++)
   {
    if(bUseWeights){wPhiToPowerP = pow(wPhi,p);}
    fQvector[h][p] += TComplex(wPhiToPowerP*TMath::Cos(h*dPhi2),wPhiToPowerP*TMath::Sin(h*dPhi2));
   } //  for(Int_t p=0;p<kMaxPower;p++)
  } // for(Int_t h=0;h<kMaxHarmonic;h++)
 } //  for(Int_t i=0;i<fParticles;i++) // loop over particles


} // void CalculateQvectors()

//=======================================================================================================================


TComplex AliAnalysisTaskStudentsMW::Q(Int_t n, Int_t p)
{
 // Using the fact that Q{-n,p} = Q{n,p}^*. 
 
 if(n>=0){return fQvector[n][p];} 
 return TComplex::Conjugate(fQvector[-n][p]);
 
} // TComplex AliAnalysisTaskStudentsMW::Q(Int_t n, Int_t p)


//========================================================================================================================


TComplex AliAnalysisTaskStudentsMW::Recursion(Int_t n, Int_t* harmonic, Int_t mult = 1, Int_t skip = 0) 
{
 // Calculate multi-particle correlators by using recursion (an improved faster version) originally developed by 
 // Kristjan Gulbrandsen (gulbrand@nbi.dk). 

  Int_t nm1 = n-1;
  TComplex c(Q(harmonic[nm1], mult));
  if (nm1 == 0) return c;
  c *= Recursion(nm1, harmonic);
  if (nm1 == skip) return c;

  Int_t multp1 = mult+1;
  Int_t nm2 = n-2;
  Int_t counter1 = 0;
  Int_t hhold = harmonic[counter1];
  harmonic[counter1] = harmonic[nm2];
  harmonic[nm2] = hhold + harmonic[nm1];
  TComplex c2(Recursion(nm1, harmonic, multp1, nm2));
  Int_t counter2 = n-3;
  while (counter2 >= skip) {
    harmonic[nm2] = harmonic[counter1];
    harmonic[counter1] = hhold;
    ++counter1;
    hhold = harmonic[counter1];
    harmonic[counter1] = harmonic[nm2];
    harmonic[nm2] = hhold + harmonic[nm1];
    c2 += Recursion(nm1, harmonic, multp1, counter2);
    --counter2;
  }
  harmonic[nm2] = harmonic[counter1];
  harmonic[counter1] = hhold;

  if (mult == 1) return c-c2;
  return c-Double_t(mult)*c2;

} // TComplex AliFlowAnalysisWithMultiparticleCorrelations::Recursion(Int_t n, Int_t* harmonic, Int_t mult = 1, Int_t skip = 0) 



//========================================================================================================================


void AliAnalysisTaskStudentsMW::Correlation(Int_t Number, Int_t h1, Int_t h2, Int_t h3, Int_t h4, Int_t h5, Int_t h6, Int_t h7, Int_t h8)
{
	
     if(h1+h2+h3+h4+h5+h6+h7+h8!=0.){return;} //protection against anisotropic correlators
	
    // Calculate n-particle correlations from Q-vectors (using recursion):	

         
        if(2==Number)
        {
         Int_t harmonicsTwoNum[2] = {h1,h2};     
         Int_t harmonicsTwoDen[2] = {0,0};       
         TComplex twoRecursion = Recursion(2,harmonicsTwoNum)/Recursion(2,harmonicsTwoDen).Re();
         Double_t wTwoRecursion = Recursion(2,harmonicsTwoDen).Re();
         fRecursion[0][0]->Fill(0.5,twoRecursion.Re()); // <<cos(h1*phi1+h2*phi2)>>
         fRecursion[0][0]->Fill(1.5,wTwoRecursion); //weight 
         fRecursion[1][0]->Fill(0.5,twoRecursion.Im()); // <<sin(h1*phi1+h2*phi2)>>
         fRecursion[1][0]->Fill(1.5,wTwoRecursion); //weight
	
         }//  2-p correlation
        
        if(3==Number)
        {
         Int_t harmonicsThreeNum[3] = {h1,h2,h3};       
         Int_t harmonicsThreeDen[3] = {0,0,0};       
         TComplex threeRecursion = Recursion(3,harmonicsThreeNum)/Recursion(3,harmonicsThreeDen).Re();
         Double_t wThreeRecursion = Recursion(3,harmonicsThreeDen).Re();
         fRecursion[0][1]->Fill(0.5,threeRecursion.Re()); // <<cos(h1*phi1+h2*phi2+h3*phi3)>>
         fRecursion[0][1]->Fill(1.5,wThreeRecursion); //weight
         fRecursion[1][1]->Fill(0.5,threeRecursion.Im()); // <<sin(h1*phi1+h2*phi2+h3*phi3)>>
         fRecursion[1][1]->Fill(1.5,wThreeRecursion); //weight

         } //  3-p correlation
        
        if(4==Number)
        {
         Int_t harmonicsFourNum[4] = {h1,h2,h3,h4};       
         Int_t harmonicsFourDen[4] = {0,0,0,0};       
         TComplex fourRecursion = Recursion(4,harmonicsFourNum)/Recursion(4,harmonicsFourDen).Re();
         Double_t wFourRecursion = Recursion(4,harmonicsFourDen).Re();
         fRecursion[0][2]->Fill(0.5,fourRecursion.Re()); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4)>>
         fRecursion[0][2]->Fill(1.5,wFourRecursion); //weight
         fRecursion[1][2]->Fill(0.5,fourRecursion.Im()); // <<<sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4)>>
         fRecursion[1][2]->Fill(1.5,wFourRecursion); //weight

         }//  4-p correlation
        
        if(5==Number)
        {
         Int_t harmonicsFiveNum[5] = {h1,h2,h3,h4,h5};       
         Int_t harmonicsFiveDen[5] = {0,0,0,0,0};       
         TComplex fiveRecursion = Recursion(5,harmonicsFiveNum)/Recursion(5,harmonicsFiveDen).Re();
         Double_t wFiveRecursion = Recursion(5,harmonicsFiveDen).Re();
         fRecursion[0][3]->Fill(0.5,fiveRecursion.Re()); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5)>>
         fRecursion[0][3]->Fill(1.5,wFiveRecursion);
         fRecursion[1][3]->Fill(0.5,fiveRecursion.Im()); // <<<sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5)>>
         fRecursion[1][3]->Fill(1.5,wFiveRecursion);
        }//  5-p correlation

        if(6==Number)
        {
         Int_t harmonicsSixNum[6] = {h1,h2,h3,h4,h5,h6};       
         Int_t harmonicsSixDen[6] = {0,0,0,0,0,0};       
         TComplex sixRecursion = Recursion(6,harmonicsSixNum)/Recursion(6,harmonicsSixDen).Re();
         Double_t wSixRecursion = Recursion(6,harmonicsSixDen).Re();
         fRecursion[0][4]->Fill(0.5,sixRecursion.Re()); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6)>>
         fRecursion[0][4]->Fill(1.5,wSixRecursion);
         fRecursion[1][4]->Fill(0.5,sixRecursion.Im()); // <<<sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6)>>
         fRecursion[1][4]->Fill(1.5,wSixRecursion);

         }//  6-p correlation
        
        
        if(7==Number)
        {
         Int_t harmonicsSevenNum[7] = {h1,h2,h3,h4,h5,h6,h7};       
         Int_t harmonicsSevenDen[7] = {0,0,0,0,0,0,0};       
         TComplex sevenRecursion = Recursion(7,harmonicsSevenNum)/Recursion(7,harmonicsSevenDen).Re();
         Double_t wSevenRecursion = Recursion(7,harmonicsSevenDen).Re();
         fRecursion[0][5]->Fill(0.5,sevenRecursion.Re()); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7)>>
         fRecursion[0][5]->Fill(1.5,wSevenRecursion);
         fRecursion[1][5]->Fill(0.5,sevenRecursion.Im()); // <<<sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7)>>
         fRecursion[1][5]->Fill(1.5,wSevenRecursion);
  

        }//  7-p correlation
        
        
        if(8==Number)
        {
         Int_t harmonicsEightNum[8] = {h1,h2,h3,h4,h5,h6,h7,h8};       
         Int_t harmonicsEightDen[8] = {0,0,0,0,0,0,0,0};       
         TComplex eightRecursion = Recursion(8,harmonicsEightNum)/Recursion(8,harmonicsEightDen).Re();
         Double_t wEightRecursion = Recursion(8,harmonicsEightDen).Re();
         fRecursion[0][6]->Fill(0.5,eightRecursion.Re()); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7+h8*phi8)>>
         fRecursion[0][6]->Fill(1.5,wEightRecursion);
         fRecursion[1][6]->Fill(0.5,eightRecursion.Im()); // <<<sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7+h8*phi8)>>
         fRecursion[1][6]->Fill(1.5,wEightRecursion);
        
        }//  8-p correlation

        if(Number!=2 && Number!=3 && Number!=4 && Number!=5 && Number!=6 && Number!=7 && Number!=8)
        {
         return;
        }
      
 }//void Correlation() 

//========================================================================================================================

