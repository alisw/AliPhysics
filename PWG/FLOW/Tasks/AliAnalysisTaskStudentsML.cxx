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
#include <TArrayF.h>
#include <vector>
#include "TMath.h"
#include "TF1.h"



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
 bUseWeights(kFALSE),
 kNumber(2),  //number of correlation
 kh1(2), kh2(-2), kh3(0), kh4(0), kh5(0), kh6(0), kh7(0), kh8(0),  //harmonics
 kSum((kh1<0?-1*kh1:kh1)+(kh2<0?-1*kh2:kh2)+(kh3<0?-1*kh3:kh3)+(kh4<0?-1*kh4:kh4)
                + (kh5<0?-1*kh5:kh5)+(kh6<0?-1*kh6:kh6)+(kh7<0?-1*kh7:kh7)+(kh8<0?-1*kh8:kh8)), // We will not go beyond 8-p correlations
 kMaxHarmonic(kSum+1),
 kMaxPower(fMaxCorrelator+1), 
 fParticles(0),
 fCentral(0),
 fMinCentrality(0.),
 fMaxCentrality(100.),
 fAngles(NULL),
 fWeights(NULL),
 fBin(NULL),
 func1(NULL),
 fCentrality(NULL),
 fCounterHistogram(NULL),
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
 bUseWeights(kFALSE),
 kNumber(2),  //number of correlation
 kh1(2), kh2(-2), kh3(0), kh4(0), kh5(0), kh6(0), kh7(0), kh8(0),  //harmonics
 kSum((kh1<0?-1*kh1:kh1)+(kh2<0?-1*kh2:kh2)+(kh3<0?-1*kh3:kh3)+(kh4<0?-1*kh4:kh4)
                + (kh5<0?-1*kh5:kh5)+(kh6<0?-1*kh6:kh6)+(kh7<0?-1*kh7:kh7)+(kh8<0?-1*kh8:kh8)), // We will not go beyond 8-p correlations
 kMaxHarmonic(kSum+1),
 kMaxPower(fMaxCorrelator+1),
 fParticles(0), 
 fCentral(0.),
 fMinCentrality(0.),
 fMaxCentrality(100.), 
 fAngles(NULL),
 fWeights(NULL),
 fBin(NULL),
 func1(NULL),
 // Final results:
 fCentrality(NULL),
 fCounterHistogram(NULL),
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

} // void AliAnalysisTaskStudentsML::UserCreateOutputObjects() 

//================================================================================================================

void AliAnalysisTaskStudentsML::UserExec(Option_t *) 
{
 // Main loop (called for each event).

 // a.0) Get pointer to AOD event;
 // a.1) Centrality; 
 // b) Start analysis over AODs;
 // c) Reset event-by-event objects;
 // d) PostData.
 
 fCounterHistogram->Fill(0.5); // counter hist 1st bin

 // a) Get pointer to AOD event:
 AliAODEvent *aAOD = dynamic_cast<AliAODEvent*>(InputEvent()); // from TaskSE
 if(!aAOD){return;}
 
 fCounterHistogram->Fill(1.5); // counter hist 2nd bin

 //a.1) Centrality;
 
 AliMultSelection *ams = (AliMultSelection*)aAOD->FindListObject("MultSelection");
 if(!ams){return;}
 
 fCounterHistogram->Fill(2.5); // counter hist 3rd bin

 //func1->SetParameter(2,gRandom->Uniform(0.,TMath::TwoPi())); //*** for testing. for each event psi_2 is different


 // b) Start analysis over AODs:
 Int_t nTracks = aAOD->GetNumberOfTracks(); // number of all tracks in current event 
 
 //fParticles=gRandom->Uniform(50.0,500.0); /*** for testing.

 Int_t nCounter=0;
 
 
 //first loop: =======================

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
   
   fParticles += 1; //one more particle passed the cuts
  
   
   // Do some analysis only with the particles which passed the cuts
   // ... your analysis code ... 

 	

  } // if( (-0.8 < eta) && (eta < 0.8) && (0.2 < pt) && (pt < 5.0) )  
 } // for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks
   
 fAngles = new TArrayD(fParticles);

 //second loop: ======================= 

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
  
   fAngles->AddAt(phi,nCounter);
   nCounter += 1;

  } // if( (-0.8 < eta) && (eta < 0.8) && (0.2 < pt) && (pt < 5.0) )  
 } // for(Int_t iTrack=0;iTrack<nTracks;iTrack++) // starting a loop over all tracks
 
 
 //b.1) analysis
 if(fParticles>0){fMultiHisto->Fill(fParticles);}
 
 //*** for testing.
 /*fAngles = new TArrayD(fParticles); 
 
 for(Int_t k=0;k<fParticles;k++) 
  {  
   fAngles->AddAt(func1->GetRandom(),nCounter);
   nCounter += 1;
  } */
 
 if(fParticles>=8) //do the correlation only if there are more than 8 particles in the event
 { 

   if(ams->GetMultiplicityPercentile("V0M") >= fMinCentrality && ams->GetMultiplicityPercentile("V0M") < fMaxCentrality)
   {
    if(0. <= (ams->GetMultiplicityPercentile("V0M")) && (ams->GetMultiplicityPercentile("V0M")) < 5. ){fCentral=2.5;}
    if(5. <= (ams->GetMultiplicityPercentile("V0M")) && (ams->GetMultiplicityPercentile("V0M")) < 10. ){fCentral=7.5;}
    if(10. <= (ams->GetMultiplicityPercentile("V0M")) && (ams->GetMultiplicityPercentile("V0M")) < 20. ){fCentral=15.5;}
    if(20. <= (ams->GetMultiplicityPercentile("V0M")) && (ams->GetMultiplicityPercentile("V0M")) < 30. ){fCentral=25.5;}
    if(30. <= (ams->GetMultiplicityPercentile("V0M")) && (ams->GetMultiplicityPercentile("V0M")) < 40. ){fCentral=35.5;}
    if(40. <= (ams->GetMultiplicityPercentile("V0M")) && (ams->GetMultiplicityPercentile("V0M")) < 50. ){fCentral=45.5;}
    if(50. <= (ams->GetMultiplicityPercentile("V0M")) && (ams->GetMultiplicityPercentile("V0M")) < 60. ){fCentral=55.5;}
    if(60. <= (ams->GetMultiplicityPercentile("V0M")) && (ams->GetMultiplicityPercentile("V0M")) < 70. ){fCentral=65.5;}
    if(70. <= (ams->GetMultiplicityPercentile("V0M")) && (ams->GetMultiplicityPercentile("V0M")) < 80. ){fCentral=75.5;}
   }
   else
   {
    return; // this event do not belong to the centrality class specified for this particular analysis
   }


    Correlation();  //do the correlation
    fCentrality->Fill(fCentral,fRecursion[0][kNumber-2]->GetBinContent(1),fRecursion[0][kNumber-2]->GetBinContent(2)); //safe output 
    
    
    fRecursion[0][kNumber-2]->Reset(); //Reset
    fRecursion[1][kNumber-2]->Reset(); //Reset
 
 
    //Testing with nested loops
    
   /* for(Int_t m=0;m<fParticles;m++) //making all possible pairs for 500 particles 
    { 
      for(Int_t k=0;k<fParticles;k++) 
     	{  
          if(m==k){continue;} //no autocorellation
          else
          { 
 	   fCentrality->Fill(55.5,TMath::Cos(2.0*( (fAngles->GetAt(m)) - (fAngles->GetAt(k)) ) )); //fill TProfile
	  }
	} //for(Int_t k=0;k<fParticles;k++)
    } //for(Int_t m=0;m<fParticles;m++) 
  */


 } //if(fParticles>=8)

 

 

 // c) Reset event-by-event objects:
 fParticles=0;
 fCentral=0.;
 nCounter=0;
 delete fAngles;
 fAngles=NULL;


 
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

   for(Int_t cs=0;cs<2;cs++) 
   {
     for(Int_t c=0;c<fMaxCorrelator;c++)
     {
   
      fRecursion[cs][c] = NULL; //! [cs]: real (0) or imaginary part (1) ....
   
     }  
    }  //for(Int_t cs=0;cs<2;cs++)

   for(Int_t js=0;js<17;js++) 
   {
     for(Int_t j=0;j<9;j++)
     {
   
      Qvector[js][j] = NULL; //! 
   
     } 
   } 
   
   

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
 // d) Book histogram to hold multiplicty distribution;
 // e) Book histogram to debug

 
 

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
 fMultiHisto = new TH1F("fMultiHisto","Multiplicity",500,0.,500.); //histogram for multiplicity
 fMultiHisto->GetXaxis()->SetTitle("Multiplicity M");
 fControlHistogramsList->Add(fMultiHisto);

 // e) Book histogram to debug
 fCounterHistogram = new TH1F("fCounterHistogram","Histogram for some checks",3,0.,3.); //histogram for multiplicity
 fControlHistogramsList->Add(fCounterHistogram);

} //void AliAnalysisTaskStudentsML::BookControlHistograms()

//=======================================================================================================================

void AliAnalysisTaskStudentsML::BookFinalResultsHistograms()
{
 // Book all histograms to hold the final results.
  
 fCentrality = new TProfile("fCentrality","Result Analysis Centrality Dependence",80,0.,80.); //histogram for multiplicity
 fCentrality->GetXaxis()->SetTitle("Centrality");
 fCentrality->GetYaxis()->SetTitle("flow");
 fFinalResultsList->Add(fCentrality);
 fCentrality->Sumw2();

 

 Cosmetics();
 
 //*** for testing 
 /*func1=new TF1("func1","(1/(TMath::TwoPi()))*(1 + 2*[1]*TMath::Cos(2*(x-[2])))",0.,TMath::TwoPi()); //defining function func1 with parameter [1] for v_2 and [2] as pis_2
 func1->SetParameter(1,0.05); //setting v_2 (aka [1]) fix to 0.05*/

} // void AliAnalysisTaskStudentsML::BookFinalResultsHistograms()

//=======================================================================================================================

void AliAnalysisTaskStudentsML::Cosmetics()
{
 // Book everything here.
  
 for(Int_t cs=0;cs<2;cs++) 
 {
  for(Int_t c=0;c<fMaxCorrelator;c++)
  {
   
   fRecursion[cs][c] = new TProfile("","",2,0.,2.); 
   fRecursion[cs][c]->Sumw2();
 
   //NOTE for fRecursion: 1.) [cs] will say if its the real (0) or imaginary part (1) 
   // 2.) [c] gives gives the kind of correlation. [n] is the (n+2)-particle correlation
   //3.) first bin holds value of correlation (single event). Second bin holds the weight! 
   
   
  } // end of for(Int_t c=0;c<fMaxCorrelator;c++) 
 } // end of for(Int_t cs=0;cs<2;cs++) 

} // void Cosmetics()

//=======================================================================================================================

void AliAnalysisTaskStudentsML::CalculateQvectors()
{
 // Calculate Q-vectors.

 // a) Make sure all Q-vectors are initially zero;
 // b) Calculate Q-vectors for available angles and weights. 

 // a) Make sure all Q-vectors are initially zero:
 for(Int_t h=0;h<kMaxHarmonic;h++)
 {
  for(Int_t p=0;p<kMaxPower;p++)
  {
   Qvector[h][p] = TComplex(0.,0.);
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
  for(Int_t h=0;h<kMaxHarmonic;h++)
  {
   for(Int_t p=0;p<kMaxPower;p++)
   {
    if(bUseWeights){wPhiToPowerP = pow(wPhi,p);}
    Qvector[h][p] += TComplex(wPhiToPowerP*TMath::Cos(h*dPhi2),wPhiToPowerP*TMath::Sin(h*dPhi2));
   } //  for(Int_t p=0;p<kMaxPower;p++)
  } // for(Int_t h=0;h<kMaxHarmonic;h++)
 } //  for(Int_t i=0;i<fParticles;i++) // loop over particles

} // void CalculateQvectors()

//=======================================================================================================================


TComplex AliAnalysisTaskStudentsML::Q(Int_t n, Int_t p)
{
 // Using the fact that Q{-n,p} = Q{n,p}^*. 
 
 if(n>=0){return Qvector[n][p];} 
 return TComplex::Conjugate(Qvector[-n][p]);
 
} // TComplex AliAnalysisTaskStudentsML::Q(Int_t n, Int_t p)




//========================================================================================================================


TComplex AliAnalysisTaskStudentsML::Recursion(Int_t n, Int_t* harmonic, Int_t mult = 1, Int_t skip = 0) 
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


void AliAnalysisTaskStudentsML::Correlation()
{
	
    
    // Calculate Q-vectors for available angles and weights;
    CalculateQvectors();
	
    // Calculate n-particle correlations from Q-vectors (using recursion):	

         
        if(kNumber==2)
        {
         Int_t harmonicsTwoNum[2] = {kh1,kh2};     
         Int_t harmonicsTwoDen[2] = {0,0};       
         TComplex twoRecursion = Recursion(2,harmonicsTwoNum)/Recursion(2,harmonicsTwoDen).Re();
         Double_t wTwoRecursion = Recursion(2,harmonicsTwoDen).Re();
         fRecursion[0][0]->Fill(0.5,twoRecursion.Re(),wTwoRecursion); // <<cos(h1*phi1+h2*phi2)>>
         fRecursion[0][0]->Fill(1.5,wTwoRecursion,1.); //weight 
         fRecursion[1][0]->Fill(0.5,twoRecursion.Im(),wTwoRecursion); // <<sin(h1*phi1+h2*phi2)>>
         fRecursion[1][0]->Fill(1.5,wTwoRecursion,1.); //weight
         }//  2-p correlation
        
        if(kNumber==3)
        {
         Int_t harmonicsThreeNum[3] = {kh1,kh2,kh3};       
         Int_t harmonicsThreeDen[3] = {0,0,0};       
         TComplex threeRecursion = Recursion(3,harmonicsThreeNum)/Recursion(3,harmonicsThreeDen).Re();
         Double_t wThreeRecursion = Recursion(3,harmonicsThreeDen).Re();
         fRecursion[0][1]->Fill(0.5,threeRecursion.Re(),wThreeRecursion); // <<cos(h1*phi1+h2*phi2+h3*phi3)>>
         fRecursion[0][1]->Fill(1.5,wThreeRecursion,1.); //weight
         fRecursion[1][1]->Fill(0.5,threeRecursion.Im(),wThreeRecursion); // <<sin(h1*phi1+h2*phi2+h3*phi3)>>
         fRecursion[1][1]->Fill(1.5,wThreeRecursion,1.); //weight
         } //  3-p correlation
        
        if(kNumber==4)
        {
         Int_t harmonicsFourNum[4] = {kh1,kh2,kh3,kh4};       
         Int_t harmonicsFourDen[4] = {0,0,0,0};       
         TComplex fourRecursion = Recursion(4,harmonicsFourNum)/Recursion(4,harmonicsFourDen).Re();
         Double_t wFourRecursion = Recursion(4,harmonicsFourDen).Re();
         fRecursion[0][2]->Fill(0.5,fourRecursion.Re(),wFourRecursion); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4)>>
         fRecursion[0][2]->Fill(1.5,wFourRecursion,1.); //weigh
         fRecursion[1][2]->Fill(0.5,fourRecursion.Im(),wFourRecursion); // <<<sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4)>>
         fRecursion[1][2]->Fill(1.5,wFourRecursion,1.); //weigh
         }//  4-p correlation
        
        if(kNumber==5)
        {
         Int_t harmonicsFiveNum[5] = {kh1,kh2,kh3,kh4,kh5};       
         Int_t harmonicsFiveDen[5] = {0,0,0,0,0};       
         TComplex fiveRecursion = Recursion(5,harmonicsFiveNum)/Recursion(5,harmonicsFiveDen).Re();
         Double_t wFiveRecursion = Recursion(5,harmonicsFiveDen).Re();
         fRecursion[0][3]->Fill(0.5,fiveRecursion.Re(),wFiveRecursion); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5)>>
         fRecursion[0][3]->Fill(1.5,wFiveRecursion,1.);
         fRecursion[1][3]->Fill(0.5,fiveRecursion.Im(),wFiveRecursion); // <<<sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5)>>
         fRecursion[1][3]->Fill(1.5,wFiveRecursion,1.);
        }//  5-p correlation

        if(kNumber==6)
        {
         Int_t harmonicsSixNum[6] = {kh1,kh2,kh3,kh4,kh5,kh6};       
         Int_t harmonicsSixDen[6] = {0,0,0,0,0,0};       
         TComplex sixRecursion = Recursion(6,harmonicsSixNum)/Recursion(6,harmonicsSixDen).Re();
         Double_t wSixRecursion = Recursion(6,harmonicsSixDen).Re();
         fRecursion[0][4]->Fill(0.5,sixRecursion.Re(),wSixRecursion); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6)>>
         fRecursion[0][4]->Fill(1.5,wSixRecursion,1.);
         fRecursion[1][4]->Fill(0.5,sixRecursion.Im(),wSixRecursion); // <<<sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6)>>
         fRecursion[1][4]->Fill(1.5,wSixRecursion,1.);
         }//  6-p correlation
        
        
        if(kNumber==7)
        {
         Int_t harmonicsSevenNum[7] = {kh1,kh2,kh3,kh4,kh5,kh6,kh7};       
         Int_t harmonicsSevenDen[7] = {0,0,0,0,0,0,0};       
         TComplex sevenRecursion = Recursion(7,harmonicsSevenNum)/Recursion(7,harmonicsSevenDen).Re();
         Double_t wSevenRecursion = Recursion(7,harmonicsSevenDen).Re();
         fRecursion[0][5]->Fill(0.5,sevenRecursion.Re(),wSevenRecursion); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7)>>
         fRecursion[0][5]->Fill(1.5,wSevenRecursion,1.);
         fRecursion[1][5]->Fill(0.5,sevenRecursion.Im(),wSevenRecursion); // <<<sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7)>>
         fRecursion[1][5]->Fill(1.5,wSevenRecursion,1.);
        }//  7-p correlation
        
        
        if(kNumber==8)
        {
         Int_t harmonicsEightNum[8] = {kh1,kh2,kh3,kh4,kh5,kh6,kh7,kh8};       
         Int_t harmonicsEightDen[8] = {0,0,0,0,0,0,0,0};       
         TComplex eightRecursion = Recursion(8,harmonicsEightNum)/Recursion(8,harmonicsEightDen).Re();
         Double_t wEightRecursion = Recursion(8,harmonicsEightDen).Re();
         fRecursion[0][6]->Fill(0.5,eightRecursion.Re(),wEightRecursion); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7+h8*phi8)>>
         fRecursion[0][6]->Fill(1.5,wEightRecursion,1.);
         fRecursion[1][6]->Fill(0.5,eightRecursion.Im(),wEightRecursion); // <<<sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7+h8*phi8)>>
         fRecursion[1][6]->Fill(1.5,wEightRecursion,1.);
        }//  8-p correlation

        if(kNumber!=2 && kNumber!=3 && kNumber!=4 && kNumber!=5 && kNumber!=6 && kNumber!=7 && kNumber!=8)
        {
         return;
        }

    
      
 }//void Correlation() 


//========================================================================================================================


