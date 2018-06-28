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

/*******************************************************************************
* Analysis task for anisotropic flow analysis of data taken by ALICE 				   *
* with different methods for two- and multiparticle correlations     				   *
*																			                                         *
* Author: Cindy Mordasini (cindy.mordasini@cern.ch)		            					   *
*******************************************************************************/

#include "Riostream.h"
#include "AliAnalysisTaskTwoMultiCorrelations.h"
#include "AliLog.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"
#include "TFile.h"
#include "TComplex.h"
#include "TMath.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskTwoMultiCorrelations)

//==============================================================================

AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights):
  AliAnalysisTaskSE(name),
  fNparticlesCorrelations(2),          // Number of m-particle correlations and harmonics (2-14)
  fHarmonicOne(2),                     // Harmonic n_1, default value for v_2{2}
  fHarmonicTwo(-2),                    // Harmonic n_2, default value for v_2{4}
  fHarmonicThree(0),                   // Harmonic n_3
  fHarmonicFour(0),                    // Harmonic n_4
  fHarmonicFive(0),                    // Harmonic n_5
  fHarmonicSix(0),                     // Harmonic n_6
  fHarmonicSeven(0),                   // Harmonic n_7
  fHarmonicEight(0),                   // Harmonic n_8
  fHarmonicNine(0),                    // Harmonic n_9
  fHarmonicTen(0),                     // Harmonic n_10
  fHarmonicEleven(0),                  // Harmonic n_11
  fHarmonicTwelve(0),                  // Harmonic n_12
  fHarmonicThirteen(0),                // Harmonic n_13
  fHarmonicFourteen(0),                // Harmonic n_14
  fMinCentrality(0.0),                 // Minimum of centrality
  fMaxCentrality(100.0),               // Maximum of centrality
  fUseParticleWeights(kFALSE),         // Use non-unit particle weights
  fDoNestedLoops(kFALSE),              // Cross-check the results with nested loops
  fOutputList(NULL),                   // Main list holding all the output objects
  fControlOutputList(NULL),            // List holding all the control objects
  fDraftOutputList(NULL),              // List holding all the intermediate objects
  fFinalOutputList(NULL),              // List holding all the final results
  fCorrelationWithQvectorsProfile(NULL),  // m-particle correlation estimated with Q-vectors
  fCorrelationWithNestedLoopsProfile(NULL),  // 2-p correlation estimated with nested loops
  fCorrelationWithQvectorsSaProfile(NULL),   // 2-p correlation estimated with stand-alone Q-vectors
  fControlPhiHisto(NULL),              // Control histogram for the azimuthal angles
  fCentralityHisto(NULL),              // Control histogram for the centrality
  fAverageMulti(NULL)                  // Control histogram for the average multiplicity
{
// Constructor of the class

    AliDebug(2, "AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights)");

    // Creation of a new main list
    fOutputList = new TList();
    fOutputList->SetName("outputAnalysis");
    fOutputList->SetOwner(kTRUE);

    // Definition of the input and output slots
    DefineOutput(1, TList::Class());

  if(useParticleWeights)
  {
    // not needed for the time being, maybe insert here the call for the file with particle weights???
  }

} // End of the constructor

//******************************************************************************

AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations():
  AliAnalysisTaskSE(),
  fNparticlesCorrelations(2),          // Number of m-particle correlations and harmonics (2-14)
  fHarmonicOne(2),                     // Harmonic n_1
  fHarmonicTwo(-2),                    // Harmonic n_2
  fHarmonicThree(0),                   // Harmonic n_3
  fHarmonicFour(0),                    // Harmonic n_4
  fHarmonicFive(0),                    // Harmonic n_5
  fHarmonicSix(0),                     // Harmonic n_6
  fHarmonicSeven(0),                   // Harmonic n_7
  fHarmonicEight(0),                   // Harmonic n_8
  fHarmonicNine(0),                    // Harmonic n_9
  fHarmonicTen(0),                     // Harmonic n_10
  fHarmonicEleven(0),                  // Harmonic n_11
  fHarmonicTwelve(0),                  // Harmonic n_12
  fHarmonicThirteen(0),                // Harmonic n_13
  fHarmonicFourteen(0),                // Harmonic n_14
  fMinCentrality(0.0),                 // Minimum of centrality
  fMaxCentrality(100.0),               // Maximum of centrality
  fUseParticleWeights(kFALSE),         // Use non-unit particle weights
  fDoNestedLoops(kFALSE),              // Cross-check the results with nested loops
  fOutputList(NULL),                   // Main list holding all the output objects
  fControlOutputList(NULL),            // List holding all the control objects
  fDraftOutputList(NULL),              // List holding all the intermediate objects
  fFinalOutputList(NULL),              // List holding all the final results
  fCorrelationWithQvectorsProfile(NULL),  // m-particle correlation estimated with Q-vectors
  fCorrelationWithNestedLoopsProfile(NULL),  // 2-p correlation estimated with nested loops
  fCorrelationWithQvectorsSaProfile(NULL),   // 2-p correlation estimated with stand-alone Q-vectors
  fControlPhiHisto(NULL),              // Control histogram for the azimuthal angles
  fCentralityHisto(NULL),              // Control histogram for the centrality
  fAverageMulti(NULL)                  // Control histogram for the average multiplicity
{
// Dummy constructor of the class

    AliDebug(2, "AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights)");

} // End of the dummy constructor

//******************************************************************************

AliAnalysisTaskTwoMultiCorrelations::~AliAnalysisTaskTwoMultiCorrelations()
{
// Destructor of the class
/// Delete the main TList => delete automatically everything holds in it

  if(fOutputList) delete fOutputList;

} // End of the destructor

//==============================================================================

void AliAnalysisTaskTwoMultiCorrelations::UserCreateOutputObjects()
{
// Method called at every worker node to initialise the lists
// Organisation of the method
  // 1.) First part of the trick to avoid name clashes
  // 2.) Booking and nesting of all the lists
  // 3.) Booking of all the objects
  // 4.) Second part of the trick to avoid name clashes

 // 1.) First part of the trick to avoid name clashes
  Bool_t oldHistAddStatus = TH1::AddDirectoryStatus(); 
  TH1::AddDirectory(kFALSE);

 // 2.) Booking and nesting of all the lists
  this->BookAndNestAllLists();

 // 3.) Booking of all the objects
  this->BookControlList();
  this->BookDraftList();
  this->BookFinalList();

 // 4.) Second part of the trick to avoid name clashes
  TH1::AddDirectory(oldHistAddStatus);

  PostData(1,fOutputList);

} // End of AliAnalysisTaskTwoMultiCorrelations::UserCreateOutputObjects()

//******************************************************************************

void AliAnalysisTaskTwoMultiCorrelations::UserExec(Option_t *)
{
// Method called for each event, contains all the calculations
// Note to self: find a way to include non-unit particle weights and select them via fUseParticleWeights
// Organisation of the method
  // 1.) Obtention of the pointer to the AOD event
  // 2.) Check if the multiplicity is higher than the considered number of correlation (to avoid the problem of dividing by zero)
  // 3.) Gestion of the centrality
  // 4.) Definition of the variables common to all correlation methods
  // 5.) Filling of the azimuthal angles and particle weights
  // 6.) Filling of some control histograms
  // 7.) Computation of the correlation with different methods for the current event
  // 8.) Release of the allocated memory
  // 9.) PostData

// 1.) Obtention of the pointer to the AOD event
  AliAODEvent *currentEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!currentEvent){return;}

// 2.) Check if the multiplicity is higher than the considered number of correlation (to avoid the problem of dividing by zero)
  Int_t nParticles = currentEvent->GetNumberOfTracks();    // Multiplicity
  fAverageMulti->Fill(0.5, nParticles);
  if (nParticles >= fNparticlesCorrelations)
  {

// 3.) Gestion of the centrality
  AliMultSelection *ams = (AliMultSelection*)currentEvent->FindListObject("MultSelection");
  if(!ams){return;}
  if(ams->GetMultiplicityPercentile("V0M") >= fMinCentrality && ams->GetMultiplicityPercentile("V0M") < fMaxCentrality)
  {
    fCentralityHisto->Fill(ams->GetMultiplicityPercentile("V0M"));
  }
  else
  {
    return; // this event does not belong to the centrality class specified for this particular analysis
  }

// 4.) Definition of the variables common to all correlation methods

  Int_t harmonics[14] = {fHarmonicOne, fHarmonicTwo, fHarmonicThree, fHarmonicFour, fHarmonicFive, fHarmonicSix, fHarmonicSeven, fHarmonicEight, fHarmonicNine, fHarmonicTen, fHarmonicEleven, fHarmonicTwelve, fHarmonicThirteen, fHarmonicFourteen};  // Harmonics n_1,... n_14 (max)
  Double_t *phi = new Double_t[nParticles]();              // Azimuthal angles
  Double_t *particleWeights = new Double_t[nParticles]();  // Particle weights

  for (Int_t iParticle = 0; iParticle < nParticles; iParticle++)
  {
    // Pointer to a particle
    AliAODTrack *currentParticle = dynamic_cast<AliAODTrack*>(currentEvent->GetTrack(iParticle));

// 5.) Filling of the azimuthal angles and particle weights
    phi[iParticle] = currentParticle->Phi();
    particleWeights[iParticle] = 1.;
    //if (fUseParticleWeights) {TO IMPLEMENT WITH ACCESS TO EXTERNAL FILE}
    //else {particleWeights[iParticle] = 1.;}

// 6.) Filling of some control histograms
    // Azimuthal angles
    fControlPhiHisto->Fill(phi[iParticle]);

  } // End of for (Int_t iParticle = 0; iParticle < nParticles; iParticle++)

// 7.) Computation of the correlation with different methods for the current event
  // Method: Q-vectors with the recursion
  ComputeCorrelationsWithQvectors(nParticles, phi, particleWeights, harmonics, fNparticlesCorrelations);
  // Method: Nested loops, used only to cross-check the results of the Q-vectors
  if ((fDoNestedLoops) && (fNparticlesCorrelations == 2)) {ComputeCorrelationsWithTwoNestedLoops(fHarmonicOne, fHarmonicTwo, nParticles, phi, particleWeights);}  // 2-p correlation
  if ((fDoNestedLoops) && (fNparticlesCorrelations == 4)) {ComputeCorrelationsWithFourNestedLoops(fHarmonicOne, fHarmonicTwo, fHarmonicThree, fHarmonicFour, nParticles, phi, particleWeights);}   // 4-p correlation
  // Method: Q-vectors with stand-alone formula pour 2-p correlation (debug for recursion)
  if (fNparticlesCorrelations == 2) {ComputeCorrelationsWithStandAloneQvectors(fHarmonicOne, 0, nParticles, phi, particleWeights);}

// 8.) Release of the allocated memory
  delete [] phi;
  delete [] particleWeights;

  } // End of if (nParticles => fNparticlesCorrelations)

// 9.) PostData
  PostData(1,fOutputList);

} // End of void AliAnalysisTaskTwoMultiCorrelations::UserExec(Option_t *)

//******************************************************************************

void AliAnalysisTaskTwoMultiCorrelations::Terminate(Option_t *)
{
// Method called at the end of the execution, once the run over the events is over
// Organisation of the method
  // 1.) Access to the merged output list
  // 2.) Filling the final profiles
  // 3.) Offline calculations
  // 4.) Creation of the output file and save of the main list in it

// 1.) Access to the merged output list
  fOutputList = (TList*)GetOutputData(1);
  if(!fOutputList){exit(1);}

// 2.) Filling the final profiles

// 3.) Offline calculations

// 4.) Creation of the output file and save of the main list in it
  TFile *outputFile = new TFile("AnalysisResults.root", "RECREATE");
  fOutputList->Write(fOutputList->GetName(),TObject::kSingleKey);
  delete outputFile;

} // End of void AliAnalysisTaskTwoMultiCorrelations::Terminate(Option_t *)

//==============================================================================

void AliAnalysisTaskTwoMultiCorrelations::BookAndNestAllLists()
{
// Method to prepare all the lists with the results in the output file
// Organisation of the method
  // 1.) Booking and nesting lists for the control objects
  // 2.) Booking and nesting lists for the intermediate objects
  // 3.) Booking and nesting lists for the final results

  TString sMethodName = "void AliAnalysisTaskTwoMultiCorrelations::BookAndNestAllLists()";
  if(!fOutputList){Fatal(sMethodName.Data(),"Main list fOutputList is NULL");}

// 1.) Booking and nesting lists for the control objects
  fControlOutputList = new TList();
  fControlOutputList->SetName("ControlOutputList");
  fControlOutputList->SetOwner(kTRUE);
  fOutputList->Add(fControlOutputList);

// 2.) Booking and nesting lists for the intermediate objects
  fDraftOutputList = new TList();
  fDraftOutputList->SetName("DraftOutputList");
  fDraftOutputList->SetOwner(kTRUE);
  fOutputList->Add(fDraftOutputList);
  
// 3.) Booking and nesting lists for the final results
  fFinalOutputList = new TList();
  fFinalOutputList->SetName("FinalOutputList");
  fFinalOutputList->SetOwner(kTRUE);
  fOutputList->Add(fFinalOutputList);

} // End of AliAnalysisTaskTwoMultiCorrelations::BookAndNestAllLists()

//==============================================================================

void AliAnalysisTaskTwoMultiCorrelations::BookControlList()
{
// Method to prepare the list with the control histograms
// Organisation of the method
  // 1.) Control histogram for the azimuthal angles

// 1.) Control histogram for the azimuthal angles
  fControlPhiHisto = new TH1F("fControlPhiHisto", "Azimuthal angles distribution", 1000, 0., 6.3);
  fControlPhiHisto->SetStats(kTRUE);
  fControlPhiHisto->GetXaxis()->SetTitle("phi");
  fControlOutputList->Add(fControlPhiHisto);

// 2.) Control histogram for the centrality
  fCentralityHisto = new TH1F("fCentralityHisto", "Centrality distribution", 9, 0., 100.);
  fCentralityHisto->SetStats(kTRUE);
  fCentralityHisto->GetXaxis()->SetTitle("Centrality percentile");
  fControlOutputList->Add(fCentralityHisto);

// 3.) Control TProfile for the average multiplicity
  fAverageMulti = new TProfile("fAverageMulti", "Average multiplicity", 1, 0., 1.);
  fAverageMulti->GetXaxis()->SetTitle("Multiplicity");
  fControlOutputList->Add(fAverageMulti);

} // End of void AliAnalysisTaskTwoMultiCorrelations::BookControlList()

//******************************************************************************

void AliAnalysisTaskTwoMultiCorrelations::BookDraftList()
{
// Method to prepare the list with the intermediate results from the computation of Q-vectors and nested loops
// Organisation of the method
  // 1.) TProfile from the method of the Q-vectors
  // 2.) TProfile from the method of the two nested loops

// 1.) TProfile from the method of the Q-vectors
  fCorrelationWithQvectorsProfile = new TProfile("fCorrelationWithQvectorsProfile", "m-particle correlation with Q-vectors", 1, 0., 1.);
  fCorrelationWithQvectorsProfile->Sumw2();
  fCorrelationWithQvectorsProfile->GetXaxis()->SetTitle("m-particle correlation");
  fDraftOutputList->Add(fCorrelationWithQvectorsProfile);

// 2.) TProfile from the method of nested loops
  fCorrelationWithNestedLoopsProfile = new TProfile("fCorrelationWithNestedLoopsProfile", "m-particle correlation with nested loops", 1, 0., 1.);
  fCorrelationWithNestedLoopsProfile->Sumw2();
  fCorrelationWithNestedLoopsProfile->GetXaxis()->SetTitle("m-particle correlation");
  fDraftOutputList->Add(fCorrelationWithNestedLoopsProfile);

// *.) TProfile from the method of the Q-vectors, stand alone formula
  fCorrelationWithQvectorsSaProfile = new TProfile("fCorrelationWithQvectorsSaProfile", "m-particle correlation with Q-vectors, stand alone", 1, 0., 1.);
  fCorrelationWithQvectorsSaProfile->Sumw2();
  fCorrelationWithQvectorsSaProfile->GetXaxis()->SetTitle("2-particle correlation");
  fDraftOutputList->Add(fCorrelationWithQvectorsSaProfile);

} // End of void AliAnalysisTaskTwoMultiCorrelations::BookDraftList()

//******************************************************************************

void AliAnalysisTaskTwoMultiCorrelations::BookFinalList()
{

} // End of void AliAnalysisTaskTwoMultiCorrelations::BookFinalList()

//==============================================================================

TComplex AliAnalysisTaskTwoMultiCorrelations::CalculateQvector(Int_t n, Int_t p, Int_t nParticles, Double_t phi[], Double_t particleWeight[])
{
// Method calculating the "general" definition of the Q-vector Q_(n,p) for arbitrary (n,p)
// Organisation of the method
  // 1.) Declaration of the local variables: Q_(n,p)
  // 2.) Computation of the Q-vector
  // 3.) Application of the property Q_(-n,p) = Q_(n,p)* and return of the result

// 1.) Declaration of the local variables: Q_(n,p)
  TComplex qVectorNp = TComplex(0,0);  // Q-vector Q_(n,p) for the current event
  Double_t pWeightPowerP = 0.;         // Particle weight to the power p

// 2.) Computation of the Q-vector
  for (Int_t iParticle = 0; iParticle < nParticles; iParticle++)
  {
    pWeightPowerP = pow(particleWeight[iParticle], p);
    qVectorNp += pWeightPowerP * TComplex::Exp((1.*n)*(TComplex::I())*phi[iParticle]);
  } // End of for (Int_t iParticle = 0; iParticle < nParticles; iParticle++)

// 3.) Application of the property Q_(-n,p) = Q_(n,p)* and return of the result
  //if (n < 0) {return TComplex::Conjugate(qVectorNp);}
  if (n < 0) {CalculateQvector(-n, p, nParticles, phi, particleWeight);}  
  return qVectorNp;

} // End of TComplex AliAnalysisTaskTwoMultiCorrelations::CalculateQvector(Int_t n, Int_t p, Int_t nParticles, Double_t phi[], Double_t particleWeight[])

//******************************************************************************

TComplex AliAnalysisTaskTwoMultiCorrelations::CalculateRecursionWithQvectors(Int_t nParticles, Double_t phi[], Double_t particleWeight[], Int_t nCorr, Int_t harmonics[], Int_t p, Int_t skip)
{
// Method calculating the recursion with the Q-vectors for the numerator of the m-particle correlation according to the generic framework
// Inspired by the method originally developped by Kristjan Gulbrandsen (gulbrand@nbi.dk)
  /// p = mult and nCorr = n in previous version
// Organisation of the method
  // 1.) Declaration of the local variables
  // 2.) Stop conditions of the recursion
  // 3.) Computation of the recursion
  // 4.) Return of the result after the recursion

// 1.) Declaration of the local variables
  Int_t nMinusOne = 0;                 // Harmonic (nCorr-1)
  TComplex stopQvector = TComplex(0,0);// Q-vector used to stop the recursion
  Int_t pPlusOne = 0;                  // (p + 1)
  Int_t nMinusTwo = 0;                 // Harmonic (n-2)
  Int_t counterOne = 0;                // First counter for the intermediate indices
  Int_t hHold = 0;                     // Temporary harmonic
  TComplex tempQvector = TComplex(0,0);// Temporary Q-vector
  Int_t counterTwo = 0;                // Second counter for the intermediate indices

// 2.) Stop conditions of the recursion
  nMinusOne = nCorr-1;
  stopQvector = CalculateQvector(harmonics[nMinusOne], p, nParticles, phi, particleWeight);
 
  if (nMinusOne == 0) {return stopQvector;}
  stopQvector *= CalculateRecursionWithQvectors(nParticles, phi, particleWeight, nMinusOne, harmonics);
  if (nMinusOne == skip) {return stopQvector;}

// 3.) Computation of the recursion
  pPlusOne = p+1;
  nMinusTwo = nCorr-2;
  hHold = harmonics[counterOne];
  harmonics[counterOne] = harmonics[nMinusTwo];
  harmonics[nMinusTwo] = hHold + harmonics[nMinusOne];
  tempQvector = CalculateRecursionWithQvectors(nParticles, phi, particleWeight, nMinusOne, harmonics, pPlusOne, nMinusTwo);
  counterTwo = nCorr-3;

  while (counterTwo >= skip)
  {
    harmonics[nMinusTwo] = harmonics[counterOne];
    harmonics[counterOne] = hHold;
    ++counterOne;
    
    hHold = harmonics[counterOne];
    harmonics[counterOne] = harmonics[nMinusTwo];
    harmonics[nMinusTwo] = hHold + harmonics[nMinusOne];
    tempQvector += CalculateRecursionWithQvectors(nParticles, phi, particleWeight, nMinusOne, harmonics, pPlusOne, counterTwo);
    --counterTwo;

  } // End of the while (counterTwo >= skip)

  harmonics[nMinusTwo] = harmonics[counterOne];
  harmonics[counterOne] = hHold;

// 4.) Return of the result after the recursion
  if (p == 1) {return stopQvector - tempQvector;}
  return stopQvector - (Double_t(p)*tempQvector);

} // End of void AliAnalysisTaskTwoMultiCorrelations::CalculateRecursionWithQvectors(Int_t nParticles, Double_t *phi[], Double_t *particleWeight[], Int_t n, Double_t *harmonics[], Int_t p = 1, Int_t skip = 0)

//******************************************************************************

void AliAnalysisTaskTwoMultiCorrelations::ComputeCorrelationsWithQvectors(Int_t nParticles, Double_t phi[], Double_t particleWeight[], Int_t harmonics[], Int_t nCorr)
{
// Method to compute the m-particle correlation using the numerator for the generic framework computed with the recursion method
// Organisation of the method
  // 1.) Declaration of the local variables
  // 2.) Filling of the array for the numerator
  // 3.) Computation of the m-particle correlation

// 1.) Declaration of the local variables
  Int_t *numeratorHarmonics = new Int_t[nCorr];   // Harmonics with values from harmonics array
  Int_t *denominatorHarmonics = new Int_t[nCorr]; // Harmonics with value 0 for the denominator

  Double_t denominator = 0.;           // Denominator of <m>_(n1...nm) and event weight for the TProfile
  TComplex mParticleCorrelation = TComplex(0,0);  // m-particle correlation
  TComplex eventWeight = TComplex(0,0);// Event weight
 
// 2.) Filling of the array for the numerator
  for (Int_t iCorr = 0; iCorr < nCorr; iCorr++)
  {
    numeratorHarmonics[iCorr] = harmonics[iCorr];
    denominatorHarmonics[iCorr] = 0;
  } // End of for (Int_t iCorr = 0; iCorr < nCorr; iCorr++)

// 3.) Computation of the m-particle correlation
  denominator = CalculateRecursionWithQvectors(nParticles, phi, particleWeight, nCorr, denominatorHarmonics).Re();
  mParticleCorrelation = CalculateRecursionWithQvectors(nParticles, phi, particleWeight, nCorr, numeratorHarmonics)/denominator;
  //eventWeight = (CalculateRecursionWithQvectors(nParticles, phi, particleWeight, nCorr, denominatorHarmonics)).Re();

// 4.) Filling of the TProfile
  fCorrelationWithQvectorsProfile->Fill(0.5, mParticleCorrelation.Re(), denominator);

  delete [] numeratorHarmonics;
  delete [] denominatorHarmonics;

} // End of void AliAnalysisTaskTwoMultiCorrelations::ComputeCorrelationsWithQvectors(Int_t n, Int_t nParticles, Double_t *phi[], Double_t *particleWeight[], Int_t *harmonics[], Int_t nCorr)

//******************************************************************************

void AliAnalysisTaskTwoMultiCorrelations::ComputeCorrelationsWithTwoNestedLoops(Int_t nOne, Int_t nTwo, Int_t nParticles, Double_t phi[], Double_t particleWeight[])
{
// Method to compute the 2-particle correlation with two nested loops (for cross-checking results)
// Organisation of the method
  // 1.) Declaration of the local variables
  // 2.) Computation of the 2-p correlation, <2>
  // 3.) Reset of the local variables before changing the event

// 1.) Declaration of the local variables
  Double_t twoParticleCorrelation = 0.;     // Single-event 2-p correlation, <2>_(n1...nm)
  Double_t totalParticleWeight = 0.;        // Particle weight for the 2-p correlation

// 2.) Computation of the 2-p correlation, <2>
  for (Int_t k = 0; k < nParticles; k++)
  {
    for (Int_t l = 0; l < nParticles; l++)
    {
      // Removal of the autocorrelations k == l 
      if (k == l) continue;
      
      // Computation of <2> for a pair of particles
      twoParticleCorrelation = cos(nOne*phi[k] + nTwo*phi[l]);
      totalParticleWeight = particleWeight[k] * particleWeight[l];

      // Filling of the TProfile
      fCorrelationWithNestedLoopsProfile->Fill(0.5, twoParticleCorrelation, totalParticleWeight);

    } // End of the loop over the second particle of the pair
  } // End of the loop over the first particle of the pair

// 3.) Reset of the local variables before changing the event
  twoParticleCorrelation = 0.;
  totalParticleWeight = 0.;

} // End of void AliAnalysisTaskTwoMultiCorrelations::ComputeCorrelationsWithTwoNestedLoops(Int_t nOne, Int nTwo, Int_t nParticles, Double_t phi[], Double_t particleWeight[])

//******************************************************************************

void AliAnalysisTaskTwoMultiCorrelations::ComputeCorrelationsWithFourNestedLoops(Int_t nOne, Int_t nTwo, Int_t nThree, Int_t nFour, Int_t nParticles, Double_t phi[], Double_t particleWeight[])
{
// Method to compute the 4-particle correlation with four nested loops (for cross-checking results)
// Organisation of the method
  // 1.) Declaration of the local variables
  // 2.) Computation of the 4-p correlation, <4>
  // 3.) Reset of the local variables before changing the event

// 1.) Declaration of the local variables
  Double_t fourParticleCorrelation = 0.;    // Single-event 4-p correlation, <4>_(n1,... nm)
  Double_t totalParticleWeight = 0.;

// 2.) Computation of the 4-p correlation, <4>
  for (Int_t k = 0; k < nParticles; k++)
  {
    for (Int_t l = 0; l < nParticles; l++)
    {
      // Removal of the autocorrelations k == l
      if (k == l) continue;

      for (Int_t m = 0; m < nParticles; m++)
      {
        // Removal of the autocorrelations k == m, l == m
        if ((k == m) || (l == m)) continue;

        for (Int_t n = 0; n < nParticles; n++)
        {
          // Removal of the autocorrelations k == n, l == n, m == n
          if ((k == n) || (l == n) || (m == n)) continue;

          // Computation of <4> for a quadruplet of particles
          fourParticleCorrelation = cos(nOne*phi[k] + nTwo*phi[l] + nThree*phi[m] + nFour*phi[n]);
          totalParticleWeight = particleWeight[k]*particleWeight[l]*particleWeight[m]*particleWeight[n];

          // Filling of the TProfile
          fCorrelationWithNestedLoopsProfile->Fill(0.5, fourParticleCorrelation, totalParticleWeight);

        } // End of the loop over the fourth particle of the quadruplet
      } // End of the loop over the third particle of the quadruplet
    } // End of the loop over the second particle of the quadruplet
  } // End of the loop over the first particle of the quadruplet

// 3.) Reset of the local variables before changing the event
  fourParticleCorrelation = 0.;
  totalParticleWeight = 0.;

} // End of void AliAnalysisTaskTwoMultiCorrelations::ComputeCorrelationsWithFourNestedLoops(Int_t nOne, Int_t nTwo, Int_t nThree, Int_t nFour, Int_t nParticles, Double_t phi[], Double_t particleWeight[])

//******************************************************************************

void AliAnalysisTaskTwoMultiCorrelations::ComputeCorrelationsWithStandAloneQvectors(Int_t n, Int_t p, Int_t nParticles, Double_t phi[], Double_t particleWeight[])
{
// Method to compute the 2-p correlation with Q-vectors, stand alone formula (second cross-check for the recursion method
// Organisation of the method
  // 1.) Declaration of the local variables
  // 2.) Computation of the 2-p correlation
  // 3.) Filling of the TProfile
  // 4.) Security: reset to 0 of the variables

// 1.) Declaration of the local variables
  Double_t twoParticleCorrelation = 0.;
  Int_t twoParticleWeight = 0;

// 2.) Computation of the 2-p correlation
  twoParticleCorrelation = (1./(2.*TMath::Binomial(nParticles, 2))) * (CalculateQvector(n, p, nParticles, phi, particleWeight).Rho2() - nParticles);
  twoParticleWeight = nParticles * (nParticles - 1);

// 3.) Filling of the TProfile
  fCorrelationWithQvectorsSaProfile->Fill(0.5, twoParticleCorrelation, twoParticleWeight);

// 4.) Security: reset to 0 of the variables
  twoParticleCorrelation = 0.;
  twoParticleWeight = 0;

} // End of void AliAnalysisTaskTwoMultiCorrelations::ComputeCorrelationsWithStandAloneQvectors()
