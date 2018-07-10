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
#include <vector>
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

// ============================================================================================== //

AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights):
  AliAnalysisTaskSE(name),
  fNparticlesCorrelations(2), // Number of m-particle correlations and harmonics (2-14)
  fHarmonicOne(2),  // Harmonic n_1, default value for 2-p correlation
  fHarmonicTwo(-2), // Harmonic n_2, default value for 2-p correlation
  fHarmonicThree(0),  // Harmonic n_3
  fHarmonicFour(0), // Harmonic n_4
  fHarmonicFive(0), // Harmonic n_5
  fHarmonicSix(0),  // Harmonic n_6
  fHarmonicSeven(0),  // Harmonic n_7
  fHarmonicEight(0),  // Harmonic n_8
  fHarmonicNine(0), // Harmonic n_9
  fHarmonicTen(0),  // Harmonic n_10
  fHarmonicEleven(0), // Harmonic n_11
  fHarmonicTwelve(0), // Harmonic n_12
  fHarmonicThirteen(0), // Harmonic n_13
  fHarmonicFourteen(0), // Harmonic n_14
  fMinCentrality(0.0),  // Minimum of centrality
  fMaxCentrality(100.0),  // Maximum of centrality
  fUseParticleWeights(kFALSE),  // Use non-unit particle weights
  fDoNestedLoops(kFALSE), // Cross-check the results with nested loops
  fOutputList(NULL),  // Main list holding all the output objects
  fControlOutputList(NULL), // List holding all the control objects
  fDraftOutputList(NULL), // List holding all the intermediate objects
  fFinalOutputList(NULL), // List holding all the final results
  fPtControlHisto(NULL),  // Control histogram for the transverse momentum
  fEtaControlHisto(NULL), // Control histogram for the pseudorapidity
  fControlPhiHisto(NULL), // Control histogram for the azimuthal angles
  fCentralityHisto(NULL), // Control histogram for the centrality
  fMultiplicityDist(NULL),  // Control histogram for the multiplicity distribution
  fCorrelationWithQvectorsProfile(NULL),  // m-particle correlation estimated with Q-vectors
  fCorrelationWithNestedLoopsProfile(NULL), // 2-p correlation estimated with nested loops
  fCorrelationWithQvectorsSaProfile(NULL)  // 2-p correlation estimated with stand-alone Q-vectors
  //fEstimatedFlowWithQcProfile(NULL) // Anisotropic flow estimated with Q-cumulants
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

// ********************************************************************************************** //

AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations():
  AliAnalysisTaskSE(),
  fNparticlesCorrelations(2), // Number of m-particle correlations and harmonics (2-14)
  fHarmonicOne(2),  // Harmonic n_1
  fHarmonicTwo(-2), // Harmonic n_2
  fHarmonicThree(0),  // Harmonic n_3
  fHarmonicFour(0), // Harmonic n_4
  fHarmonicFive(0), // Harmonic n_5
  fHarmonicSix(0),  // Harmonic n_6
  fHarmonicSeven(0),  // Harmonic n_7
  fHarmonicEight(0),  // Harmonic n_8
  fHarmonicNine(0), // Harmonic n_9
  fHarmonicTen(0),  // Harmonic n_10
  fHarmonicEleven(0), // Harmonic n_11
  fHarmonicTwelve(0), // Harmonic n_12
  fHarmonicThirteen(0), // Harmonic n_13
  fHarmonicFourteen(0), // Harmonic n_14
  fMinCentrality(0.0),  // Minimum of centrality
  fMaxCentrality(100.0),  // Maximum of centrality
  fUseParticleWeights(kFALSE),  // Use non-unit particle weights
  fDoNestedLoops(kFALSE), // Cross-check the results with nested loops
  fOutputList(NULL),  // Main list holding all the output objects
  fControlOutputList(NULL), // List holding all the control objects
  fDraftOutputList(NULL), // List holding all the intermediate objects
  fFinalOutputList(NULL), // List holding all the final results
  fPtControlHisto(NULL),  // Control histogram for the transverse momentum
  fEtaControlHisto(NULL), // Control histogram for the pseudorapidity
  fControlPhiHisto(NULL), // Control histogram for the azimuthal angles
  fCentralityHisto(NULL), // Control histogram for the centrality
  fMultiplicityDist(NULL),  // Control histogram for the multiplicity distribution
  fCorrelationWithQvectorsProfile(NULL),  // m-particle correlation estimated with Q-vectors
  fCorrelationWithNestedLoopsProfile(NULL), // 2-p correlation estimated with nested loops
  fCorrelationWithQvectorsSaProfile(NULL)  // 2-p correlation estimated with stand-alone Q-vectors
  //fEstimatedFlowWithQcProfile(NULL) // Anisotropic flow estimated with Q-cumulants
{
// Dummy constructor of the class

    AliDebug(2, "AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights)");

} // End of the dummy constructor

// ********************************************************************************************** //

AliAnalysisTaskTwoMultiCorrelations::~AliAnalysisTaskTwoMultiCorrelations()
{
// Destructor of the class
  // Delete the main TList => delete automatically everything holds in it

  if(fOutputList) delete fOutputList;

} // End of the destructor

// ============================================================================================== //

void AliAnalysisTaskTwoMultiCorrelations::UserCreateOutputObjects()
{
// Method called at every worker node to initialise the lists
// Organisation of the method
  // First part of the trick to avoid name clashes
  // Booking and nesting of all the lists
  // Booking of all the objects
  // Second part of the trick to avoid name clashes

// First part of the trick to avoid name clashes
  Bool_t oldHistAddStatus = TH1::AddDirectoryStatus(); 
  TH1::AddDirectory(kFALSE);

// Booking and nesting of all the lists
  this->BookAndNestAllLists();

// Booking of all the objects
  this->BookControlList();
  this->BookDraftList();
  this->BookFinalList();

// Second part of the trick to avoid name clashes
  TH1::AddDirectory(oldHistAddStatus);

  PostData(1,fOutputList);

} // End of AliAnalysisTaskTwoMultiCorrelations::UserCreateOutputObjects()

// ********************************************************************************************** //

void AliAnalysisTaskTwoMultiCorrelations::UserExec(Option_t *)
{
// Method called for each event, contains all the calculations
  // Note to self: find a way to include non-unit particle weights and select them via fUseParticleWeights
// Organisation of the method
  // Obtention of the pointer to the AOD event
  // Gestion of the centrality
  // Loop over the multiplicity of the event
    // Obtention of the pointer to the current particle
    // Definition and control histograms of some variables for the cuts
    // Determination of the number of particles that pass the cuts
  // Definition of the varibles for the analysis
  // Do the analysis only if the number of particles after the cuts is higher than the one for the m-particle correlation
    // Gestion of the centrality
    // Loop over the number of remaining particles
      // Obtention of the pointer to the current particle
      // Filling of the azimuthal angles and particle weights and the control histograms
    // Computation of the correlations with different possible methods
    // Release of the allocated memory
  // PostData

//cout << "TEST" << endl;
//exit(0);


// Obtention of the pointer to the AOD event (from TaskSE)
  AliAODEvent *currentEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!currentEvent){return;}  // Protection against NULL pointer

// Gestion of the centrality
  AliMultSelection *ams = (AliMultSelection*)currentEvent->FindListObject("MultSelection");
  if(!ams){return;}

// Loop over the multiplicity of the event
  Int_t nParticles = currentEvent->GetNumberOfTracks(); // Number of particles of the event
  Int_t nParticlesAfterCuts = 0;  // Number of particles that pass the various cuts (pT, eta, ...)
  Int_t *goodIndices = new Int_t[nParticles](); // Vector where the index of a particle which passes the cut is saved with 1 and if it does not pass, the index corresponds to 0

  for (Int_t iParticle = 0; iParticle < nParticles; iParticle++)
  {
  // Obtention of the pointer to the current particle
    AliAODTrack *currentParticle = dynamic_cast<AliAODTrack*>(currentEvent->GetTrack(iParticle));
    if(!currentParticle){continue;} // Protection against NULL pointer
    if(!currentParticle->TestFilterBit(128)){continue;} // Filter bit 128 denotes TPC-only tracks

  // Definition and control histograms of some variables for the cuts
    Double_t pT = currentParticle->Pt();  // Transverse momentum
    Double_t eta = currentParticle->Eta();  // Pseudorapidity

    fPtControlHisto->Fill(pT);
    fEtaControlHisto->Fill(eta);

  // Determination of the number of particles which pass the cuts and set the flag 1/0 of goodIndices[iParticle]
    if ( (-0.8 < eta) && (eta < 0.8) && (0.2 < pT) && (pT < 5.0) )
    {
      nParticlesAfterCuts++;
      goodIndices[iParticle] = 1;
    }
    else {goodIndices[iParticle] = 0;}

  } // End of for (Int_t iParticle = 0; iParticle < nParticles; iParticle++)

// Filling of the control histogram for the multiplicity distribution
  fMultiplicityDist->Fill(nParticlesAfterCuts);

// Definition of the varibles for the analysis
  Double_t *phi = new Double_t[nParticlesAfterCuts]();  // Azimuthal angles
  Double_t *particleWeights = new Double_t[nParticlesAfterCuts]();  // Particle weights
  Int_t harmonics[14] = {fHarmonicOne, fHarmonicTwo, fHarmonicThree, fHarmonicFour, fHarmonicFive, fHarmonicSix, fHarmonicSeven, fHarmonicEight, fHarmonicNine, fHarmonicTen, fHarmonicEleven, fHarmonicTwelve, fHarmonicThirteen, fHarmonicFourteen};  // Harmonics (n_1,... n_14)
  Int_t index = 0;  // Index of the "good index", increased when the loop reaches a 1 in goodIndices

// Do the analysis only if the number of particles after the cuts is higher than the one for the m-particle correlation
  // (in order to avoid the division by zero) 
  if (nParticlesAfterCuts >= fNparticlesCorrelations)
  {
  // Gestion of the centrality
    if(ams->GetMultiplicityPercentile("V0M") >= fMinCentrality && ams->GetMultiplicityPercentile("V0M") < fMaxCentrality)
    {
      fCentralityHisto->Fill(ams->GetMultiplicityPercentile("V0M"));
    }
    else {return;} // this event does not belong to the centrality class specified for this particular analysis

  // Loop over the number of remaining particles
    for (Int_t iiParticle = 0; iiParticle < nParticles; iiParticle++)
    {
    // Obtention of the pointer to the current particle
      AliAODTrack *keptParticle = dynamic_cast<AliAODTrack*>(currentEvent->GetTrack(iiParticle));
      if(!keptParticle){continue;} // Protection against NULL pointer
      if(!keptParticle->TestFilterBit(128)){continue;} // Filter bit 128 denotes TPC-only tracks

    // Filling of the azimuthal angles and particle weights and the control histograms only if goodIndices[iiParticle] == 1
      if (goodIndices[iiParticle] == 1)
      {
        phi[index] = keptParticle->Phi();
        particleWeights[index] = 1.;
        //if (fUseParticleWeights) {continue;}  // Note to self: add the gestion of non-unit weight from external file
        //else {particleWeights[iiParticle] = 1.;}

        fControlPhiHisto->Fill(phi[index]);
        index++;
      }
      else {continue;} // End of if ( (-0.8 < eta) && (eta < 0.8) && (0.2 < pT) && (pT < 5.0) )
    } // End of for (Int_t iiParticle = 0; iiParticle < nParticlesAfterCuts; iiParticle++)

  // Computation of the correlations with different possible methods
    // Q-vectors with recursion formula
      for (Int_t iParticleCorr = 2; iParticleCorr <= fNparticlesCorrelations; iParticleCorr = iParticleCorr + 2)
      {
        //ComputeCorrelationsWithQvectors(nParticles, phi, particleWeights, harmonics, fNparticlesCorrelations);
        ComputeCorrelationsWithQvectors(nParticlesAfterCuts, phi, particleWeights, harmonics, iParticleCorr);
      } // End of for (Int_t iCorr = 2; iCorr <= fNparticlesCorrelations; iCorr = iCorr + 2)

    // Nested loops (for cross-checking the results from the Q-vectors)
      if ((fDoNestedLoops) && (fNparticlesCorrelations == 2)) {ComputeCorrelationsWithTwoNestedLoops(fHarmonicOne, fHarmonicTwo, nParticlesAfterCuts, phi, particleWeights);}
      if ((fDoNestedLoops) && (fNparticlesCorrelations == 4)) {ComputeCorrelationsWithFourNestedLoops(fHarmonicOne, fHarmonicTwo, fHarmonicThree, fHarmonicFour, nParticlesAfterCuts, phi, particleWeights);}

    // Q-vectors with stand-alone formula pour 2-p correlation (debug for recursion)
      if (fNparticlesCorrelations == 2) {ComputeCorrelationsWithStandAloneQvectors(fHarmonicOne, 0, nParticlesAfterCuts, phi, particleWeights);}

  // Release of the allocated memory
    delete [] goodIndices;
    delete [] phi;
    delete [] particleWeights;

  } // End of if (nParticlesAfterCuts >= fNparticlesCorrelations)

// PostData
  PostData(1,fOutputList);

} // End of void AliAnalysisTaskTwoMultiCorrelations::UserExec(Option_t *)

//******************************************************************************

void AliAnalysisTaskTwoMultiCorrelations::Terminate(Option_t *)
{
// Method called at the end of the execution, once the run over the events is over
// Organisation of the method
  // 1.) Access to the merged output list
  // 2.) Estimation of anisotropic flow with Q-cumulants
  // 3.) Creation of the output file and save of the main list in it

// 1.) Access to the merged output list
  fOutputList = (TList*)GetOutputData(1);
  if(!fOutputList){exit(1);}

// 2.) Estimation of anisotropic flow with Q-cumulants (CANNOT BE MERGED)
  //EstimateFlowWithCumulants(fNparticlesCorrelations);

// 3.) Creation of the output file and save of the main list in it
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
  // Control histogram for the transverse momentum
  // Control histogram for the pseudorapidity
  // 1.) Control histogram for the azimuthal angles

// Control histrogram for the transverse momentum
  fPtControlHisto = new TH1F("fPtControlHisto", "Transverse momentum distribution", 1000, 0., 20.);
  fPtControlHisto->SetStats(kTRUE);
  fPtControlHisto->GetXaxis()->SetTitle("p_{T}");
  fControlOutputList->Add(fPtControlHisto);

// Control histogram for the pseudorapidity
  fEtaControlHisto = new TH1F("fEtaControlHisto", "Pseudorapidity distribution", 1000, -1., 1.);
  fEtaControlHisto->SetStats(kTRUE);
  fEtaControlHisto->GetXaxis()->SetTitle("eta");
  fControlOutputList->Add(fEtaControlHisto);

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
  fMultiplicityDist = new TH1F("fMultiplicityDist", "Multiplicity distribution", 3000, 0., 3000.);
  fMultiplicityDist->GetXaxis()->SetTitle("Multiplicity");
  fControlOutputList->Add(fMultiplicityDist);

} // End of void AliAnalysisTaskTwoMultiCorrelations::BookControlList()

//******************************************************************************

void AliAnalysisTaskTwoMultiCorrelations::BookDraftList()
{
// Method to prepare the list with the intermediate results from the computation of Q-vectors and nested loops
// Organisation of the method
  // 1.) TProfile from the method of the Q-vectors
  // 2.) TProfile from the method of the two nested loops
  // 3.) TProfile from the method of the Q-vectors, stand alone formula

// 1.) TProfile from the method of the Q-vectors
  fCorrelationWithQvectorsProfile = new TProfile("fCorrelationWithQvectorsProfile", "m-particle correlation with Q-vectors", static_cast<int>(fNparticlesCorrelations)/2, 0., static_cast<double>(fNparticlesCorrelations)/2.);
  fCorrelationWithQvectorsProfile->SetStats(kFALSE);
  fCorrelationWithQvectorsProfile->Sumw2();
  //fCorrelationWithQvectorsProfile->GetXaxis()->SetTitle("m-particle correlation");
  fDraftOutputList->Add(fCorrelationWithQvectorsProfile);

// 2.) TProfile from the method of nested loops
  fCorrelationWithNestedLoopsProfile = new TProfile("fCorrelationWithNestedLoopsProfile", "m-particle correlation with nested loops", 1, 0., 1.);
  fCorrelationWithNestedLoopsProfile->SetStats(kFALSE);
  fCorrelationWithNestedLoopsProfile->Sumw2();
  //fCorrelationWithNestedLoopsProfile->GetXaxis()->SetTitle("m-particle correlation");
  fDraftOutputList->Add(fCorrelationWithNestedLoopsProfile);

// 3.) TProfile from the method of the Q-vectors, stand alone formula
  fCorrelationWithQvectorsSaProfile = new TProfile("fCorrelationWithQvectorsSaProfile", "m-particle correlation with Q-vectors, stand alone", 1, 0., 1.);
  fCorrelationWithQvectorsSaProfile->SetStats(kTRUE);
  fCorrelationWithQvectorsSaProfile->Sumw2();
  //fCorrelationWithQvectorsSaProfile->GetXaxis()->SetTitle("2-particle correlation");
  fDraftOutputList->Add(fCorrelationWithQvectorsSaProfile);

} // End of void AliAnalysisTaskTwoMultiCorrelations::BookDraftList()

//******************************************************************************

void AliAnalysisTaskTwoMultiCorrelations::BookFinalList()
{
// Method to prepare the list with the final results for the anisotropic flow
// Organisation of the method
  // 1.) TProfile of the anisotropic flow estimated with Q-cumulants
  // 2.) TProfile of the errors on the anisotropic flow estimated with Q-cumulants

// 1.) TProfile of the anisotropic flow estimated with Q-cumulants
  /*fEstimatedFlowWithQcProfile = new TProfile("fEstimatedFlowWithQcProfile", "Flow v_n estimated with Q-cumulants", 1, 0., 1.);
  fEstimatedFlowWithQcProfile->SetStats(kTRUE);
  //fEstimatedFlowWithQcProfile->Sumw2();
  fEstimatedFlowWithQcProfile->GetXaxis()->SetTitle("Estimated v_n");
  fFinalOutputList->Add(fEstimatedFlowWithQcProfile);*/

// 2.) TProfile of the errors on the anisotropic flow estimated with Q-cumulants
  /*fEstimatedFlowWithQcErrorsProfile = new TProfile("fEstimatedFlowWithQcErrorsProfile", "Error on v_n estimated with QC", 1, 0., 1.);
  fEstimatedFlowWithQcErrorsProfile->SetStats(kTRUE);
  fEstimatedFlowWithQcErrorsProfile->GetXaxis()->SetTitle("Error on estimated v_n");
  fFinalOutputList->Add(fEstimatedFlowWithQcErrorsProfile);*/

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
    qVectorNp += pWeightPowerP * TComplex::Exp((static_cast<double>(n))*(TComplex::I())*phi[iParticle]);
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
  std::vector<Int_t> numeratorHarmonics(nCorr);   // Harmonics with values from harmonics array
  std::vector<Int_t> denominatorHarmonics(nCorr); // Harmonics with value 0 for the denominator
  memset(numeratorHarmonics.data(), 0, sizeof(Int_t) * numeratorHarmonics.size());
  memset(denominatorHarmonics.data(), 0, sizeof(Int_t) * denominatorHarmonics.size());

  Double_t denominator = 0.;           // Denominator of <m>_(n1...nm) and event weight for the TProfile
  TComplex mParticleCorrelation = TComplex(0,0);  // m-particle correlation
  TComplex eventWeight = TComplex(0,0);// Event weight

  Double_t middleBin = 0.;             // Compute the middle value of the bin corresponding to <<m>> in the TProfile
 
// 2.) Filling of the array for the numerator
  for (Int_t iCorr = 0; iCorr < nCorr; iCorr++)
  {
    numeratorHarmonics[iCorr] = harmonics[iCorr];
    denominatorHarmonics[iCorr] = 0;
  } // End of for (Int_t iCorr = 0; iCorr < nCorr; iCorr++)

// 3.) Computation of the m-particle correlation
  denominator = (CalculateRecursionWithQvectors(nParticles, phi, particleWeight, nCorr, denominatorHarmonics.data())).Re();
  cout << Form("Denominator: %f", denominator) << endl;
  mParticleCorrelation = CalculateRecursionWithQvectors(nParticles, phi, particleWeight, nCorr, numeratorHarmonics.data())/denominator;
  //eventWeight = (CalculateRecursionWithQvectors(nParticles, phi, particleWeight, nCorr, denominatorHarmonics)).Re();

// 4.) Filling of the TProfile
  middleBin = (static_cast<double>(nCorr) - 1.)/2.;
  fCorrelationWithQvectorsProfile->Fill(middleBin, mParticleCorrelation.Re(), denominator);

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

// 4.) Security: reset of the variables to 0.
  twoParticleCorrelation = 0.;
  twoParticleWeight = 0;

} // End of void AliAnalysisTaskTwoMultiCorrelations::ComputeCorrelationsWithStandAloneQvectors()

//******************************************************************************
/*
void AliAnalysisTaskTwoMultiCorrelations::EstimateFlowWithCumulants(Int_t maxMcorr)
{
// Offline method to compute the estimation of v_n{np} using the m-p correlations
// obtained with the Q-vectors
// Organisation of the method
  // 1.) Declaration of the local variables
  // 2.) Filling of the values of the correlations from the TProfile
  // *.) Computation of the Q-cumulant and the anisotropic flow for all possible cases of maxMcorr
  // *.) Security: reset of the variables to 0.

// 1.) Declaration of the local variables
  Double_t cNforMaxMcorr = 0.;         // Q-cumulant c_n{MaxMcorr}
  Double_t vNforMaxMcorr = 0.;         // Estimated anisotropic flow v_n{MaxMcorr}

  std::vector<Double_t> mParticleCorrelation(static_cast<int>(maxMcorr/2)); // Vector with the m-particle correlations: <<2>>, <<4>>, <<6>>, <<8>>
  memset(mParticleCorrelation.data(), 0, sizeof(Double_t) * mParticleCorrelation.size());

// 2.) Filling of the values of the correlations from the TProfile
  for (Int_t i = 0; i < static_cast<int>(maxMcorr/2); i++)
  {
    mParticleCorrelation[i] = fCorrelationWithQvectorsProfile->GetBinContent(i+1);
  } // End of for (Int_t i = 0; i < static_cast<int>(maxMcorr/2); i++)

// *.) Computation of the Q-cumulant and the anisotropic flow for all possible cases of maxMcorr
  switch(maxMcorr)
  {
    case 2 :  // Computation of c_n{2} and v_n{2}
      cNforMaxMcorr = mParticleCorrelation[0];
      vNforMaxMcorr = pow(cNforMaxMcorr, 0.5);
      break;
    case 4 :
      cNforMaxMcorr = mParticleCorrelation[1] - 2.*pow(mParticleCorrelation[0], 2.);
      vNforMaxMcorr = pow((-1.*cNforMaxMcorr), 0.25);
      break;
    case 6 :
      cNforMaxMcorr = mParticleCorrelation[2] - 9.*mParticleCorrelation[0]*mParticleCorrelation[1] + 12.*pow(mParticleCorrelation[0], 3.);
      vNforMaxMcorr = pow((cNforMaxMcorr/4.), (1./6.));
      break;
    case 8 :
      cNforMaxMcorr = mParticleCorrelation[3] - 16.*mParticleCorrelation[0]*mParticleCorrelation[2] - 18.*pow(mParticleCorrelation[1], 2.) + 144.*mParticleCorrelation[1]*pow(mParticleCorrelation[2], 2.) - 144.*pow(mParticleCorrelation[0], 4.);
      vNforMaxMcorr = pow((-1.*cNforMaxMcorr/33.), 0.125);
      break;
  } // End of switch

// *.) Filling of the TProfile
  fEstimatedFlowWithQcProfile->Fill(0.5, vNforMaxMcorr);
  cout << Form("cN: %f", cNforMaxMcorr) << endl;

// *.) Security_ reset of the variables to 0.
  cNforMaxMcorr = 0.;
  vNforMaxMcorr = 0.;
  for (Int_t i = 0; i < static_cast<int>(maxMcorr/2); i++)
  {
    mParticleCorrelation[i] = 0.;
  } // End of for (Int_t i = 0; i < static_cast<int>(maxMcorr/2); i++)

  //return vNforMaxMcorr;

}*/ // End of Double_t AliAnalysisTaskTwoMultiCorrelations::ComputeCorrelationsWithStandAloneQvectors()
