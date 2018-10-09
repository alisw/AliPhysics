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

//--------------------------------------------------------------------------------------//
// Analysis task for the data analysis of the correlations between the anisotropic flow //
// harmonics v_n with the Pb-Pb data taken by the ALICE experiment.                     //
// The current script computes the multiparticle correlators using the method of the    //
// Q-vectors for a maximum of 6 different harmonics and 8 particles).                   //
//                                                                                      //
// Author: Cindy Mordasini (cindy.mordasini@cern.ch)                                    //
//--------------------------------------------------------------------------------------//

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

//======================================================================================//

AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights):
  AliAnalysisTaskSE(name),
  fMaxNumberCorrelations(8), // General parameters
  fHighHarmonic(6),
  fNumberDifferentHarmonics(2),
  fUseParticleWeights(kFALSE),
  fComputeNestedLoops(kFALSE),
  fComputeSine(kFALSE),
  fHarmonicOne(2),  // Harmonics
  fHarmonicTwo(-2),
  fHarmonicThree(4),
  fHarmonicFour(-4),
  fHarmonicFive(0),
  fHarmonicSix(0),
  fHarmonicSeven(0),
  fHarmonicEight(0),
  fMinCentrality(0.0),  // Ranges
  fMaxCentrality(100.0),
  fMinPtCut(0.2),
  fMaxPtCut(5.0),
  fMinEtaCut(-0.8),
  fMaxEtaCut(0.8),
  fOutputMainList(NULL),  // Structure of the output file
  fPreCutControlList(NULL),
  fPostCutControlList(NULL),
  fCorrelationResultsList(NULL),
  fCentralityPreCutHisto(NULL), // Control histograms before the cuts
  fMultiplicityPreCutHisto(NULL),
  fPtPreCutControlHisto(NULL),
  fEtaPreCutControlHisto(NULL),
  fPhiPreCutHisto(NULL),
  fMultiplicityPostCutHisto(NULL),  // Control histograms after the cuts
  fPtPostCutControlHisto(NULL),
  fEtaPostCutControlHisto(NULL),
  fPhiPostCutHisto(NULL)
{
// Constructor of the class
  AliDebug(2, "AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights)");

// Creation of a new main list
  fOutputMainList = new TList();
  fOutputMainList->SetName("outputAnalysis");
  fOutputMainList->SetOwner(kTRUE);

// Definition of the input and output slots
  DefineOutput(1, TList::Class());

// Initialisation of the fQvectors to zero
  InitialiseArraysOfQvectors();

// Initialisation of the pointers for the TProfiles to NULL
  InitialiseArraysOfTProfiles();

// Use of the particle weights
  if(useParticleWeights)
  {
    // not needed for the time being, maybe insert here the call for the file with particle weights???
  }
} // End of the constructor

//======================================================================================//

AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations():
  AliAnalysisTaskSE(),
  fMaxNumberCorrelations(8), // General parameters
  fHighHarmonic(6),
  fNumberDifferentHarmonics(2),
  fUseParticleWeights(kFALSE),
  fComputeNestedLoops(kFALSE),
  fComputeSine(kFALSE),
  fHarmonicOne(2),  // Harmonics
  fHarmonicTwo(-2),
  fHarmonicThree(4),
  fHarmonicFour(-4),
  fHarmonicFive(0),
  fHarmonicSix(0),
  fHarmonicSeven(0),
  fHarmonicEight(0),
  fMinCentrality(0.0),  // Ranges
  fMaxCentrality(100.0),
  fMinPtCut(0.2),
  fMaxPtCut(5.0),
  fMinEtaCut(-0.8),
  fMaxEtaCut(0.8),
  fOutputMainList(NULL),  // Structure of the output file
  fPreCutControlList(NULL),
  fPostCutControlList(NULL),
  fCorrelationResultsList(NULL),
  fCentralityPreCutHisto(NULL), // Control histograms before the cuts
  fMultiplicityPreCutHisto(NULL),
  fPtPreCutControlHisto(NULL),
  fEtaPreCutControlHisto(NULL),
  fPhiPreCutHisto(NULL),
  fMultiplicityPostCutHisto(NULL),  // Control histograms after the cuts
  fPtPostCutControlHisto(NULL),
  fEtaPostCutControlHisto(NULL),
  fPhiPostCutHisto(NULL)
  //fEstimatedFlowWithQcProfile(NULL) // Anisotropic flow estimated with Q-cumulants
{
// Dummy constructor of the class

    AliDebug(2, "AliAnalysisTaskTwoMultiCorrelations::AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights)");

// Initialisation of the ffQvectors to zero
  InitialiseArraysOfQvectors();

// Initialisation of the pointers for the TProfiles to NULL
  InitialiseArraysOfTProfiles();
} // End of the dummy constructor

//======================================================================================//

AliAnalysisTaskTwoMultiCorrelations::~AliAnalysisTaskTwoMultiCorrelations()
{
// Destructor of the class
  // Delete the main TList => delete automatically everything holds in it

  if(fOutputMainList) {delete fOutputMainList;}
} // End of the destructor

//======================================================================================//

void AliAnalysisTaskTwoMultiCorrelations::UserCreateOutputObjects()
{
// Method called at every worker node to initialise the lists

// First part of the trick to avoid name clashes
  Bool_t oldHistAddStatus = TH1::AddDirectoryStatus(); 
  TH1::AddDirectory(kFALSE);

// Booking and nesting of all the lists
  this->BookAndNestAllLists();

// Booking of all the secondary lists
  this->BookPreCutControlList();
  this->BookPostCutControlList();
  this->BookCorrelationResultsList();

// Second part of the trick to avoid name clashes
  TH1::AddDirectory(oldHistAddStatus);

  PostData(1,fOutputMainList);
} // End of AliAnalysisTaskTwoMultiCorrelations::UserCreateOutputObjects()

//======================================================================================//

void AliAnalysisTaskTwoMultiCorrelations::UserExec(Option_t *)
{
// Method called for each event, contains all the calculations
/// Note to self: find a way to include non-unit particle weights and select them via fUseParticleWeights

// Obtention of the pointer to the AOD event (from TaskSE)
  AliAODEvent *currentEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!currentEvent){return;}  // Protection against NULL pointer

// Gestion of the centrality
  AliMultSelection *ams = (AliMultSelection*)currentEvent->FindListObject("MultSelection");
  if(!ams){return;}
  if(ams->GetMultiplicityPercentile("V0M") >= fMinCentrality && ams->GetMultiplicityPercentile("V0M") < fMaxCentrality)
  {
    fCentralityPreCutHisto->Fill(ams->GetMultiplicityPercentile("V0M"));
  }
  else {return;} // This event does not belong to the centrality class specified for this particular analysis

// Loop over the multiplicity of the event
  Int_t nParticles = currentEvent->GetNumberOfTracks(); // Number of particles in the event
  Int_t nParticlesAfterCuts = 0;  // Number of particles remaining after the cuts
  Int_t *goodIndices = new Int_t[nParticles](); // Array where the index of a particle which passes the cuts is saved with 1 and if not with 0

  // Filling of the control histogram for the multiplicity before the cuts
  fMultiplicityPreCutHisto->Fill(nParticles);

  for (Int_t iParticle = 0; iParticle < nParticles; iParticle++)
  {
  // Obtention of the pointer to the current particle
    AliAODTrack *currentParticle = dynamic_cast<AliAODTrack*>(currentEvent->GetTrack(iParticle));
    if(!currentParticle){continue;} // Protection against NULL pointer
    if(!currentParticle->TestFilterBit(128)){continue;} // Filter bit 128 denotes TPC-only tracks

  // Filling of some control histograms before the application of the cuts
    Double_t preCutPt = currentParticle->Pt();  // Transverse momentum
    Double_t preCutEta = currentParticle->Eta();  // Pseudorapidity
    Double_t preCutPhi = currentParticle->Phi(); // Azimuthal angles

    fPtPreCutControlHisto->Fill(preCutPt);
    fEtaPreCutControlHisto->Fill(preCutEta);
    fPhiPreCutHisto->Fill(preCutPhi);

  // Determination of the number of particles which pass the cuts and set the flag 1/0 of goodIndices[iParticle]
    if ( (fMinEtaCut < preCutEta) && (preCutEta < fMaxEtaCut) && (fMinPtCut < preCutPt) && (preCutPt < fMaxPtCut) )
    {
      nParticlesAfterCuts++;
      goodIndices[iParticle] = 1;
    }
    else {goodIndices[iParticle] = 0;}

  } // End: for (Int_t iParticle = 0; iParticle < nParticles; iParticle++)

// Filling of the control histogram for the multiplicity after the cuts
  fMultiplicityPostCutHisto->Fill(nParticlesAfterCuts);

// Definition of the varibles for the analysis
  Double_t *pt = new Double_t[nParticlesAfterCuts](); // Array of transverse momenta
  Double_t *eta = new Double_t[nParticlesAfterCuts](); // Array of pseudorapidity
  Double_t *phi = new Double_t[nParticlesAfterCuts]();  // Array of azimuthal angles
  Double_t *particleWeights = new Double_t[nParticlesAfterCuts]();  // Array of üarticle weights
  Int_t fullHarmonicsArray[8] = {fHarmonicOne, fHarmonicTwo, fHarmonicThree, fHarmonicFour, fHarmonicFive, fHarmonicSix, fHarmonicSeven, fHarmonicEight};  // Harmonics (n_1,... n_8)
  Int_t *harmonics = new Int_t[2*fNumberDifferentHarmonics]();
  Int_t index = 0;  // Index of the "good index", increased when the loop reaches a 1 in goodIndices

// Filling of the harmonics' array
  for (Int_t iArray = 0; iArray < 2*fNumberDifferentHarmonics; iArray++)
  {
    harmonics[iArray] = fullHarmonicsArray[iArray];
  } // End: for (Int_t iArray = 0; iArray < 2*fNumberDifferentHarmonics; iArray++)

// Do the analysis only if the number of particles after the cuts is higher than the one for the m-particle correlation (in order to avoid the division by zero) 
  if (nParticlesAfterCuts > (2*fNumberDifferentHarmonics))
  {
  // Loop over the initial number of particles, index takes care of keeping only the ones which pass the cuts
    for (Int_t iiParticle = 0; iiParticle < nParticles; iiParticle++)
    {
    // Obtention of the pointer to the current particle
      AliAODTrack *keptParticle = dynamic_cast<AliAODTrack*>(currentEvent->GetTrack(iiParticle));
      if(!keptParticle){continue;} // Protection against NULL pointer
      if(!keptParticle->TestFilterBit(128)){continue;} // Filter bit 128 denotes TPC-only tracks

    // Filling of the azimuthal angles and particle weights and the control histograms only if goodIndices[iiParticle] == 1
      if (goodIndices[iiParticle] == 1)
      {
        pt[index] = keptParticle->Pt();
        eta[index] = keptParticle->Eta(); 
        phi[index] = keptParticle->Phi();
        particleWeights[index] = 1.;
        //if (fUseParticleWeights) {continue;}  // Note to self: add the gestion of non-unit weight from external file
        //else {particleWeights[iiParticle] = 1.;}

        fPtPostCutControlHisto->Fill(pt[index]);
        fEtaPostCutControlHisto->Fill(eta[index]);
        fPhiPostCutHisto->Fill(phi[index]);

        index++;
      }
      else {continue;} // End of if ( (-0.8 < eta) && (eta < 0.8) && (0.2 < pT) && (pT < 5.0) )
    } // End: for (Int_t iiParticle = 0; iiParticle < nParticlesAfterCuts; iiParticle++)

  // Computation of the Q-vectors for the current event
    CalculateQvectors(nParticlesAfterCuts, phi, particleWeights);

  // Computation of the multiparticle correlations
  /// 2-particle correlation
    for (Int_t iDiffHarmo = 0; iDiffHarmo < fNumberDifferentHarmonics; iDiffHarmo++)
    {
      Int_t numeratorFirstTwoParticle[2] = {harmonics[2*iDiffHarmo],harmonics[(2*iDiffHarmo)+1]};  // iDiffHarmo-th pair of harmonics for the numerator
      Int_t numeratorLastTwoParticle[2] = {harmonics[2],harmonics[3]};  // the alternative expression is computed only for 2 different harmonics, which implies harmonics has a length of 4
      ComputeTwoParticleCorrelationWithQvectors(numeratorFirstTwoParticle, numeratorLastTwoParticle, iDiffHarmo);

      if (fComputeNestedLoops) {ComputeCorrelationsWithTwoNestedLoops(numeratorFirstTwoParticle, numeratorLastTwoParticle, iDiffHarmo, nParticles, phi, particleWeights);}
    } // End: for (Int_t iDiffHarmo = 0; iDiffHarmo < fNumberDifferentHarmonics; iDiffHarmo++)

  /// 4-particle correlation
    if (fNumberDifferentHarmonics >= 2)
    {
      Int_t iCombi = 0; // i-th combination of 4 ordonnated elements without repetition out of fNumberDifferentHarmonics*2
      for (Int_t i = 0; i < fNumberDifferentHarmonics-1; i++)
      {
      // Loop for the first two particles' harmonics
        for (Int_t j = i+1; j < fNumberDifferentHarmonics; j++)
        {
          Int_t numeratorFourParticle[4] = {harmonics[2*i],harmonics[(2*i)+1], harmonics[2*j],harmonics[(2*j)+1]};
          ComputeFourParticleCorrelationWithQvectors(numeratorFourParticle, iCombi);

          if (fComputeNestedLoops) {ComputeCorrelationsWithFourNestedLoops(numeratorFourParticle, iCombi, nParticles, phi, particleWeights);}
          iCombi++;
        } // End: for (Int_t j = i+1; j <= fkNumberDifferentHarmonics; j++)
      } // End: for (Int_t i = 0; i < fkNumberDifferentHarmonics, i++)
    } // End: if (fNumberDifferentHarmonics >= 2)

  /// 6-particle correlation
    if (fNumberDifferentHarmonics >= 3)
    {
      Int_t jCombi = 0; // (i-1)-th combination of 6 ordonnated elements without repetition out of fkNumberDifferentHarmonics*2
      for (Int_t i = 0; i < fNumberDifferentHarmonics-2; i++)
      {
      // Loop over the first two particles' harmonics
        for (Int_t j = i+1; j < fNumberDifferentHarmonics-1; j++)
        {
        // Loop over the second particle pair's harmonics
          for (Int_t k = j+1; k < fNumberDifferentHarmonics; k++)
          {
          // Loop over the third particle pair's harmonics
            Int_t numeratorSixParticle[6] = {harmonics[2*i],harmonics[(2*i)+1], harmonics[2*j],harmonics[(2*j)+1],harmonics[2*k],harmonics[(2*k)+1]};
            ComputeSixParticleCorrelationWithQvectors(numeratorSixParticle, jCombi);
            jCombi++;
          } // End: for (Int_t k = j+1; k < fNumberDifferentHarmonics; k++)
        } // End: for (Int_t j = i+1; j < fNumberDifferentHarmonics-1; j++)
      } // End: for (Int_t i = 0; i < fNumberDifferentHarmonics-2; i++)
    } // End: if (fNumberDifferentHarmonics >= 3)

  /// 8-particle correlation
    if (fNumberDifferentHarmonics >= 4)
    {
      ComputeEightParticleCorrelationWithQvectors(harmonics);
    } // End: if (fNumberDifferentHarmonics >= 4)   

  // Release of the allocated memory
    delete [] harmonics;
    delete [] goodIndices;
    delete [] pt;
    delete [] eta;
    delete [] phi;
    delete [] particleWeights;

  } // End of if (nParticlesAfterCuts >= fNparticlesCorrelations)

// PostData
  PostData(1,fOutputMainList);
} // End of void AliAnalysisTaskTwoMultiCorrelations::UserExec(Option_t *)

//======================================================================================//

void AliAnalysisTaskTwoMultiCorrelations::Terminate(Option_t *)
{
// Method called at the end of the execution, once the run over the events is over
// Organisation of the method
  // 1.) Access to the merged output list
  // 2.) Estimation of anisotropic flow with Q-cumulants
  // 3.) Creation of the output file and save of the main list in it

// 1.) Access to the merged output list
  fOutputMainList = (TList*)GetOutputData(1);
  if(!fOutputMainList){exit(1);}

// 3.) Creation of the output file and save of the main list in it
  TFile *outputFile = new TFile("AnalysisResults.root", "RECREATE");
  fOutputMainList->Write(fOutputMainList->GetName(),TObject::kSingleKey);
  delete outputFile;
} // End of void AliAnalysisTaskTwoMultiCorrelations::Terminate(Option_t *)

//======================================================================================//

void AliAnalysisTaskTwoMultiCorrelations::InitialiseArraysOfQvectors()
{
// Method to initialise the Q-vectors to zero
  for (Int_t iHarmo = 0; iHarmo < (fHighHarmonic*fMaxNumberCorrelations)+1; iHarmo++)
  {
    for (Int_t iPower = 0; iPower < fMaxNumberCorrelations+1; iPower++)
    {
      fQvectors[iHarmo][iPower] = TComplex(0.,0.);
    } // End: for (Int_t iPower = 0; iPower < fMaxNumberCorrelations+1; iPower++)
  } // End: for (Int_t iHarmo = 0; iHarmo < maxHarmo; iHarmo++)

} // End: void InitialiseArraysOffQvectors() 

//======================================================================================//

void AliAnalysisTaskTwoMultiCorrelations::InitialiseArraysOfTProfiles()
{
// Method to initialise the 1-d and 2-d arrays of TProfiles containing the results
  for (Int_t cs = 0; cs < 2; cs++)
  {
    for (Int_t t = 0; t < 4; t++) // 2-p and 6-p correlations
    {
      fTwoParticleCorrelationProfile[cs][t] = NULL;
      fSixParticleCorrelationProfile[cs][t] = NULL;
      fTwoCosineAverageProfile[cs][t] = NULL;
      fTwoNestedCorrelationProfile[cs][t] = NULL;  
      fTwoCosineAverageNestedProfile[cs][t] = NULL;
    } // End: for (Int_t t = 0; t < 4; t++)

    for (Int_t f = 0; f < 6; f++) // 4-p correlations
    {
      fFourParticleCorrelationProfile[cs][f] = NULL;
      fFourNestedCorrelationProfile[cs][f] = NULL;
    } // End: for (Int_t f = 0; f < 6; f++)

    fEightParticleCorrelationProfile[cs] = NULL;  // 8-p correlations

  } // End: for (Int_t cs = 0; cs < endCs; cs++)
} // End: void InitialiseArraysOfTProfiles()

//======================================================================================//

void AliAnalysisTaskTwoMultiCorrelations::BookAndNestAllLists()
{
// Method to book and nest all the lists where the results are kept in the output file

  TString sMethodName = "void AliAnalysisTaskTwoMultiCorrelations::BookAndNestAllLists()";
  if(!fOutputMainList){Fatal(sMethodName.Data(),"Main list fOutputMainList is NULL");}

// Control histograms before the cuts
  fPreCutControlList = new TList();
  fPreCutControlList->SetName("PreCutControlList");
  fPreCutControlList->SetOwner(kTRUE);
  fOutputMainList->Add(fPreCutControlList);

// Control histograms after the cuts
  fPostCutControlList = new TList();
  fPostCutControlList->SetName("PostCutControlList");
  fPostCutControlList->SetOwner(kTRUE);
  fOutputMainList->Add(fPostCutControlList);

// Results for the multiparticle correlations
  fCorrelationResultsList = new TList();
  fCorrelationResultsList->SetName("CorrelationResultsList");
  fCorrelationResultsList->SetOwner(kTRUE);
  fOutputMainList->Add(fCorrelationResultsList);
} // End: AliAnalysisTaskTwoMultiCorrelations::BookAndNestAllLists()

//======================================================================================//

void AliAnalysisTaskTwoMultiCorrelations::BookPreCutControlList()
{
// Method to prepare the list with the control histograms before the application of the cuts
/// Centrality distribution
/// Multiplicity distribution
/// Transverse momentum distribution
/// Pseudorapidity distribution
/// Azimuthal angles

// Centrality distribution
  fCentralityPreCutHisto = new TH1D("fCentralityPreCutHisto", "Centrality distribution before the cuts", 100, 0., 100.);
  fCentralityPreCutHisto->SetStats(kTRUE);
  fCentralityPreCutHisto->GetXaxis()->SetTitle("Centrality percentile");
  fPreCutControlList->Add(fCentralityPreCutHisto);

// Multiplicity distribution
  fMultiplicityPreCutHisto = new TH1D("fMultiplicityPreCutHisto", "Multiplicity distribution before the cuts", 10000000, 0., 10000000.);
  fMultiplicityPreCutHisto->SetStats(kTRUE);
  fMultiplicityPreCutHisto->GetXaxis()->SetTitle("Multiplicity");
  fPreCutControlList->Add(fMultiplicityPreCutHisto);

// Transverse momentum distribution
  fPtPreCutControlHisto = new TH1D("fPtPreCutControlHisto", "Transverse momentum distribution before the cuts", 1000, 0., 20.);
  fPtPreCutControlHisto->SetStats(kTRUE);
  fPtPreCutControlHisto->GetXaxis()->SetTitle("p_{T}");
  fPreCutControlList->Add(fPtPreCutControlHisto);

// Pseudorapidity distribution
  fEtaPreCutControlHisto = new TH1D("fEtaPreCutControlHisto", "Pseudorapidity distribution before the cuts", 1000, -1., 1.);
  fEtaPreCutControlHisto->SetStats(kTRUE);
  fEtaPreCutControlHisto->GetXaxis()->SetTitle("eta");
  fPreCutControlList->Add(fEtaPreCutControlHisto);

// Azimuthal angles
  fPhiPreCutHisto = new TH1D("fPhiPreCutHisto", "Azimuthal angles distribution before the cuts", 1000, 0., 6.3);
  fPhiPreCutHisto->SetStats(kTRUE);
  fPhiPreCutHisto->GetXaxis()->SetTitle("phi");
  fPreCutControlList->Add(fPhiPreCutHisto);
} // End: void AliAnalysisTaskTwoMultiCorrelations::BookPreCutControlList()

//======================================================================================//

void AliAnalysisTaskTwoMultiCorrelations::BookPostCutControlList()
{
// Method to prepare the list with the control histograms after the application of the cuts
/// Multiplicity distribution
/// Transverse momentum distribution
/// Pseudorapidity distribution
/// Azimuthal angles

// Multiplicity distribution
  fMultiplicityPostCutHisto = new TH1D("fMultiplicityPostCutHisto", "Multiplicity distribution after the cuts", 10000000, 0., 10000000.);
  fMultiplicityPostCutHisto->SetStats(kTRUE);
  fMultiplicityPostCutHisto->GetXaxis()->SetTitle("Multiplicity");
  fPostCutControlList->Add(fMultiplicityPostCutHisto);

// Transverse momentum distribution
  fPtPostCutControlHisto = new TH1D("fPtPostCutControlHisto", "Transverse momentum distribution after the cuts", 1000, 0., 20.);
  fPtPostCutControlHisto->SetStats(kTRUE);
  fPtPostCutControlHisto->GetXaxis()->SetTitle("p_{T}");
  fPostCutControlList->Add(fPtPostCutControlHisto);

// Pseudorapidity distribution
  fEtaPostCutControlHisto = new TH1D("fEtaPostCutControlHisto", "Pseudorapidity distribution after the cuts", 1000, -1., 1.);
  fEtaPostCutControlHisto->SetStats(kTRUE);
  fEtaPostCutControlHisto->GetXaxis()->SetTitle("eta");
  fPostCutControlList->Add(fEtaPostCutControlHisto);

// Azimuthal angles
  fPhiPostCutHisto = new TH1D("fPhiPostCutHisto", "Azimuthal angles distribution after the cuts", 1000, 0., 6.3);
  fPhiPostCutHisto->SetStats(kTRUE);
  fPhiPostCutHisto->GetXaxis()->SetTitle("phi");
  fPostCutControlList->Add(fPhiPostCutHisto);
} // End: void AliAnalysisTaskTwoMultiCorrelations::BookPreCutControlList()

//======================================================================================//

void AliAnalysisTaskTwoMultiCorrelations::BookCorrelationResultsList()
{
// Method to prepare the list with the results from the computation of Q-vectors and nested loops
/// TProfiles for 2-p and 6-p correlations
/// TProfiles for 4-p correlations
/// TProfiles for 8-p correlations

  Int_t endCs = 1;  // Number of components to save (1: only Cosine, 2: Cosine and Sine)
  if (fComputeSine) {endCs = 2;}

  for (Int_t cs = 0; cs < endCs; cs++)
  {
    for (Int_t t = 0; t < 4; t++)
    {
    // TProfiles for 2-p and 6-p correlations
    // cs: {cos, sin}, t: {m,n,p,q} for 2-p and t: {mnp, mnq, mpq, npq} for 6-p
      fTwoParticleCorrelationProfile[cs][t] = new TProfile("", "", 1, 0., 1.);
      fTwoParticleCorrelationProfile[cs][t]->SetStats(kTRUE);
      fTwoParticleCorrelationProfile[cs][t]->SetName(Form("fTwoParticleCorrelationProfile_cs%d_combi%d", cs, t));
      fTwoParticleCorrelationProfile[cs][t]->SetTitle(Form("Q-vectors, #LT #LT 2 #GT #GT, cs = %d, combi = %d", cs, t));
      fTwoParticleCorrelationProfile[cs][t]->Sumw2();
      fTwoParticleCorrelationProfile[cs][t]->GetXaxis()->SetTitle("#LT #LT 2 #GT #GT");
      fCorrelationResultsList->Add(fTwoParticleCorrelationProfile[cs][t]);

      fSixParticleCorrelationProfile[cs][t] = new TProfile("", "", 1, 0., 1.);
      fSixParticleCorrelationProfile[cs][t]->SetStats(kTRUE);
      fSixParticleCorrelationProfile[cs][t]->SetName(Form("fSixParticleCorrelationProfile_cs%d_combi%d", cs, t));
      fSixParticleCorrelationProfile[cs][t]->SetTitle(Form("Q-vectors, #LT #LT 6 #GT #GT, cs = %d, combi = %d", cs, t));
      fSixParticleCorrelationProfile[cs][t]->Sumw2();
      fSixParticleCorrelationProfile[cs][t]->GetXaxis()->SetTitle("#LT #LT 6 #GT #GT");
      fCorrelationResultsList->Add(fSixParticleCorrelationProfile[cs][t]);

      fTwoCosineAverageProfile[cs][t] = new TProfile("", "", 2, 0., 2.);
      fTwoCosineAverageProfile[cs][t]->SetStats(kTRUE);
      fTwoCosineAverageProfile[cs][t]->SetName(Form("fTwoCosineAverageProfile_cs%d_combi%d", cs, t));
      fTwoCosineAverageProfile[cs][t]->SetTitle(Form("Q-vectors, #LT #LT cos(m (phi1-phi2)) #GT #LT cos(n (phi1-phi2)) #GT #GT, cs = %d, combi = %d", cs, t)); // only [cos,sin][0] is useful)
      fTwoCosineAverageProfile[cs][t]->Sumw2();
      fTwoCosineAverageProfile[cs][t]->GetXaxis()->SetTitle("#LT #LT cos(m (phi1-phi2)) #GT #LT cos(n (phi1-phi2)) #GT #GT");
      fCorrelationResultsList->Add(fTwoCosineAverageProfile[cs][t]);

      if (fComputeNestedLoops)
      {
        fTwoNestedCorrelationProfile[cs][t] = new TProfile("", "", 1, 0., 1.);
        fTwoNestedCorrelationProfile[cs][t]->SetStats(kTRUE);
        fTwoNestedCorrelationProfile[cs][t]->SetName(Form("fTwoNestedCorrelationProfile_cs%d_combi%d", cs, t));
        fTwoNestedCorrelationProfile[cs][t]->SetTitle(Form("Nested loops, #LT #LT 2 #GT #GT, cs = %d, combi = %d", cs, t));
        fTwoNestedCorrelationProfile[cs][t]->Sumw2();
        fTwoNestedCorrelationProfile[cs][t]->GetXaxis()->SetTitle("#LT #LT 2 #GT #GT");
        fCorrelationResultsList->Add(fTwoNestedCorrelationProfile[cs][t]);

        fTwoCosineAverageNestedProfile[cs][t] = new TProfile("", "", 1, 0., 1.);
        fTwoCosineAverageNestedProfile[cs][t]->SetStats(kTRUE);
        fTwoCosineAverageNestedProfile[cs][t]->SetName(Form("fTwoCosineAverageNestedProfile_cs%d_combi%d", cs, t));
        fTwoCosineAverageNestedProfile[cs][t]->SetTitle(Form("Nested loops, #LT #LT cos(m (phi1-phi2)) #GT #LT cos(n (phi1-phi2)) #GT #GT, cs = %d, combi = %d", cs, t));
        fTwoCosineAverageNestedProfile[cs][t]->Sumw2();
        fTwoCosineAverageNestedProfile[cs][t]->GetXaxis()->SetTitle("#LT #LT cos(m (phi1-phi2)) #GT #LT cos(n (phi1-phi2)) #GT #GT");
        fCorrelationResultsList->Add(fTwoCosineAverageNestedProfile[cs][t]);
      } // End: if (fComputeNestedLoops)
    } // End: for (Int_t t = 0; t < 4; t++)

    for (Int_t t = 0; t < 6; t++)
    {
    // TProfiles for 4-p correlations
    // cs: {cos, sin}, f: {mn,mp,mq,np,nq,pq}
      fFourParticleCorrelationProfile[cs][t] = new TProfile("", "", 1, 0., 1.);
      fFourParticleCorrelationProfile[cs][t]->SetStats(kTRUE);
      fFourParticleCorrelationProfile[cs][t]->SetName(Form("fFourParticleCorrelationProfile_cs%d_combi%d", cs, t));
      fFourParticleCorrelationProfile[cs][t]->SetTitle(Form("Q-vectors, #LT #LT 4 #GT #GT, cs = %d, combi = %d", cs, t));
      fFourParticleCorrelationProfile[cs][t]->Sumw2();
      fFourParticleCorrelationProfile[cs][t]->GetXaxis()->SetTitle("#LT #LT 4 #GT #GT");
      fCorrelationResultsList->Add(fFourParticleCorrelationProfile[cs][t]);

      if (fComputeNestedLoops)
      {
        fFourNestedCorrelationProfile[cs][t] = new TProfile("", "", 1, 0., 1.);
        fFourNestedCorrelationProfile[cs][t]->SetStats(kTRUE);
        fFourNestedCorrelationProfile[cs][t]->SetName(Form("fFourNestedCorrelationProfile_cs%d_combi%d", cs, t));
        fFourNestedCorrelationProfile[cs][t]->SetTitle(Form("Nested loops, #LT #LT 4 #GT #GT, cs = %d, combi = %d", cs, t));
        fFourNestedCorrelationProfile[cs][t]->Sumw2();
        fFourNestedCorrelationProfile[cs][t]->GetXaxis()->SetTitle("#LT #LT 4 #GT #GT");
        fCorrelationResultsList->Add(fFourNestedCorrelationProfile[cs][t]);
      } // End: if (fComputeNestedLoops)
    } // End: for (Int_t t = 0; t < 6; t++)

    // TProfiles for 8-p correlations
    fEightParticleCorrelationProfile[cs] = new TProfile("", "", 1, 0., 1.);
    fEightParticleCorrelationProfile[cs]->SetStats(kTRUE);
    fEightParticleCorrelationProfile[cs]->SetName(Form("fEightParticleCorrelationProfile_cs%d", cs));
    fEightParticleCorrelationProfile[cs]->SetTitle(Form("Q-vectors, #LT #LT 8 #GT #GT, cs = %d", cs));
    fEightParticleCorrelationProfile[cs]->Sumw2();
    fEightParticleCorrelationProfile[cs]->GetXaxis()->SetTitle("#LT #LT 8 #GT #GT");
    fCorrelationResultsList->Add(fEightParticleCorrelationProfile[cs]);
  } // End: for (Int_t cs = 0; cs < endCs; cs++)
} // End: void AliAnalysisTaskTwoMultiCorrelations::BookCorrelationResultsList()

//======================================================================================//

void AliAnalysisTaskTwoMultiCorrelations::CalculateQvectors(Int_t nParticles, Double_t phi[], Double_t particleWeight[])
{
// Method calculating the "general" definition of the Q-vector Q_(n,p) for arbitrary (n,p)
  Double_t pWeightPowerP = 0.;  // Particle weight rised to the power p

// Filling of the array of Q-vectors
  for (Int_t iParticle = 0; iParticle < nParticles; iParticle++)
  {
    for (Int_t iHarmo = 0; iHarmo < (fHighHarmonic*fMaxNumberCorrelations)+1; iHarmo++)
    {
      for (Int_t iPower = 0; iPower < fMaxNumberCorrelations+1; iPower++)
      {
        pWeightPowerP = pow(particleWeight[iParticle], iPower);
        fQvectors[iHarmo][iPower] += TComplex(pWeightPowerP*TMath::Cos(iHarmo*phi[iParticle]),pWeightPowerP*TMath::Sin(iHarmo*phi[iParticle]));
      } // End: for (Int_t iPower = 0; iPower < maxPower; iPower++)
    } // End: for (Int_t iHarmo = 0; iHarmo < maxHarmo; iHarmo++)
  } // End: for (Int_t iParticle = 0; iParticle < nParticles; iParticle++)
} // End: void AliAnalysisTaskTwoMultiCorrelations::CalculateQvectors(Int_t nParticles, Double_t phi[], Double_t particleWeight[])

//======================================================================================//

TComplex AliAnalysisTaskTwoMultiCorrelations::Q(Int_t n, Int_t p)
{
// Method used to simply the usage of the Q-vectors with the fact that Q(-n,p)=Q(n,p)*
  if (n >= 0) {return fQvectors[n][p];}
  return TComplex::Conjugate(fQvectors[-n][p]);
} // End: TComplex AliAnalysisTaskTwoMultiCorrelations::Q(Int_t n, Int_p) 

//======================================================================================//

TComplex AliAnalysisTaskTwoMultiCorrelations::CalculateRecursionWithQvectors(Int_t nPartCorr, Int_t harmonics[], Int_t p, Int_t skip)
{
// Calculate the recursion for the numerator of the m-particle correlations using the Q-vectors according to the generic framework
// The recursion method was originally developped by Kristjan Gulbrandsen (gulbrand@nbi.dk)
  Int_t nMinusOne = 0;  // Harmonic n-1
  TComplex stopQvector = TComplex(0.,0.); // Stop condition of the recusion
  Int_t pPlusOne = 0; // Power p+1
  Int_t nMinusTwo = 0;  // Harmonic n-2
  Int_t counterOne = 0; // First counter for the intermediate indices
  Int_t hHold = 0;  // Temporary harmonic
  TComplex tempQvector = TComplex(0,0); // Temporary Q-vector
  Int_t counterTwo = 0; // Second counter for the intermediate indices

// Stop conditions of the recursion
  nMinusOne = nPartCorr-1;
  stopQvector = Q(harmonics[nMinusOne],p);

  if (nMinusOne == 0) {return stopQvector;}
  stopQvector *= CalculateRecursionWithQvectors(nMinusOne, harmonics);
  if (nMinusOne == skip) {return stopQvector;}

// Computation of the recursion itself
  pPlusOne = p+1;
  nMinusTwo = nPartCorr-2;
  hHold = harmonics[counterOne];
  harmonics[counterOne] = harmonics[nMinusTwo];
  harmonics[nMinusTwo] = hHold + harmonics[nMinusOne];
  tempQvector = CalculateRecursionWithQvectors(nMinusOne, harmonics, pPlusOne, nMinusTwo);
  counterTwo = nPartCorr-3;

  while (counterTwo >= skip)
  {
    harmonics[nMinusTwo] = harmonics[counterOne];
    harmonics[counterOne] = hHold;
    ++counterOne;
    
    hHold = harmonics[counterOne];
    harmonics[counterOne] = harmonics[nMinusTwo];
    harmonics[nMinusTwo] = hHold + harmonics[nMinusOne];
    tempQvector += CalculateRecursionWithQvectors(nMinusOne, harmonics, pPlusOne, counterTwo);
    --counterTwo;
  } // End: while (counterTwo >= skip)

  harmonics[nMinusTwo] = harmonics[counterOne];
  harmonics[counterOne] = hHold;

// Return of the result after the recursion
  if (p == 1) {return stopQvector - tempQvector;}
  return stopQvector - (Double_t(p)*tempQvector);

// Reset of the variables to zero
  nMinusOne = 0;
  stopQvector = TComplex(0.,0.);
  pPlusOne = 0;
  nMinusTwo = 0;
  counterOne = 0;
  hHold = 0;
  tempQvector = TComplex(0,0);
  counterTwo = 0;

} // End: TComplex AliAnalysisTaskTwoMultiCorrelations::CalculateRecursionWithQvectors(TComplex Qvectors[][9], Int_t nPartCorr, Int_t harmonics[], Int_t p, Int_t skip)

//======================================================================================//

void AliAnalysisTaskTwoMultiCorrelations::ComputeTwoParticleCorrelationWithQvectors(Int_t numeratorFirstTwoParticleHarmonics[], Int_t numeratorLastTwoParticleHarmonics[], Int_t indexTProfile)
{
// Method to compute the 2-p correlation with the Q-vectors (<<cos(h1*phi1+h2*phi2)>>,<<sin(h1*phi1+h2*phi2)>>)
  Int_t denominatorTwoParticleHarmonics[2] = {0,0}; // Harmonics with value 0 for the denominator
  Double_t denominator = 0.;  // Denominator (and event weight of the TProfile)
  TComplex firstTwoParticleCorrelation = TComplex(0.,0.); // 2-p correlation for the two first particle, used in the "normal" expression
  TComplex lastTwoParticleCorrelation = TComplex(0.,0.);  // 2-p correlation to compute <<cos(m(phi1-phi2))><cos(n(phi1-phi2))>>

// Computation of the 2-p correlations
  denominator = (CalculateRecursionWithQvectors(2, denominatorTwoParticleHarmonics)).Re();
  firstTwoParticleCorrelation = (CalculateRecursionWithQvectors(2, numeratorFirstTwoParticleHarmonics))/denominator;
  lastTwoParticleCorrelation = (CalculateRecursionWithQvectors(2, numeratorLastTwoParticleHarmonics))/denominator;

// Filling of the TProfiles
  fTwoParticleCorrelationProfile[0][indexTProfile]->Fill(0.5,firstTwoParticleCorrelation.Re(),denominator);  // Cosine component
  fTwoCosineAverageProfile[0][indexTProfile]->Fill(1.5,lastTwoParticleCorrelation.Re());
  fTwoCosineAverageProfile[0][indexTProfile]->Fill(0.5,(fTwoParticleCorrelationProfile[0][indexTProfile]->GetBinContent(1))*(fTwoCosineAverageProfile[0][indexTProfile]->GetBinContent(2)));
  if (fComputeSine) {fTwoParticleCorrelationProfile[1][indexTProfile]->Fill(0.5,firstTwoParticleCorrelation.Im(),denominator);} // Sine component

// Reset of the variables to zero
  denominator = 0;
  firstTwoParticleCorrelation = TComplex(0.,0.);
  lastTwoParticleCorrelation = TComplex(0.,0.);
} // End: void AliAnalysisTaskTwoMultiCorrelations::ComputeTwoParticleCorrelationWithQvectors(Int_t numeratorFirstTwoParticleHarmonics[], Int_t numeratorLastTwoParticleHarmonics[], Int_t indexTProfile)

//======================================================================================//

void AliAnalysisTaskTwoMultiCorrelations::ComputeFourParticleCorrelationWithQvectors(Int_t numeratorFourParticleHarmonics[], Int_t indexTProfile)
{
// Method to compute the 4-p correlation with the Q-vectors (<<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4)>>,<<sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4)>>)
  Int_t denominatorFourParticleHarmonics[4] = {0,0,0,0};  // Harmonics with value 0 for the denominator
  Double_t denominator = 0.;  // Denominator (and event weight of the TProfile)
  TComplex fourParticleCorrelation = TComplex(0.,0.); // 4-p correlation

// Computation of the 4-p correlation
  denominator = (CalculateRecursionWithQvectors(4, denominatorFourParticleHarmonics)).Re();
  fourParticleCorrelation = CalculateRecursionWithQvectors(4, numeratorFourParticleHarmonics)/denominator;

// Filling of the TProfiles
  fFourParticleCorrelationProfile[0][indexTProfile]->Fill(0.5,fourParticleCorrelation.Re(),denominator); // Cosine component 
  if (fComputeSine) {fFourParticleCorrelationProfile[1][indexTProfile]->Fill(0.5,fourParticleCorrelation.Im(),denominator);} // Sine component
} // End: void AliAnalysisTaskTwoMultiCorrelations::ComputeFourParticleCorrelationWithQvectors(Int_t numeratorFourParticleHarmonics[], Int_t indexTProfile)

//======================================================================================//

void AliAnalysisTaskTwoMultiCorrelations::ComputeSixParticleCorrelationWithQvectors(Int_t numeratorSixParticleHarmonics[], Int_t indexTProfile)
{
// Method to compute the 6-p correlation with the Q-vectors (<<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6)>>,<<sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6)>>)
  Int_t denominatorSixParticleHarmonics[6] = {0,0,0,0,0,0}; // Harmonics with value 0 for the denominator
  Double_t denominator = 0.;  // Denominator = event weight for the TProfile
  TComplex sixParticleCorrelation = TComplex(0.,0.);  // 6-p correlation

// Computation of the 6-p correlation
  denominator = (CalculateRecursionWithQvectors(6, denominatorSixParticleHarmonics)).Re();
  sixParticleCorrelation = CalculateRecursionWithQvectors(6, numeratorSixParticleHarmonics)/denominator;

// Filling of the TProfiles
  fSixParticleCorrelationProfile[0][indexTProfile]->Fill(0.5, sixParticleCorrelation.Re(), denominator); // Cosine component
  if (fComputeSine) {fSixParticleCorrelationProfile[1][indexTProfile]->Fill(0.5, sixParticleCorrelation.Im(), denominator);} // Sine component
} // End: void AliAnalysisTaskTwoMultiCorrelations::ComputeSixParticleCorrelationWithQvectors(Int_t numeratorSixParticleHarmonics[], Int_t indexTProfile)

//======================================================================================//

void AliAnalysisTaskTwoMultiCorrelations::ComputeEightParticleCorrelationWithQvectors(Int_t numeratorEightParticleHarmonics[])
{
// Method to compute the 8-p correlation with the Q-vectors (<<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7+h8*phi8)>>,<<sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7+h8*phi8)>>)
  Int_t denominatorEightParticleHarmonics[8] = {0,0,0,0,0,0,0,0}; // Harmonics with value 0 for the denominator
  Double_t denominator = 0.;  // Denominator = event weight for the TProfile
  TComplex eightParticleCorrelation = TComplex(0.,0.);  // (<<cos>>,<<sin>>)

// Computation of the 8-p correlation
  denominator = (CalculateRecursionWithQvectors(8, denominatorEightParticleHarmonics)).Re();
  eightParticleCorrelation = CalculateRecursionWithQvectors(8, numeratorEightParticleHarmonics)/denominator;

// Filling of the TProfiles
  fEightParticleCorrelationProfile[0]->Fill(0.5, eightParticleCorrelation.Re(), denominator); // Cosine component
  if(fComputeSine) {fEightParticleCorrelationProfile[1]->Fill(0.5, eightParticleCorrelation.Im(), denominator);} // Sine component
} // End: void AliAnalysisTaskTwoMultiCorrelations::ComputeEightParticleCorrelationWithQvectors(Int_t numeratorEightParticleHarmonics[], Int_t indexTProfile)

//======================================================================================//

void AliAnalysisTaskTwoMultiCorrelations::ComputeCorrelationsWithTwoNestedLoops(Int_t numeratorTwoParticleHarmonics[], Int_t numeratorTwoLastParticleHarmonics[], Int_t indexTProfile, Int_t nParticles, Double_t phi[], Double_t particleWeight[])
{
// Method to compute the 2-particle correlation with two nested loops (for cross-checking results)
  Double_t twoParticleCos = 0.;  // Cos component of the correlation = single-event 2-p correlation, <2>_(n1...nm)
  Double_t twoParticleSecondCos = 0.; // Second cosine average in the average with two cosine
  Double_t twoParticleSin = 0.; // Sin component of the correlation
  Double_t totalParticleWeight = 0.;  // Particle weight for the 2-p correlation

// Computation of the 2-p single-event average, <2>
  for (Int_t k = 0; k < nParticles; k++)
  {
    for (Int_t l = 0; l < nParticles; l++)
    {
    // Removal of the autocorrelations k == l 
      if (k == l) {continue;}
      
    // Computation of <2> for a pair of particles
      twoParticleCos = TMath::Cos(numeratorTwoParticleHarmonics[0]*phi[k] + numeratorTwoParticleHarmonics[1]*phi[l]);
      twoParticleSecondCos = TMath::Cos(numeratorTwoLastParticleHarmonics[0]*phi[k] + numeratorTwoLastParticleHarmonics[1]*phi[l]);

      if (fComputeSine) {twoParticleSin = TMath::Sin(numeratorTwoParticleHarmonics[0]*phi[k] + numeratorTwoParticleHarmonics[1]*phi[l]);}
      totalParticleWeight = particleWeight[k] * particleWeight[l];

    // Filling of the TProfile
      fTwoNestedCorrelationProfile[0][indexTProfile]->Fill(0.5, twoParticleCos,totalParticleWeight);
      fTwoCosineAverageNestedProfile[0][indexTProfile]->Fill(1.5,twoParticleSecondCos);
      if (fComputeSine) {fTwoNestedCorrelationProfile[1][indexTProfile]->Fill(0.5, twoParticleSin,totalParticleWeight);}
    } // End of the loop over the second particle of the pair
  } // End of the loop over the first particle of the pair

// Filling of the TProfile bin for <<cos><cos>>
  fTwoCosineAverageNestedProfile[0][indexTProfile]->Fill(0.5,(fTwoNestedCorrelationProfile[0][indexTProfile]->GetBinContent(1))*(fTwoCosineAverageNestedProfile[0][indexTProfile]->GetBinContent(2)));

// Reset of the local variables before changing the event
  twoParticleCos = 0.;
  twoParticleSecondCos = 0.;
  twoParticleSin = 0.;
  totalParticleWeight = 0.;

} // End: void AliAnalysisTaskTwoMultiCorrelations::ComputeCorrelationsWithTwoNestedLoops(Int_t numeratorTwoParticleHarmonics[], Int_t numeratorTwoLastParticleHarmonics[], Int_t indexTProfile, Int_t nParticles, Double_t phi[], Double_t particleWeight[])

//======================================================================================//

void AliAnalysisTaskTwoMultiCorrelations::ComputeCorrelationsWithFourNestedLoops(Int_t numeratorFourParticleHarmonics[], Int_t indexTProfile, Int_t nParticles, Double_t phi[], Double_t particleWeight[])
{
// Method to compute the 4-particle correlation with four nested loops (for cross-checking results)
  Double_t fourParticleCos = 0.;  // Cos component of the correlation = single-event 4-p correlation, <4>_(n1,... nm)
  Double_t fourParticleSin = 0.;  // Sin component of the correlation
  Double_t totalParticleWeight = 0.;  // Particle weight for the 4-p correlation

// Computation of the 4-p correlation, <4>
  for (Int_t k = 0; k < nParticles; k++)
  {
    for (Int_t l = 0; l < nParticles; l++)
    {
      // Removal of the autocorrelations k == l
      if (k == l) {continue;}

      for (Int_t m = 0; m < nParticles; m++)
      {
        // Removal of the autocorrelations k == m, l == m
        if ((k == m) || (l == m)) {continue;}

        for (Int_t n = 0; n < nParticles; n++)
        {
          // Removal of the autocorrelations k == n, l == n, m == n
          if ((k == n) || (l == n) || (m == n)) {continue;}

          // Computation of <4> for a quadruplet of particles
          fourParticleCos = TMath::Cos(numeratorFourParticleHarmonics[0]*phi[k] + numeratorFourParticleHarmonics[1]*phi[l] + numeratorFourParticleHarmonics[2]*phi[m] + numeratorFourParticleHarmonics[3]*phi[n]);
          if (fComputeSine) {fourParticleSin = TMath::Sin(numeratorFourParticleHarmonics[0]*phi[k] + numeratorFourParticleHarmonics[1]*phi[l] + numeratorFourParticleHarmonics[2]*phi[m] + numeratorFourParticleHarmonics[3]*phi[n]);}
          totalParticleWeight = particleWeight[k]*particleWeight[l]*particleWeight[m]*particleWeight[n];

          // Filling of the TProfile
          fFourNestedCorrelationProfile[0][indexTProfile]->Fill(0.5, fourParticleCos,totalParticleWeight);
          if (fComputeSine) {fFourNestedCorrelationProfile[1][indexTProfile]->Fill(0.5, fourParticleSin,totalParticleWeight);}

        } // End of the loop over the fourth particle of the quadruplet
      } // End of the loop over the third particle of the quadruplet
    } // End of the loop over the second particle of the quadruplet
  } // End of the loop over the first particle of the quadruplet

// 3.) Reset of the local variables before changing the event
  fourParticleCos = 0.;
  fourParticleSin = 0.;
  totalParticleWeight = 0.;

} // End: void AliAnalysisTaskTwoMultiCorrelations::ComputeCorrelationsWithFourNestedLoops(Int_t numeratorFourParticleHarmonics[], Int_t indexTProfile, Int_t nParticles, Double_t phi[], Double_t particleWeight[])

//======================================================================================//
