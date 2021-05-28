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

#include "AliAnalysisTaskHOCFA.h"
#include "Riostream.h"
#include "TList.h"
#include "TComplex.h"
#include "TFile.h"
#include "TH2I.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TProfile.h"
#include "TMath.h"
#include <vector>
#include <TExMap.h>
#include <TClonesArray.h>
#include "TDirectoryFile.h"
#include "TSystem.h"
#include "AliJBaseTrack.h"
#include <sstream>

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskHOCFA)

/* -------------------------------------------------------------------------- /
/ Methods inherited from AliAnalysisTaskSE.                                   /
/ -------------------------------------------------------------------------- */
AliAnalysisTaskHOCFA::AliAnalysisTaskHOCFA():
  fInputList(0),
  fDebugLevel(0),
  fMainList(0),
  fNCentralityBins(9), fCentrality(-1.), fCentralityBin(-1),
  fMultiplicity(0), fMultiplicityMin(6),
  fUseWeights(kTRUE),
  fNCombi(6)
{
// Dummy constructor of the class.
  printf("AliAnalysisTaskHOCFA Dummy constructor\n");
  InitialiseArrayMembers();
}

// ------------------------------------------------------------------------- //
AliAnalysisTaskHOCFA::AliAnalysisTaskHOCFA(const char *name, Bool_t useWeights):
  fInputList(0),
  fDebugLevel(0),
  fMainList(0),
  fNCentralityBins(9), fCentrality(-1.), fCentralityBin(-1),
  fMultiplicity(0), fMultiplicityMin(6),
  fUseWeights(kTRUE),
  fNCombi(6)
{
// Constructor of the class.
  printf("AliAnalysisTaskHOCFA Constructor\n");
  InitialiseArrayMembers();

// Define the main TList.
  fMainList = new TList();
  fMainList->SetName("analysisCFA");
  fMainList->SetOwner(kTRUE);
}

// ------------------------------------------------------------------------- //
AliAnalysisTaskHOCFA::~AliAnalysisTaskHOCFA()
{
// Destructor of the class.
  if (fMainList) {delete fMainList;}
}

// ------------------------------------------------------------------------- //
void AliAnalysisTaskHOCFA::UserCreateOutputObjects()
{
// Declare the outputs of the task at the beginning of the analysis.
  if (fDebugLevel > 1) {printf("AliAnalysisTaskHOCFA::UserCreateOutputObjects() reached\n");}

  this->BookFinalResults();
}

// ------------------------------------------------------------------------- //
void AliAnalysisTaskHOCFA::UserExec(Option_t *option)
{
// Execute the analysis for each event in the input sample.

// Get the centrality and multiplicity of the trimmed event.
  fCentralityBin = GetCentralityBin(fCentrality);
  fHistoCent[fCentralityBin]->Fill(fCentrality);
  fMultiplicity = fInputList->GetEntriesFast();

// Reject all events with not enough tracks to calculate the correlators.
  if (fMultiplicity < fMultiplicityMin) {return;}
  fHistoMulti[fCentralityBin]->Fill(fMultiplicity);

// Get the information on the selected tracks.
  Double_t *iPhi = new Double_t[fMultiplicity](); // Azimuthal angles of all particles.
  Double_t *iWeights = new Double_t[fMultiplicity]();  // Particle weights for all particles.
  Double_t iEffCorr = 1.; // Efficiency (Inverse gives the pT-weight).
  Double_t iPhiModuleCorr = 1.; // (phi, eta, PVz)-weight.

  for (Int_t iTrack = 0; iTrack < fMultiplicity; iTrack++) {
    AliJBaseTrack *aTrack = (AliJBaseTrack*)fInputList->At(iTrack);
    if (!aTrack) {continue;}

    fHistoPt[fCentralityBin]->Fill(aTrack->Pt());
    fHistoEta[fCentralityBin]->Fill(aTrack->Eta());
    iPhi[iTrack] = aTrack->Phi();
    fHistoPhi[fCentralityBin]->Fill(iPhi[iTrack]);
    fHistoCharge[fCentralityBin]->Fill(aTrack->GetCharge());

    if (fUseWeights) {  // Get the NUE/NUA corrections.
      iEffCorr = aTrack->GetTrackEff();
      iPhiModuleCorr = aTrack->GetWeight();
    }
    iWeights[iTrack] = (1.0/iEffCorr)/iPhiModuleCorr;
  }

// Compute the Q-vectors and multiparticle correlations.
  CalculateQvectors(fMultiplicity, iPhi, iWeights);
  ComputeAllTerms();

// Reset the variables for the next event.
  fMultiplicity = 0;
  delete [] iPhi;
  delete [] iWeights;
  iEffCorr = 1.;
  iPhiModuleCorr = 1.;
}

// ------------------------------------------------------------------------- //
void AliAnalysisTaskHOCFA::Terminate(Option_t *)
{
// All the additional steps after the loop over the events.
  printf("AliAnalysisTaskHOCFA is done! \n");
}

/* -------------------------------------------------------------------------- /
/ Methods specific for this class.                                            /
/ -------------------------------------------------------------------------- */
void AliAnalysisTaskHOCFA::InitialiseArrayMembers()
{
// Initialise the arrays in the data members.
  for (Int_t i = 0; i < 16; i++) {
    fCentralityList[i] = NULL;
    fCentralityArray[i] = 0.;
    fHistoCent[i] = NULL;
    fHistoMulti[i] = NULL;
    fHistoPt[i] = NULL;
    fHistoEta[i] = NULL;
    fHistoPhi[i] = NULL;
    fHistoCharge[i] = NULL;
    for (Int_t j = 0; j < 6; j++) {fCorrelTerms[j][i] = NULL;}
  }
  fCentralityArray[16] = 0.;

  for (Int_t i = 0; i < 6; i++) {
    for (Int_t j = 0; j < 3; j++) {fHarmoArray[i][j] = 0;}
  }

  Int_t powers[15][3] = { {1,0,0}, {0,1,0}, {0,0,1}, {2,0,0}, {1,1,0}, {1,0,1}, {0,1,1}, {1,1,1}, {2,1,0}, {2,0,1}, {3,0,0}, {2,1,1}, {3,1,0}, {4,0,0}, {4,1,0} };
  for (Int_t i = 0; i < 15; i++){
    for (Int_t j = 0; j < 3; j++) {fPowers[i][j] = powers[i][j];}
  }

  for (Int_t iHarmo = 0; iHarmo < 81; iHarmo++) {
    for (Int_t iPower = 0; iPower < 11; iPower++) {
      fQvectors[iHarmo][iPower] = TComplex(0., 0.);
    }
  }
}

// ------------------------------------------------------------------------- //
void AliAnalysisTaskHOCFA::SetCentralityArray(TString values)
{
// Insert the elements of the string into the centrality array.
  printf("Set the centrality binning.\n");

  Float_t value = 0;  // Current centrality value in the string.
  Int_t index = 0;  // Current index in the array.
  std::stringstream sString;  // Streamer for the string.
  sString << values;

// Get the values in the string one by one.
  while (sString >> value) {
    if (index > fNCentralityBins+1) {return;}
    fCentralityArray[index] = value;
    index++;
  }
}

// ------------------------------------------------------------------------- //
Int_t AliAnalysisTaskHOCFA::GetCentralityBin(Float_t cent)
{
// Fetch the bin corresponding to the given centrality.
  for (Int_t i = 0; i < fNCentralityBins+1; i++) {
    if (cent >= fCentralityArray[i]) {continue;}
    else {return i-1;}
  }
  return -1;  // Not a valid centrality value.
}

// ------------------------------------------------------------------------- //
void AliAnalysisTaskHOCFA::BookFinalResults()
{
// Book the needed lists and TProfiles for the results of this analysis.
  if (fDebugLevel > 1) {printf("AliAnalysisTaskHOCFA::BookFinalResults() reached\n");}
  if (!fMainList) {Fatal("AliAnalysisTaskHOCFA::BookFinalResults()", "FATAL: fMainList not found.");}

  for (Int_t i = 0; i < fNCentralityBins; i++) {
  // Book the needed centrality lists.
    fCentralityList[i] = new TList();
    fCentralityList[i]->SetName(Form("Centrality_%.1f-%.1f", fCentralityArray[i], fCentralityArray[i+1]));
    fCentralityList[i]->SetOwner(kTRUE);
    fMainList->Add(fCentralityList[i]);

  // Book the histograms in the correct lists.
    fHistoCent[i] = new TH1F("", "", 100, 0., 100.);
    fHistoCent[i]->SetName(Form("fHistoCent_Bin%d", i));
    fHistoCent[i]->SetTitle(Form("Centrality distribution, bin%d", i));
    fHistoCent[i]->SetStats(kTRUE);
    fCentralityList[i]->Add(fHistoCent[i]);

    fHistoMulti[i] = new TH1I("", "", 5000, 0., 5000.);
    fHistoMulti[i]->SetName(Form("fHistoMulti_Bin%d", i));
    fHistoMulti[i]->SetTitle(Form("Multiplicity distribution, bin%d", i));
    fHistoMulti[i]->SetStats(kTRUE);
    fCentralityList[i]->Add(fHistoMulti[i]);

    fHistoPt[i] = new TH1D("", "", 500, 0., 5.);
    fHistoPt[i]->SetName(Form("fHistoPt_Bin%d", i));
    fHistoPt[i]->SetTitle(Form("Transverse momentum distribution, bin%d", i));
    fHistoPt[i]->SetStats(kTRUE);
    fCentralityList[i]->Add(fHistoPt[i]);

    fHistoEta[i] = new TH1D("", "", 200, -1., 1.);
    fHistoEta[i]->SetName(Form("fHistoEta_Bin%d", i));
    fHistoEta[i]->SetTitle(Form("Pseudorapidity distribution, bin%d", i));
    fHistoEta[i]->SetStats(kTRUE);
    fCentralityList[i]->Add(fHistoEta[i]);

    fHistoPhi[i] = new TH1D("", "", 630, -TMath::Pi(),TMath::Pi());
    fHistoPhi[i]->SetName(Form("fHistoPhi_Bin%d", i));
    fHistoPhi[i]->SetTitle(Form("Azimuthal angle distribution, bin%d", i));
    fHistoPhi[i]->SetStats(kTRUE);
    fCentralityList[i]->Add(fHistoPhi[i]);

    fHistoCharge[i] = new TH1I("", "", 6, -3., 3.);
    fHistoCharge[i]->SetName(Form("fHistoCharge_Bin%d", i));
    fHistoCharge[i]->SetTitle(Form("Charge distribution, bin%d", i));
    fHistoCharge[i]->SetStats(kTRUE);
    fCentralityList[i]->Add(fHistoCharge[i]);

    for (Int_t j = 0; j < fNCombi; j++) {
      fCorrelTerms[j][i] = new TProfile("", "", 15, 0., 15.);
      fCorrelTerms[j][i]->SetName(Form("fCorrelTerms_Combi%d%d%d_Bin%d", fHarmoArray[j][0], fHarmoArray[j][1], fHarmoArray[j][2], i));
      fCorrelTerms[j][i]->Sumw2();
      fCentralityList[i]->Add(fCorrelTerms[j][i]);
    }
  }
}

// ------------------------------------------------------------------------- //
void AliAnalysisTaskHOCFA::SetHarmoArray(TString combiString)
{
// Divide the string of harmonics into the different combinations of 3 harmonics.
  printf("Select the arrays of combinations of harmonics.\n");

  Int_t value = 0;  // Current harmonic in the string.
  Int_t row = 0;  // Current row in the array.
  Int_t col = 0;  // Current column in the array.
  std::stringstream sString;  // Streamer for the string.
  sString << combiString;

// Get the values in the string one by one.
  while (sString >> value) {
    fHarmoArray[row][col] = value;
    if (col == 2) {col = 0; row++;}
    else {col++;}
  }
}

// ------------------------------------------------------------------------- //
void AliAnalysisTaskHOCFA::CalculateQvectors(Long64_t multiplicity, Double_t angles[], Double_t pWeights[])
{
// Calculate the needed Q-vectors.
  Double_t iAngle = 0.; // Azimuthal angle of the current track.
  Double_t iWeight = 0.;  // Particle weight of the current track.
  Double_t iWeightToPowerP = 0.;  // Particle weight rised to the power p.

// Ensure all the Q-vectors are initially zero.
  for (Int_t iHarmo = 0; iHarmo < 81; iHarmo++){
    for (Int_t iPower = 0; iPower < 11; iPower++){
      fQvectors[iHarmo][iPower] = TComplex(0., 0.);
    }
  }

// Compute the Q-vectors.
  for (Long64_t iTrack = 0; iTrack < multiplicity; iTrack++){
    iAngle = angles[iTrack];
    iWeight = pWeights[iTrack];
    for (Int_t iHarmo = 0; iHarmo < 81; iHarmo++){
      for (Int_t iPower = 0; iPower < 11; iPower++){
        iWeightToPowerP = TMath::Power(iWeight, iPower);
        fQvectors[iHarmo][iPower] += TComplex(iWeightToPowerP*TMath::Cos(iHarmo*iAngle), iWeightToPowerP*TMath::Sin(iHarmo*iAngle));
      }
    }
  }

// Reset the variables.
  iAngle = 0.;
  iWeight = 0.;
  iWeightToPowerP = 0.;
} 

// ------------------------------------------------------------------------- //
TComplex AliAnalysisTaskHOCFA::Q(Int_t n, Int_t p)
{
// Alias for fQvectors to make it more easy to use.
  if (n >= 0) {return fQvectors[n][p];}
  return TComplex::Conjugate(fQvectors[-n][p]); // Use that Q*(n,p) = Q(-n,p).
}

// ------------------------------------------------------------------------- //
TComplex AliAnalysisTaskHOCFA::CalculateRecursion(Int_t n, Int_t *harmonic, Int_t mult, Int_t skip)
{
// Calculate the multi-particle correlators by using the recursion method.
// Improved, faster version originally developed by Kristjan Gulbrandsen (gulbrand@nbi.dk).
  Int_t nm1 = n-1;
  TComplex c(Q(harmonic[nm1], mult));
  if (nm1 == 0) return c;
  c *= CalculateRecursion(nm1, harmonic);
  if (nm1 == skip) return c;

  Int_t multp1 = mult+1;
  Int_t nm2 = n-2;
  Int_t counter1 = 0;
  Int_t hhold = harmonic[counter1];
  harmonic[counter1] = harmonic[nm2];
  harmonic[nm2] = hhold + harmonic[nm1];
  TComplex c2(CalculateRecursion(nm1, harmonic, multp1, nm2));
  Int_t counter2 = n-3;

  while (counter2 >= skip) {
    harmonic[nm2] = harmonic[counter1];
    harmonic[counter1] = hhold;
    ++counter1;
    hhold = harmonic[counter1];
    harmonic[counter1] = harmonic[nm2];
    harmonic[nm2] = hhold + harmonic[nm1];
    c2 += CalculateRecursion(nm1, harmonic, multp1, counter2);
    --counter2;
  }
  harmonic[nm2] = harmonic[counter1];
  harmonic[counter1] = hhold;

  if (mult == 1) return c-c2;
  return c-Double_t(mult)*c2;
}

// ------------------------------------------------------------------------- //
void AliAnalysisTaskHOCFA::ComputeAllTerms()
{
// Compute all the terms needed for the ACs/SCs for all the observables.
  if (fDebugLevel > 5) {printf("Compute all the needed correlators.\n");}
  Int_t nTerms = 15;  // Number of terms to calculate for the ACs/SCs.
  Int_t nHarmos = 3;  // Total number of harmonics in an observable.
  //Int_t powers[nTerms][nHarmos] = { {2,0,0}, {0,2,0}, {0,0,2}, {4,0,0}, {2,2,0}, {2,0,2}, {0,2,2}, {2,2,2}, {4,2,0}, {4,0,2}, {6,0,0}, {4,2,2}, {6,2,0}, {8,0,0}, {8,2,0} };

// Loop over the number of combinations of harmonics to analyse.
  for (Int_t iCombi = 0; iCombi < fNCombi; iCombi++){
  // Customise the corresponding TProfile.
    fCorrelTerms[iCombi][fCentralityBin]->SetTitle(Form("Harmonics: (%d,%d,%d), centrality bin: %d", fHarmoArray[iCombi][0], fHarmoArray[iCombi][1], fHarmoArray[iCombi][2], fCentralityBin));

  // Fill the list of harmonics for the correlator itself.
    for (Int_t jTerm = 0; jTerm < nTerms; jTerm++) {
      Int_t jHarmonics[7] = {0}; // Array of harmonics for the recursion.
      Int_t counter = 0;
      Int_t nPart = 0;  // Number of particles in the correlator.

      for (Int_t iHarmo = 0; iHarmo < nHarmos; iHarmo++) {
        if (fPowers[jTerm][iHarmo] == 0) {continue;} // Skip the unneeded harmonics.

        for (Int_t jPower = 0; jPower < fPowers[jTerm][iHarmo]; jPower++) {
          jHarmonics[counter] = fHarmoArray[iCombi][iHarmo];
          counter++;
        }

        nPart += 2*fPowers[jTerm][iHarmo];
      }

    // Calculate the multiparticle correlator itself using the recursion method.
      CalculateCorrelator(iCombi, jTerm, nPart, jHarmonics);
    }
  }
}

// ------------------------------------------------------------------------- //
void AliAnalysisTaskHOCFA::CalculateCorrelator(Int_t combi, Int_t bin, Int_t nParticles, Int_t harmonics[])
{
// Calculate the multiparticle correlator corresponding to harmonics[].
  if (fDebugLevel > 5) {printf("Calculate the correlator for the provided harmonics[].\n");}
  TComplex cCorrel = TComplex(0., 0.);
  Double_t eventWeight = 0.;  // Event weight for this correlator.
  Double_t rCorrel = 0.;

  Int_t twoHarmoNum[2] = {harmonics[0], -1*harmonics[0]};
  Int_t twoHarmoDen[2] = {0};
  Int_t fourHarmoNum[4] = {harmonics[0], harmonics[1], -1*harmonics[0], -1*harmonics[1]};
  Int_t fourHarmoDen[4] = {0};
  Int_t sixHarmoNum[6] = {harmonics[0], harmonics[1], harmonics[2], -1*harmonics[0], -1*harmonics[1], -1*harmonics[2]};
  Int_t sixHarmoDen[6] = {0};
  Int_t eightHarmoNum[8] = {harmonics[0], harmonics[1], harmonics[2], harmonics[3], -1*harmonics[0], -1*harmonics[1], -1*harmonics[2], -1*harmonics[3]};
  Int_t eightHarmoDen[8] = {0};
  Int_t tenHarmoNum[10] = {harmonics[0], harmonics[1], harmonics[2], harmonics[3], harmonics[4], -1*harmonics[0], -1*harmonics[1], -1*harmonics[2], -1*harmonics[3], -1*harmonics[4]};
  Int_t tenHarmoDen[10] = {0};

// Compute the needed quantities.
  switch(nParticles) {
    case 2 :
      eventWeight = ( CalculateRecursion(2, twoHarmoDen) ).Re();
      cCorrel = ( CalculateRecursion(2, twoHarmoNum) )/eventWeight;
      rCorrel = cCorrel.Re();
      break;
    case 4 :
      eventWeight = ( CalculateRecursion(4, fourHarmoDen) ).Re();
      cCorrel = ( CalculateRecursion(4, fourHarmoNum) )/eventWeight;
      rCorrel = cCorrel.Re();
      break;
    case 6 :    
      eventWeight = ( CalculateRecursion(6, sixHarmoDen) ).Re();
      cCorrel = ( CalculateRecursion(6, sixHarmoNum) )/eventWeight;
      rCorrel = cCorrel.Re();
      break;
    case 8 :
      eventWeight = ( CalculateRecursion(8, eightHarmoDen) ).Re();
      cCorrel = ( CalculateRecursion(8, eightHarmoNum) )/eventWeight;
      rCorrel = cCorrel.Re();
      break;
    case 10 :
      eventWeight = ( CalculateRecursion(10, tenHarmoDen) ).Re();
      cCorrel = ( CalculateRecursion(10, tenHarmoNum) )/eventWeight;
      rCorrel = cCorrel.Re();
      break;
    default :
      printf("Error: invalid number of particles (odd or larger than ten).\n");
      break;
  }

// Fill the corresponding bin in the right TProfile.
  fCorrelTerms[combi][fCentralityBin]->Fill( (Float_t)bin + 0.5, rCorrel, eventWeight );
  fCorrelTerms[combi][fCentralityBin]->GetXaxis()->SetBinLabel(bin, Form("{%d,%d,%d}", fPowers[bin][0], fPowers[bin][1], fPowers[bin][2]));

// Reset the local variables for the next call.
  cCorrel = TComplex(0., 0.);
  eventWeight = 0.;
  rCorrel = 0.;
}

// ------------------------------------------------------------------------- //
