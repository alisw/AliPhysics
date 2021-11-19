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
#include "TSystem.h"
#include "Riostream.h"
#include "TClonesArray.h"
#include "TList.h"
#include "TString.h"
#include "TComplex.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "AliJBaseTrack.h"
#include "TMath.h"

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
  fMultiplicity(0), fMultiplicityMin(10),
  fPtMin(0.2), fPtMax(5.),
  fUseWeightsNUE(kTRUE), fUseWeightsNUA(kFALSE),
  fGetSC3h(kTRUE),
  fGetEtaGap(kFALSE), fEtaGap(0.),
  fNCombi(6),
  fHistoConfig(0)
{
// Dummy constructor of the class.
  printf("AliAnalysisTaskHOCFA Dummy constructor\n");
  InitialiseArrayMembers();
}

// ------------------------------------------------------------------------- //
AliAnalysisTaskHOCFA::AliAnalysisTaskHOCFA(const char *name):
  fInputList(0),
  fDebugLevel(0),
  fMainList(0),
  fNCentralityBins(9), fCentrality(-1.), fCentralityBin(-1),
  fMultiplicity(0), fMultiplicityMin(10),
  fPtMin(0.2), fPtMax(5.),
  fUseWeightsNUE(kTRUE), fUseWeightsNUA(kFALSE),
  fGetSC3h(kTRUE),
  fGetEtaGap(kFALSE), fEtaGap(0.),
  fNCombi(6),
  fHistoConfig(0)
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
// Declare the output of the task at the beginning of the analysis.
  if (fDebugLevel > 1) {
    printf("AliAnalysisTaskHOCFA::UserCreateOutputObjects() reached.\n");
  }

  BookFinalResults();

  // Save the configuration of the analysis.
  fHistoConfig->Fill(0.5, fPtMin);
  fHistoConfig->Fill(1.5, fPtMax);
  fHistoConfig->Fill(2.5, fMultiplicityMin);
  if (fUseWeightsNUE) {fHistoConfig->Fill(3.5, 1.);}
  if (fUseWeightsNUA) {fHistoConfig->Fill(4.5, 1.);}
  if (fGetEtaGap) {fHistoConfig->Fill(5.5, 1.);}
}

// ------------------------------------------------------------------------- //
void AliAnalysisTaskHOCFA::UserExec(Option_t *option)
{
// Execute the analysis for each event in the input sample.

// Get the centrality and multiplicity of the trimmed event.
  fMultiplicity = fInputList->GetEntriesFast();
  fCentralityBin = GetCentralityBin(fCentrality);

// Reject all events with not enough tracks to calculate the correlators.
  if (fMultiplicity < fMultiplicityMin) {return;}
  fHistoCent[fCentralityBin]->Fill(fCentrality);
  fHistoMulti[fCentralityBin]->Fill(fMultiplicity);

// Get the information on the selected tracks.
  double *iEta = new double[fMultiplicity]();
  double *iPhi = new double[fMultiplicity]();
  double *iWeights = new double[fMultiplicity]();
  double iEffCorr = 1.;         // Efficiency (pT-weight = 1/efficiency).
  double iPhiModuleCorr = 1.;   // (phi, eta, PVz)-weight.

  for (Long64_t iTrack = 0; iTrack < fMultiplicity; iTrack++) {
    AliJBaseTrack *aTrack = (AliJBaseTrack*)fInputList->At(iTrack);
    if (!aTrack) {continue;}
    iEta[iTrack] = aTrack->Eta();
    iPhi[iTrack] = aTrack->Phi();

    fHistoPt[fCentralityBin]->Fill(aTrack->Pt());
    fHistoEta[fCentralityBin]->Fill(iEta[iTrack]);
    fHistoPhi[fCentralityBin]->Fill(iPhi[iTrack]);
    fHistoCharge[fCentralityBin]->Fill(aTrack->GetCharge());

    if (fUseWeightsNUE) {iEffCorr = aTrack->GetTrackEff();}
    if(fUseWeightsNUA) {iPhiModuleCorr = aTrack->GetWeight();}
    iWeights[iTrack] = (1.0/iEffCorr)/iPhiModuleCorr;
  }

// Compute the Q-vectors and multiparticle correlations.
  CalculateQvectors(fMultiplicity, iPhi, iWeights);
  ComputeAllTerms();

// If true, calculate the 2-particle correlators for v1 to v8 using eta gaps.
  if (fGetEtaGap) {ComputeEtaGaps(fMultiplicity, iPhi, iWeights, iEta);}

// Reset the variables for the next event.
  fMultiplicity = 0;
  delete [] iEta; delete [] iPhi; delete [] iWeights;
  iEffCorr = 1.; iPhiModuleCorr = 1.;
}

// ------------------------------------------------------------------------- //
void AliAnalysisTaskHOCFA::Terminate(Option_t *)
{
// All the additional steps after the loop over the events.
  printf("AliAnalysisTaskHOCFA is done!\n");
}

/* -------------------------------------------------------------------------- /
/ Methods specific to this class.                                             /
/ -------------------------------------------------------------------------- */
void AliAnalysisTaskHOCFA::InitialiseArrayMembers()
{
// Initialise the arrays in the data members.
  for (int i = 0; i < 16; i++) {
    fCentralityList[i] = NULL;
    fCentralityArray[i] = 0.;
    fHistoCent[i] = NULL;
    fHistoMulti[i] = NULL;
    fHistoPt[i] = NULL;
    fHistoEta[i] = NULL;
    fHistoPhi[i] = NULL;
    fHistoCharge[i] = NULL;
    fCorrelEtaGap[i] = NULL;
    for (int j = 0; j < 6; j++) {
      fCorrelTerms[j][i] = NULL;
      fErrorTermsSC3h[j][i] = NULL;
      fErrorTermsAC41[j][i] = NULL;
    }
  }
  fCentralityArray[16] = 0.;

  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 3; j++) {fHarmoArray[i][j] = 0;}
  }

  int powers[15][3] = { {1,0,0}, {0,1,0}, {0,0,1}, {2,0,0}, {1,1,0}, {1,0,1}, {0,1,1}, {1,1,1}, {2,1,0}, {2,0,1}, {3,0,0}, {2,1,1}, {3,1,0}, {4,0,0}, {4,1,0} };
  for (int i = 0; i < 15; i++){
    for (int j = 0; j < 3; j++) {fPowers[i][j] = powers[i][j];}
  }

  for (int iHarmo = 0; iHarmo < 81; iHarmo++) {
    for (int iPower = 0; iPower < 11; iPower++) {
      fQvectors[iHarmo][iPower] = TComplex(0., 0.);
    }
  }

  if (fDebugLevel > 5) {printf("'ArrayMembers' initialised.\n");}
}

// ------------------------------------------------------------------------- //
void AliAnalysisTaskHOCFA::SetCentralityArray(TString values)
{
// Insert the elements of the string into the centrality array.
  printf("Set the centrality binning.\n");

  float value = 0;              // Current string element.
  int index = 0;                // Current index in the array.
  std::stringstream sString;    // Streamer for the string.
  sString << values;

// Get the values in the string one by one.
  while (sString >> value) {
    if (index > fNCentralityBins+1) {return;}
    fCentralityArray[index] = value;
    index++;
  }
}

// ------------------------------------------------------------------------- //
int AliAnalysisTaskHOCFA::GetCentralityBin(float cent)
{
// Fetch the bin corresponding to the given centrality.
  for (int i = 0; i < fNCentralityBins+1; i++) {
    if (cent >= fCentralityArray[i]) {continue;}
    else {return i-1;}
  }
  return -1;  // Not a valid centrality value.
}

// ------------------------------------------------------------------------- //
void AliAnalysisTaskHOCFA::BookFinalResults()
{
// Book the needed objects for the analysis results.
  if (fDebugLevel > 1) {printf("AliAnalysisTaskHOCFA::BookFinalResults() reached.\n");}
  if (!fMainList) {Fatal("AliAnalysisTaskHOCFA::BookFinalResults()", "FATAL: fMainList not found.");}

// Book the configuration histogram.
  fHistoConfig = new TProfile("", "", 6, 0., 6.);
  fHistoConfig->SetName("fHistoConfig");
  fHistoConfig->SetTitle("Configuration of the analysis");
  fHistoConfig->SetStats(kFALSE);
  fHistoConfig->GetXaxis()->SetBinLabel(1, "Min pT");
  fHistoConfig->GetXaxis()->SetBinLabel(2, "Max pT");
  fHistoConfig->GetXaxis()->SetBinLabel(3, "Min M");
  fHistoConfig->GetXaxis()->SetBinLabel(4, "NUE?");
  fHistoConfig->GetXaxis()->SetBinLabel(5, "NUA?");
  fHistoConfig->GetXaxis()->SetBinLabel(6, "Eta gap?");
  fMainList->Add(fHistoConfig);

  for (int i = 0; i < fNCentralityBins; i++) {
  // Book the centrality lists only for the chosen range.
    fCentralityList[i] = new TList();
    fCentralityList[i]->SetName(Form(
      "Centrality_%.1f-%.1f", fCentralityArray[i], fCentralityArray[i+1]));
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

    fCorrelEtaGap[i] = new TProfile("", "", 8, 0., 8.);
    fCorrelEtaGap[i]->SetName(Form("fCorrelEtaGap_Bin%d", i));
    fCorrelEtaGap[i]->SetTitle(Form("|#it{#Delta #eta}| > %.1f, bin%d", fEtaGap, i));
    fCorrelEtaGap[i]->SetStats(kTRUE);
    if (fGetEtaGap) {fCentralityList[i]->Add(fCorrelEtaGap[i]);}

    for (Int_t j = 0; j < fNCombi; j++) {
      fCorrelTerms[j][i] = new TProfile("", "", 15, 0., 15.);
      fCorrelTerms[j][i]->SetName(Form("fCorrelTerms_Combi%d%d%d_Bin%d", fHarmoArray[j][0], fHarmoArray[j][1], fHarmoArray[j][2], i));
      fCorrelTerms[j][i]->Sumw2();
      fCentralityList[i]->Add(fCorrelTerms[j][i]);

      fErrorTermsSC3h[j][i] = new TProfile("", "", 35, 0., 35.);
      fErrorTermsSC3h[j][i]->SetName(Form("fErrorTermsSC3h_Combi%d%d%d_Bin%d", fHarmoArray[j][0], fHarmoArray[j][1], fHarmoArray[j][2], i));
      fErrorTermsSC3h[j][i]->Sumw2();
      if (fGetSC3h) {fCentralityList[i]->Add(fErrorTermsSC3h[j][i]);}

      fErrorTermsAC41[j][i] = new TProfile("", "", 54, 0., 54.);
      fErrorTermsAC41[j][i]->SetName(Form("fErrorTermsAC41_Combi%d%d_Bin%d", fHarmoArray[j][0], fHarmoArray[j][1], i));
      fErrorTermsAC41[j][i]->Sumw2();
      if (!fGetSC3h) {fCentralityList[i]->Add(fErrorTermsAC41[j][i]);}
    }
  }
}

// ------------------------------------------------------------------------- //
void AliAnalysisTaskHOCFA::SetHarmoArray(TString combiString)
{
// Convert the string of harmonics into the corresponding combinations.
  printf("Convert the arrays of combinations of harmonics.\n");

  int value = 0;    // Current harmonic in the string.
  int row = 0;      // Current row in the array.
  int col = 0;      // Current column in the array.
  std::stringstream sString;
  sString << combiString;

// Get the values in the string one by one.
  while (sString >> value) {
    fHarmoArray[row][col] = value;
    if (col == 2) {col = 0; row++;}
    else {col++;}
  }
}

// ------------------------------------------------------------------------- //
void AliAnalysisTaskHOCFA::CalculateQvectors(Long64_t multiplicity,
  double angles[], double pWeights[])
{
// Calculate the needed Q-vectors.
  double iAngle = 0.;           // Azimuthal angle of the current track.
  double iWeight = 0.;          // Particle weight of the current track.
  double iWeightToPowerP = 0.;  // Particle weight rised to the power p.

// Ensure all the Q-vectors are initially zero.
  for (int iHarmo = 0; iHarmo < 81; iHarmo++) {
    for (int iPower = 0; iPower < 11; iPower++) {
      fQvectors[iHarmo][iPower] = TComplex(0., 0.);
    }
  }

// Compute the Q-vectors.
  for (Long64_t iTrack = 0; iTrack < multiplicity; iTrack++) {
    iAngle = angles[iTrack];
    iWeight = pWeights[iTrack];

    for (int iHarmo = 0; iHarmo < 81; iHarmo++) {
      for (int iPower = 0; iPower < 11; iPower++) {
        iWeightToPowerP = TMath::Power(iWeight, iPower);
        fQvectors[iHarmo][iPower] += TComplex(iWeightToPowerP*TMath::Cos(iHarmo*iAngle), iWeightToPowerP*TMath::Sin(iHarmo*iAngle));
      }
    }
  }

// Reset the variables.
  iAngle = 0.; iWeight = 0.; iWeightToPowerP = 0.;
} 

// ------------------------------------------------------------------------- //
TComplex AliAnalysisTaskHOCFA::Q(int n, int p)
{
// Alias for fQvectors to make it more easy to use.
  if (n >= 0) {return fQvectors[n][p];}
  return TComplex::Conjugate(fQvectors[-n][p]); // Use Q*(n,p) = Q(-n,p).
}

// ------------------------------------------------------------------------- //
TComplex AliAnalysisTaskHOCFA::CalculateRecursion(int n, int *harmonic,
  int mult, int skip)
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
  const int nTerms = 9;    // Number of terms in the CFA (9 for AC(m,n).)
  int nHarmos = 3;    // Total number of harmonics in an observable.
  //if (fGetSC3h) {nHarmos = 3;}
  //Int_t powers[nTerms][nHarmos] = { {2,0,0}, {0,2,0}, {0,0,2}, {4,0,0}, {2,2,0}, {2,0,2}, {0,2,2}, {2,2,2}, {4,2,0}, {4,0,2}, {6,0,0}, {4,2,2}, {6,2,0}, {8,0,0}, {8,2,0} };

// Loop over the number of combinations to analyse.
  for (int iCombi = 0; iCombi < fNCombi; iCombi++) {
  // Customise the corresponding TProfile.
    fCorrelTerms[iCombi][fCentralityBin]->SetTitle(Form("Harmonics: (%d,%d,%d), centrality bin: %d", fHarmoArray[iCombi][0], fHarmoArray[iCombi][1], fHarmoArray[iCombi][2], fCentralityBin));
    if (fGetSC3h) {
      fErrorTermsSC3h[iCombi][fCentralityBin]->SetTitle(Form("Errors for SC(k,l,m), Harmonics: (%d,%d,%d), centrality bin: %d", fHarmoArray[iCombi][0], fHarmoArray[iCombi][1], fHarmoArray[iCombi][2], fCentralityBin));
    } else {
      fErrorTermsAC41[iCombi][fCentralityBin]->SetTitle(Form("Errors for AC_41, Harmonics: (%d,%d), centrality bin: %d", fHarmoArray[iCombi][0], fHarmoArray[iCombi][1], fCentralityBin));
    }

  // Define the variables needed for the error propagation.
    int i = 0;  // Index for the input array.
    double inputArray[nTerms] = {0.};     // Single-event correlators.
    double inputWeights[nTerms] = {0.};   // Corresponding event weights.
    double *tempData = new double[2];

  // Fill the list of harmonics for the correlator itself.
    for (int jTerm = 0; jTerm < 15; jTerm++) {
      if (!fGetSC3h && (fPowers[jTerm][2] != 0)) {continue;}    // Skip terms with non-zero 3rd harmonics.
      if (fGetSC3h) {
        if (jTerm == 3) {continue;}   // Skip v^4 in SC.
        if (jTerm == 8) {break;}      // All powers for SC(k,l,m) have been read.
      }
      int jHarmonics[7] = {0};   // Harmonics for the recursion.
      int counter = 0;
      int nPart = 0;    // Number of particles in the correlator.

      for (int iHarmo = 0; iHarmo < nHarmos; iHarmo++) {
        if (fPowers[jTerm][iHarmo] == 0) {continue;} // Skip the unneeded harmonics.
        //if (!fGetSC3h && (fPowers[jTerm][iHarmo] != 0))

        for (int jPower = 0; jPower < fPowers[jTerm][iHarmo]; jPower++) {
          jHarmonics[counter] = fHarmoArray[iCombi][iHarmo];
          counter++;
        }

        nPart += 2*fPowers[jTerm][iHarmo];
      }
    if (fDebugLevel > 5) {printf("Harmo (%d,%d,%d). nPart %d\n", jHarmonics[0], jHarmonics[1],jHarmonics[2], nPart);}

    // Calculate the multiparticle correlator itself using the recursion method.
      CalculateCorrelator(iCombi, jTerm, nPart, jHarmonics, tempData);

    // Fill the input array to calculate all the terms needed for the error propagation.
      inputArray[i] = tempData[0];
      inputWeights[i] = tempData[1];
      i++;
    }

    if (fDebugLevel > 5) {printf("Fill the output arrays.\n");}
    // Fill the output array and profiles (always for non-normalised.)
    // The lists for 'normalised' contain permutations of the same elements.
    const int outputSize = 54;
    double outputArray[outputSize] = {0.};
    double outputWeights[outputSize] = {0.};

    int last_element_nTerm = (fGetSC3h) ? 7 : nTerms;
    int nOut = nTerms;
    for (int iOut = 0; iOut < last_element_nTerm; iOut++) {
      // Fill the nTerms first elements.
      outputArray[iOut] = inputArray[iOut];
      outputWeights[iOut] = inputWeights[iOut];
    }
    for (int iOut = 0; iOut < last_element_nTerm; iOut++) {
      // Fill the combinations of single-event correlators.
      for (int jOut = iOut; jOut < last_element_nTerm; jOut++) {
        outputArray[nOut] = inputArray[iOut]*inputArray[jOut];
        outputWeights[nOut] = inputWeights[iOut]*inputWeights[jOut];
        nOut++;
      }
    }

    if (fDebugLevel > 5) {printf("Fill the output profiles.\n");}
    int last_element_output = (fGetSC3h) ? 35 : outputSize;
    for (int iO = 0; iO < last_element_output; iO++) {
      if (fGetSC3h) {
        fErrorTermsSC3h[iCombi][fCentralityBin]->Fill((float)iO + 0.5, outputArray[iO], outputWeights[iO]);
        fErrorTermsSC3h[iCombi][fCentralityBin]->GetXaxis()->SetBinLabel(iO+1, Form("corr%d", iO+1));
      } else {
        fErrorTermsAC41[iCombi][fCentralityBin]->Fill((float)iO + 0.5, outputArray[iO], outputWeights[iO]);
        fErrorTermsAC41[iCombi][fCentralityBin]->GetXaxis()->SetBinLabel(iO+1, Form("corr%d", iO+1));      
      }
    }

    delete [] tempData;
  }
}

// ------------------------------------------------------------------------- //
void AliAnalysisTaskHOCFA::CalculateCorrelator(int combi, int bin,
  int nParticles, int harmonics[], double *errorTerms)
{
// Calculate the multiparticle correlator corresponding to harmonics[].
  if (fDebugLevel > 5) {printf("Calculate correlator for 'harmonics[]'.\n");}
  TComplex cCorrel = TComplex(0., 0.);
  double eventWeight = 0.;  // Event weight for this correlator.
  double rCorrel = 0.;

  int twoHarmoNum[2] = {harmonics[0], -1*harmonics[0]};
  int twoHarmoDen[2] = {0};
  int fourHarmoNum[4] = {harmonics[0], harmonics[1], -1*harmonics[0], -1*harmonics[1]};
  int fourHarmoDen[4] = {0};
  int sixHarmoNum[6] = {harmonics[0], harmonics[1], harmonics[2], -1*harmonics[0], -1*harmonics[1], -1*harmonics[2]};
  int sixHarmoDen[6] = {0};
  int eightHarmoNum[8] = {harmonics[0], harmonics[1], harmonics[2], harmonics[3], -1*harmonics[0], -1*harmonics[1], -1*harmonics[2], -1*harmonics[3]};
  int eightHarmoDen[8] = {0};
  int tenHarmoNum[10] = {harmonics[0], harmonics[1], harmonics[2], harmonics[3], harmonics[4], -1*harmonics[0], -1*harmonics[1], -1*harmonics[2], -1*harmonics[3], -1*harmonics[4]};
  int tenHarmoDen[10] = {0};

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
    printf("Error: invalid number of particles.\n");
    break;
  }

// Fill the corresponding bin in the right TProfile.
  fCorrelTerms[combi][fCentralityBin]->Fill( (float)bin + 0.5, rCorrel, eventWeight );
  fCorrelTerms[combi][fCentralityBin]->GetXaxis()->SetBinLabel(bin+1, Form("{%d,%d,%d}", fPowers[bin][0], fPowers[bin][1], fPowers[bin][2]));

// Fill the temporary informations.
  errorTerms[0] = rCorrel;
  errorTerms[1] = eventWeight;

// Reset the local variables for the next call.
  cCorrel = TComplex(0., 0.);
  rCorrel = 0.;
  eventWeight = 0.;
}

// ------------------------------------------------------------------------- //
void AliAnalysisTaskHOCFA::ComputeEtaGaps(Long64_t multiplicity,
  double angles[], double pWeights[], double pseudorapidity[])
{
// Calculate the two-particle correlators (v1..v8) with a given eta gap.
  if (fDebugLevel > 5) {printf("Calculate 2-particle correlators with eta gap'.\n");}
  TComplex Qminus[8] = {TComplex(0., 0.)};  // Q-vectors for the negative subevent.
  float Mminus[8] = {0.};                   // Multiplicity of the negative subevent.
  TComplex Qplus[8] = {TComplex(0., 0.)};
  float Mplus[8] = {0.};

  float iAngle = 0.;      // Azimuthal angle.
  float iEta = 0.;        // Pseudorapidity.
  float iWeight = 1.;     // Particle weight.
  float iWeightToP = 1.;  // Particle weight raised to the power p.

  TComplex cCorrel = TComplex(0., 0.);
  double rCorrel = 0.;

  for (Long64_t iPart = 0; iPart < multiplicity; iPart++) {
    // Read the information for the current track.
    iAngle  = angles[iPart];
    iWeight = pWeights[iPart];
    iEta    = pseudorapidity[iPart];
    if (fUseWeightsNUE || fUseWeightsNUA) {iWeightToP = iWeight;}

    // Compute the Q-vectors for each subevent.
    if (iEta < 0.) {
      for (int iHarmo = 0; iHarmo < 8; iHarmo++) {
        // Compute only if the particle belongs to the subevent.
        if (iEta < ((-0.5)*fEtaGap)) {
          Qminus[iHarmo] += TComplex(iWeightToP*TMath::Cos((iHarmo+1)*iAngle),
            iWeightToP*TMath::Sin((iHarmo+1)*iAngle));
          Mminus[iHarmo] += iWeightToP;
        }
      }

    } else if (iEta > 0.) {
      for (int iHarmo = 0; iHarmo < 8; iHarmo++) {
        if (iEta > 0.5*fEtaGap) {
          Qplus[iHarmo] += TComplex(iWeightToP*TMath::Cos((iHarmo+1)*iAngle),
            iWeightToP*TMath::Sin((iHarmo+1)*iAngle));
          Mplus[iHarmo] += iWeightToP;
        }
      }

    } else {continue;}  // Case iEta = 0.
  }   // All the particles have been treated.

  // Compute the two-particle correlators themselves.
  for (int iHarmo = 0; iHarmo < 8; iHarmo++) {
    if (!( (Qminus[iHarmo].TComplex::Rho() > 0.) && (Qplus[iHarmo].TComplex::Rho() > 0.) )) {continue;}
    if (!( (Mminus[iHarmo] > 0.) && (Mplus[iHarmo] > 0.) )) {continue;}

    cCorrel = Qminus[iHarmo]*TComplex::Conjugate(Qplus[iHarmo]);
    rCorrel = (cCorrel.Re())/(Mminus[iHarmo]*Mplus[iHarmo]);
    fCorrelEtaGap[fCentralityBin]->Fill((float)iHarmo + 0.5, rCorrel, Mminus[iHarmo]*Mplus[iHarmo]);
    fCorrelEtaGap[fCentralityBin]->GetXaxis()->SetBinLabel(iHarmo + 1, Form("v_{%d}", iHarmo+1));

    // Reset the correlators to prevent mixing between harmonics.
    cCorrel = TComplex(0., 0.);
    rCorrel = 0.;
  }
// All the 2-particle correlators with eta gaps have been calculated.
}