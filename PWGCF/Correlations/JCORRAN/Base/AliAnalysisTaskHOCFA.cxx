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

#include <sstream>

#include "Riostream.h"
#include "TClonesArray.h"
#include "TComplex.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TList.h"
#include "TMath.h"
#include "TProfile.h"
#include "TString.h"
#include "TSystem.h"

#include "AliJBaseTrack.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskHOCFA)

/* -------------------------------------------------------------------------- /
/ Methods inherited from AliAnalysisTaskSE.                                   /
/ -------------------------------------------------------------------------- */
AliAnalysisTaskHOCFA::AliAnalysisTaskHOCFA():
  fInputList(0),
  fMainList(0),
  fDebugLevel(0),
  fMultiplicity(0), fCentrality(-1.),
  fNCentralityBins(9), fCentralityBin(-1), fMultiplicityMin(10),
  fPtMin(0.2), fPtMax(5.),
  fEtaGap(0.), fApplyEtaGap(kFALSE), 
  fUseWeightsNUE(kTRUE), fUseWeightsNUA(kFALSE), fUseWeightsCent(kFALSE),
  fNCombi(7), fGetSC(kTRUE), fGetLowerHarmos(kTRUE),
  fHistoConfig(0)
{
// Dummy constructor of the class.
  printf("AliAnalysisTaskHOCFA Dummy constructor\n");
  InitialiseArrayMembers();
}

// ------------------------------------------------------------------------- //
AliAnalysisTaskHOCFA::AliAnalysisTaskHOCFA(const char *name):
  fInputList(0),
  fMainList(0),
  fDebugLevel(0),
  fMultiplicity(0), fCentrality(-1.),
  fNCentralityBins(9), fCentralityBin(-1), fMultiplicityMin(10),
  fPtMin(0.2), fPtMax(5.),
  fEtaGap(0.), fApplyEtaGap(kFALSE), 
  fUseWeightsNUE(kTRUE), fUseWeightsNUA(kFALSE), fUseWeightsCent(kFALSE),
  fNCombi(7), fGetSC(kTRUE), fGetLowerHarmos(kTRUE),
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
  if (fDebugLevel > 5) {printf("AliAnalysisTaskHOCFA::UserCreateOutputObjects().\n");}
  BookFinalResults();

  // Save the configuration of the analysis.
  fHistoConfig->Fill(0., fMultiplicityMin);
  fHistoConfig->Fill(1., fPtMin);
  fHistoConfig->Fill(2., fPtMax);
  if (fApplyEtaGap) {fHistoConfig->Fill(3., fEtaGap);}
  if (fUseWeightsNUE) {fHistoConfig->Fill(4., 1.);}
  if (fUseWeightsNUA) {fHistoConfig->Fill(5., 1.);}
  if (fGetSC) {fHistoConfig->Fill(6., 1.);}
  if (fGetLowerHarmos) {fHistoConfig->Fill(7., 1.);}
}

// ------------------------------------------------------------------------- //
void AliAnalysisTaskHOCFA::UserExec(Option_t *option)
{
// Execute the analysis for each event in the input sample.

  // Check if the event has enough tracks for the largest correlator.
  // If not, we reject it, else we fill the QA for the multiplicity and centrality.
  fMultiplicity = fInputList->GetEntriesFast();
  if (fMultiplicity < fMultiplicityMin) {return;}

  fCentralityBin = GetCentralityBin(fCentrality);
  fHistoMulti[fCentralityBin]->Fill(fMultiplicity);
  fHistoCent[fCentralityBin]->Fill(fCentrality);

// Get the information on all the tracks preselected by the catalyst.
  double *iEta = new double[fMultiplicity]();
  double *iPhi = new double[fMultiplicity]();
  double *iWeights = new double[fMultiplicity]();
  double iEffCorr = 1.;         // Efficiency (pT-weight = 1/efficiency).
  double iPhiModuleCorr = 1.;   // (phi, eta, PVz)-weight.
  float iCentWeight = 1.;       // Centrality weight for LHC15o.

  for (Long64_t iTrack = 0; iTrack < fMultiplicity; iTrack++) {
    AliJBaseTrack *aTrack = (AliJBaseTrack*)fInputList->At(iTrack);
    if (!aTrack) {continue;}
    iEta[iTrack] = aTrack->Eta();
    iPhi[iTrack] = aTrack->Phi();

    if (fUseWeightsNUE) {iEffCorr = aTrack->GetTrackEff();}
    if (fUseWeightsNUA) {iPhiModuleCorr = aTrack->GetWeight();}
    if (fDebugLevel > 10) printf("iEffCorr: %.6f iPhiModuleCorr: %.6f \n", iEffCorr, iPhiModuleCorr);
    if (fUseWeightsCent) {iCentWeight = aTrack->GetCentWeight();} // Same for all tracks in an event.
    iWeights[iTrack] = (1.0/iEffCorr)/iPhiModuleCorr;

    fHistoPt[fCentralityBin]->Fill(aTrack->Pt(), (1./iEffCorr));
    fHistoEta[fCentralityBin]->Fill(iEta[iTrack]);
    fHistoPhi[fCentralityBin]->Fill(iPhi[iTrack], (1./iPhiModuleCorr));
    fHistoCharge[fCentralityBin]->Fill(aTrack->GetCharge());
  } // Go to the next track.

// Fill the QA for the centrality*weight for this event.
  fHistoCentCorrect[fCentralityBin]->Fill(fCentrality, iCentWeight);

// Compute the Q-vectors and multiparticle correlations.
  CalculateQvectors(fMultiplicity, iPhi, iWeights);
  ComputeAllTerms(iCentWeight);

// If true, calculate the 2-particle correlators for v1 to v8 using eta gaps.
  if (fApplyEtaGap) {ComputeEtaGaps(fMultiplicity, iPhi, iWeights, iEta, iCentWeight);}

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
    fHistoCentCorrect[i] = NULL;
    fHistoMulti[i] = NULL;
    fHistoPt[i] = NULL;
    fHistoEta[i] = NULL;
    fHistoPhi[i] = NULL;
    fHistoCharge[i] = NULL;

    fCorrel2p[i] = NULL;
    fCorrel2pEtaGap[i] = NULL;
    for (int j = 0; j < 7; j++) {
      fCorrel2h[j][i] = NULL;
      //fErrorTermsSC3h[j][i] = NULL;
      //fErrorTermsAC41[j][i] = NULL;
    }
    fCorrel3h[i] = NULL;
  }
  fCentralityArray[16] = 0.;

  int harmo2h[13][2] = {{3,2}, {4,2}, {4,3}, {5,2}, {5,3}, {5,4}, {6,2},  // Lower: 7 first.
                        {6,3}, {6,4}, {7,2}, {7,3}, {8,2}, {8,3}};
  int powers[13][2] = { {1,1}, {2,0}, {2,1}, {3,0}, {3,1}, {4,0}, {4,1},  // For the 2-h AC.
                               {0,2}, {1,2}, {0,3}, {1,3}, {0,4}, {1,4}};
  for (int i = 0; i < 13; i++) {
    for (int j = 0; j < 2; j++) {
      fHarmoArray2h[i][j] = harmo2h[i][j];
      fPowers[i][j] = powers[i][j];
    }
  }

  int harmo3h[9][3] = { {2,3,4}, {2,3,5}, {2,4,5}, {2,3,6}, {2,3,7},  // Lower: 3 first.
                        {2,3,8}, {2,4,6}, {3,4,5}, {3,4,6}};
  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 3; j++) {fHarmoArray3h[i][j] = harmo3h[i][j];}
  }

  for (int iHarmo = 0; iHarmo < 81; iHarmo++) {
    for (int iPower = 0; iPower < 11; iPower++) {
      fQvectors[iHarmo][iPower] = TComplex(0., 0.);
    }
  }
  if (fDebugLevel > 5) {printf("'ArrayMembers' have been initialised.\n");}
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
int AliAnalysisTaskHOCFA::GetCentralityBin(float myCent)
{
// Get the bin corresponding to the passed centrality value.
  for (int i = 0; i < fNCentralityBins+1; i++) {
    if (myCent >= fCentralityArray[i]) {continue;}
    else {return i-1;}
  }
  return -1;  // Not a valid centrality value.
}

// ------------------------------------------------------------------------- //
void AliAnalysisTaskHOCFA::BookFinalResults()
{
// Book the needed objects for the analysis results.
  if (fDebugLevel > 5) {printf("AliAnalysisTaskHOCFA::BookFinalResults() reached.\n");}
  if (!fMainList) {Fatal("AliAnalysisTaskHOCFA::BookFinalResults()", "FATAL: fMainList not found.");}

  // Book the configuration profile.
  fHistoConfig = new TProfile("", "", 8, -0.5, 7.5);
  fHistoConfig->SetName("fHistoConfig");
  fHistoConfig->SetTitle("Configuration of the analysis");
  fHistoConfig->SetStats(kFALSE);
  fHistoConfig->GetXaxis()->SetBinLabel(1, "Min M");
  fHistoConfig->GetXaxis()->SetBinLabel(2, "Min pT");
  fHistoConfig->GetXaxis()->SetBinLabel(3, "Max pT");
  fHistoConfig->GetXaxis()->SetBinLabel(4, "#eta gap?");
  fHistoConfig->GetXaxis()->SetBinLabel(5, "NUE?");
  fHistoConfig->GetXaxis()->SetBinLabel(6, "NUA?");
  fHistoConfig->GetXaxis()->SetBinLabel(7, "SC?");
  fHistoConfig->GetXaxis()->SetBinLabel(8, "Lower 2-h?");
  fMainList->Add(fHistoConfig);

  for (int i = 0; i < fNCentralityBins; i++) {
    // Book the centrality lists only for the chosen range.
    fCentralityList[i] = new TList();
    fCentralityList[i]->SetName(Form("Centrality_%.1f-%.1f",
      fCentralityArray[i], fCentralityArray[i+1]));
    fCentralityList[i]->SetOwner(kTRUE);
    fMainList->Add(fCentralityList[i]);

    // Book the histograms in the correct lists.
    fHistoCent[i] = new TH1F("", "", 100, 0., 100.);
    fHistoCent[i]->SetName(Form("fHistoCent_Bin%d", i));
    fHistoCent[i]->SetTitle(Form("Centrality distribution, bin%d", i));
    fHistoCent[i]->SetStats(kTRUE);
    fCentralityList[i]->Add(fHistoCent[i]);

    fHistoCentCorrect[i] = new TH1F("", "", 100, 0., 100.);
    fHistoCentCorrect[i]->SetName(Form("fHistoCentCorrect_Bin%d", i));
    fHistoCentCorrect[i]->SetTitle(Form("Corrected centrality distribution, bin%d", i));
    fHistoCentCorrect[i]->SetStats(kTRUE);
    fCentralityList[i]->Add(fHistoCentCorrect[i]);

    fHistoMulti[i] = new TH1I("", "", 5000, 0., 5000.);
    fHistoMulti[i]->SetName(Form("fHistoMulti_Bin%d", i));
    fHistoMulti[i]->SetTitle(Form("Multiplicity distribution, bin%d", i));
    fHistoMulti[i]->SetStats(kTRUE);
    fCentralityList[i]->Add(fHistoMulti[i]);

    fHistoPt[i] = new TH1D("", "", 500, 0., 5.);
    fHistoPt[i]->SetName(Form("fHistoPt_Bin%d", i));
    fHistoPt[i]->SetTitle(Form("Corrected transverse momentum distribution, bin%d", i));
    fHistoPt[i]->SetStats(kTRUE);
    fCentralityList[i]->Add(fHistoPt[i]);

    fHistoEta[i] = new TH1D("", "", 200, -1., 1.);
    fHistoEta[i]->SetName(Form("fHistoEta_Bin%d", i));
    fHistoEta[i]->SetTitle(Form("Pseudorapidity distribution, bin%d", i));
    fHistoEta[i]->SetStats(kTRUE);
    fCentralityList[i]->Add(fHistoEta[i]);

    fHistoPhi[i] = new TH1D("", "", 630, -TMath::Pi(), TMath::Pi());
    fHistoPhi[i]->SetName(Form("fHistoPhi_Bin%d", i));
    fHistoPhi[i]->SetTitle(Form("Corrected azimuthal angle distribution, bin%d", i));
    fHistoPhi[i]->SetStats(kTRUE);
    fCentralityList[i]->Add(fHistoPhi[i]);

    fHistoCharge[i] = new TH1I("", "", 6, -3., 3.);
    fHistoCharge[i]->SetName(Form("fHistoCharge_Bin%d", i));
    fHistoCharge[i]->SetTitle(Form("Charge distribution, bin%d", i));
    fHistoCharge[i]->SetStats(kTRUE);
    fCentralityList[i]->Add(fHistoCharge[i]);

    // Book the profiles for the resulting correlators.
    fCorrel2p[i] = new TProfile("", "", 8, 0., 8.);
    fCorrel2p[i]->SetName(Form("fCorrel2p_Bin%d", i));
    fCorrel2p[i]->SetTitle(Form("2-p terms, no #Delta #eta, bin%d", i));
    fCorrel2p[i]->SetStats(kTRUE);
    fCorrel2p[i]->Sumw2();  // NEW (seems not needed but safer to have it)
    fCentralityList[i]->Add(fCorrel2p[i]);

    fCorrel2pEtaGap[i] = new TProfile("", "", 8, 0., 8.);
    fCorrel2pEtaGap[i]->SetName(Form("ffCorrel2pEtaGap_Bin%d", i));
    fCorrel2pEtaGap[i]->SetTitle(Form("2-p terms, |#it{#Delta #eta}| > %.1f, bin%d", fEtaGap, i));
    fCorrel2pEtaGap[i]->Sumw2();  // NEW (seems not needed but safer to have it)
    fCorrel2pEtaGap[i]->SetStats(kTRUE);
    if (fApplyEtaGap) {fCentralityList[i]->Add(fCorrel2pEtaGap[i]);}

    for (int j = 0; j < fNCombi; j++) {   // fNCombi = 1 (SC), 7 (lower AC) or 6 (higher AC).
      fCorrel2h[j][i] = new TProfile("", "", 13, 0., 13.);
      if (fGetSC) {
        fCorrel2h[j][i]->SetName(Form("fCorrel2h_allSC_Bin%d", i));
        fCorrel2h[j][i]->SetTitle(Form("2-h terms for SC, bin%d", i));
      } else {
        if (fGetLowerHarmos) {  // We measure only the 7 first 2-h terms.
          fCorrel2h[j][i]->SetName(Form("fCorrel2h_AC%d%d_Bin%d",
            fHarmoArray2h[j][0], fHarmoArray2h[j][1], i));
           fCorrel2h[j][i]->SetTitle(Form("2-h terms for AC(%d,%d), bin%d",
            fHarmoArray2h[j][0], fHarmoArray2h[j][1], i));
        } else {  // Else we measure only the 6 last 2-h terms.
          fCorrel2h[j][i]->SetName(Form("fCorrel2h_AC%d%d_Bin%d",
            fHarmoArray2h[j+7][0], fHarmoArray2h[j+7][1], i));
          fCorrel2h[j][i]->SetTitle(Form("2-h terms for AC(%d,%d), bin%d",
            fHarmoArray2h[j+7][0], fHarmoArray2h[j+7][1], i));
        }
      }   // The names and titles for the 2-h profile have been set.
      fCorrel2h[j][i]->Sumw2();
      fCentralityList[i]->Add(fCorrel2h[j][i]);
    } // Go to the next combination in the list.

    fCorrel3h[i] = new TProfile("", "", 9, 0., 9.);
    fCorrel3h[i]->SetName(Form("fCorrel3h_Bin%d", i));
    fCorrel3h[i]->SetTitle(Form("3-h terms for SC, bin%d", i));
    fCorrel3h[i]->Sumw2();  // NEW    
    fCorrel3h[i]->SetStats(kTRUE);
    if (fGetSC) {fCentralityList[i]->Add(fCorrel3h[i]);}

  /*
    fErrorTermsSC3h[j][i] = new TProfile("", "", 35, 0., 35.);
    fErrorTermsSC3h[j][i]->SetName(Form("fErrorTermsSC3h_Combi%d%d%d_Bin%d", fHarmoArray[j][0], fHarmoArray[j][1], fHarmoArray[j][2], i));
    fErrorTermsSC3h[j][i]->Sumw2();
    if (fGetSC) {fCentralityList[i]->Add(fErrorTermsSC3h[j][i]);}

    fErrorTermsAC41[j][i] = new TProfile("", "", 54, 0., 54.);
    fErrorTermsAC41[j][i]->SetName(Form("fErrorTermsAC41_Combi%d%d_Bin%d", fHarmoArray[j][0], fHarmoArray[j][1], i));
    fErrorTermsAC41[j][i]->Sumw2();
    if (!fGetSC) {fCentralityList[i]->Add(fErrorTermsAC41[j][i]);}
  */
    
  } // Go to the next centrality bin in the list.
}

// ------------------------------------------------------------------------- //
void AliAnalysisTaskHOCFA::CalculateQvectors(Long64_t myMulti,
  double myAngles[], double myWeights[])
{
// Calculate the Q-vectors needed for the analysis.
  double iAngle = 0.;           // Azimuthal angle of the current track.
  double iWeight = 0.;          // Particle weight of the current track.
  double iWeightToPowerP = 0.;  // Particle weight rised to the power p.

  // Ensure all the Q-vectors are initially zero.
  for (int iHarmo = 0; iHarmo < 81; iHarmo++) {
    for (int iPower = 0; iPower < 11; iPower++) {
      fQvectors[iHarmo][iPower] = TComplex(0., 0.);
    }
  }

  // Compute the Q-vectors themselves.
  for (Long64_t iTrack = 0; iTrack < myMulti; iTrack++) {
    iAngle = myAngles[iTrack];
    iWeight = myWeights[iTrack];

    for (int iHarmo = 0; iHarmo < 81; iHarmo++) {
      for (int iPower = 0; iPower < 11; iPower++) {
        iWeightToPowerP = TMath::Power(iWeight, iPower);
        fQvectors[iHarmo][iPower] += TComplex(iWeightToPowerP*TMath::Cos(iHarmo*iAngle),
          iWeightToPowerP*TMath::Sin(iHarmo*iAngle));
      }
    }
  } // Go to the next track.

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
void AliAnalysisTaskHOCFA::ComputeAllTerms(float myCentWeight)
{
// Compute all the terms needed for the ACs/SCs for all the observables.
  if (fDebugLevel > 5) {printf("Compute all the needed correlators.\n");}

  int nPart = 2;        // We start with the case of 2-particle correlations.
  int lHarmo[7] = {0};  // We work with "symmetric" correlators of max 14 particles > max 7 harmonics.
  int lPower[2] = {0};  // Powers of {vm2,vn2} for the AC case.

  // Fill the profile for the 2-p as it is common to all analysis cases.
  for (int iBin = 0; iBin < 8; iBin++){
    lHarmo[0] = iBin+1; // Only the first element needs to be filled in the 2-p case.
    CalculateCorrelator(nPart, lHarmo, fCorrel2p[fCentralityBin], iBin, lPower, myCentWeight);
  }

  // Fill the 2-harmonic profile according to the analysis configuration.
  // fNCombi = 1 (for SC, 13 combis in 1 profile), 7 (for lower AC) or 6 (for higher AC).
  for (int iProf = 0; iProf < fNCombi; iProf++) {
    for (int iBin = 0; iBin < 13; iBin++) {
      // Define the number of particles and the harmonics for AC/SC individually.
      if (fGetSC) {
        nPart = 4;  // 2-h terms in SC are always 4-particle correlators.
        lHarmo[0] = fHarmoArray2h[iBin][0];
        lHarmo[1] = fHarmoArray2h[iBin][1];
      } else {
        nPart = 0;  // Reset the number of particles for this bin.
        int cHarmo = 0; // Current index in the harmonic array.
        for (int iHarmo = 0; iHarmo < 2; iHarmo++) {
          if (fPowers[iBin][iHarmo] == 0) {continue;}   // Skip the unneeded harmonics.
          lPower[iHarmo] = fPowers[iBin][iHarmo];   // We change the power element only if non-zero.

          for (int jPower = 0; jPower < fPowers[iBin][iHarmo]; jPower++) {
            lHarmo[cHarmo] = fHarmoArray2h[iProf][iHarmo];   // Write the harmo as many times as its power is.
            cHarmo++;
          }

          nPart += 2*fPowers[iBin][iHarmo]; // 2-h terms in AC have twice as many particles as their cumulant order.
        } // Go to the next harmonic of the pair.
      } // We have the harmonic array and the number of particles at this point.
      if (fDebugLevel > 5) {printf("Harmo: (%d,%d), nPart: %d\n", lHarmo[0], lHarmo[1], nPart);}

      // Calculate the multiparticle correlator itself using the recursion method.
      CalculateCorrelator(nPart, lHarmo, fCorrel2h[iProf][fCentralityBin], iBin, lPower, myCentWeight);

      // Reset the variables for safety purposes.
      lPower[0] = 0.; lPower[1] = 0.;
    } // Go to the next bin in the current profile.
  } // Go to the next profile in the array.

  // Fill the 3-harmonic profile only if the analysis is configured for SC.
  if (fGetSC) {
    nPart = 6;  // 3-h terms in SC are always 6-particle correlators.
    for (int iBin = 0; iBin < 9; iBin++) {
      // Define the harmonic array for this bin.
      for (int iH = 0; iH < 3; iH++) {lHarmo[iH] = fHarmoArray3h[iBin][iH];}

      // Calculate the multiparticle correlator for this combination of harmonics.
      CalculateCorrelator(nPart, lHarmo, fCorrel3h[fCentralityBin], iBin, lPower, myCentWeight);
    } // Go to the next combination of 3 harmonics.
  }

  // TODO: Calculate the terms needed for the error propagation formulas.
}

// ------------------------------------------------------------------------- //
void AliAnalysisTaskHOCFA::CalculateCorrelator(int myMulti, int myHarmos[],
  TProfile *myProfile, int myBin, int myPowers[], float myCentWeight)
{
// Calculate the multiparticle correlator corresponding to harmonics[].
  if (fDebugLevel > 5) {printf("Calculate correlators for the provided harmonics.\n");}
  TComplex cCorrel = TComplex(0., 0.);
  double eventWeight = 0.;  // Event weight for this correlator.
  double rCorrel = 0.;      // Final value to fill in the profile.

  int numer2h[2] = {myHarmos[0], -1*myHarmos[0]};
  int denom2h[2] = {0};
  int numer4h[4] = {myHarmos[0], myHarmos[1],
                    -1*myHarmos[0], -1*myHarmos[1]};
  int denom4h[4] = {0};
  int numer6h[6] = {myHarmos[0], myHarmos[1], myHarmos[2],
                    -1*myHarmos[0], -1*myHarmos[1], -1*myHarmos[2]};
  int denom6h[6] = {0};
  int numer8h[8] = {myHarmos[0], myHarmos[1], myHarmos[2], myHarmos[3],
                    -1*myHarmos[0], -1*myHarmos[1], -1*myHarmos[2], -1*myHarmos[3]};
  int denom8h[8] = {0};
  int numer10h[10] = {myHarmos[0], myHarmos[1], myHarmos[2], myHarmos[3], myHarmos[4],
                    -1*myHarmos[0], -1*myHarmos[1], -1*myHarmos[2], -1*myHarmos[3], -1*myHarmos[4]};
  int denom10h[10] = {0};

  // Compute the denominator (= event weight) and numerator with the provided harmonics.
  switch(myMulti) {
  case 2 :
    eventWeight = ( CalculateRecursion(2, denom2h) ).Re();
    cCorrel = ( CalculateRecursion(2, numer2h) )/eventWeight;
    rCorrel = cCorrel.Re();
    break;
  case 4 :
    eventWeight = ( CalculateRecursion(4, denom4h) ).Re();
    cCorrel = ( CalculateRecursion(4, numer4h) )/eventWeight;
    rCorrel = cCorrel.Re();
    break;
  case 6 :    
    eventWeight = ( CalculateRecursion(6, denom6h) ).Re();
    cCorrel = ( CalculateRecursion(6, numer6h) )/eventWeight;
    rCorrel = cCorrel.Re();
    break;
  case 8 :
    eventWeight = ( CalculateRecursion(8, denom8h) ).Re();
    cCorrel = ( CalculateRecursion(8, numer8h) )/eventWeight;
    rCorrel = cCorrel.Re();
    break;
  case 10 :
    eventWeight = ( CalculateRecursion(10, denom10h) ).Re();
    cCorrel = ( CalculateRecursion(10, numer10h) )/eventWeight;
    rCorrel = cCorrel.Re();
    break;
  default :
    printf("Error: invalid number of particles.\n");
    break;
  }

  // Fill the corresponding bin in the right TProfile.
  myProfile->Fill( (float)myBin + 0.5, rCorrel, eventWeight*myCentWeight );

  if (myMulti == 2) {   // Bins are filled with v1-v8.
    myProfile->GetXaxis()->SetBinLabel(myBin+1, Form("v_{%d}", myBin+1));
  } else if (fGetSC) {
    if (myMulti == 4) {  // Bins are filled with 2-h combinations.
      myProfile->GetXaxis()->SetBinLabel(myBin+1, Form("(%d,%d)", myHarmos[0], myHarmos[1]));
    } else {
      myProfile->GetXaxis()->SetBinLabel(myBin+1, Form("(%d,%d,%d)", myHarmos[0], myHarmos[1], myHarmos[2]));
    }
  } else {  // Bins are filled with the power patterns.
    myProfile->GetXaxis()->SetBinLabel(myBin+1, Form("{%d,%d}", myPowers[0], myPowers[1]));
  }

/*
// Fill the temporary informations.
  errorTerms[0] = rCorrel;
  errorTerms[1] = eventWeight;
*/

// Reset the local variables for the next call.
  cCorrel = TComplex(0., 0.);
  rCorrel = 0.;
  eventWeight = 0.;
}

// ------------------------------------------------------------------------- //
void AliAnalysisTaskHOCFA::ComputeEtaGaps(Long64_t multiplicity,
  double angles[], double pWeights[], double pseudorapidity[], float myCentWeight)
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
    fCorrel2pEtaGap[fCentralityBin]->Fill((float)iHarmo + 0.5, rCorrel, Mminus[iHarmo]*Mplus[iHarmo]*myCentWeight);
    fCorrel2pEtaGap[fCentralityBin]->GetXaxis()->SetBinLabel(iHarmo + 1, Form("v_{%d}", iHarmo+1));

    // Reset the correlators to prevent mixing between harmonics.
    cCorrel = TComplex(0., 0.);
    rCorrel = 0.;
  }
// All the 2-particle correlators with eta gaps have been calculated.
}