
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

#include "AliAnalysisAnaTwoMultiCorrelations.h"
#include "Riostream.h"
#include "TList.h"
#include "TComplex.h"
#include "TFile.h"
#include "TH2I.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TProfile.h"
#include "TMath.h"
#include <vector>
#include <TExMap.h>
#include <TClonesArray.h>
#include "TDirectoryFile.h"
#include "TSystem.h"
#include "AliJBaseTrack.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisAnaTwoMultiCorrelations)

/* ========================================================================================== /
/ Mandatory methods needed for AliAnalysisTaskSE.                                             /
/ ========================================================================================== */
AliAnalysisAnaTwoMultiCorrelations::AliAnalysisAnaTwoMultiCorrelations() :
  fInputList(0),
// 1. General parameters for the configuration of the analysis.
  fACHarmoOne(-1), fACHarmoTwo(-1), fACHarmoThree(-1),
  fComputeACs(kFALSE), fWriteMinimum(kFALSE),
// 2. Parameters related to the centrality.
  fNumberBinsMulti(30000), fTotalCentralityBin(9),
  fCentrality(-1),
  fCentralityBin(-1), fInitialMultiplicity(0),
  fCentralityMin(0.), fCentralityMax(100.),
  fMultiplicityMin(6),
// 6. Parameters related to the efficiency and acceptance weights.
  fUseParticleWeights(kFALSE),
  fUsePtWeights(kFALSE), fUsePhiWeights(kFALSE), fUseEtaWeights(kFALSE),
  fEfficiency(NULL),
  fFirstEvent(kTRUE),
  fUseJEfficiency(kTRUE),
  fFilterbitIndex(0),
// 7. Parameters related to the multi-particle correlations.
  fHighestHarmonic(8), fLargestCorrelators(10),
  fReducedQPower(0),
  fMPCList(NULL),
// 8. Parameters related to the 2-particle correlations with eta gaps.
  fComputeEtaGaps(kFALSE),
  fTPCEtaList(NULL)
{
/* Dummy constructor of the class............................................................ /
/ 1. Initialise the arrays of data members.                                                  */
// 1. Initialise the arrays of data members.
  InitialiseArraysOfDataMembers();

} // End: AliAnalysisAnaTwoMultiCorrelations().

/* ----------------------------------------------------------------------------------------- */
AliAnalysisAnaTwoMultiCorrelations::AliAnalysisAnaTwoMultiCorrelations(const char *name, Bool_t useParticleWeights) :
// 1. General parameters for the configuration of the analysis.
  fACHarmoOne(-1), fACHarmoTwo(-1), fACHarmoThree(-1),
// 2. Parameters related to the centrality.
  fNumberBinsMulti(30000), fTotalCentralityBin(9),
  fCentrality(-1),
  fCentralityBin(-1), fInitialMultiplicity(0),
  fCentralityMin(0.), fCentralityMax(100.),
  fMultiplicityMin(6),
// 6. Parameters related to the efficiency and acceptance weights.
  fUseParticleWeights(kFALSE),
  fUsePtWeights(kFALSE), fUsePhiWeights(kFALSE), fUseEtaWeights(kFALSE),
  fEfficiency(NULL),
  fFirstEvent(kTRUE),
  fUseJEfficiency(kTRUE),
  fFilterbitIndex(0),
// 7. Parameters related to the multi-particle correlations.
  fHighestHarmonic(8), fLargestCorrelators(10),
  fReducedQPower(0),
  fMPCList(NULL),
// 8. Parameters related to the 2-particle correlations with eta gaps.
  fComputeEtaGaps(kFALSE),
  fTPCEtaList(NULL)
{
/* Constructor of the class.................................................................. /
/ 1. Create the mother list. (The rights on everything it holds are given to the list.)       /
/ 2. Define the input and output slots.                                                       /
/ 3. Initialise the arrays of data members.                                                  */
  fMPCList = new TList();
  fMPCList->SetName("MPC");
  fMPCList->SetOwner(kTRUE);

// 3. Initialise the arrays of data members.
  InitialiseArraysOfDataMembers();

} // End: AliAnalysisAnaTwoMultiCorrelations(const char *, Bool_t)

/* ----------------------------------------------------------------------------------------- */
AliAnalysisAnaTwoMultiCorrelations::~AliAnalysisAnaTwoMultiCorrelations()
{
/* Destructor of the class. */
  if (fMPCList) {delete fMPCList;}
} // End: ~AliAnalysisAnaTwoMultiCorrelations()

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisAnaTwoMultiCorrelations::UserCreateOutputObjects()
{
/* Define the outputs of the task at the beginning of the analysis........................... /
/ 1. Avoid name clashes. (Part 1.)                                                            /
/ 2. JEfficiency for NUA correction : DongJo (If chosen).                                     /
/ 3. Book the lists and their content.                                                        /
/ 4. 
/ 5. Avoid name clashes. (Part 2.)                                                           */
  TString sMethodName = "void Ali-TwoMultiCorrelations::UserCreateOutputObjects()";

// 1. Avoid name clashes.
  //Bool_t oldHistAddStatus = TH1::AddDirectoryStatus();
  //TH1::AddDirectory(kFALSE);

// 2. JEfficiency for NUA correction : DongJo
  //  Now it gets from ALiJCorrectionMap Task a.l.a AliJCatalyst

// 3. Book the lists and their content.
  this->BookAllLists();
// 5. Avoid name clashes.
  //TH1::AddDirectory(oldHistAddStatus);
  fFirstEvent = kTRUE;

}// End: void UserCreateOutputObjects().

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisAnaTwoMultiCorrelations::UserExec(Option_t *)
{
/* Execute the chosen analysis for each event of the input dataset........................... /                                                                     */
  TString sMethod = "void Ali-TwoMultiCorrelations::UserExec(Option_t *)";
  AnalyseRecoEvent();

} // End: void UserExec(Option_t *).

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisAnaTwoMultiCorrelations::Terminate(Option_t *)
{
/* Save the outputs at the end of the execution of the script................................ /
/ 1. Access the mother list.                                                                  /
/ 2. Create the output file and save inside the mother list.                                 */
  TString sMethod = "void Ali-TwoMultiCorrelations::Terminate(Option_t *)";

} // End: void Terminate(Option_t *).

/* ========================================================================================== /
/ Methods called in the constructors and setters.                                             /
/ ========================================================================================== */
/* ----------------------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisAnaTwoMultiCorrelations::InitialiseArraysOfDataMembers()
{
/* Initialise to zero all the elements belonging to arrays of data members. ---------------- */
// Q-vectors.
  for (Int_t iHarmo = 0; iHarmo < 81; iHarmo++)
  {
    for (Int_t iPower = 0; iPower < 11; iPower++)
    {
      fQvectors[iHarmo][iPower] = TComplex(0., 0.);
    } // End: iPower.
  } // End: iHarmo.

} // End: void InitialiseArraysOfDataMembers().

/* ========================================================================================== /
/ Methods called in "UserExec".                                                               /
/ ========================================================================================== */
void AliAnalysisAnaTwoMultiCorrelations::AnalyseRecoEvent()
{
/* Do the analysis at reco level to get the multiparticle correlations....................... /
/ 1. Check the presence of the pointer to an AOD event.                                       /
/ 2. Check if the event belongs to the current centrality bin.                                /
/ 3. Apply the event selection with the HMOs criteria.                                        /
/ 4. If chosen: get the run number and load the correct NUE correction.                       /
/ 5. Identify which tracks pass the selection.                                                /
/ 6. Keep the events with enough tracks to have a meaningful event weight.                    /
/ 7. Get the arrays needed to compute the multiparticle correlations.                         /
/ 8. Get the particle weights if chosen.                                                      /
/ 9. Compute all the possible multi-particle correlations.                                    /
/ 10. Reset the variables for the next event.                                                */
  TString sMethod = "void Ali-TwoMultiCorrelations::AnalyseRecoEvent()";
// 2. Get the centrality of this event.
  //fHistoNumberEvents[fCentralityBin]->Fill(0.5);

// 4. If chosen: get the run number and load the correct NUE correction.
  if (fFirstEvent)
  {
    fFirstEvent = kFALSE;
  } // End: if (fUseJEfficiency && fFirstEvent).

// 5. Identify which tracks pass the selection.
  long long finalMultiplicity = 0;    // Multiplicity after the track selection.

// 5.1 Get the current number of events.
  finalMultiplicity = fInputList->GetEntriesFast();
  fCentralityBin = GetCentralityBin(fCentrality); // from the task fTwoMultiAna->SetEventCentrality( fCent );

// 5.3 Save the number of events after the track selection.
  //fHistoNumberEvents[fCentralityBin]->Fill(2.5);

// 6. Keep the events with enough tracks to have a meaningful event weight.
//    Save the final number of tracks and final number of events.
  if (finalMultiplicity < fMultiplicityMin) {return;}
  //if (!fWriteMinimum) {fHistoMultiplicity[fCentralityBin][1]->Fill(finalMultiplicity);}
  //HistoNumberEvents[fCentralityBin]->Fill(3.5);

// 7. Get the arrays needed to compute the multiparticle correlations.
  Int_t iIndex = 0;           // Index of the selected track in the final arrays.
  Float_t iEffCorr = 1.;      // Efficiency correction.
  Float_t iEffInverse = 1.;   // Inverse of the efficiency correction.
  Float_t *iPt  = new Float_t[finalMultiplicity](); // Transverse momentum.
  Float_t *iEta = new Float_t[finalMultiplicity](); // Pseudorapidity for the eta gaps.
  Float_t *iPhi = new Float_t[finalMultiplicity](); // Azimuthal angles.
  Float_t *iParticleWeights = new Float_t[finalMultiplicity](); // Particle weights.

  for (Int_t iiTrack = 0; iiTrack < finalMultiplicity; iiTrack++)
  {
  // 7.2 Get a pointer to the selected track.
    AliJBaseTrack *aaTrack = (AliJBaseTrack*)fInputList->At(iiTrack); // load track
    if (!aaTrack) {continue;}
  // 7.3 Get all the needed variables.
    iPt[iIndex] = aaTrack->Pt();
    iEta[iIndex] = aaTrack->Eta();
    iPhi[iIndex] = aaTrack->Phi();
    iParticleWeights[iIndex] = 1.;
    if (fUseJEfficiency)
    {
      iEffCorr = aaTrack->GetTrackEff();//fEfficiency->GetCorrection( pt, fEffFilterBit, fCent);
      iEffInverse = 1.0/iEffCorr;
      Double_t phi_module_corr = aaTrack->GetWeight();// doing it in AliJCatalyst while filling track information.
      iParticleWeights[iIndex] = iEffInverse/phi_module_corr;
    } // End: if (fUseJEfficiency).

  // 7.4 Increase the index in the observables' arrays.
    iIndex++;
  } // End: for (Int_t iiTrack = 0; iiTrack < multiplicity; iiTrack++).

// 9. Compute all the possible multi-particle correlations.
  CalculateQvectors(finalMultiplicity, iPhi, iParticleWeights);
  ComputeReducedQvectors(finalMultiplicity);
  if (!fComputeACs) // Compute the correlators needed for SCs.
  ComputeSCsCorrelators(finalMultiplicity, iPhi, iParticleWeights);
  if (fComputeACs)  // Compute the correlators needed for ACs.
  {ComputeACsCorrelators(finalMultiplicity, iPhi, iParticleWeights);}
  if (fComputeEtaGaps)
  {ComputeTPCWithEtaGaps(finalMultiplicity, iPhi, iParticleWeights, iEta);}

// 10. Reset the variables for the next event.
  fInitialMultiplicity = 0;
  finalMultiplicity = 0;
  iIndex = 0;
  iEffCorr = 0.;
  iEffInverse = 0.;
  delete [] iPt;
  delete [] iEta;
  delete [] iPhi;
  delete [] iParticleWeights;

} // End: void AnalyseRecoEvent().

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisAnaTwoMultiCorrelations::CalculateQvectors(long long numberOfParticles, Float_t angles[], Float_t pWeights[])
{
/* Calculate all the Q-vectors for the given arrays of azimuthal angles and particle weights. */
  Float_t   iAngle            = 0.;   // Azimuthal angle of the current track.
  Float_t   iWeight           = 0.;   // Particle weight of the current track.
  Float_t   iWeightToPowerP   = 0.;   // Particle weight rised to the power p.

// Ensure all the Q-vectors are initially zero.
  for (Int_t iHarmo = 0; iHarmo < 81; iHarmo++)
  {
    for (Int_t iPower = 0; iPower < 11; iPower++)
    {
      fQvectors[iHarmo][iPower] = TComplex(0., 0.);
    }
  }

// Compute the Q-vectors.
  for (long long iTrack = 0; iTrack < numberOfParticles; iTrack++)
  {
    iAngle = angles[iTrack];
    iWeight = pWeights[iTrack];
    for (Int_t iHarmo = 0; iHarmo < 81; iHarmo++)
    {
      for (Int_t iPower = 0; iPower < 11; iPower++)
      {
        iWeightToPowerP = TMath::Power(iWeight, iPower);
        fQvectors[iHarmo][iPower] += TComplex(iWeightToPowerP*TMath::Cos(iHarmo*iAngle), iWeightToPowerP*TMath::Sin(iHarmo*iAngle));
      }   // End of the loop over the powers.
    }   // End of the loop over the harmonics.
  }   // End of the loop over the tracks.

// Reset the variables.
  iAngle            = 0.;
  iWeight           = 0.;
  iWeightToPowerP   = 0.;
}   // End of "void CalculateQvectors()".

/* ----------------------------------------------------------------------------------------- */
TComplex AliAnalysisAnaTwoMultiCorrelations::Q(Int_t n, Int_t p)
{
/* Simplify the use of the Q-vectors in the next methods. */
  if (n >= 0) {return fQvectors[n][p];}
  return TComplex::Conjugate(fQvectors[-n][p]);   // Use that Q*(n,p) = Q(-n,p).
}   // End of "TComplex Q()".

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisAnaTwoMultiCorrelations::ComputeReducedQvectors(long long numberOfParticles)
{
/* Compute the modulus of reduced Q-vector for the current event and save it in the correct histogram. */
  Double_t  reducedQ    = 0.;   // Modulus of reduced Q-vector.
  Int_t     iHarmo      = 0;    // Harmonic corresponding to the current histogram.

// Loop over the harmonics to get each q_n.
  for (Int_t i = 0; i < 8; i++)
  {
    iHarmo = i + 1;
    reducedQ = Q(iHarmo, fReducedQPower).Rho()/(TMath::Sqrt(1.*numberOfParticles));
    fHistoReducedQvectors[fCentralityBin][i]->Fill(reducedQ);
  }

// Reset the variables.
  reducedQ  = 0.;
  iHarmo    = 0;
}   // End of "void ComputeReducedQvectors()".

/* ----------------------------------------------------------------------------------------- */
TComplex AliAnalysisAnaTwoMultiCorrelations::CalculateRecursion(Int_t n, Int_t *harmonic, Int_t mult, Int_t skip)
{
/* Calculate the multi-particle correlators by using the recursion method.                  /
/ Improved, faster version originally developed by Kristjan Gulbrandsen (gulbrand@nbi.dk). */
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
}   // End of "TComplex CalculateRecursion()".

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisAnaTwoMultiCorrelations::ComputeSCsCorrelators(long long numberOfParticles, Float_t angles[], Float_t pWeights[])
{
/* Compute and save all the multiparticle correlators for the SCs for the current event. */
  Int_t twoZerosArray[2] = {0};     // 2d array with zeros for the denominator.
  Int_t twoPartArray[2] = {0};      // 2d array for the current harmonic.
  Int_t fourZerosArray[4] = {0};    // 4d array with zeros for the denominator.
  Int_t fourPartArray[4] = {0};     // 4d array for the current couple of harmonics.
  Int_t sixZerosArray[6] = {0};     // 6d array with zeros for the denominator.
  Int_t sixPartArray[6] = {0};      // 6d array for the current triplet of harmonics.
  Double_t eventWeight = 0.;        // Event weight and denominator of the correlator.
  TComplex complexCorrelator = TComplex(0., 0.);   // Complex form of the correlator.
  Double_t realCorrelator = 0.;     // Real part of the correlator.
  Int_t iBin = 1;                   // Bin for the current correlator.
  Double_t iMiddleBin = 0.;         // Middle of the bin for the current correlator.

// Compute the 2-particle correlators for v_1 to v_8.
  eventWeight = (CalculateRecursion(2, twoZerosArray)).Re();
  for (Int_t n = 1; n < 9; n++)   // Loop over the harmonics.
  {
  // Compute the correlator.
    twoPartArray[0]     = n;
    twoPartArray[1]     = -1*n;
    complexCorrelator   = (CalculateRecursion(2, twoPartArray))/eventWeight;
    realCorrelator      = complexCorrelator.Re();

  // Fill the corresponding bin in the profile.
    iMiddleBin          = (1.*n) - 0.5;
    fProfileTwoPartCorrel[fCentralityBin]->Fill(iMiddleBin, realCorrelator, eventWeight);
    fProfileTwoPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(n, Form("%d", n));

  // Reset the variables for the next harmonic.
    complexCorrelator   = TComplex(0., 0.);
    realCorrelator      = 0.;
    iMiddleBin          = 0.;
  }   // End of the loop over the harmonics for the 2-p correlators.

  eventWeight           = 0.;   // Reset of the denominator.

// Compute the 4-particle correlators.
  eventWeight           = (CalculateRecursion(4, fourZerosArray)).Re();
  for(Int_t m = 1; m < 9; m++)
  {
    for (Int_t n = m; n < 9; n++)
    {
    // Compute the correlator.
      fourPartArray[0]  = m;
      fourPartArray[1]  = n;
      fourPartArray[2]  = -1*m;
      fourPartArray[3]  = -1*n;
      complexCorrelator = (CalculateRecursion(4, fourPartArray))/eventWeight;
      realCorrelator    = complexCorrelator.Re();

    // Fill the corresponding bin in the profile.
      iMiddleBin        = (1.*iBin) - 0.5;
      fProfileFourPartCorrel[fCentralityBin]->Fill(iMiddleBin, realCorrelator, eventWeight);
      fProfileFourPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(iBin, Form("(%d,%d)", m, n));

    // Reset of the variables for the next pair of harmonics.
      iBin++;
      complexCorrelator = TComplex(0., 0.);
      realCorrelator    = 0.;
      iMiddleBin        = 0.;
    }   // End of the loop over the second harmonic for the 4-p correlators.
  }   // End of the loop over the first harmonic for the 4-p correlators.

  iBin                  = 1;
  if (!fWriteMinimum) // Compute the cross-check histograms only in full writing mode.
  {
    for(Int_t m = 2; m < 9; m++)
    {
      for (Int_t n = 1; n < m; n++)
      {
      // Compute the correlator.
        fourPartArray[0]  = m;
        fourPartArray[1]  = n;
        fourPartArray[2]  = -1*m;
        fourPartArray[3]  = -1*n;
        complexCorrelator = (CalculateRecursion(4, fourPartArray))/eventWeight;
        realCorrelator    = complexCorrelator.Re();

      // Fill the corresponding bin in the profile.
        iMiddleBin        = (1.*iBin) - 0.5;
        fProfileFourPartCorrelCheck[fCentralityBin]->Fill(iMiddleBin, realCorrelator, eventWeight);
        fProfileFourPartCorrelCheck[fCentralityBin]->GetXaxis()->SetBinLabel(iBin, Form("(%d,%d)", m, n));

      // Reset of the variables for the next pair of harmonics.
       iBin++;
       complexCorrelator  = TComplex(0., 0.);
       realCorrelator     = 0.;
       iMiddleBin         = 0.;
      }   // End of the loop over the second harmonic for cross-check.
    }   // End of the loop over the first harmonic for cross-check.
  } // End: if (!fWriteMinimum).
  eventWeight           = 0.;
  iBin                  = 1;

// Compute the needed 6-particle correlators.
  eventWeight           = (CalculateRecursion(6, sixZerosArray)).Re();
  for (Int_t l = 1; l < 7; l++)
  {
    for (Int_t m = 2; m < 8; m++)
    {
      if (l >= m) {continue;}   // Remove the cases in (l,l,n).
      for (Int_t n = 3; n < 9; n++)
      {
        if ((l >= n) || (m >= n)) {continue;}   // Remove the cases in (l,m,l) or (l,m,m).

      // Compute the correlator.
        sixPartArray[0] = l;
        sixPartArray[1] = m;
        sixPartArray[2] = n;
        sixPartArray[3] = -1*l;
        sixPartArray[4] = -1*m;
        sixPartArray[5] = -1*n;
        complexCorrelator = (CalculateRecursion(6, sixPartArray))/eventWeight;
        realCorrelator  = complexCorrelator.Re();

      // Fill the corresponding bin in the TProfile.
        iMiddleBin = (1.*iBin) - 0.5;
        fProfileSixPartCorrel[fCentralityBin]->Fill(iMiddleBin, realCorrelator, eventWeight);
        fProfileSixPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(iBin, Form("(%d,%d,%d)", l, m, n));

      // Reset of the variables for the next pair of harmonics.
        iBin++;
        complexCorrelator = TComplex(0., 0.);
        realCorrelator  = 0.;
        iMiddleBin      = 0.; 
      }   // End of the loop over the third harmonic of the triplet.
    }   // End of the loop over the second harmonic of the triplet.
  }   // End of the loop over for the first harmonic of the triplet.

// Reset the variables for the next event.
  eventWeight           = 0.;
  complexCorrelator     = TComplex(0., 0.);
  realCorrelator        = 0.;
  iBin                  = 1;
  iMiddleBin            = 0.;

} // End of "void ComputeSCsCorrelators()".

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisAnaTwoMultiCorrelations::ComputeACsCorrelators(long long numberOfParticles, Float_t angles[], Float_t pWeights[])
{
/* Compute and save all the multiparticle correlators for the ACs for the current event. */

// Compute the 2-particle correlations.
  Int_t twoZerosArray[2] = {0}; // Null 2d array for the event weight and denominator.
  Double_t eventWeight = (CalculateRecursion(2, twoZerosArray)).Re();
    // Event weight and denominator of the correlator.

/// <<2>>_{k,-k}, first bin.
  Int_t twoPartArray[2] = {fACHarmoOne, -1*fACHarmoOne};  // Array of harmonics.
  TComplex complexCorrelator = (CalculateRecursion(2, twoPartArray))/eventWeight;
    // Complex form of the multiparticle correlator.
  Double_t realCorrelator = complexCorrelator.Re();  // Real part of the multiparticle correlator.
  fProfileTwoPartCorrel[fCentralityBin]->Fill(0.5, realCorrelator, eventWeight);
  fProfileTwoPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(1, Form("%d", twoPartArray[0]));

  realCorrelator = 0.;  // Partial reset for the next 2-particle correlator.
  complexCorrelator = TComplex(0., 0.);

/// <<2>>_{l,-l}, second bin.
  twoPartArray[0] = fACHarmoTwo; twoPartArray[1] = -1*fACHarmoTwo;
  complexCorrelator = (CalculateRecursion(2, twoPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileTwoPartCorrel[fCentralityBin]->Fill(1.5, realCorrelator, eventWeight);
  fProfileTwoPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(2, Form("%d", twoPartArray[0]));

  realCorrelator = 0.;  // Partial reset for the next 2-particle correlator.
  complexCorrelator = TComplex(0., 0.);

/// <<2>>_{m,-m}, third bin.
  twoPartArray[0] = fACHarmoThree; twoPartArray[1] = -1*fACHarmoThree;
  complexCorrelator = (CalculateRecursion(2, twoPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileTwoPartCorrel[fCentralityBin]->Fill(2.5, realCorrelator, eventWeight);
  fProfileTwoPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(3, Form("%d", twoPartArray[0]));

// Full reset for the 4-particle correlators.
  eventWeight = 0.;
  realCorrelator = 0.;
  complexCorrelator = TComplex(0., 0.);

// Compute the 4-particle correlations.
  Int_t fourZerosArray[4] = {0};  // Null 4d array for the event weight and denominator.
  eventWeight = (CalculateRecursion(4, fourZerosArray)).Re();

/// <<4>>_{k,k,-k,-k}., first bin.
  Int_t fourPartArray[4] = {fACHarmoOne, fACHarmoOne, -1*fACHarmoOne, -1*fACHarmoOne};
  complexCorrelator = (CalculateRecursion(4, fourPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileFourPartCorrel[fCentralityBin]->Fill(0.5, realCorrelator, eventWeight);
  fProfileFourPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(1, Form("(%d,%d)", fourPartArray[0], fourPartArray[1]));

  realCorrelator = 0.;  // Partial reset for the next 4-particle correlator.
  complexCorrelator = TComplex(0., 0.);

/// <<4>>_{k,l,-k,-l}, second bin.
  fourPartArray[1] = fACHarmoTwo; fourPartArray[3] = -1*fACHarmoTwo;
  complexCorrelator = (CalculateRecursion(4, fourPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileFourPartCorrel[fCentralityBin]->Fill(1.5, realCorrelator, eventWeight);
  fProfileFourPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(2, Form("(%d,%d)", fourPartArray[0], fourPartArray[1]));

  realCorrelator = 0.;  // Partial reset for the next 4-particle correlator.
  complexCorrelator = TComplex(0., 0.);

/// <<4>>_{k,m,-k,-m}, third bin.
  fourPartArray[1] = fACHarmoThree; fourPartArray[3] = -1*fACHarmoThree;
  complexCorrelator = (CalculateRecursion(4, fourPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileFourPartCorrel[fCentralityBin]->Fill(2.5, realCorrelator, eventWeight);
  fProfileFourPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(3, Form("(%d,%d)", fourPartArray[0], fourPartArray[1]));

  realCorrelator = 0.;  // Partial reset for the next 4-particle correlator.
  complexCorrelator = TComplex(0., 0.);

/// <<4>>_{l,m,-l,-m}, fourth bin.
  fourPartArray[0] = fACHarmoTwo; fourPartArray[2] = -1*fACHarmoTwo;
  complexCorrelator = (CalculateRecursion(4, fourPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileFourPartCorrel[fCentralityBin]->Fill(3.5, realCorrelator, eventWeight);
  fProfileFourPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(4, Form("(%d,%d)", fourPartArray[0], fourPartArray[1]));

// Full reset for the 6-particle correlators.
  eventWeight = 0.;
  realCorrelator = 0.;
  complexCorrelator = TComplex(0., 0.);

// Compute the 6-particle correlations.
  Int_t sixZerosArray[6] = {0};   // Null 6d array for the event weight and denominator.
  eventWeight = (CalculateRecursion(6, sixZerosArray)).Re();

/// <<6>>_{k,k,k,-k,-k,-k}, first bin.
  Int_t sixPartArray[6] = {fACHarmoOne, fACHarmoOne, fACHarmoOne,
    -1*fACHarmoOne, -1*fACHarmoOne, -1*fACHarmoOne};
  complexCorrelator = (CalculateRecursion(6, sixPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileSixPartCorrel[fCentralityBin]->Fill(0.5, realCorrelator, eventWeight);
  fProfileSixPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(1, Form("(%d,%d,%d)", sixPartArray[0], sixPartArray[1], sixPartArray[2]));

  realCorrelator = 0.;  // Partial reset for the next 6-particle correlator.
  complexCorrelator = TComplex(0., 0.);

/// <<6>>_{k,k,l,-k,-k,-l}, second bin.
  sixPartArray[2] = fACHarmoTwo; sixPartArray[5] = -1*fACHarmoTwo;
  complexCorrelator = (CalculateRecursion(6, sixPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileSixPartCorrel[fCentralityBin]->Fill(1.5, realCorrelator, eventWeight);
  fProfileSixPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(2, Form("(%d,%d,%d)", sixPartArray[0], sixPartArray[1], sixPartArray[2]));

  realCorrelator = 0.;  // Partial reset for the next 6-particle correlator.
  complexCorrelator = TComplex(0., 0.);

/// <<6>>_{k,k,m,-k,-k,-m}, third bin.
  sixPartArray[2] = fACHarmoThree; sixPartArray[5] = -1*fACHarmoThree;
  complexCorrelator = (CalculateRecursion(6, sixPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileSixPartCorrel[fCentralityBin]->Fill(2.5, realCorrelator, eventWeight);
  fProfileSixPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(3, Form("(%d,%d,%d)", sixPartArray[0], sixPartArray[1], sixPartArray[2]));

  realCorrelator = 0.;  // Partial reset for the next 6-particle correlator.
  complexCorrelator = TComplex(0., 0.);

/// <<6>>_{k,l,m,-k,-l,-m}, fourth bin.
  sixPartArray[1] = fACHarmoTwo; sixPartArray[4] = -1*fACHarmoTwo;
  complexCorrelator = (CalculateRecursion(6, sixPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileSixPartCorrel[fCentralityBin]->Fill(3.5, realCorrelator, eventWeight);
  fProfileSixPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(4, Form("(%d,%d,%d)", sixPartArray[0], sixPartArray[1], sixPartArray[2]));

// Full reset for the 8-particle correlators.
  eventWeight = 0.;
  realCorrelator = 0.;
  complexCorrelator = TComplex(0., 0.);

// Compute the 8-particle correlations.
  Int_t eightZerosArray[8] = {0};   // Null 8d array for the event weight and denominator.
  eventWeight = (CalculateRecursion(8, eightZerosArray)).Re();

/// <<8>>_{k,k,k,k,-k,-k,-k,-k}, first bin.
  Int_t eightPartArray[8] = {fACHarmoOne, fACHarmoOne, fACHarmoOne, fACHarmoOne,
    -1*fACHarmoOne, -1*fACHarmoOne, -1*fACHarmoOne, -1*fACHarmoOne};
  complexCorrelator = (CalculateRecursion(8, eightPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileEightPartCorrel[fCentralityBin]->Fill(0.5, realCorrelator, eventWeight);
  fProfileEightPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(1, Form("(%d,%d,%d,%d)", eightPartArray[0], eightPartArray[1], eightPartArray[2], eightPartArray[3]));

  realCorrelator = 0.;  // Partial reset for the next 8-particle correlator.
  complexCorrelator = TComplex(0., 0.);

/// <<8>>_{k,k,k,l,-k,-k,-k,-l}, second bin.
  eightPartArray[3] = fACHarmoTwo; eightPartArray[7] = -1*fACHarmoTwo;
  complexCorrelator = (CalculateRecursion(8, eightPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileEightPartCorrel[fCentralityBin]->Fill(1.5, realCorrelator, eventWeight);
  fProfileEightPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(2, Form("(%d,%d,%d,%d)", eightPartArray[0], eightPartArray[1], eightPartArray[2], eightPartArray[3]));

  realCorrelator = 0.;  // Partial reset for the next 8-particle correlator.
  complexCorrelator = TComplex(0., 0.);

/// <<8>>_{k,k,l,m,-k,-k,-l,-m}, third bin.
  eightPartArray[2] = fACHarmoTwo; eightPartArray[6] = -1*fACHarmoTwo;
  eightPartArray[3] = fACHarmoThree; eightPartArray[7] = -1*fACHarmoThree;
  complexCorrelator = (CalculateRecursion(8, eightPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileEightPartCorrel[fCentralityBin]->Fill(2.5, realCorrelator, eventWeight);
  fProfileEightPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(3, Form("(%d,%d,%d,%d)", eightPartArray[0], eightPartArray[1], eightPartArray[2], eightPartArray[3]));

/// Full reset for the 10-particle correlators.
  eventWeight = 0.;
  realCorrelator = 0.;
  complexCorrelator = TComplex(0., 0.);

// Compute the 10-particle correlations.
  Int_t tenZerosArray[10] = {0};  // Null 10d array for the event weight and denominator.
  eventWeight = (CalculateRecursion(10, tenZerosArray)).Re();

/// <<10>>_{k,k,k,k,l,-k,-k,-k,-k,-l}, first bin.
  Int_t tenPartArray[10] = {fACHarmoOne, fACHarmoOne, fACHarmoOne, fACHarmoOne, fACHarmoTwo,
    -1*fACHarmoOne, -1*fACHarmoOne, -1*fACHarmoOne, -1*fACHarmoOne, -1*fACHarmoTwo};
  complexCorrelator = (CalculateRecursion(10, tenPartArray))/eventWeight;
  realCorrelator = complexCorrelator.Re();
  fProfileTenPartCorrel[fCentralityBin]->Fill(0.5, realCorrelator, eventWeight);
  fProfileTenPartCorrel[fCentralityBin]->GetXaxis()->SetBinLabel(1, Form("(%d,%d,%d,%d,%d)", tenPartArray[0], tenPartArray[1], tenPartArray[2], tenPartArray[3], tenPartArray[4]));

// Reset all variables for the next event.
  eventWeight = 0.;
  realCorrelator = 0.;
  complexCorrelator = TComplex(0., 0.);

} // End of "void ComputeACsCorrelators()".

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisAnaTwoMultiCorrelations::ComputeTPCWithEtaGaps(long long numberOfParticles, Float_t angles[], Float_t pWeights[], Float_t pseudorapidity[])
{
/* Compute the 2-particle correlators using eta gaps for the current event. */
  TComplex  Qminus[11][8]   = {{TComplex(0., 0.)}};   // Q-vectors for the negative subset of the eta range, for v_1 to v_8.
  TComplex  Qplus[11][8]    = {{TComplex(0., 0.)}};   // Q-vectors for the positive subset of the eta range, for v_1 to v_8.
  Float_t   Mminus[11][8]   = {{0.}};                 // Multiplicity in the negative subset of the eta range.
  Float_t   Mplus[11][8]    = {{0.}};                 // Multiplicity in the positive subset of the eta range.
  Float_t   etaGaps[11]     = {1., 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.};  // Possible values for the gaps.
  Float_t   iAngle          = 0.;                     // Azimuthal angle of the current particle.
  Float_t   iWeight         = 1.;                     // Particle weight of the current particle (default: unit weight).
  Float_t   iEta            = 0.;                     // Pseudorapidity of the current particle.
  Float_t   iWeightToP      = 1.;                     // Particle weight rised to the power p.
  TComplex  complexCorrel   = TComplex(0., 0.);       // Complex value of the 2-p correlator.
  Double_t  realCorrel      = 0.;                     // Real value of the 2-p correlator.

// Compute the Q-vectors for the negative and positive subsets of the eta range.
  for (Int_t iPart = 0; iPart < numberOfParticles; iPart++)
  {
  // Read the right elements in the provided arrays.
    iAngle  = angles[iPart];
    iWeight = pWeights[iPart];
    iEta    = pseudorapidity[iPart];
    if (fUseParticleWeights) {iWeightToP = iWeight;}   // All weights are multiplied to get the final one.

  // Compute the Q-vectors.
    if (iEta < 0.)    // Negative subset of the eta range.
    {
      for (Int_t iHarmo = 0; iHarmo < 8; iHarmo++)
      {
        for (Int_t iGap = 0; iGap < 11; iGap++)
        {
          if (iEta < ((-0.5)*etaGaps[iGap]))    // Compute only if the particle is in the range.
          {
            Qminus[iGap][iHarmo] += TComplex(iWeightToP*TMath::Cos((iHarmo+1)*iAngle), iWeightToP*TMath::Sin((iHarmo+1)*iAngle));
            Mminus[iGap][iHarmo] += iWeightToP;
          }
          else {continue;}
        }   // End of the loop over the gaps.
      }   // End of the loop over the harmonics.
    }   // End of the condition "negative subset".
    else if (iEta > 0.)   // Positive subset of the eta range.
    {
      for (Int_t iHarmo = 0; iHarmo < 8; iHarmo++)
      {
        for (Int_t iGap = 0; iGap < 11; iGap++)
        {
          if (iEta > (0.5*etaGaps[iGap]))   // Compute only if the particle is in the range.
          {
            Qplus[iGap][iHarmo] += TComplex(iWeightToP*TMath::Cos((iHarmo+1)*iAngle), iWeightToP*TMath::Sin((iHarmo+1)*iAngle));
            Mplus[iGap][iHarmo] += iWeightToP;
          }
        }   // End of the loop over the gaps.
      }   // End of the loop over the harmonics.
    }   // End of the condition "positive subset".
    else {continue;}    // Particle with iEta = 0.
  }   // End of the loop over the particles for the Q-vectors.

// Compute the 2-p correlators using Qminus and Qplus.
  for (Int_t iHarmo = 0; iHarmo < 8; iHarmo++)
  {
    for (Int_t iGap = 0; iGap < 11; iGap++)
    {
      if (!( (Qminus[iGap][iHarmo].TComplex::Rho() > 0.) && (Qplus[iGap][iHarmo].TComplex::Rho() > 0.) )) {continue;}
      if (!( (Mminus[iGap][iHarmo] > 0.) && (Mplus[iGap][iHarmo] > 0.) )) {continue;}

      complexCorrel = Qminus[iGap][iHarmo]*TComplex::Conjugate(Qplus[iGap][iHarmo]);
      realCorrel    = (complexCorrel.Re())/(Mminus[iGap][iHarmo]*Mplus[iGap][iHarmo]);
      fProfileTPCEta[fCentralityBin][iGap]->Fill(iHarmo + 0.5, realCorrel, Mminus[iGap][iHarmo]*Mplus[iGap][iHarmo]);

    // Reset the 2-particle correlator.
      complexCorrel = TComplex(0.,0.);
      realCorrel    = 0.;
    }   // End of the loop over the gaps.
  }   // End of the loop over the harmonics.

}   // End of "void ComputeTPCWithEtaGaps()".

/* ========================================================================================== /
/ Methods called in "Terminate".                                                              /
/ ========================================================================================== */
/* ----------------------------------------------------------------------------------------- */
/* ========================================================================================== /
/ Methods called in "UserCreateOutputObjects".                                                /
/ ========================================================================================== */
void AliAnalysisAnaTwoMultiCorrelations::BookAllLists()
{
/* Book all the lists in fMainList........................................................... /
// 1. Check if the mother list exists.                                                        /
// 2. Book the daughter list for the minimum QA histograms.                                   /
// 3. Book the daughter list for the multiplicity histograms.                                 /
// 4. Book the daughter list for the event QA histograms.                                     /
// 5. Book the daughter list for the track QA histograms.                                     /
// 6. Book the daughter list for the MPC histograms.                                          /
// 7. Book the daughter list for the eta gaps profiles.                                      */
  TString sMethodName = "void Ali-TwoMultiCorrelations::BookAllLists()";

// 1. Check if the mother list exists.
  if (!fMPCList) {Fatal(sMethodName.Data(), "ERROR: 'fMainList' does not exist.");}
  BookMPCList();
} // End: void BookAllLists().

/* ----------------------------------------------------------------------------------------- */
void AliAnalysisAnaTwoMultiCorrelations::BookMPCList()
{
/* Book the histograms and profiles for the MPC............................................. /
// 1. */
  TString sMethodName = "void Ali-TwoMultiCorrelations::BookMPCList()";

// Define the numbers of bins depending on if the code computes SCs or ACs.
// The default binnings are for SCs.
  Int_t numberBinsTwo = 8;
  Int_t numberBinsFour = 36;
  Int_t numberBinsSix = 56;
  if (fComputeACs) {numberBinsTwo = 3; numberBinsFour = 4; numberBinsSix = 4;}

// Declare all the profiles with their settings.
  for (Int_t iCent = 0; iCent < 9; iCent++)
  {
    if (iCent >= fTotalCentralityBin) {break;}  // Histograms outside of the analysis range.
    if (!fWriteMinimum) // Reduced q-vectors part of QA.
    {
    // Modulus of reduced Q-vectors.
      for (Int_t i = 0; i < 8; i++)
      {
        fHistoReducedQvectors[iCent][i] = new TH1F("", "", 1000, 0., 10.);
        fHistoReducedQvectors[iCent][i]->SetName(Form("fHistoReducedQvectors%d_Bin%d",
            i, iCent));
        fHistoReducedQvectors[iCent][i]->SetTitle(Form("Distribution of q_{%d}, Centrality bin %d", i+1, iCent));
        fHistoReducedQvectors[iCent][i]->SetStats(kTRUE);
        fHistoReducedQvectors[iCent][i]->GetXaxis()->SetTitle(Form("q_{%d}", i+1));
        fHistoReducedQvectors[iCent][i]->GetYaxis()->SetTitle(Form("dN/dq_{%d}", i+1));
        fMPCList->Add(fHistoReducedQvectors[iCent][i]);
      } // End: for (Int_t i = 0; i < 8; i++).
    } // End: if (!fWriteMinimum).

  // Profile with the 2-particle correlators.
    fProfileTwoPartCorrel[iCent] = new TProfile(Form("fProfileTwoPartCorrel_Bin%d", iCent),
      Form("#LT#LT2#GT#GT_{n,-n}, Centrality bin %d (no #eta gap)", iCent),
      numberBinsTwo, 0., numberBinsTwo);
    fProfileTwoPartCorrel[iCent]->Sumw2();
    fProfileTwoPartCorrel[iCent]->SetStats(kTRUE);
    fProfileTwoPartCorrel[iCent]->GetXaxis()->SetTitle("n");
    fMPCList->Add(fProfileTwoPartCorrel[iCent]);

  // Profiles with the 4-particle correlators.
    fProfileFourPartCorrel[iCent] = new TProfile(Form("fProfileFourPartCorrel_Bin%d", iCent),
      Form("#LT#LT4#GT#GT_{m,n,-m,-n}, Centrality bin %d", iCent), numberBinsFour, 0., numberBinsFour);
    fProfileFourPartCorrel[iCent]->Sumw2();
    fProfileFourPartCorrel[iCent]->SetStats(kTRUE);
    fProfileFourPartCorrel[iCent]->GetXaxis()->SetTitle("(m,n)");
    fMPCList->Add(fProfileFourPartCorrel[iCent]);

    if (!fWriteMinimum) // Cross-check is part of QA.
    {
      fProfileFourPartCorrelCheck[iCent] = new TProfile(
        Form("fProfileFourPartCorrelCheck_Bin%d", iCent),
        Form("#LT#LT4#GT#GT_{m,n,-m,-n}, Centrality bin %d", iCent),
          28, 0., 28.);
      fProfileFourPartCorrelCheck[iCent]->Sumw2();
      fProfileFourPartCorrelCheck[iCent]->SetStats(kTRUE);
      fProfileFourPartCorrelCheck[iCent]->GetXaxis()->SetTitle("(m,n)");
      fMPCList->Add(fProfileFourPartCorrelCheck[iCent]);
    } // End: if (!fWriteMinimum).

  // Profiles with the needed 6-particle correlators.
    fProfileSixPartCorrel[iCent] = new TProfile(Form("fProfileSixPartCorrel_Bin%d" ,iCent),
      Form("#LT#LT6#GT#GT_{m,n,o,-m,-n,-o}, Centrality bin %d", iCent),
      numberBinsSix, 0., numberBinsSix);
    fProfileSixPartCorrel[iCent]->Sumw2();
    fProfileSixPartCorrel[iCent]->SetStats(kTRUE);
    fMPCList->Add(fProfileSixPartCorrel[iCent]);

    if (fComputeACs)
    {
    // Profiles with the needed 8-particle correlators for ACs.
      fProfileEightPartCorrel[iCent] = new TProfile(Form("fProfileEightPartCorrel_Bin%d" ,iCent),
        Form("#LT#LT8#GT#GT_{m,n,o,p,-m,-n,-o,-p}, Centrality bin %d", iCent),
        3, 0., 3.);
      fProfileEightPartCorrel[iCent]->Sumw2();
      fProfileEightPartCorrel[iCent]->SetStats(kTRUE);
      fMPCList->Add(fProfileEightPartCorrel[iCent]);

    // Profiles with the needed 10-particle correlators for ACs.
      fProfileTenPartCorrel[iCent] = new TProfile(Form("fProfileTenPartCorrel_Bin%d" ,iCent),
        Form("#LT#LT10#GT#GT_{m,n,o,p,q,-m,-n,-o,-p,-q}, Centrality bin %d", iCent),
        1, 0., 1.);
      fProfileTenPartCorrel[iCent]->Sumw2();
      fProfileTenPartCorrel[iCent]->SetStats(kTRUE);
      fMPCList->Add(fProfileTenPartCorrel[iCent]);
    } // End: if (fComputeACs).
  } // End: iCent.
} // End: void BookMPCList().

// 5. Determine the index of the centrality bin for the current event.
Int_t AliAnalysisAnaTwoMultiCorrelations::GetCentralityBin(Float_t fFstCentrality) {
//    The analysis always uses the same division of centralities.
  if ( (fFstCentrality >= 0.) && (fFstCentrality < 5.) ) {fCentralityBin = 0;}
  if ( (fFstCentrality >= 5.) && (fFstCentrality < 10.) ) {fCentralityBin = 1;}
  if ( (fFstCentrality >= 10.) && (fFstCentrality < 20.) ) {fCentralityBin = 2;}
  if ( (fFstCentrality >= 20.) && (fFstCentrality < 30.) ) {fCentralityBin = 3;}
  if ( (fFstCentrality >= 30.) && (fFstCentrality < 40.) ) {fCentralityBin = 4;}
  if ( (fFstCentrality >= 40.) && (fFstCentrality < 50.) ) {fCentralityBin = 5;}
  if ( (fFstCentrality >= 50.) && (fFstCentrality < 60.) ) {fCentralityBin = 6;}
  if ( (fFstCentrality >= 60.) && (fFstCentrality < 70.) ) {fCentralityBin = 7;}
  if ( (fFstCentrality >= 70.) && (fFstCentrality < 80.) ) {fCentralityBin = 8;}

  return fCentralityBin;
}
