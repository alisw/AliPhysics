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

#include "AliAnalysisSPCRun2.h"

#include "Riostream.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisSPCRun2)

// ------------------------------------------------------------------------- //
AliAnalysisSPCRun2::AliAnalysisSPCRun2():
  fInputList(NULL),
  fHistList(NULL),
  fDebugLevel(0),
  fCentrality(0.), fcent_0(0.), fcent_1(0.), fcent_2(0.), fcent_3(0.), fcent_4(0.),
  fcent_5(0.), fcent_6(0.), fcent_7(0.), fcent_8(0.), fcent_9(0.),
  fcent_10(0.), fcent_11(0.), fcent_12(0.), fcent_13(0.), fcent_14(0.),
  fcent_15(0.), fcent_16(0.),
  fCentralityBins(16), fMinNumberPart(14),
  bUseWeightsNUE(kTRUE), bUseWeightsNUA(kTRUE),
  bComputeEtaGap(kFALSE), fEtaGap(0.8),
  bSaveAllQA(kTRUE),
  fCounterHistogram(NULL), fProfileTrackCuts(NULL)
{
// Dummy constructor of the class.
  printf("AliAnalysisSPCRun2::AliAnalysisSPCRun2()\n");
  InitializeArrays();
}

// ------------------------------------------------------------------------- //
AliAnalysisSPCRun2::AliAnalysisSPCRun2(const char *name):
  fInputList(NULL),
  fHistList(NULL),
  fDebugLevel(0),
  fCentrality(0.), fcent_0(0.), fcent_1(0.), fcent_2(0.), fcent_3(0.), fcent_4(0.),
  fcent_5(0.), fcent_6(0.), fcent_7(0.), fcent_8(0.), fcent_9(0.),
  fcent_10(0.), fcent_11(0.), fcent_12(0.), fcent_13(0.), fcent_14(0.),
  fcent_15(0.), fcent_16(0.),
  fCentralityBins(16), fMinNumberPart(14),
  bUseWeightsNUE(kTRUE), bUseWeightsNUA(kTRUE),
  bComputeEtaGap(kFALSE), fEtaGap(0.8),
  bSaveAllQA(kTRUE),
  fCounterHistogram(NULL), fProfileTrackCuts(NULL)
{
// Constructor of the class.
  printf("AliAnalysisSPCRun2::AliAnalysisSPCRun2(const char *name)\n");

  // Define the base list.
  fHistList = new TList();
  fHistList->SetName("outputStudentAnalysis");
  fHistList->SetOwner(kTRUE);

  // Initialize all arrays within the data members.
  InitializeArrays();

  //DefineOutput(1, TList::Class()); 
}


// ------------------------------------------------------------------------- //
AliAnalysisSPCRun2::~AliAnalysisSPCRun2()
{
// Destructor of the class.
  if (fHistList) {delete fHistList;}
}

// ------------------------------------------------------------------------- //
void AliAnalysisSPCRun2::UserCreateOutputObjects()
{
// Called at every worker node to declare the output objects.
  // Trick to avoid name clashes, part 1.
  // Not needed anymore, but kept in case.
  ///Bool_t oldHistAddStatus = TH1::AddDirectoryStatus(); 
  ///TH1::AddDirectory(kFALSE);

  // Book all the output lists and objects.
  BookAndNestAllLists();
  if (bSaveAllQA) {BookControlHistograms();}
  BookFinalResultsHistograms();

  // Fill the corresponding profile with the values of each cut. 
  fProfileTrackCuts->Fill(0.5, fMinNumberPart);
  if (bUseWeightsNUE) {fProfileTrackCuts->Fill(1.5, 1);} 
  if (bUseWeightsNUA) {fProfileTrackCuts->Fill(2.5, 1);} 

  // *) Trick to avoid name clashes, part 2.
  ///TH1::AddDirectory(oldHistAddStatus);
  ///PostData(1,fHistList);
}

// ------------------------------------------------------------------------- //
void AliAnalysisSPCRun2::UserExec(Option_t *option)
{
// Execute the analysis for each provided event in the input sample.
  // Fill the first bin of the counter histogram, to check if UserExec is entered properly.
  fCounterHistogram->Fill(0.5);  

  // Get the class bin corresponding to the centrality of the current event.
  Int_t centralityBin = SelectCentrality(fCentrality);

  // Start the analysis over the AODs.
  // 'DoMixed' and Fiser-Yates from the Run1 analysis task not included as not needed here.

  // First check if there are enough tracks, and fill the QA accordingly.
  Int_t nTracks = fInputList->GetEntriesFast();
  if (bSaveAllQA) {fMultHistogram[centralityBin]->Fill(nTracks);}
  if (nTracks < fMinNumberPart) {return;}
  fCentralityHistogram[centralityBin]->Fill(fCentrality);

  // Get the selected tracks information and store them into arrays for further use.
  Double_t* angles = new Double_t[nTracks];   // Azimuthal angles.
  Double_t* pt = new Double_t[nTracks];       // Transverse momenta.
  Double_t* eta = new Double_t[nTracks];      // Pseudorapidity.
  Double_t* weights = new Double_t[nTracks];  // Particle weights.

  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++)
  {
    Double_t iPhi = 0.;
    Double_t iPt = 0.;
    Double_t iEta = 0.;
    Double_t iCharge = 0.;
    Double_t iWeight = 1.;  // Default value: unit particle weight.

    // Parse the selected track and get their information into the arrays.
    AliJBaseTrack *aTrack = (AliJBaseTrack*)fInputList->At(iTrack);
    if (!aTrack) {continue;}

    iPhi = aTrack->Phi();
    iPt = aTrack->Pt();
    iEta = aTrack->Eta();
    iCharge = aTrack->GetCharge();

    Double_t iEffCorr = 1.;
    Double_t iEffInverse = 1.;
    Double_t phi_module_corr = 1.; // Doing it in AliJCatalyst while filling track information.

    if (bUseWeightsNUE) {
      iEffCorr = aTrack->GetTrackEff();
      iEffInverse = 1.0/iEffCorr;
    }
    if (bUseWeightsNUA) {phi_module_corr = aTrack->GetWeight();}
    //printf("iEffCorr: %.6f iPhiModuleCorr: %.6f \n", iEffCorr, phi_module_corr);
    iWeight = iEffInverse/phi_module_corr;

    angles[iTrack] = iPhi; 
    pt[iTrack] = iPt; 
    eta[iTrack] = iEta; 
    weights[iTrack] = iWeight;

    // Fill the track QA if needed.
    if (bSaveAllQA) {
      fPhiHistogram[centralityBin][0].Fill(iPhi, (1./phi_module_corr)); 
      fPhiWeightProfile[centralityBin]->Fill(iPhi,(1./phi_module_corr));
      fEtaHistogram[centralityBin][0].Fill(iEta);
      fPTHistogram[centralityBin][0].Fill(iPt, (1./iEffCorr));
      fChargeHistogram[centralityBin]->Fill(iCharge); 
    }

  } // Go to the next iTrack.

  // Calculate the correlators needed for the SPC.
  MainTask(centralityBin, nTracks, angles, weights);
  if (bComputeEtaGap) {ComputeTPCWithEtaGaps(centralityBin, nTracks, angles, weights, eta);}

  // Reset event-by-event objects.
  nTracks = 0;
  delete [] angles; 
  delete [] pt;
  delete [] eta;
  delete [] weights;
}

// ------------------------------------------------------------------------- //
void AliAnalysisSPCRun2::Terminate(Option_t *option)
{
// Execute the additional steps after the loop over the events is over.
}

// ------------------------------------------------------------------------- //
// ANALYSIS-SPECIFIC METHODS.                                                //
// ------------------------------------------------------------------------- //
void AliAnalysisSPCRun2::InitializeArrays()
{
// Initialise the arrays present in the data members.
  // Q-vectors.
  for (Int_t js = 0; js < 113; js++) 
  {
    for (Int_t j = 0; j < 15; j++) {fQvector[js][j] = TComplex(0.,0.);} 
  }

  // Harmonics combinations.
  for (Int_t i = 0; i < 12; i++)
  {
    for (Int_t j = 0; j < 8; j++) {fHarmosArray[i][j] = 0;}
  }

  //Centrality-dependent objects.
  for (Int_t icent = 0; icent < fCentralityBins; icent++)
  {
    // Lists.
    fCentralityList[icent] = NULL;
    fControlHistogramsList[icent] = NULL;
    fFinalResultsList[icent] = NULL;

    // QC histograms.
    fCentralityHistogram[icent] = NULL;
    fMultHistogram[icent] = NULL;
    fPTHistogram[icent] = NULL;
    fPhiHistogram[icent] = NULL;
    fEtaHistogram[icent] = NULL;
    fChargeHistogram[icent] = NULL;
    fPhiWeightProfile[icent] = NULL;

    //Output histograms.
    fResults[icent] = NULL;
    fResultsAlternativeError[icent] = NULL;
    fCovResults[icent] = NULL;
    fJoinedCovResults[icent] = NULL;
    fProfileTPCEta[icent] = NULL;
  } // Go to the next centrality bin.
}

// ------------------------------------------------------------------------- //
void AliAnalysisSPCRun2::BookAndNestAllLists()
{
// Book all the lists nested in the base list 'fHistList'.
  TString sMethodName = "void AliAnalysisSPCRun2::BookAndNestAllLists()";
  if (!fHistList) {Fatal(sMethodName.Data(),"fHistList is NULL");}

  for (Int_t icent = 0; icent < fCentralityBins; icent++)
  {
    // Check if the value of this centrality bin is negative. If yes: break.
    if (fCentralityArray[icent+1] < 0) {
      break; // The next edge is a breaking point -> This bin does not exist.
    }

    fCentralityList[icent] = new TList();
    fCentralityList[icent]->SetName(Form("MultCut_task=>%.1f-%.1f",
      fCentralityArray[icent],fCentralityArray[icent+1]));
    fCentralityList[icent]->SetOwner(kTRUE);
    fHistList->Add(fCentralityList[icent]);

    fControlHistogramsList[icent] = new TList();
    fControlHistogramsList[icent]->SetName("ControlHistograms");
    fControlHistogramsList[icent]->SetOwner(kTRUE);
    if (bSaveAllQA) {fCentralityList[icent]->Add(fControlHistogramsList[icent]);}
    
    fFinalResultsList[icent] = new TList();
    fFinalResultsList[icent]->SetName("FinalResults");
    fFinalResultsList[icent]->SetOwner(kTRUE);
    fCentralityList[icent]->Add(fFinalResultsList[icent]);
 }  // Go to the next centrality bin.

}

// ------------------------------------------------------------------------- //
void AliAnalysisSPCRun2::BookControlHistograms()
{
// Book all control histograms.
 for (Int_t icent = 0; icent < fCentralityBins; icent++)
 {
  // Check if the value of this centrality bin is negative. If yes: break.
  if (fCentralityArray[icent+1] < 0) {break;}

  // Transverse momentum spectrum.
  fPTHistogram[icent] = new TH1F("fPTHistAfterTrackSelection","Pt Distribution", 1000, 0., 10.);
  fPTHistogram[icent]->GetXaxis()->SetTitle("P_t");
  fPTHistogram[icent]->SetLineColor(4);
  fControlHistogramsList[icent]->Add(fPTHistogram[icent]); 
   
  // Azimuthal angle spectrum.
  fPhiHistogram[icent] = new TH1F("fPhiHistAfterTrackSelection","Phi Distribution",
    1000, -TMath::Pi(), TMath::Pi());
  fPhiHistogram[icent]->GetXaxis()->SetTitle("Phi");
  fPhiHistogram[icent]->SetLineColor(4);
  fControlHistogramsList[icent]->Add(fPhiHistogram[icent]);

  fPhiWeightProfile[icent] = new TProfile("fPhiWeightProfile","Phi Weights",
    100, -TMath::Pi(), TMath::Pi());
  fPhiWeightProfile[icent]->GetXaxis()->SetTitle("#varphi");
  fPhiWeightProfile[icent]->GetYaxis()->SetTitle("weight");
  fControlHistogramsList[icent]->Add(fPhiWeightProfile[icent]);

  // Pseudorapidity spectrum.
  fEtaHistogram[icent] = new TH1F("fEtaHistAfterTrackSelection","Eta Distribution", 1000,-1.,1.);
  fEtaHistogram[icent]->GetXaxis()->SetTitle("Eta");
  fEtaHistogram[icent]->SetLineColor(4);
  fControlHistogramsList[icent]->Add(fEtaHistogram[icent]);

  // Multiplicity spectrum.
  fMultHistogram[icent] = new TH1F("fMultiHistoAfterTrackSelection","Multiplicity", 30000,0.,30000.);
  fMultHistogram[icent]->GetXaxis()->SetTitle("Multiplicity M");
  fControlHistogramsList[icent]->Add(fMultHistogram[icent]);

  // Electric charge spectrum.
  fChargeHistogram[icent] = new TH1I("ChargeAfterCut","ChargeAfterCut",11, -5.5, 5.5); 
  fControlHistogramsList[icent]->Add(fChargeHistogram[icent]); 

  // Centrality distribution.
  fCentralityHistogram[icent]= new TH1F("fCentralityHistogramAfter","CentralityHistogramAfter",
    22, 0., 110.);
  fCentralityHistogram[icent]->GetXaxis()->SetTitle("Centrality");
  fCentralityHistogram[icent]->SetLineColor(4);
  fControlHistogramsList[icent]->Add(fCentralityHistogram[icent]);
  } // Go to the next centrality bin.

}

// ------------------------------------------------------------------------- //
void AliAnalysisSPCRun2::BookFinalResultsHistograms()
{
// Book all histograms to hold the final results.
  for (Int_t icent = 0; icent < fCentralityBins; icent++)
  {
    // Check if the value of this centrality bin is negative. If yes: break.
    if (fCentralityArray[icent+1] < 0) {break;}

    fResults[icent] = new TProfile("fResults", "Numerators and denominators for SPC",
      24, 0., 24.,"s");
    fResults[icent]->GetXaxis()->SetTitle("");
    fResults[icent]->GetYaxis()->SetTitle("");
    fResults[icent]->Sumw2();
    fFinalResultsList[icent]->Add(fResults[icent]);

    fResultsAlternativeError[icent] = new TProfile("fResultsAlternativeError",
      "Numerators and denominators for SPC (not 's' option)", 24, 0.,24.);
    fResultsAlternativeError[icent]->GetXaxis()->SetTitle("");
    fResultsAlternativeError[icent]->GetYaxis()->SetTitle("");
    fResultsAlternativeError[icent]->Sumw2();
    fFinalResultsList[icent]->Add(fResultsAlternativeError[icent]);

    fCovResults[icent] = new TProfile("fCovResults",
      "Results for the covariance terms", 48, 0.,48., "s");
    fCovResults[icent]->GetXaxis()->SetTitle("");
    fCovResults[icent]->GetYaxis()->SetTitle("");
    fCovResults[icent]->Sumw2();
    fFinalResultsList[icent]->Add(fCovResults[icent]); 

    fJoinedCovResults[icent] = new TProfile("fJoinedCovResults",
      "Joined Covariance calculated as one correlator <z> and not as <x*y>", 24, 0., 24.);
    fJoinedCovResults[icent]->GetXaxis()->SetTitle("");
    fJoinedCovResults[icent]->GetYaxis()->SetTitle("");
    fJoinedCovResults[icent]->Sumw2();
    fFinalResultsList[icent]->Add(fJoinedCovResults[icent]);

    fProfileTPCEta[icent] = new TProfile("fProfileTPCEta","2-particle correlators with eta gap",
      9, 0., 9.,"s"); 
    fProfileTPCEta[icent]->GetXaxis()->SetTitle("");
    fProfileTPCEta[icent]->GetYaxis()->SetTitle("");
    fProfileTPCEta[icent]->Sumw2();
    if (bComputeEtaGap) {fFinalResultsList[icent]->Add(fProfileTPCEta[icent]);}        
  } // Go to the next centrality bin.

  // Debugging histograms.
  fCounterHistogram = new TH1F("fCounterHistogram","Histogram for some checks", 3,0.,3.);
  fHistList->Add(fCounterHistogram);

  //Profile to save the cut values for track selection.
  fProfileTrackCuts = new TProfile("", "", 3, 0., 3.);
  fProfileTrackCuts->SetName("fProfileTrackCuts");
  fProfileTrackCuts->SetTitle("Configuration of the track selection");
  fProfileTrackCuts->SetStats(kFALSE);
  fProfileTrackCuts->GetXaxis()->SetBinLabel(1, "Multiplicity min");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(2, "NUE-Weights");
  fProfileTrackCuts->GetXaxis()->SetBinLabel(3, "NUA-Weights");
  fHistList->Add(fProfileTrackCuts);
}

// ------------------------------------------------------------------------- //
void AliAnalysisSPCRun2::CalculateQvectors(Int_t c_nParticles,
  Double_t* c_angles, Double_t* c_weights)
{
// Calculate Q-vectors. 
  // Ensure first that all Q-vectors are initially zero.
  for (Int_t h = 0; h < 113; h++)
  {
    for (Int_t p = 0; p < 15; p++) {fQvector[h][p] = TComplex(0.,0.);}
  } // Go to the next harmonic.

  // Calculate the Q-vectors for the given angles and weights. 
  Double_t dPhi2 = 0.;         // Particle azimuthal angle.
  Double_t wPhi = 1.;          // Particle weight.
  Double_t wPhiToPowerP = 1.;  // Particle weight raised to power p.

  for (Int_t i = 0; i < c_nParticles; i++)
  {
    dPhi2 = c_angles[i];
    if (bUseWeightsNUE || bUseWeightsNUA) {wPhi = c_weights[i];}
    for (Int_t h = 0; h < 113; h++)
    {
      for (Int_t p = 0; p < 15; p++)
      {
        if (bUseWeightsNUE || bUseWeightsNUA) {wPhiToPowerP = pow(wPhi,p);}
        fQvector[h][p] += TComplex(wPhiToPowerP*TMath::Cos(h*dPhi2),
          wPhiToPowerP*TMath::Sin(h*dPhi2));
      } // Go to the next power.
    } // Go to the next harmonic.
  }  // Go to the next particle.
}

// ------------------------------------------------------------------------- //
TComplex AliAnalysisSPCRun2::Q(Int_t n, Int_t p)
{
// Simplify the calculations using the fact that Q{-n,p} = Q{n,p}^*. 
  if ( n >= 0) {return fQvector[n][p];} 
  return TComplex::Conjugate(fQvector[-n][p]);
}

// ------------------------------------------------------------------------- //
TComplex AliAnalysisSPCRun2::Recursion(Int_t n, Int_t* harmonic, Int_t mult = 1, Int_t skip = 0)
{
// Calculate multi-particle correlators by using recursion (an improved faster version)
// originally developed by Kristjan Gulbrandsen (gulbrand@nbi.dk). 

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
}

// ------------------------------------------------------------------------- //
Int_t AliAnalysisSPCRun2::SelectCentrality(Double_t centValue)
{
// Get the centrality bin corresponding to the given centrality value.
  for (Int_t iCent = 0; iCent < fCentralityBins+1; iCent++)
  {
    if (fCentralityArray[iCent] < 0) {return -1;}
    if (centValue >= fCentralityArray[iCent]) {continue;} 
    else {return iCent-1;} 
  } 

  // We went through all centrality edges without returning. This means that the measured value
  // is bigger than the maximum centrality we want for our analysis.
  return -1;
}

// ------------------------------------------------------------------------- //
void AliAnalysisSPCRun2::SetInitializeCentralityArray()
{
// Fill the centrality array based on the content of the addTask.
  TString sMethodName = "void AliAnalysisSPC::BookAndNestAllLists()";

  Float_t listCentralities[17] = {
    fcent_0, fcent_1, fcent_2, fcent_3, fcent_4, fcent_5, fcent_6, fcent_7,
    fcent_8, fcent_9, fcent_10, fcent_11, fcent_12, fcent_13, fcent_14, fcent_15,
    fcent_16};
  for (Int_t i = 0; i < 17; i++) {fCentralityArray[i] = listCentralities[i];}

  // Protections.
  /// We need at least one well defined centrality bin.
  if (fCentralityArray[0] < 0 || fCentralityArray[1] < 0) {
    Fatal(sMethodName.Data(),"First Centrality bin not defined");
  }

  /// The next bin should be a valid boundary, i.e. > 0, but also smaller than the previous boundary
  /// -> Wrong ordering
  for (Int_t icent = 0; icent < fCentralityBins; icent++)
  {
    if (fCentralityArray[icent+1] > 0. && fCentralityArray[icent+1] < fCentralityArray[icent]) {
      Fatal(sMethodName.Data(),"Wrong ordering of centrality bounderies");
    }
  }

}

// ------------------------------------------------------------------------- //
void AliAnalysisSPCRun2::MainTask(Int_t centBin, Int_t mult,
  Double_t* m_angles, Double_t* m_weights)
{
// Calculate the SPC and error components for all the combinations of harmonics.
  if (fDebugLevel > 5) {printf("\n\nEntering now MainTask()\n");}

  if (mult >= fMinNumberPart) {
    // Calculate all the Q-vectors for this current event.
    CalculateQvectors(mult, m_angles, m_weights);

    // Declare the arrays to later fill all the needed bins for the correlators
    // and the error terms.
    Double_t *dataCorrelation = new Double_t[3];  // cosine, weight, sine.
    Double_t correlationNum;
    Double_t weightCorrelationNum;
    Double_t correlationDenom;
    Double_t weightCorrelationDenom;

    Double_t correlationJoinedCov;
    Double_t weightCorrelationJoinedCov;

    // Loop over the combinations of harmonics and calculate the corresponding SPC num and den.
    for (Int_t j = 0; j < 12; j++) {
      if (fHarmosArray[j][0] == 0) {continue;}  // Skip null correlator list.

      // Calculate the numerator.
      if (fDebugLevel > 5) {printf("Calculating numerator.\n");}
      Int_t hArrayNum[7] = {0};
      for (int iH = 0; iH < 7; iH++) {hArrayNum[iH] = fHarmosArray[j][iH+1];}
      Correlation(fHarmosArray[j][0], 7, hArrayNum, dataCorrelation);
      correlationNum = dataCorrelation[0];
      weightCorrelationNum = dataCorrelation[1];

      // Calculate the denominator.
      if (fDebugLevel > 5) {printf("Calculating denominator.\n");}
      Int_t nPartDen = 2*fHarmosArray[j][0];
      Int_t hArrayDen[14] = {0};
      for (int iH = 0; iH < 7; iH++) {
        hArrayDen[2*iH] = hArrayNum[iH];
        hArrayDen[2*iH+1] = -1*hArrayNum[iH];
      }
      Correlation(nPartDen, 14, hArrayDen, dataCorrelation);
      correlationDenom = dataCorrelation[0];
      weightCorrelationDenom = dataCorrelation[1];

      // Calculate the joined product of num*den as one single term.
      // Works only if the order of the correlator is <= 4
      if (fHarmosArray[j][0] <= 4) {
        if (fDebugLevel > 5) {printf("Calculating joined term.\n");}
        Int_t nPartJoinedCov = 3*fHarmosArray[j][0];
        Int_t hArrayJoinedCov[12] = {0};
        for (int o = 0; o < 4; o++) {hArrayJoinedCov[o] = fHarmosArray[j][o+1];}
        hArrayJoinedCov[4] = fHarmosArray[j][1];
        hArrayJoinedCov[5] = -1*fHarmosArray[j][1];
        hArrayJoinedCov[6] = fHarmosArray[j][2];
        hArrayJoinedCov[7] = -1*fHarmosArray[j][2];
        hArrayJoinedCov[8] = fHarmosArray[j][3];
        hArrayJoinedCov[9] = -1*fHarmosArray[j][3];
        hArrayJoinedCov[10] = fHarmosArray[j][4];
        hArrayJoinedCov[11] = -1*fHarmosArray[j][4];
        Correlation(nPartJoinedCov, 12, hArrayJoinedCov, dataCorrelation);  

        correlationJoinedCov = dataCorrelation[0];
        weightCorrelationJoinedCov = dataCorrelation[1];
      } // End: if (fHarmosArray[j][0] <= 4)

      // Fill the corresponding bins in the results profiles.
      fResults[centBin]->Fill(2.*(Float_t)(j)+0.5, correlationNum, weightCorrelationNum);
      fResults[centBin]->Fill(2.*(Float_t)(j)+1.5, correlationDenom, weightCorrelationDenom);
      fResultsAlternativeError[centBin]->Fill(2.*(Float_t)(j)+0.5,
        correlationNum, weightCorrelationNum);
      fResultsAlternativeError[centBin]->Fill(2.*(Float_t)(j)+1.5,
        correlationDenom, weightCorrelationDenom);

      fCovResults[centBin]->Fill(4.*(Float_t)(j)+0.5, correlationNum*correlationDenom,
        weightCorrelationNum*weightCorrelationDenom); // w_D*N*w_D*D
      fCovResults[centBin]->Fill(4.*(Float_t)(j)+1.5,
        weightCorrelationNum*weightCorrelationDenom, 1.); //w_N*w_D
      fCovResults[centBin]->Fill(4.*(Float_t)(j)+2.5, weightCorrelationNum, 1.); // w_N
      fCovResults[centBin]->Fill(4.*(Float_t)(j)+3.5, weightCorrelationDenom, 1.); // w_D 

      fJoinedCovResults[centBin]->Fill(2.*(Float_t)(j)+0.5, correlationJoinedCov,
        weightCorrelationJoinedCov); // Joined Cov Term z 
      fJoinedCovResults[centBin]->Fill(2.*(Float_t)(j)+1.5, weightCorrelationJoinedCov, 1.);

      correlationNum = 0;
      weightCorrelationNum = 0;
      correlationDenom = 0;
      weightCorrelationDenom = 0;
      correlationJoinedCov = 0;
      weightCorrelationJoinedCov = 0;

    } // End of the loop over the harmonic combinations (j).

    delete [] dataCorrelation;
  } // End: if (mult >= fMinNumberPart)

}

// ------------------------------------------------------------------------- //
void AliAnalysisSPCRun2::Correlation(Int_t c_nPart, Int_t c_nHarmo, Int_t* harmo, Double_t *correlData)
{
// Calculate the correlators for the provided set of harmonics using Q-vectors.
  if (fDebugLevel > 5) {for (Int_t i = 0; i < c_nHarmo; i++) {printf("%d ", harmo[i]);}}
  if (fDebugLevel > 5) {printf("\n");}

  // Protection against anisotropic correlators.
  Int_t sumHarmo = 0;
  for (Int_t i = 0; i < c_nHarmo; i++) {sumHarmo += harmo[i];}
  if (sumHarmo != 0) {printf("\nOups, this correlator is not isotropic(sum = %d). Bye\n", sumHarmo); return;}
  //printf("\nYeah, this correlator is isotropic(sum = %d).\n", sumHarmo);

  switch (c_nPart) {
  case 2:{
        Int_t harmonicsTwoNum[2] = {harmo[0], harmo[1]};     
        Int_t harmonicsTwoDen[2] = {0,0};      
        Double_t wTwoRecursion = Recursion(2,harmonicsTwoDen).Re(); 
        TComplex twoRecursion = Recursion(2,harmonicsTwoNum)/wTwoRecursion;
    
        correlData[0] = twoRecursion.Re();  // <cos(h1*phi1+h2*phi2)>
        correlData[1] = wTwoRecursion;      // weight
        correlData[2] = twoRecursion.Im();  // <sin(h1*phi1+h2*phi2)>
      }
    break;
  case 3:{
      Int_t harmonicsThreeNum[3] = {harmo[0],harmo[1],harmo[2]};       
      Int_t harmonicsThreeDen[3] = {0,0,0}; 
      Double_t wThreeRecursion = Recursion(3,harmonicsThreeDen).Re();      
      TComplex threeRecursion = Recursion(3,harmonicsThreeNum)/wThreeRecursion;
  
      correlData[0] = threeRecursion.Re(); // <cos(h1*phi1+h2*phi2+h3*phi3)>
      correlData[1] = wThreeRecursion;     // weight
      correlData[2] = threeRecursion.Im(); // <sin(h1*phi1+h2*phi2+h3*phi3)>
    }
    break;
  case 4:{
      Int_t harmonicsFourNum[4] = {harmo[0],harmo[1],harmo[2],harmo[3]};       
      Int_t harmonicsFourDen[4] = {0,0,0,0}; 
      Double_t wFourRecursion = Recursion(4,harmonicsFourDen).Re();      
      TComplex fourRecursion = Recursion(4,harmonicsFourNum)/wFourRecursion;
  
      correlData[0] = fourRecursion.Re(); // <cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4)>
      correlData[1] = wFourRecursion;     // weight
      correlData[2] = fourRecursion.Im(); // <sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4)>
    }
    break;
  case 5:{
      Int_t harmonicsFiveNum[5] = {harmo[0],harmo[1],harmo[2],harmo[3],harmo[4]};       
      Int_t harmonicsFiveDen[5] = {0,0,0,0,0};       
      Double_t wFiveRecursion = Recursion(5,harmonicsFiveDen).Re();
      TComplex fiveRecursion = Recursion(5,harmonicsFiveNum)/wFiveRecursion;
  
      correlData[0] = fiveRecursion.Re(); // <cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5)>
      correlData[1] = wFiveRecursion;     // weight
      correlData[2] = fiveRecursion.Im(); // <sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5)>
    }
    break;
  case 6:{
      Int_t harmonicsSixNum[6] = {harmo[0],harmo[1],harmo[2],harmo[3],harmo[4],harmo[5]};       
      Int_t harmonicsSixDen[6] = {0,0,0,0,0,0}; 
      Double_t wSixRecursion = Recursion(6,harmonicsSixDen).Re();      
      TComplex sixRecursion = Recursion(6,harmonicsSixNum)/wSixRecursion;
  
      correlData[0] = sixRecursion.Re(); // <cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6)>
      correlData[1] = wSixRecursion;     // weight
      correlData[2] = sixRecursion.Im(); // <sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6)>
    }
    break;
  case 7:{
      Int_t harmonicsSevenNum[7] = {harmo[0],harmo[1],harmo[2],harmo[3],harmo[4],harmo[5],harmo[6]};       
      Int_t harmonicsSevenDen[7] = {0,0,0,0,0,0,0};  
      Double_t wSevenRecursion = Recursion(7,harmonicsSevenDen).Re();     
      TComplex sevenRecursion = Recursion(7,harmonicsSevenNum)/wSevenRecursion;
  
      correlData[0] = sevenRecursion.Re(); // <cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7)>
      correlData[1] = wSevenRecursion;     // weight
      correlData[2] = sevenRecursion.Im(); // <sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7)>
    }
    break;
  case 8:{
      Int_t harmonicsEightNum[8] = {harmo[0],harmo[1],harmo[2],harmo[3],
        harmo[4],harmo[5],harmo[6],harmo[7]};       
      Int_t harmonicsEightDen[8] = {0,0,0,0,0,0,0,0};       
      Double_t wEightRecursion = Recursion(8,harmonicsEightDen).Re();
      TComplex eightRecursion = Recursion(8,harmonicsEightNum)/wEightRecursion;
  
      correlData[0] = eightRecursion.Re(); // <cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7+h8*phi8)>
      correlData[1] = wEightRecursion;     // weight
      correlData[2] = eightRecursion.Im(); // <sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6+h7*phi7+h8*phi8)>
    }
    break;    
  case 9:{
      Int_t harmonicsNineNum[9] = {harmo[0],harmo[1],harmo[2],harmo[3],harmo[4],
        harmo[5],harmo[6],harmo[7],harmo[8]};       
      Int_t harmonicsNineDen[9] = {0,0,0,0,0,0,0,0,0};    
      Double_t wNineRecursion = Recursion(9,harmonicsNineDen).Re();   
      TComplex nineRecursion = Recursion(9,harmonicsNineNum)/wNineRecursion;
  
      correlData[0] = nineRecursion.Re(); 
      correlData[1] = wNineRecursion;
      correlData[2] = nineRecursion.Im(); } 
    break;
  case 10:{
      Int_t harmonicsTenNum[10] = {harmo[0],harmo[1],harmo[2],harmo[3],harmo[4],
        harmo[5],harmo[6],harmo[7],harmo[8],harmo[9]};       
      Int_t harmonicsTenDen[10] = {0,0,0,0,0,0,0,0,0,0};  
      Double_t wTenRecursion = Recursion(10,harmonicsTenDen).Re();     
      TComplex tenRecursion = Recursion(10,harmonicsTenNum)/wTenRecursion;
  
      correlData[0] = tenRecursion.Re(); 
      correlData[1] = wTenRecursion;
      correlData[2] = tenRecursion.Im();}
    break;
  case 12:{
      Int_t harmonicsTwelveNum[12] = {harmo[0],harmo[1],harmo[2],harmo[3],harmo[4],harmo[5],
        harmo[6],harmo[7],harmo[8],harmo[9],harmo[10],harmo[11]};       
      Int_t harmonicsTwelveDen[12] = {0,0,0,0,0,0,0,0,0,0,0,0};     
      Double_t wTwelveRecursion = Recursion(12,harmonicsTwelveDen).Re();  
      TComplex twelveRecursion = Recursion(12,harmonicsTwelveNum)/wTwelveRecursion;
  
      correlData[0] = twelveRecursion.Re(); 
      correlData[1] = wTwelveRecursion;
      correlData[2] = twelveRecursion.Im();}
    break;
  case 14:{
      Int_t harmonicsFourteenNum[14] = {harmo[0],harmo[1],harmo[2],harmo[3],harmo[4],harmo[5],harmo[6],
        harmo[7],harmo[8],harmo[9],harmo[10],harmo[11],harmo[12],harmo[13]};       
      Int_t harmonicsFourteenDen[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0}; 
      Double_t wFourteenRecursion = Recursion(14,harmonicsFourteenDen).Re();      
      TComplex fourteenRecursion = Recursion(14,harmonicsFourteenNum)/wFourteenRecursion;
  
      correlData[0] = fourteenRecursion.Re(); 
      correlData[1] = wFourteenRecursion;
      correlData[2] = fourteenRecursion.Im();}
    break;
  default:
    return;
  }

}

// ------------------------------------------------------------------------- //
void AliAnalysisSPCRun2::ComputeTPCWithEtaGaps(Int_t centBin, Int_t mult, Double_t* m_angles,
  Double_t* m_weights, Double_t* m_pseudo)
{
// Compute the 2-particle correlators using eta gaps for the current event.
  TComplex  Qminus[8]   = {TComplex(0., 0.)};   // Q-vectors for the negative subset of the eta range, for v_1 to v_8.
  TComplex  Qplus[8]    = {TComplex(0., 0.)};   // Q-vectors for the positive subset of the eta range, for v_1 to v_8.
  Float_t   Mminus[8]   = {0.};                 // Multiplicity in the negative subset of the eta range.
  Float_t   Mplus[8]    = {0.};                 // Multiplicity in the positive subset of the eta range.
  Float_t   iAngle          = 0.;                     // Azimuthal angle of the current particle.
  Float_t   iWeight         = 1.;                     // Particle weight of the current particle (default: unit weight).
  Float_t   iEta            = 0.;                     // Pseudorapidity of the current particle.
  Float_t   iWeightToP      = 1.;                     // Particle weight rised to the power p.
  TComplex  complexCorrel   = TComplex(0., 0.);       // Complex value of the 2-p correlator.
  Double_t  realCorrel      = 0.;                     // Real value of the 2-p correlator.

  fProfileTPCEta[centBin]->Fill(8.5, fEtaGap,1.); //Fill Eta Gap for saving purpose

// Compute the Q-vectors for the negative and positive subsets of the eta range.
  for (Int_t iPart = 0; iPart < mult; iPart++)
  {
  // Read the right elements in the provided arrays.
    iAngle  = m_angles[iPart];
    iWeight = m_weights[iPart];
    iEta    = m_pseudo[iPart];
    if (bUseWeightsNUE || bUseWeightsNUA) {iWeightToP = iWeight;}   // All weights are multiplied to get the final one.

  // Compute the Q-vectors.
    if (iEta < 0.)    // Negative subset of the eta range.
    {
      for (Int_t iHarmo = 0; iHarmo < 8; iHarmo++)
      {
        if (iEta < ((-0.5)*fEtaGap))    // Compute only if the particle is in the range.
        {
          Qminus[iHarmo] += TComplex(iWeightToP*TMath::Cos((iHarmo+1)*iAngle), iWeightToP*TMath::Sin((iHarmo+1)*iAngle));
          Mminus[iHarmo] += iWeightToP;
        }
        else {continue;}
      }   // End of the loop over the harmonics.
    }   // End of the condition "negative subset".
    else if (iEta > 0.)   // Positive subset of the eta range.
    {
      for (Int_t iHarmo = 0; iHarmo < 8; iHarmo++)
      {
        if (iEta > (0.5*fEtaGap))   // Compute only if the particle is in the range.
        {
          Qplus[iHarmo] += TComplex(iWeightToP*TMath::Cos((iHarmo+1)*iAngle), iWeightToP*TMath::Sin((iHarmo+1)*iAngle));
          Mplus[iHarmo] += iWeightToP;
        } 
      }   // End of the loop over the harmonics.
    }   // End of the condition "positive subset".
    else {continue;}    // Particle with iEta = 0.
  }   // End of the loop over the particles for the Q-vectors.

// Compute the 2-p correlators using Qminus and Qplus.
  for (Int_t iHarmo = 0; iHarmo < 8; iHarmo++)
  {
    if (!( (Qminus[iHarmo].TComplex::Rho() > 0.) && (Qplus[iHarmo].TComplex::Rho() > 0.) )) {continue;}
    if (!( (Mminus[iHarmo] > 0.) && (Mplus[iHarmo] > 0.) )) {continue;}

    complexCorrel = Qminus[iHarmo]*TComplex::Conjugate(Qplus[iHarmo]);
    realCorrel    = (complexCorrel.Re())/(Mminus[iHarmo]*Mplus[iHarmo]);
    fProfileTPCEta[centBin]->Fill(iHarmo + 0.5, realCorrel, Mminus[iHarmo]*Mplus[iHarmo]); //GANESHA declare

    // Reset the 2-particle correlator.
    complexCorrel = TComplex(0.,0.);
    realCorrel    = 0.;    
  }   // End of the loop over the harmonics.

}
