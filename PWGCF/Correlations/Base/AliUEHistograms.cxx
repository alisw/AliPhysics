/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id: AliUEHistograms.cxx 20164 2007-08-14 15:31:50Z morsch $ */

//
//
// encapsulates several AliUEHist objects for a full UE analysis plus additional control histograms
//
//
// Author: Jan Fiete Grosse-Oetringhaus, Sara Vallero

#include "AliUEHistograms.h"

#include "AliCFContainer.h"
#include "AliBasicParticle.h"
#include "AliVParticle.h"
#include "AliAODTrack.h"

#include "TList.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TMath.h"
#include "TLorentzVector.h"

ClassImp(AliUEHistograms)

const Int_t AliUEHistograms::fgkUEHists = 3;

AliUEHistograms::AliUEHistograms(const char* name, const char* histograms, const char* binning) : 
  TNamed(name, name),
  fNumberDensitypT(0),
  fSumpT(0),
  fNumberDensityPhi(0),
  fCorrelationpT(0),
  fCorrelationEta(0),
  fCorrelationPhi(0),
  fCorrelationR(0),
  fCorrelationLeading2Phi(0),
  fCorrelationMultiplicity(0),
  fYields(0),
  fYieldsEtaPhiPT(0),
  fInvYield2(0),
  fEventCount(0),
  fEventCountDifferential(0),
  fVertexContributors(0),
  fCentralityDistribution(0),
  fCentralityCorrelation(0),
  fITSClusterMap(0),
  fControlConvResoncances(0),
  fEfficiencyCorrectionTriggers(0),
  fEfficiencyCorrectionAssociated(0),
  fSelectCharge(0),
  fTriggerSelectCharge(0),
  fAssociatedSelectCharge(0),
  fTriggerRestrictEta(-1),
  fEtaOrdering(kFALSE),
  fCutConversionsV(-1),
  fCutK0sV(-1),
  fCutLambdaV(-1),
  fCutPhiV(-1),
  fCutRhoV(-1),
  fCutCustomMass(-1),
  fCutCustomFirst(-1),
  fCutCustomSecond(-1),
  fCutCustomV(-1),
  fRejectResonanceDaughters(-1),
  fOnlyOneEtaSide(0),
  fOnlyOneAssocEtaSide(0),
  fWeightPerEvent(kFALSE),
  fPtOrder(kTRUE),
  fTwoTrackCutMinRadius(0.8),
  fCheckEventNumberInCorrelation(kFALSE),
  fRunNumber(0),
  fMergeCount(1)
{
  // Constructor
  //
  // the string histograms defines which histograms are created:
  //    1 = NumberDensitypT
  //    2 = SumpT
  //    3 = NumberDensityPhi
  //    4 = NumberDensityPhiCentrality (other multiplicity for Pb)
  
  AliLog::SetClassDebugLevel("AliCFContainer", -1);
  AliLog::SetClassDebugLevel("AliCFGridSparse", -3);

  fTwoTrackDistancePt[0] = 0;
  fTwoTrackDistancePt[1] = 0;
  
  TString histogramsStr(histograms);
  
  TString defaultBinningStr;
  defaultBinningStr = "eta: -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0\n"
    "p_t_assoc: 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 8.0\n"
    "p_t_leading: 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 24.0, 24.5, 25.0, 25.5, 26.0, 26.5, 27.0, 27.5, 28.0, 28.5, 29.0, 29.5, 30.0, 30.5, 31.0, 31.5, 32.0, 32.5, 33.0, 33.5, 34.0, 34.5, 35.0, 35.5, 36.0, 36.5, 37.0, 37.5, 38.0, 38.5, 39.0, 39.5, 40.0, 40.5, 41.0, 41.5, 42.0, 42.5, 43.0, 43.5, 44.0, 44.5, 45.0, 45.5, 46.0, 46.5, 47.0, 47.5, 48.0, 48.5, 49.0, 49.5, 50.0\n"
    "p_t_leading_course: 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0\n"
    "p_t_eff: 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0\n"
    "vertex_eff: -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10\n"
  ;

  if (histogramsStr.Contains("4") || histogramsStr.Contains("5") || histogramsStr.Contains("6")) // Dphi Corr
  {
    if (histogramsStr.Contains("C"))
      defaultBinningStr += "multiplicity: 0, 20, 40, 60, 80, 100.1\n"; // course
    else
      defaultBinningStr += "multiplicity: 0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100.1\n";

    if (histogramsStr.Contains("5"))
      defaultBinningStr += "vertex: -7, -5, -3, -1, 1, 3, 5, 7\n";
    else
      defaultBinningStr += "vertex: -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10\n";
    
    if (histogramsStr.Contains("R"))
      defaultBinningStr += "delta_phi: -1.570796, -1.483530, -1.396263, -1.308997, -1.221730, -1.134464, -1.047198, -0.959931, -0.872665, -0.785398, -0.698132, -0.610865, -0.523599, -0.436332, -0.349066, -0.261799, -0.174533, -0.087266, 0.0, 0.087266, 0.174533, 0.261799, 0.349066, 0.436332, 0.523599, 0.610865, 0.698132, 0.785398, 0.872665, 0.959931, 1.047198, 1.134464, 1.221730, 1.308997, 1.396263, 1.483530, 1.570796, 1.658063, 1.745329, 1.832596, 1.919862, 2.007129, 2.094395, 2.181662, 2.268928, 2.356194, 2.443461, 2.530727, 2.617994, 2.705260, 2.792527, 2.879793, 2.967060, 3.054326, 3.141593, 3.228859, 3.316126, 3.403392, 3.490659, 3.577925, 3.665191, 3.752458, 3.839724, 3.926991, 4.014257, 4.101524, 4.188790, 4.276057, 4.363323, 4.450590, 4.537856, 4.625123, 4.712389\n" // this binning starts at -pi/2 and is modulo 3 
	"delta_eta: -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,2.1, 2.2, 2.3, 2.4\n"
      ;
    else // for TTR studies
      defaultBinningStr += "delta_phi: -1.570796, -1.483530, -1.396263, -1.308997, -1.221730, -1.134464, -1.047198, -0.959931, -0.872665, -0.785398, -0.698132, -0.610865, -0.523599, -0.436332, -0.349066, -0.261799, -0.174533, -0.087266, -0.043633, -0.021817, 0.0, 0.021817, 0.043633, 0.087266, 0.174533, 0.261799, 0.349066, 0.436332, 0.523599, 0.610865, 0.698132, 0.785398, 0.872665, 0.959931, 1.047198, 1.134464, 1.221730, 1.308997, 1.396263, 1.483530, 1.570796, 1.658063, 1.745329, 1.832596, 1.919862, 2.007129, 2.094395, 2.181662, 2.268928, 2.356194, 2.443461, 2.530727, 2.617994, 2.705260, 2.792527, 2.879793, 2.967060, 3.054326, 3.141593, 3.228859, 3.316126, 3.403392, 3.490659, 3.577925, 3.665191, 3.752458, 3.839724, 3.926991, 4.014257, 4.101524, 4.188790, 4.276057, 4.363323, 4.450590, 4.537856, 4.625123, 4.712389\n"
	"delta_eta: -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, -0.05, -0.025, 0, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0\n"
      ;
  }
  else // UE
    defaultBinningStr += "multiplicity: -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 200.0\n"
      "delta_phi: -1.570796, -1.483530, -1.396263, -1.308997, -1.221730, -1.134464, -1.047198, -0.959931, -0.872665, -0.785398, -0.698132, -0.610865, -0.523599, -0.436332, -0.349066, -0.261799, -0.174533, -0.087266, 0.0, 0.087266, 0.174533, 0.261799, 0.349066, 0.436332, 0.523599, 0.610865, 0.698132, 0.785398, 0.872665, 0.959931, 1.047198, 1.134464, 1.221730, 1.308997, 1.396263, 1.483530, 1.570796, 1.658063, 1.745329, 1.832596, 1.919862, 2.007129, 2.094395, 2.181662, 2.268928, 2.356194, 2.443461, 2.530727, 2.617994, 2.705260, 2.792527, 2.879793, 2.967060, 3.054326, 3.141593, 3.228859, 3.316126, 3.403392, 3.490659, 3.577925, 3.665191, 3.752458, 3.839724, 3.926991, 4.014257, 4.101524, 4.188790, 4.276057, 4.363323, 4.450590, 4.537856, 4.625123, 4.712389\n"
      "vertex: -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10\n"
      "delta_eta: 0\n"
    ;

  // combine customBinning with defaultBinningStr -> use customBinning where available and otherwise defaultBinningStr
  TString binningStr = AliUEHist::CombineBinning(defaultBinningStr, TString(binning));

  if (histogramsStr.Contains("1"))
    fNumberDensitypT = new AliUEHist("NumberDensitypT", binningStr);
  if (histogramsStr.Contains("2"))
    fSumpT = new AliUEHist("SumpT", binningStr);
  
  if (histogramsStr.Contains("3"))
    fNumberDensityPhi = new AliUEHist("NumberDensityPhi", binningStr);
  else if (histogramsStr.Contains("4") || histogramsStr.Contains("5") || histogramsStr.Contains("6")) {
    TString configStr = "NumberDensityPhiCentrality";
    
    if (histogramsStr.Contains("5") || histogramsStr.Contains("6"))
      configStr += "Vtx";

    if (histogramsStr.Contains("S"))
      configStr += "Sparse";
    else if (histogramsStr.Contains("D"))
      configStr += "Double";
    
    fNumberDensityPhi = new AliUEHist(configStr, binningStr);
  }
  
  // do not add this hists to the directory
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  if (!histogramsStr.Contains("4") && !histogramsStr.Contains("5") && !histogramsStr.Contains("6"))
  {
    fCorrelationpT  = new TH2F("fCorrelationpT", ";p_{T,lead} (MC);p_{T,lead} (RECO)", 200, 0, 50, 200, 0, 50);
    fCorrelationEta = new TH2F("fCorrelationEta", ";#eta_{lead} (MC);#eta_{T,lead} (RECO)", 200, -1, 1, 200, -1, 1);
    fCorrelationPhi = new TH2F("fCorrelationPhi", ";#varphi_{lead} (MC);#varphi_{T,lead} (RECO)", 200, 0, TMath::TwoPi(), 200, 0, TMath::TwoPi());
  }
  else
  {
    fCorrelationpT  = new TH2F("fCorrelationpT", ";Centrality;p_{T} (RECO)", 100, 0, 100.001, 200, 0, 50);
    fCorrelationEta = new TH2F("fCorrelationEta", ";Centrality;#eta (RECO)", 100, 0, 100.001, 200, -5, 5);
    fCorrelationPhi = new TH2F("fCorrelationPhi", ";Centrality;#varphi (RECO)", 100, 0, 100.001, 200, 0, TMath::TwoPi());
  }
  
  fCorrelationR =   new TH2F("fCorrelationR", ";R;p_{T,lead} (MC)", 200, 0, 2, 200, 0, 50);
  fCorrelationLeading2Phi = new TH2F("fCorrelationLeading2Phi", ";#Delta #varphi;p_{T,lead} (MC)", 200, -TMath::Pi(), TMath::Pi(), 200, 0, 50);
  fCorrelationMultiplicity = new TH2F("fCorrelationMultiplicity", ";MC tracks;Reco tracks", 100, -0.5, 99.5, 100, -0.5, 99.5);
  fYields = new TH3F("fYields", ";centrality;pT;eta", 100, 0, 100, 40, 0, 20, 200, -10, 10);
  fYieldsEtaPhiPT = new TH3F("fYieldsEtaPhiPT", ";pT;eta;phi", 40,0,20,60,-3,3,200,0,TMath::TwoPi());
  fInvYield2 = new TH2F("fInvYield2", ";centrality;pT;1/pT dNch/dpT", 100, 0, 100, 80, 0, 20);

  if (!histogramsStr.Contains("4") && !histogramsStr.Contains("5") && !histogramsStr.Contains("6"))
  {
    fEventCount = new TH2F("fEventCount", ";step;event type;count", AliUEHist::fgkCFSteps+2, -2.5, -0.5 + AliUEHist::fgkCFSteps, 3, -0.5, 2.5);
    fEventCount->GetYaxis()->SetBinLabel(1, "ND");
    fEventCount->GetYaxis()->SetBinLabel(2, "SD");
    fEventCount->GetYaxis()->SetBinLabel(3, "DD");
  }
  else
  {
    fEventCount = new TH2F("fEventCount", ";step;centrality;count", AliUEHist::fgkCFSteps+2, -2.5, -0.5 + AliUEHist::fgkCFSteps, fNumberDensityPhi->GetEventHist()->GetNBins(1), fNumberDensityPhi->GetEventHist()->GetAxis(1, 0)->GetXbins()->GetArray());
  }
  
  fEventCountDifferential = new TH3F("fEventCountDifferential", ";p_{T,lead};step;event type", 100, 0, 50, AliUEHist::fgkCFSteps, -0.5, -0.5 + AliUEHist::fgkCFSteps, 3, -0.5, 2.5);
  fEventCountDifferential->GetZaxis()->SetBinLabel(1, "ND");
  fEventCountDifferential->GetZaxis()->SetBinLabel(2, "SD");
  fEventCountDifferential->GetZaxis()->SetBinLabel(3, "DD");
  
  fVertexContributors = new TH1F("fVertexContributors", ";contributors;count", 100, -0.5, 99.5);
  
  if (fNumberDensityPhi)
  {
    fCentralityDistribution = new TH1F("fCentralityDistribution", ";centrality;count", fNumberDensityPhi->GetEventHist()->GetNBins(1), fNumberDensityPhi->GetEventHist()->GetAxis(1, 0)->GetXbins()->GetArray());
    fCentralityCorrelation = new TH2F("fCentralityCorrelation", ";centrality;multiplicity", 404, 0, 101, 200, 0, 4000);
  }
  
  fITSClusterMap = new TH3F("fITSClusterMap", "; its cluster map; centrality; pT", 256, -0.5, 255.5, 20, 0, 100.001, 100, 0, 20);
  
  fControlConvResoncances = new TH2F("fControlConvResoncances", ";id;delta mass", 6, -0.5, 5.5, 500, -0.5, 0.5);
  
  TH1::AddDirectory(oldStatus);
}

//_____________________________________________________________________________
AliUEHistograms::AliUEHistograms(const AliUEHistograms &c) :
  TNamed(c.fName, c.fTitle),
  fNumberDensitypT(0),
  fSumpT(0),
  fNumberDensityPhi(0),
  fCorrelationpT(0),
  fCorrelationEta(0),
  fCorrelationPhi(0),
  fCorrelationR(0),
  fCorrelationLeading2Phi(0),
  fCorrelationMultiplicity(0),
  fYields(0),
  fYieldsEtaPhiPT(0),
  fInvYield2(0),
  fEventCount(0),
  fEventCountDifferential(0),
  fVertexContributors(0),
  fCentralityDistribution(0),
  fCentralityCorrelation(0),
  fITSClusterMap(0),
  fControlConvResoncances(0),
  fEfficiencyCorrectionTriggers(0),
  fEfficiencyCorrectionAssociated(0),
  fSelectCharge(0),
  fTriggerSelectCharge(0),
  fAssociatedSelectCharge(0),
  fTriggerRestrictEta(-1),
  fEtaOrdering(kFALSE),
  fCutConversionsV(-1),
  fCutK0sV(-1),
  fCutLambdaV(-1),
  fCutPhiV(-1),
  fCutRhoV(-1),
  fCutCustomMass(-1),
  fCutCustomFirst(-1),
  fCutCustomSecond(-1),
  fCutCustomV(-1),
  fRejectResonanceDaughters(-1),
  fOnlyOneEtaSide(0),
  fOnlyOneAssocEtaSide(0),
  fWeightPerEvent(kFALSE),
  fPtOrder(kTRUE),
  fTwoTrackCutMinRadius(0.8),
  fCheckEventNumberInCorrelation(kFALSE),
  fRunNumber(0),
  fMergeCount(1)
{
  //
  // AliUEHistograms copy constructor
  //

  fTwoTrackDistancePt[0] = 0;
  fTwoTrackDistancePt[1] = 0;

  ((AliUEHistograms &) c).Copy(*this);
}

//____________________________________________________________________
AliUEHistograms::~AliUEHistograms()
{
  // Destructor
  
  DeleteContainers();
}

void AliUEHistograms::DeleteContainers()
{
  if (fNumberDensitypT)
  {
    delete fNumberDensitypT;
    fNumberDensitypT = 0;
  }
  
  if (fSumpT)
  {
    delete fSumpT;
    fSumpT = 0;
  }
  
  if (fNumberDensityPhi)
  {
    delete fNumberDensityPhi;
    fNumberDensityPhi = 0;
  }
  
  if (fCorrelationpT)
  {
    delete fCorrelationpT;
    fCorrelationpT = 0;
  }
  
  if (fCorrelationEta)
  {
    delete fCorrelationEta;
    fCorrelationEta = 0;
  }
  
  if (fCorrelationPhi)
  {
    delete fCorrelationPhi;
    fCorrelationPhi = 0;
  }
  
  if (fCorrelationR)
  {
    delete fCorrelationR;
    fCorrelationR = 0;
  }

  if (fCorrelationLeading2Phi)
  {
    delete fCorrelationLeading2Phi;
    fCorrelationLeading2Phi = 0;
  }
  
  if (fCorrelationMultiplicity)
  {
    delete fCorrelationMultiplicity;
    fCorrelationMultiplicity = 0;
  }
  
  if (fYields)
  {
    delete fYields;
    fYields = 0;
  }

  if (fYieldsEtaPhiPT)
  {
    delete fYieldsEtaPhiPT;
    fYieldsEtaPhiPT = 0;
  }
  
  if (fInvYield2)
  {
    delete fInvYield2;
    fInvYield2 = 0;
  }
  
  if (fEventCount)
  {
    delete fEventCount;
    fEventCount = 0;
  }

  if (fEventCountDifferential)
  {
    delete fEventCountDifferential;
    fEventCountDifferential = 0;
  }
  
  if (fVertexContributors)
  {
    delete fVertexContributors;
    fVertexContributors = 0;
  }
  
  if (fCentralityDistribution)
  {
    delete fCentralityDistribution;
    fCentralityDistribution = 0;
  }
  
  if (fCentralityCorrelation)
  {
    delete fCentralityCorrelation;
    fCentralityCorrelation = 0;
  }
  
  if (fITSClusterMap)
  {
    delete fITSClusterMap;
    fITSClusterMap = 0;
  }
  
  for (Int_t i=0; i<2; i++)
    if (fTwoTrackDistancePt[i])
    {
      delete fTwoTrackDistancePt[i];
      fTwoTrackDistancePt[i] = 0;
    }
    
  if (fControlConvResoncances)
  {
    delete fControlConvResoncances;
    fControlConvResoncances = 0;
  }
    
  if (fEfficiencyCorrectionTriggers)
  {
    if (fEfficiencyCorrectionTriggers == fEfficiencyCorrectionAssociated)
      fEfficiencyCorrectionAssociated = 0;
    delete fEfficiencyCorrectionTriggers;
    fEfficiencyCorrectionTriggers = 0;
  }
  
  if (fEfficiencyCorrectionAssociated)
  {
    delete fEfficiencyCorrectionAssociated;
    fEfficiencyCorrectionAssociated = 0;
  }
}

AliUEHist* AliUEHistograms::GetUEHist(Int_t id)
{
  // returns AliUEHist object, useful for loops
  
  switch (id)
  {
    case 0: return fNumberDensitypT; break;
    case 1: return fSumpT; break;
    case 2: return fNumberDensityPhi; break;
  }
  
  return 0;
}

//____________________________________________________________________
Int_t AliUEHistograms::CountParticles(TList* list, Float_t ptMin)
{
  // counts the number of particles in the list with a pT above ptMin
  // TODO eta cut needed here?
  
  Int_t count = 0;
  for (Int_t j=0; j<list->GetEntries(); j++)
    if (((AliVParticle*) list->At(j))->Pt() > ptMin)
      count++;
      
  return count;
}
  
//____________________________________________________________________
void AliUEHistograms::Fill(Int_t eventType, Float_t zVtx, AliUEHist::CFStep step, AliVParticle* leading, TList* toward, TList* away, TList* min, TList* max)
{
  // fills the UE event histograms
  //
  // this function needs the leading (track or jet or ...) and four lists of AliVParticles which contain the particles/tracks to be filled in the four regions
  
  // if leading is not set, just fill event statistics
  if (leading)
  {
    Int_t multiplicity = 0;
    
    // TODO configurable?
    Float_t ptMin = 0.15;
    if (leading->Pt() > ptMin)
      multiplicity++;
    
    multiplicity += CountParticles(toward, ptMin);
    multiplicity += CountParticles(away, ptMin);
    multiplicity += CountParticles(min, ptMin);
    multiplicity += CountParticles(max, ptMin);
     
    FillRegion(AliUEHist::kToward, zVtx, step, leading, toward, multiplicity);
    FillRegion(AliUEHist::kAway,   zVtx, step, leading, away, multiplicity);
    FillRegion(AliUEHist::kMin,    zVtx, step, leading, min, multiplicity);
    FillRegion(AliUEHist::kMax,    zVtx, step, leading, max, multiplicity);
 
    Double_t vars[3];
    vars[0] = leading->Pt();
    vars[1] = multiplicity;
    vars[2] = zVtx;
    for (Int_t i=0; i<fgkUEHists; i++)
      if (GetUEHist(i))
        GetUEHist(i)->GetEventHist()->Fill(vars, step);
  
    fEventCountDifferential->Fill(leading->Pt(), step, eventType);
  }
  
  FillEvent(eventType, step);
}
  
//____________________________________________________________________
void AliUEHistograms::FillRegion(AliUEHist::Region region, Float_t zVtx, AliUEHist::CFStep step, AliVParticle* leading, TList* list, Int_t multiplicity)
{
  // loops over AliVParticles in list and fills the given region at the given step
  //
  // See also Fill(...)

  for (Int_t i=0; i<list->GetEntries(); i++)
  {
    AliVParticle* particle = (AliVParticle*) list->At(i);
    
    Double_t vars[6];
    vars[0] = particle->Eta();
    vars[1] = particle->Pt();
    vars[2] = leading->Pt();
    vars[3] = multiplicity;
    vars[4] = leading->Phi() - particle->Phi();
    if (vars[4] > 1.5 * TMath::Pi()) 
      vars[4] -= TMath::TwoPi();
    if (vars[4] < -0.5 * TMath::Pi())
      vars[4] += TMath::TwoPi();
    vars[5] = zVtx;
    
    if (fNumberDensitypT)
      fNumberDensitypT->GetTrackHist(region)->Fill(vars, step);
      
    if (fSumpT)
      fSumpT->GetTrackHist(region)->Fill(vars, step, particle->Pt());
    
    // fill all in toward region (is anyway as function of delta phi!)
    if (fNumberDensityPhi)
      fNumberDensityPhi->GetTrackHist(AliUEHist::kToward)->Fill(vars, step);
  }
}

//____________________________________________________________________
void AliUEHistograms::Fill(AliVParticle* leadingMC, AliVParticle* leadingReco)
{
  // fills the correlation histograms
  
  if (leadingMC)
  {
    fCorrelationpT->Fill(leadingMC->Pt(), leadingReco->Pt());
    if (leadingMC->Pt() > 0.5)
    {
      fCorrelationEta->Fill(leadingMC->Eta(), leadingReco->Eta());
      fCorrelationPhi->Fill(leadingMC->Phi(), leadingReco->Phi());
    }
    
    Float_t phiDiff = leadingMC->Phi() - leadingReco->Phi();
    if (phiDiff > TMath::Pi())
      phiDiff -= TMath::TwoPi();
    if (phiDiff < -TMath::Pi())
      phiDiff += TMath::TwoPi();
      
    Float_t etaDiff = leadingMC->Eta() - leadingReco->Eta();
    
    fCorrelationR->Fill(TMath::Sqrt(phiDiff * phiDiff + etaDiff * etaDiff), leadingMC->Pt());
    fCorrelationLeading2Phi->Fill(phiDiff, leadingMC->Pt());
  }
  else
  {
    fCorrelationpT->Fill(1.0, leadingReco->Pt());
    if (leadingReco->Pt() > 0.5)
    {
      fCorrelationEta->Fill(0.0, leadingReco->Eta());
      fCorrelationPhi->Fill(0.0, leadingReco->Phi());
    }
  }
}

//____________________________________________________________________
void AliUEHistograms::FillCorrelations(Double_t centrality, Float_t zVtx, AliUEHist::CFStep step, TObjArray* particles, TObjArray* mixed, Float_t weight, Bool_t firstTime, Bool_t twoTrackCuts, Float_t bSign, Float_t twoTrackEfficiencyCutValue, Bool_t applyEfficiency)
{
  // fills the fNumberDensityPhi histogram
  //
  // this function need a list of AliVParticles which contain the particles/tracks to be filled
  //
  // if mixed is non-0, mixed events are filled, the trigger particle is from particles, the associated from mixed
  // if weight < 0, then the pt of the associated particle is filled as weight
  
  Bool_t fillpT = kFALSE;
  if (weight < 0)
    fillpT = kTRUE;
  
  if (twoTrackCuts && (twoTrackEfficiencyCutValue > 0) && !fTwoTrackDistancePt[0])
  {
    // do not add this hists to the directory
    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);

    fTwoTrackDistancePt[0] = new TH3F("fTwoTrackDistancePt[0]", ";#Delta#eta;#Delta#varphi^{*}_{min};#Delta p_{T}", 100, -0.15, 0.15, 100, -0.05, 0.05, 20, 0, 10);
    fTwoTrackDistancePt[1] = (TH3F*) fTwoTrackDistancePt[0]->Clone("fTwoTrackDistancePt[1]");

    TH1::AddDirectory(oldStatus);
  }

  // Eta() is extremely time consuming, therefore cache it for the inner loop here:
  TObjArray* input = (mixed) ? mixed : particles;
  TArrayF eta(input->GetEntriesFast());
  for (Int_t i=0; i<input->GetEntriesFast(); i++)
    eta[i] = ((AliVParticle*) input->UncheckedAt(i))->Eta();
  
  // if particles is not set, just fill event statistics
  if (particles)
  {
    Int_t jMax = particles->GetEntriesFast();
    if (mixed)
      jMax = mixed->GetEntriesFast();
    
    TH1* triggerWeighting = 0;
    if (fWeightPerEvent)
    {
      TAxis* axis = fNumberDensityPhi->GetTrackHist(AliUEHist::kToward)->GetGrid(0)->GetGrid()->GetAxis(2);
      triggerWeighting = new TH1F("triggerWeighting", "", axis->GetNbins(), axis->GetXbins()->GetArray());
    
      for (Int_t i=0; i<particles->GetEntriesFast(); i++)
      {
	AliVParticle* triggerParticle = (AliVParticle*) particles->UncheckedAt(i);
	
	// some optimization
	Float_t triggerEta = triggerParticle->Eta();

	if (fTriggerRestrictEta > 0 && TMath::Abs(triggerEta) > fTriggerRestrictEta)
	  continue;

	if (fOnlyOneEtaSide != 0)
	{
	  if (fOnlyOneEtaSide * triggerEta < 0)
	    continue;
	}
	
	if (fTriggerSelectCharge != 0)
	  if (triggerParticle->Charge() * fTriggerSelectCharge < 0)
	    continue;
	
	triggerWeighting->Fill(triggerParticle->Pt());
      }
    }
    
    // identify K, Lambda candidates and flag those particles
    // a TObject bit is used for this
    const UInt_t kResonanceDaughterFlag = 1 << 14;
    if (fRejectResonanceDaughters > 0)
    {
      Double_t resonanceMass = -1;
      Double_t massDaughter1 = -1;
      Double_t massDaughter2 = -1;
      const Double_t interval = 0.02;
      
      switch (fRejectResonanceDaughters)
      {
	case 1: resonanceMass = 1.2; massDaughter1 = 0.1396; massDaughter2 = 0.9383; break; // method test
	case 2: resonanceMass = 0.4976; massDaughter1 = 0.1396; massDaughter2 = massDaughter1; break; // k0
	case 3: resonanceMass = 1.115; massDaughter1 = 0.1396; massDaughter2 = 0.9383; break; // lambda
	default: AliFatal(Form("Invalid setting %d", fRejectResonanceDaughters));
      }

      for (Int_t i=0; i<particles->GetEntriesFast(); i++)
	particles->UncheckedAt(i)->ResetBit(kResonanceDaughterFlag);
      if (mixed)
	for (Int_t i=0; i<jMax; i++)
	  mixed->UncheckedAt(i)->ResetBit(kResonanceDaughterFlag);
      
      for (Int_t i=0; i<particles->GetEntriesFast(); i++)
      {
	AliVParticle* triggerParticle = (AliVParticle*) particles->UncheckedAt(i);
	
	for (Int_t j=0; j<jMax; j++)
	{
	  if (!mixed && i == j)
	    continue;
	
	  AliVParticle* particle = 0;
	  if (!mixed)
	    particle = (AliVParticle*) particles->UncheckedAt(j);
	  else
	    particle = (AliVParticle*) mixed->UncheckedAt(j);
	  
	  // check if both particles point to the same element (does not occur for mixed events, but if subsets are mixed within the same event)
	  if (fCheckEventNumberInCorrelation)
	  {
	    AliBasicParticle* triggerParticleBasic = dynamic_cast<AliBasicParticle*>(triggerParticle);
	    AliBasicParticle* particleBasic        = dynamic_cast<AliBasicParticle*>(particle);
	    if(!triggerParticleBasic || !particleBasic)
	    {
	      AliFatal("If fCheckEventNumberInCorrelation is set, particle must be derived from AliBasicParticle");
	      continue;
	    }
	
	    if(triggerParticleBasic->IsInSameEvent(particleBasic))
	      continue;
	  }
	  else if (mixed && triggerParticle->IsEqual(particle))
	    continue;
	  
	  if (triggerParticle->Charge() * particle->Charge() > 0)
	    continue;
      
	  Float_t mass = GetInvMassSquaredCheap(triggerParticle->Pt(), triggerParticle->Eta(), triggerParticle->Phi(), particle->Pt(), particle->Eta(), particle->Phi(), massDaughter1, massDaughter2);
	      
	  if (TMath::Abs(mass - resonanceMass*resonanceMass) < interval*5)
	  {
	    mass = GetInvMassSquared(triggerParticle->Pt(), triggerParticle->Eta(), triggerParticle->Phi(), particle->Pt(), particle->Eta(), particle->Phi(), massDaughter1, massDaughter2);

	    if (mass > (resonanceMass-interval)*(resonanceMass-interval) && mass < (resonanceMass+interval)*(resonanceMass+interval))
	    {
	      triggerParticle->SetBit(kResonanceDaughterFlag);
	      particle->SetBit(kResonanceDaughterFlag);
	      
// 	      Printf("Flagged %d %d %f", i, j, TMath::Sqrt(mass));
	    }
	  }
	}
      }
    }
    
    for (Int_t i=0; i<particles->GetEntriesFast(); i++)
    {
      AliVParticle* triggerParticle = (AliVParticle*) particles->UncheckedAt(i);
      
      // some optimization
      Float_t triggerEta = triggerParticle->Eta();
      
      if (fTriggerRestrictEta > 0 && TMath::Abs(triggerEta) > fTriggerRestrictEta)
	continue;

      if (fOnlyOneEtaSide != 0)
      {
	if (fOnlyOneEtaSide * triggerEta < 0)
	  continue;
      }
      
      if (fTriggerSelectCharge != 0)
	if (triggerParticle->Charge() * fTriggerSelectCharge < 0)
	  continue;
	
      if (fRejectResonanceDaughters > 0)
	if (triggerParticle->TestBit(kResonanceDaughterFlag))
	{
// 	  Printf("Skipped i=%d", i);
	  continue;
	}
	
      for (Int_t j=0; j<jMax; j++)
      {
        if (!mixed && i == j)
          continue;
      
        AliVParticle* particle = 0;
        if (!mixed)
          particle = (AliVParticle*) particles->UncheckedAt(j);
        else
          particle = (AliVParticle*) mixed->UncheckedAt(j);
        
        // check if both particles point to the same element (does not occur for mixed events, but if subsets are mixed within the same event)
        if (fCheckEventNumberInCorrelation)
        {
          AliBasicParticle* triggerParticleBasic = dynamic_cast<AliBasicParticle*>(triggerParticle);
          AliBasicParticle* particleBasic        = dynamic_cast<AliBasicParticle*>(particle);
          if(!triggerParticleBasic || !particleBasic)
            AliFatal("If fCheckEventNumberInCorrelation is set, particle must be derived from AliBasicParticle");
      
          if(triggerParticleBasic->IsInSameEvent(particleBasic))
            continue;
        }
        else if (mixed && triggerParticle->IsEqual(particle))
          continue;
        
        if (fPtOrder)
	  if (particle->Pt() >= triggerParticle->Pt())
	    continue;
	
	if (fAssociatedSelectCharge != 0)
	  if (particle->Charge() * fAssociatedSelectCharge < 0)
	    continue;

        if (fSelectCharge > 0)
        {
          // skip like sign
          if (fSelectCharge == 1 && particle->Charge() * triggerParticle->Charge() > 0)
            continue;
            
          // skip unlike sign
          if (fSelectCharge == 2 && particle->Charge() * triggerParticle->Charge() < 0)
            continue;
        }
        
        if (fOnlyOneAssocEtaSide != 0)
          if (fOnlyOneAssocEtaSide * eta[j] < 0)
            continue;

	if (fEtaOrdering)
	{
	  if (triggerEta < 0 && eta[j] < triggerEta)
	    continue;
	  if (triggerEta > 0 && eta[j] > triggerEta)
	    continue;
	}

	if (fRejectResonanceDaughters > 0)
	  if (particle->TestBit(kResonanceDaughterFlag))
	  {
// 	    Printf("Skipped j=%d", j);
	    continue;
	  }

	// conversions
	if (twoTrackCuts && fCutConversionsV > 0 && particle->Charge() * triggerParticle->Charge() < 0)
	{
	  Float_t mass = GetInvMassSquaredCheap(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), particle->Pt(), eta[j], particle->Phi(), 0.510e-3, 0.510e-3);
	  
	  if (mass < fCutConversionsV * 5)
	  {
	    mass = GetInvMassSquared(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), particle->Pt(), eta[j], particle->Phi(), 0.510e-3, 0.510e-3);
	    
	    fControlConvResoncances->Fill(0.0, mass);

	    if (mass < fCutConversionsV*fCutConversionsV) 
	      continue;
	  }
	}
	
	// K0s
	if (twoTrackCuts && fCutK0sV > 0 && particle->Charge() * triggerParticle->Charge() < 0)
	{
	  Float_t mass = GetInvMassSquaredCheap(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), particle->Pt(), eta[j], particle->Phi(), 0.1396, 0.1396);
	  
	  const Float_t kK0smass = 0.4976;
	  
	  if (TMath::Abs(mass - kK0smass*kK0smass) < fCutK0sV * 5)
	  {
	    mass = GetInvMassSquared(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), particle->Pt(), eta[j], particle->Phi(), 0.1396, 0.1396);
	    
	    fControlConvResoncances->Fill(1, mass - kK0smass*kK0smass);

	    if (mass > (kK0smass-fCutK0sV)*(kK0smass-fCutK0sV) && mass < (kK0smass+fCutK0sV)*(kK0smass+fCutK0sV))
	      continue;
	  }
	}

	// Lambda
	if (twoTrackCuts && fCutLambdaV > 0 && particle->Charge() * triggerParticle->Charge() < 0)
	{
	  Float_t mass1 = GetInvMassSquaredCheap(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), particle->Pt(), eta[j], particle->Phi(), 0.1396, 0.9383);
	  Float_t mass2 = GetInvMassSquaredCheap(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), particle->Pt(), eta[j], particle->Phi(), 0.9383, 0.1396);
	  
	  const Float_t kLambdaMass = 1.115;

	  if (TMath::Abs(mass1 - kLambdaMass*kLambdaMass) < fCutLambdaV * 5)
	  {
	    mass1 = GetInvMassSquared(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), particle->Pt(), eta[j], particle->Phi(), 0.1396, 0.9383);

	    fControlConvResoncances->Fill(2, mass1 - kLambdaMass*kLambdaMass);
	    
	    if (mass1 > (kLambdaMass-fCutLambdaV)*(kLambdaMass-fCutLambdaV) && mass1 < (kLambdaMass+fCutLambdaV)*(kLambdaMass+fCutLambdaV))
	      continue;
	  }
	  if (TMath::Abs(mass2 - kLambdaMass*kLambdaMass) < fCutLambdaV * 5)
	  {
	    mass2 = GetInvMassSquared(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), particle->Pt(), eta[j], particle->Phi(), 0.9383, 0.1396);

	    fControlConvResoncances->Fill(2, mass2 - kLambdaMass*kLambdaMass);

	    if (mass2 > (kLambdaMass-fCutLambdaV)*(kLambdaMass-fCutLambdaV) && mass2 < (kLambdaMass+fCutLambdaV)*(kLambdaMass+fCutLambdaV))
	      continue;
	  }
	}

        // Phi
	if (twoTrackCuts && fCutPhiV > 0 && particle->Charge() * triggerParticle->Charge() < 0)
	{
	  Float_t mass = GetInvMassSquaredCheap(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), particle->Pt(), eta[j], particle->Phi(), 0.4937, 0.4937);
	  
	  const Float_t kPhimass = 1.019;
	  
	  if (TMath::Abs(mass - kPhimass*kPhimass) < fCutPhiV * 5)
	  {
	    mass = GetInvMassSquared(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), particle->Pt(), eta[j], particle->Phi(), 0.4937, 0.4937);
	    
	    fControlConvResoncances->Fill(3, mass - kPhimass*kPhimass);
	    
	    if (mass > (kPhimass-fCutPhiV)*(kPhimass-fCutPhiV) && mass < (kPhimass+fCutPhiV)*(kPhimass+fCutPhiV))
	      continue;
	  }
	}	

        // Rho
	if (twoTrackCuts && fCutRhoV > 0 && particle->Charge() * triggerParticle->Charge() < 0)
        {
	  Float_t mass = GetInvMassSquaredCheap(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), particle->Pt(), eta[j], particle->Phi(), 0.1396, 0.1396);
	  
	  const Float_t kRhomass = 0.770;
	  
	  if (TMath::Abs(mass - kRhomass*kRhomass) < fCutRhoV * 5)
          {
	    mass = GetInvMassSquared(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), particle->Pt(), eta[j], particle->Phi(), 0.1396, 0.1396);
	    
	    fControlConvResoncances->Fill(4, mass - kRhomass*kRhomass);
	    
	    if (mass > (kRhomass-fCutRhoV)*(kRhomass-fCutRhoV) && mass < (kRhomass+fCutRhoV)*(kRhomass+fCutRhoV))
	      continue;
	  }
	}

        // User-defined cut
	if (twoTrackCuts && fCutCustomMass > 0 && fCutCustomFirst > 0 && fCutCustomSecond > 0 && fCutCustomV > 0 && particle->Charge() * triggerParticle->Charge() < 0)
        {
	  Float_t mass = GetInvMassSquaredCheap(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), particle->Pt(), eta[j], particle->Phi(), fCutCustomFirst, fCutCustomSecond);
	  
	  if (TMath::Abs(mass - fCutCustomMass*fCutCustomMass) < fCutCustomV * 5)
          {
	    mass = GetInvMassSquared(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), particle->Pt(), eta[j], particle->Phi(), fCutCustomFirst, fCutCustomSecond);
	    
	    fControlConvResoncances->Fill(5, mass - fCutCustomMass*fCutCustomMass);
	    
	    if (mass > (fCutCustomMass-fCutCustomV)*(fCutCustomMass-fCutCustomV) && mass < (fCutCustomMass+fCutCustomV)*(fCutCustomMass+fCutCustomV))
	      continue;
	  }
	}

	if (twoTrackCuts && twoTrackEfficiencyCutValue > 0)
	{
	  // the variables & cuthave been developed by the HBT group 
	  // see e.g. https://indico.cern.ch/materialDisplay.py?contribId=36&sessionId=6&materialId=slides&confId=142700

	  Float_t phi1 = triggerParticle->Phi();
	  Float_t pt1 = triggerParticle->Pt();
	  Float_t charge1 = triggerParticle->Charge();
	    
	  Float_t phi2 = particle->Phi();
	  Float_t pt2 = particle->Pt();
	  Float_t charge2 = particle->Charge();
	      
	  Float_t deta = triggerEta - eta[j];
	      
	  // optimization
	  if (TMath::Abs(deta) < twoTrackEfficiencyCutValue * 2.5 * 3)
	  {
	    // check first boundaries to see if is worth to loop and find the minimum
	    Float_t dphistar1 = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, fTwoTrackCutMinRadius, bSign);
	    Float_t dphistar2 = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, 2.5, bSign);
	    
	    const Float_t kLimit = twoTrackEfficiencyCutValue * 3;

	    Float_t dphistarminabs = 1e5;
	    Float_t dphistarmin = 1e5;
	    if (TMath::Abs(dphistar1) < kLimit || TMath::Abs(dphistar2) < kLimit || dphistar1 * dphistar2 < 0)
	    {
	      for (Double_t rad=fTwoTrackCutMinRadius; rad<2.51; rad+=0.01) 
	      {
		Float_t dphistar = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, rad, bSign);

		Float_t dphistarabs = TMath::Abs(dphistar);
		
		if (dphistarabs < dphistarminabs)
		{
		  dphistarmin = dphistar;
		  dphistarminabs = dphistarabs;
		}
	      }
	      
	      fTwoTrackDistancePt[0]->Fill(deta, dphistarmin, TMath::Abs(pt1 - pt2));
	      
	      if (dphistarminabs < twoTrackEfficiencyCutValue && TMath::Abs(deta) < twoTrackEfficiencyCutValue)
	      {
// 		Printf("Removed track pair %d %d with %f %f %f %f %f %f %f %f %f", i, j, deta, dphistarminabs, phi1, pt1, charge1, phi2, pt2, charge2, bSign);
		continue;
	      }

    	      fTwoTrackDistancePt[1]->Fill(deta, dphistarmin, TMath::Abs(pt1 - pt2));
	    }
	  }
	}
        
        Double_t vars[6];
        vars[0] = triggerEta - eta[j];
        vars[1] = particle->Pt();
        vars[2] = triggerParticle->Pt();
        vars[3] = centrality;
        vars[4] = triggerParticle->Phi() - particle->Phi();
        if (vars[4] > 1.5 * TMath::Pi()) 
          vars[4] -= TMath::TwoPi();
        if (vars[4] < -0.5 * TMath::Pi())
          vars[4] += TMath::TwoPi();
	vars[5] = zVtx;
	
	if (fillpT)
	  weight = particle->Pt();
	
	Double_t useWeight = weight;
	if (applyEfficiency)
	{
	  if (fEfficiencyCorrectionAssociated)
	  {
	    Int_t effVars[4];
	    // associated particle
	    effVars[0] = fEfficiencyCorrectionAssociated->GetAxis(0)->FindBin(eta[j]);
	    effVars[1] = fEfficiencyCorrectionAssociated->GetAxis(1)->FindBin(vars[1]); //pt
	    effVars[2] = fEfficiencyCorrectionAssociated->GetAxis(2)->FindBin(vars[3]); //centrality
	    effVars[3] = fEfficiencyCorrectionAssociated->GetAxis(3)->FindBin(vars[5]); //zVtx
	    
	    // 	  Printf("%d %d %d %d %f", effVars[0], effVars[1], effVars[2], effVars[3], fEfficiencyCorrectionAssociated->GetBinContent(effVars));
	  
	    useWeight *= fEfficiencyCorrectionAssociated->GetBinContent(effVars);
	  }
	  if (fEfficiencyCorrectionTriggers)
	  {
	    Int_t effVars[4];

	    effVars[0] = fEfficiencyCorrectionTriggers->GetAxis(0)->FindBin(triggerEta);
	    effVars[1] = fEfficiencyCorrectionTriggers->GetAxis(1)->FindBin(vars[2]); //pt
	    effVars[2] = fEfficiencyCorrectionTriggers->GetAxis(2)->FindBin(vars[3]); //centrality
	    effVars[3] = fEfficiencyCorrectionTriggers->GetAxis(3)->FindBin(vars[5]); //zVtx
	    useWeight *= fEfficiencyCorrectionTriggers->GetBinContent(effVars);
	  }
	}

	if (fWeightPerEvent)
	{
	  Int_t weightBin = triggerWeighting->GetXaxis()->FindBin(vars[2]);
// 	  Printf("Using weight %f", triggerWeighting->GetBinContent(weightBin));
	  useWeight /= triggerWeighting->GetBinContent(weightBin);
	}
    
        // fill all in toward region and do not use the other regions
	fNumberDensityPhi->GetTrackHist(AliUEHist::kToward)->Fill(vars, step, useWeight);

// 	Printf("%.2f %.2f --> %.2f", triggerEta, eta[j], vars[0]);
      }
 
      if (firstTime)
      {
        // once per trigger particle
        Double_t vars[3];
        vars[0] = triggerParticle->Pt();
        vars[1] = centrality;
	vars[2] = zVtx;

	Double_t useWeight = 1;
	if (fEfficiencyCorrectionTriggers && applyEfficiency)
	{
	  Int_t effVars[4];
	  
	  // trigger particle
	  effVars[0] = fEfficiencyCorrectionTriggers->GetAxis(0)->FindBin(triggerEta);
	  effVars[1] = fEfficiencyCorrectionTriggers->GetAxis(1)->FindBin(vars[0]); //pt
	  effVars[2] = fEfficiencyCorrectionTriggers->GetAxis(2)->FindBin(vars[1]); //centrality
	  effVars[3] = fEfficiencyCorrectionTriggers->GetAxis(3)->FindBin(vars[2]); //zVtx
	  useWeight *= fEfficiencyCorrectionTriggers->GetBinContent(effVars);
	}

	if (TMath::Abs(triggerEta) < 0.8 && triggerParticle->Pt() > 0)
	  fInvYield2->Fill(centrality, triggerParticle->Pt(), useWeight / triggerParticle->Pt());

	if (fWeightPerEvent)
	{
	  // leads effectively to a filling of one entry per filled trigger particle pT bin
	  Int_t weightBin = triggerWeighting->GetXaxis()->FindBin(vars[0]);
// 	  Printf("Using weight %f", triggerWeighting->GetBinContent(weightBin));
	  useWeight /= triggerWeighting->GetBinContent(weightBin);
	}
	
        fNumberDensityPhi->GetEventHist()->Fill(vars, step, useWeight);

	// QA
        fCorrelationpT->Fill(centrality, triggerParticle->Pt());
        fCorrelationEta->Fill(centrality, triggerEta);
        fCorrelationPhi->Fill(centrality, triggerParticle->Phi());
	fYields->Fill(centrality, triggerParticle->Pt(), triggerEta);
	fYieldsEtaPhiPT->Fill(triggerParticle->Pt(), triggerEta, triggerParticle->Phi());
	
/*        if (dynamic_cast<AliAODTrack*>(triggerParticle))
          fITSClusterMap->Fill(((AliAODTrack*) triggerParticle)->GetITSClusterMap(), centrality, triggerParticle->Pt());*/
      }
    }
    
    if (triggerWeighting)
    {
      delete triggerWeighting;
      triggerWeighting = 0;
    }
  }
  
  fCentralityDistribution->Fill(centrality);
  fCentralityCorrelation->Fill(centrality, particles->GetEntriesFast());
  FillEvent(centrality, step);
}
  
//____________________________________________________________________
void AliUEHistograms::FillTrackingEfficiency(TObjArray* mc, TObjArray* recoPrim, TObjArray* recoAll, TObjArray* recoPrimPID, TObjArray* recoAllPID, TObjArray* fake, Int_t particleType, Double_t centrality, Double_t zVtx)
{
  // fills the tracking efficiency objects
  //
  // mc: all primary MC particles
  // recoPrim: reconstructed primaries (again MC particles)
  // recoAll: reconstructed (again MC particles)
  // recoPrim: reconstructed primaries with checks on PID (again MC particles)
  // recoAll: reconstructed with checks on PID (again MC particles)
  // particleType is: 0 for pion, 1 for kaon, 2 for proton, 3 for others
 
  for (Int_t step=0; step<6; step++)
  {
    TObjArray* list = mc;
    if (step == 1)
      list = recoPrim;
    else if (step == 2)
      list = recoAll;
    else if (step == 3)
      list = recoPrimPID;
    else if (step == 4)
      list = recoAllPID;
    else if (step == 5)
      list = fake;
    
    if (!list)
      continue;

    for (Int_t i=0; i<list->GetEntriesFast(); i++)
    {
      AliVParticle* particle = (AliVParticle*) list->UncheckedAt(i);
      Double_t vars[5];
      vars[0] = particle->Eta();
      vars[1] = particle->Pt();
      vars[2] = particleType;
      vars[3] = centrality;
      vars[4] = zVtx;
      
      for (Int_t j=0; j<fgkUEHists; j++)
        if (GetUEHist(j))
          GetUEHist(j)->GetTrackHistEfficiency()->Fill(vars, step);
    }
  }
}

//____________________________________________________________________
void AliUEHistograms::FillFakePt(TObjArray* fake, Double_t centrality)
{
  TObjArray* tracksReco = (TObjArray*) fake->At(0);
  TObjArray* tracksMC = (TObjArray*) fake->At(1);
  
  if (tracksReco->GetEntriesFast() != tracksMC->GetEntriesFast())
    AliFatal(Form("Inconsistent arrays: %d vs %d", tracksReco->GetEntriesFast(), tracksMC->GetEntriesFast()));
  
  for (Int_t i=0; i<tracksReco->GetEntriesFast(); i++)
  {
    AliVParticle* particle1 = (AliVParticle*) tracksReco->At(i);
    AliVParticle* particle2 = (AliVParticle*) tracksMC->At(i);
    Double_t vars[3];
    vars[0] = particle1->Pt();
    vars[1] = particle2->Pt();
    vars[2] = centrality;
    for (Int_t j=0; j<fgkUEHists; j++)
      if (GetUEHist(j))
        GetUEHist(j)->GetMCRecoPtCorrelation()->Fill(vars[0],vars[1],vars[2]);
  }
}

//____________________________________________________________________
void AliUEHistograms::FillEvent(Int_t eventType, Int_t step)
{
  // fills the number of events at the given step and the given enty type
  //
  // WARNING: This function is called from Fill, so only call it for steps where Fill is not called
  
  fEventCount->Fill(step, eventType);
}

//____________________________________________________________________
void AliUEHistograms::FillEvent(Double_t centrality, Int_t step)
{
  // fills the number of events at the given step and the given centrality
  //
  // WARNING: This function is called from Fill, so only call it for steps where Fill is not called
  
  fEventCount->Fill(step, centrality);
}

//____________________________________________________________________
void AliUEHistograms::SetEtaRange(Float_t etaMin, Float_t etaMax)
{
  // sets eta min and max for all contained AliUEHist classes
  
  for (Int_t i=0; i<fgkUEHists; i++)
    if (GetUEHist(i))
      GetUEHist(i)->SetEtaRange(etaMin, etaMax);
}

//____________________________________________________________________
void AliUEHistograms::SetPtRange(Float_t ptMin, Float_t ptMax)
{
  // sets pT min and max for all contained AliUEHist classes
  
  for (Int_t i=0; i<fgkUEHists; i++)
    if (GetUEHist(i))
      GetUEHist(i)->SetPtRange(ptMin, ptMax);
}

//____________________________________________________________________
void AliUEHistograms::SetPartSpecies(Int_t species)
{
  // sets PartSpecie for all contained AliUEHist classes
  
  for (Int_t i=0; i<fgkUEHists; i++)
    if (GetUEHist(i))
      GetUEHist(i)->SetPartSpecies(species);
}

//____________________________________________________________________
void AliUEHistograms::SetZVtxRange(Float_t min, Float_t max)
{
  // sets pT min and max for all contained AliUEHist classes
  
  for (Int_t i=0; i<fgkUEHists; i++)
    if (GetUEHist(i))
      GetUEHist(i)->SetZVtxRange(min, max);
}

//____________________________________________________________________
void AliUEHistograms::SetContaminationEnhancement(TH1F* hist)
{
  // sets the contamination enhancement histogram in all contained AliUEHist classes
  
  for (Int_t i=0; i<fgkUEHists; i++)
    if (GetUEHist(i))
      GetUEHist(i)->SetContaminationEnhancement(hist);
}  

//____________________________________________________________________
void AliUEHistograms::SetCombineMinMax(Bool_t flag)
{
  // sets pT min and max for all contained AliUEHist classes
  
  for (Int_t i=0; i<fgkUEHists; i++)
    if (GetUEHist(i))
      GetUEHist(i)->SetCombineMinMax(flag);
}

//____________________________________________________________________
void AliUEHistograms::SetTrackEtaCut(Float_t value)
{
  // sets track eta cut for all contained AliUEHist classes
  
  for (Int_t i=0; i<fgkUEHists; i++)
    if (GetUEHist(i))
      GetUEHist(i)->SetTrackEtaCut(value);
}

//____________________________________________________________________
void AliUEHistograms::SetWeightPerEvent(Bool_t flag)
{
  // sets fWeightPerEvent for all contained AliUEHist classes
  
  fWeightPerEvent = flag;
  
  for (Int_t i=0; i<fgkUEHists; i++)
    if (GetUEHist(i))
      GetUEHist(i)->SetWeightPerEvent(fWeightPerEvent);
}

//____________________________________________________________________
void AliUEHistograms::Correct(AliUEHistograms* corrections)
{
  // corrects the contained histograms by calling AliUEHist::Correct
  
  for (Int_t i=0; i<fgkUEHists; i++)
    if (GetUEHist(i))
      GetUEHist(i)->Correct(corrections->GetUEHist(i));
}

//____________________________________________________________________
AliUEHistograms &AliUEHistograms::operator=(const AliUEHistograms &c)
{
  // assigment operator

  DeleteContainers();

  if (this != &c)
    ((AliUEHistograms &) c).Copy(*this);

  return *this;
}

//____________________________________________________________________
void AliUEHistograms::Copy(TObject& c) const
{
  // copy function

  AliUEHistograms& target = (AliUEHistograms &) c;

  if (fNumberDensitypT)
    target.fNumberDensitypT = dynamic_cast<AliUEHist*> (fNumberDensitypT->Clone());

  if (fSumpT)
    target.fSumpT = dynamic_cast<AliUEHist*> (fSumpT->Clone());

  if (fNumberDensityPhi)
    target.fNumberDensityPhi = dynamic_cast<AliUEHist*> (fNumberDensityPhi->Clone());

  if (fCorrelationpT)
    target.fCorrelationpT = dynamic_cast<TH2F*> (fCorrelationpT->Clone());

  if (fCorrelationEta)
    target.fCorrelationEta = dynamic_cast<TH2F*> (fCorrelationEta->Clone());

  if (fCorrelationPhi)
    target.fCorrelationPhi = dynamic_cast<TH2F*> (fCorrelationPhi->Clone());

  if (fCorrelationR)
    target.fCorrelationR = dynamic_cast<TH2F*> (fCorrelationR->Clone());

  if (fCorrelationLeading2Phi)
    target.fCorrelationLeading2Phi = dynamic_cast<TH2F*> (fCorrelationLeading2Phi->Clone());

  if (fCorrelationMultiplicity)
    target.fCorrelationMultiplicity = dynamic_cast<TH2F*> (fCorrelationMultiplicity->Clone());
  
  if (fYields)
    target.fYields = dynamic_cast<TH3F*> (fYields->Clone());

  if (fYieldsEtaPhiPT)
    target.fYieldsEtaPhiPT = dynamic_cast<TH3F*> (fYieldsEtaPhiPT->Clone());
  
  if (fInvYield2)
    target.fInvYield2 = dynamic_cast<TH2F*> (fInvYield2->Clone());
  
  if (fEventCount)
    target.fEventCount = dynamic_cast<TH2F*> (fEventCount->Clone());

  if (fEventCountDifferential)
    target.fEventCountDifferential = dynamic_cast<TH3F*> (fEventCountDifferential->Clone());
    
  if (fVertexContributors)
    target.fVertexContributors = dynamic_cast<TH1F*> (fVertexContributors->Clone());

  if (fCentralityDistribution)
    target.fCentralityDistribution = dynamic_cast<TH1F*> (fCentralityDistribution->Clone());
    
  if (fCentralityCorrelation)
    target.fCentralityCorrelation = dynamic_cast<TH2F*> (fCentralityCorrelation->Clone());

  if (fITSClusterMap)
    target.fITSClusterMap = dynamic_cast<TH3F*> (fITSClusterMap->Clone());
  
  if (fControlConvResoncances)
    target.fControlConvResoncances = dynamic_cast<TH2F*> (fControlConvResoncances->Clone());
  
  for (Int_t i=0; i<2; i++)
    if (fTwoTrackDistancePt[i])
      target.fTwoTrackDistancePt[i] = dynamic_cast<TH3F*> (fTwoTrackDistancePt[i]->Clone());

  if (fEfficiencyCorrectionTriggers)
    target.fEfficiencyCorrectionTriggers = dynamic_cast<THnF*> (fEfficiencyCorrectionTriggers->Clone());
 
 if (fEfficiencyCorrectionAssociated)
    target.fEfficiencyCorrectionAssociated = dynamic_cast<THnF*> (fEfficiencyCorrectionAssociated->Clone());
    
  target.fSelectCharge = fSelectCharge;
  target.fTriggerSelectCharge = fTriggerSelectCharge;
  target.fAssociatedSelectCharge = fAssociatedSelectCharge;
  target.fTriggerRestrictEta = fTriggerRestrictEta;
  target.fEtaOrdering = fEtaOrdering;
  target.fCutConversionsV = fCutConversionsV;
  target.fCutK0sV = fCutK0sV;
  target.fCutLambdaV = fCutLambdaV;
  target.fCutPhiV = fCutPhiV;
  target.fCutRhoV = fCutRhoV;
  target.fCutCustomMass = fCutCustomMass;
  target.fCutCustomFirst = fCutCustomFirst;
  target.fCutCustomSecond = fCutCustomSecond;
  target.fCutCustomV = fCutCustomV;
  target.fOnlyOneEtaSide = fOnlyOneEtaSide;
  target.fOnlyOneAssocEtaSide = fOnlyOneAssocEtaSide;
  target.fWeightPerEvent = fWeightPerEvent;
  target.fRunNumber = fRunNumber;
  target.fMergeCount = fMergeCount;
  target.fWeightPerEvent = fWeightPerEvent;
  target.fPtOrder = fPtOrder;
  target.fTwoTrackCutMinRadius = fTwoTrackCutMinRadius;
  target.fCheckEventNumberInCorrelation = fCheckEventNumberInCorrelation;
}

//____________________________________________________________________
Long64_t AliUEHistograms::Merge(TCollection* list)
{
  // Merge a list of AliUEHistograms objects with this (needed for
  // PROOF). 
  // Returns the number of merged objects (including this).

  if (!list)
    return 0;
  
  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;

  // collections of objects
  const Int_t kMaxLists = 21;
  TList* lists[kMaxLists];
  
  for (Int_t i=0; i<kMaxLists; i++)
    lists[i] = new TList;
  
  Int_t count = 0;
  while ((obj = iter->Next())) {
    
    AliUEHistograms* entry = dynamic_cast<AliUEHistograms*> (obj);
    if (entry == 0) 
      continue;

    if (entry->fNumberDensitypT)
      lists[0]->Add(entry->fNumberDensitypT);
    if (entry->fSumpT)
      lists[1]->Add(entry->fSumpT);
    if (entry->fNumberDensityPhi)
      lists[2]->Add(entry->fNumberDensityPhi);
    lists[3]->Add(entry->fCorrelationpT);
    lists[4]->Add(entry->fCorrelationEta);
    lists[5]->Add(entry->fCorrelationPhi);
    lists[6]->Add(entry->fCorrelationR);
    lists[7]->Add(entry->fCorrelationLeading2Phi);
    lists[8]->Add(entry->fCorrelationMultiplicity);
    lists[9]->Add(entry->fEventCount);
    lists[10]->Add(entry->fEventCountDifferential);
    lists[11]->Add(entry->fVertexContributors);
    lists[12]->Add(entry->fCentralityDistribution);
    lists[13]->Add(entry->fITSClusterMap);
    if (entry->fTwoTrackDistancePt[0])
      lists[14]->Add(entry->fTwoTrackDistancePt[0]);
    if (entry->fTwoTrackDistancePt[1])
      lists[15]->Add(entry->fTwoTrackDistancePt[1]);
    if (entry->fCentralityCorrelation)
      lists[16]->Add(entry->fCentralityCorrelation);
    if (entry->fYields)
      lists[17]->Add(entry->fYields);
    if (entry->fInvYield2)
      lists[18]->Add(entry->fInvYield2);
    if (entry->fControlConvResoncances)
      lists[19]->Add(entry->fControlConvResoncances);
    if (entry->fYieldsEtaPhiPT)
      lists[20]->Add(entry->fYieldsEtaPhiPT);

    fMergeCount += entry->fMergeCount;

    count++;
  }
    
  if (fNumberDensitypT)
    fNumberDensitypT->Merge(lists[0]);
  if (fSumpT)
    fSumpT->Merge(lists[1]);
  if (fNumberDensityPhi)
    fNumberDensityPhi->Merge(lists[2]);
  fCorrelationpT->Merge(lists[3]);
  fCorrelationEta->Merge(lists[4]);
  fCorrelationPhi->Merge(lists[5]);
  fCorrelationR->Merge(lists[6]);
  fCorrelationLeading2Phi->Merge(lists[7]);
  fCorrelationMultiplicity->Merge(lists[8]);
  fEventCount->Merge(lists[9]);
  fEventCountDifferential->Merge(lists[10]);
  fVertexContributors->Merge(lists[11]);
  fCentralityDistribution->Merge(lists[12]);
  fITSClusterMap->Merge(lists[13]);
  if (fTwoTrackDistancePt[0] && lists[14]->GetEntries() > 0)
    fTwoTrackDistancePt[0]->Merge(lists[14]);
  if (fTwoTrackDistancePt[1] && lists[15]->GetEntries() > 0)
    fTwoTrackDistancePt[1]->Merge(lists[15]);
  if (fCentralityCorrelation)
    fCentralityCorrelation->Merge(lists[16]);
  if (fYields && lists[17]->GetEntries() > 0)
    fYields->Merge(lists[17]);
  if (fInvYield2 && lists[18]->GetEntries() > 0)
    fInvYield2->Merge(lists[18]);
  if (fControlConvResoncances && lists[19]->GetEntries() > 0)
    fControlConvResoncances->Merge(lists[19]);
  if (fYieldsEtaPhiPT && lists[20]->GetEntries() > 0)
    fYieldsEtaPhiPT->Merge(lists[20]);
  
  for (Int_t i=0; i<kMaxLists; i++)
    delete lists[i];
  
  return count+1;
}

void AliUEHistograms::CopyReconstructedData(AliUEHistograms* from)
{
  // copies those histograms extracted from ESD to this object
  
  for (Int_t i=0; i<fgkUEHists; i++)
    if (GetUEHist(i))
      GetUEHist(i)->CopyReconstructedData(from->GetUEHist(i));
}

void AliUEHistograms::DeepCopy(AliUEHistograms* from)
{
  // copies the entries of this object's members from the object <from> to this object

  for (Int_t i=0; i<fgkUEHists; i++)
    if (GetUEHist(i) && from->GetUEHist(i))
      GetUEHist(i)->DeepCopy(from->GetUEHist(i));
}

void AliUEHistograms::ExtendTrackingEfficiency(Bool_t verbose)
{
  // delegates to AliUEHists

  for (Int_t i=0; i<fgkUEHists; i++)
    if (GetUEHist(i))
      GetUEHist(i)->ExtendTrackingEfficiency(verbose);
}

void AliUEHistograms::Scale(Double_t factor)
{
  // scales all contained histograms by the given factor
  
  for (Int_t i=0; i<fgkUEHists; i++)
    if (GetUEHist(i))
      GetUEHist(i)->Scale(factor);
      
  TList list;
  list.Add(fCorrelationpT);
  list.Add(fCorrelationEta);
  list.Add(fCorrelationPhi);
  list.Add(fCorrelationR);
  list.Add(fCorrelationLeading2Phi);
  list.Add(fCorrelationMultiplicity);
  list.Add(fYields);
  list.Add(fYieldsEtaPhiPT);
  list.Add(fInvYield2);
  list.Add(fEventCount);
  list.Add(fEventCountDifferential);
  list.Add(fVertexContributors);
  list.Add(fCentralityDistribution);
  list.Add(fCentralityCorrelation);
  list.Add(fITSClusterMap);
  list.Add(fTwoTrackDistancePt[0]);
  list.Add(fTwoTrackDistancePt[1]);
  list.Add(fControlConvResoncances);
  
  for (Int_t i=0; i<list.GetEntries(); i++)
    ((TH1*) list.At(i))->Scale(factor);
}

void AliUEHistograms::Reset()
{
  // delegates to AliUEHists

  for (Int_t i=0; i<fgkUEHists; i++)
    if (GetUEHist(i))
      GetUEHist(i)->Reset();
}
