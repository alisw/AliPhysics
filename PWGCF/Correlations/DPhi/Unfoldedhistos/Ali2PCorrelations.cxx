/**************************************************************************
* Copyright(c) 1998-2021, ALICE Experiment at CERN, All rights reserved.  *
*                                                                         *
* Permission to use, copy, modify and distribute this software and its    *
* documentation strictly for non-commercial purposes is hereby granted    *
* without fee, provided that the above copyright notice appears in all    *
* copies and that both the copyright notice and this permission notice    *
* appear in the supporting documentation. The authors make no claims      *
* about the suitability of this software for any purpose. It is           *
* provided "as is" without express or implied warranty.                   *
**************************************************************************/

/// \file Ali2PCorrelations.cxx
/// \brief Implementation of the Ali2PCorrelations class

#include <TObjString.h>
#include <TList.h>
#include <TParameter.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THn.h>
#include <TProfile.h>
#include <TMath.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include "AliLog.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisManager.h"
#include "Ali2PCorrelations.h"

const char *Ali2PCorrelations::TrackPairsNames[] = {"OO","OT","TO","TT"};
Int_t Ali2PCorrelations::fgkNoOfResonances = 4; ///< four resonances / conversions for the time being
Double_t Ali2PCorrelations::fgkMass[16] = {/* photon */ 0.0, /* k0 */ 0.4976, /* lambda */ 1.115, /* rho */ 0.775};
Double_t Ali2PCorrelations::fgkChildMass[2][16] = {
    {0.510e-3, 0.1396, 0.1396, 0.1396},
    {0.510e-3, 0.1396, 0.9383, 0.1396}
};
Double_t Ali2PCorrelations::fgkMassThreshold[16] = {0.04,0.01,0.05,0.04};

/// Default constructor for object serialization
Ali2PCorrelations::Ali2PCorrelations() :
    TNamed(),
    fOutput(nullptr),
    fVertexZ(0.0),
    fIxVertexZ(0),
    fCentrality(0.0),
    fSinglesOnly(true),
    fUseWeights(false),
    fUseSimulation(false),
    /* the arrays with tracks 1 and 2 information */
    fArraySize(5*1024),
    fNoOfTracks1(0),
    fId_1(nullptr),
    fCharge_1(nullptr),
    fIxEta_1(nullptr),
    fIxPhi_1(nullptr),
    fIxPt_1(nullptr),
    fFlags_1(nullptr),
    fPt_1(nullptr),
    fEta_1(nullptr),
    fPhi_1(nullptr),
    fCorrection_1(nullptr),
    fNoOfTracks2(0),
    fId_2(nullptr),
    fCharge_2(nullptr),
    fIxEta_2(nullptr),
    fIxPhi_2(nullptr),
    fIxPt_2(nullptr),
    fFlags_2(nullptr),
    fPt_2(nullptr),
    fEta_2(nullptr),
    fPhi_2(nullptr),
    fCorrection_2(nullptr),
    /* correction weights */
    fCorrectionWeights_1(nullptr),
    fCorrectionWeights_2(nullptr),
    /* efficiency correction */
    fEfficiencyCorrection_1(nullptr),
    fEfficiencyCorrection_2(nullptr),
    fPairsEfficiency_PP(nullptr),
    fPairsEfficiency_PM(nullptr),
    fPairsEfficiency_MM(nullptr),
    fPairsEfficiency_MP(nullptr),
    /* simulation pdfs */
    fPositiveTrackPdfs(nullptr),
    fNegativeTrackPdfs(nullptr),
    fPositiveTrackCurrentPdf(nullptr),
    fNegativeTrackCurrentPdf(nullptr),
    fSimEventsPerEvent(1),
    /* vertex bins */
    fNBins_vertexZ(40),
    fMin_vertexZ(-10.0),
    fMax_vertexZ(10.0),
    fWidth_vertexZ(0.5),
    /* phi origin shift */
    fNBinsPhiShift(0.0),
    /* pT1 bins */
    fNBins_pt_1(18),
    fMin_pt_1(0.2),
    fMax_pt_1(2.0),
    fWidth_pt_1(0.1),
    /* phi1 bins */
    fNBins_phi_1(72),
    fMin_phi_1(0.0),
    fMax_phi_1(TMath::Pi()*2.0),
    fWidth_phi_1(TMath::Pi()*2.0/72.0),
    /* eta1 bins */
    fNBins_eta_1(20),
    fMin_eta_1(-1.0),
    fMax_eta_1(1.0),
    fWidth_eta_1(0.1),
    fNBins_etaPhi_1(1440),
    fNBins_etaPhiPt_1(25920),
    fNBins_zEtaPhiPt_1(1036800),
    /* pT2 bins */
    fNBins_pt_2(18),
    fMin_pt_2(0.2),
    fMax_pt_2(2.0),
    fWidth_pt_2(0.1),
    /* phi2 bins */
    fNBins_phi_2(72),
    fMin_phi_2(0.0),
    fMax_phi_2(TMath::Pi()*2.0),
    fWidth_phi_2(TMath::Pi()*2.0/72.0),
    /* eta2 bins */
    fNBins_eta_2(20),
    fMin_eta_2(-1.0),
    fMax_eta_2(1.0),
    fWidth_eta_2(0.1),
    fNBins_etaPhi_2(1440),
    fNBins_etaPhiPt_2(25920),
    fNBins_zEtaPhiPt_2(1036800),
    /* pair bins */
    fNBins_deltaphi(72),
    fMin_deltaphi(0.0),
    fMax_deltaphi(TMath::Pi()*2.0),
    fWidth_deltaphi(TMath::Pi()*2.0/72.0),
    fNBins_deltaeta(20*2-1),
    fMin_deltaeta(-2.0),
    fMax_deltaeta(2.0),
    fWidth_deltaeta(4.0/(20*2-1)),

    /* accumulated values */
    fN1_1(0.0),
    fN1_2(0.0),
    fNnw1_1(0.0),
    fNnw1_2(0.0),
    fSum1Pt_1(0.0),
    fSum1Pt_2(0.0),
    fSum1Ptnw_1(0.0),
    fSum1Ptnw_2(0.0),
    /* histograms */
    fhN1_1_vsPt(nullptr),
    fhN1_1_vsEtaPhi(nullptr),
    fhSum1Pt_1_vsEtaPhi(nullptr),
    fhN1_1_vsZEtaPhiPt(nullptr),
    fhN1_2_vsPt(nullptr),
    fhN1_2_vsEtaPhi(nullptr),
    fhSum1Pt_2_vsEtaPhi(nullptr),
    fhN1_2_vsZEtaPhiPt(nullptr),
    fhN2_12_vsPtPt{nullptr,nullptr,nullptr,nullptr},
    fhN2_12_vsDEtaDPhi{nullptr,nullptr,nullptr,nullptr},
    fhSum2PtPt_12_vsDEtaDPhi{nullptr,nullptr,nullptr,nullptr},
    fhSum2DptDpt_12_vsDEtaDPhi{nullptr,nullptr,nullptr,nullptr},
    /* vs centrality profiles */
    fhN1_1_vsC(nullptr),
    fhSum1Pt_1_vsC(nullptr),
    fhN1nw_1_vsC(nullptr),
    fhSum1Ptnw_1_vsC(nullptr),
    fhN1_2_vsC(nullptr),
    fhSum1Pt_2_vsC(nullptr),
    fhN1nw_2_vsC(nullptr),
    fhSum1Ptnw_2_vsC(nullptr),
    fhN2_12_vsC{nullptr,nullptr,nullptr,nullptr},
    fhSum2PtPt_12_vsC{nullptr,nullptr,nullptr,nullptr},
    fhSum2DptDpt_12_vsC{nullptr,nullptr,nullptr,nullptr},
    fhN2nw_12_vsC{nullptr,nullptr,nullptr,nullptr},
    fhSum2PtPtnw_12_vsC{nullptr,nullptr,nullptr,nullptr},
    fhSum2DptDptnw_12_vsC{nullptr,nullptr,nullptr,nullptr},
    fhResonanceRoughMasses(nullptr),
    fhResonanceMasses(nullptr),
    fhDiscardedResonanceMasses(nullptr)
{
  for (Int_t ires = 0; ires < 16; ires++)
    fThresholdMult[ires] = 0x0;
}

/// Normal constructor
/// \param name the name for the object instance
Ali2PCorrelations::Ali2PCorrelations(const char *name) :
    TNamed(name,name),
    fOutput(nullptr),
    fVertexZ(0.0),
    fIxVertexZ(0),
    fCentrality(0.0),
    fSinglesOnly(true),
    fUseWeights(false),
    fUseSimulation(false),
    /* the arrays with tracks 1 and 2 information */
    fArraySize(5*1024),
    fNoOfTracks1(0),
    fId_1(nullptr),
    fCharge_1(nullptr),
    fIxEta_1(nullptr),
    fIxPhi_1(nullptr),
    fIxPt_1(nullptr),
    fFlags_1(nullptr),
    fPt_1(nullptr),
    fEta_1(nullptr),
    fPhi_1(nullptr),
    fCorrection_1(nullptr),
    fNoOfTracks2(0),
    fId_2(nullptr),
    fCharge_2(nullptr),
    fIxEta_2(nullptr),
    fIxPhi_2(nullptr),
    fIxPt_2(nullptr),
    fFlags_2(nullptr),
    fPt_2(nullptr),
    fEta_2(nullptr),
    fPhi_2(nullptr),
    fCorrection_2(nullptr),
    /* correction weights */
    fCorrectionWeights_1(nullptr),
    fCorrectionWeights_2(nullptr),
    /* efficiency correction */
    fEfficiencyCorrection_1(nullptr),
    fEfficiencyCorrection_2(nullptr),
    fPairsEfficiency_PP(nullptr),
    fPairsEfficiency_PM(nullptr),
    fPairsEfficiency_MM(nullptr),
    fPairsEfficiency_MP(nullptr),
    /* simulation pdfs */
    fPositiveTrackPdfs(nullptr),
    fNegativeTrackPdfs(nullptr),
    fPositiveTrackCurrentPdf(nullptr),
    fNegativeTrackCurrentPdf(nullptr),
    fSimEventsPerEvent(1),
    /* vertex bins */
    fNBins_vertexZ(40),
    fMin_vertexZ(-10.0),
    fMax_vertexZ(10.0),
    fWidth_vertexZ(0.5),
    /* phi origin shift */
    fNBinsPhiShift(0.0),
    /* pT1 bins */
    fNBins_pt_1(18),
    fMin_pt_1(0.2),
    fMax_pt_1(2.0),
    fWidth_pt_1(0.1),
    /* phi1 bins */
    fNBins_phi_1(72),
    fMin_phi_1(0.0),
    fMax_phi_1(TMath::Pi()*2.0),
    fWidth_phi_1(TMath::Pi()*2.0/72.0),
    /* eta1 bins */
    fNBins_eta_1(20),
    fMin_eta_1(-1.0),
    fMax_eta_1(1.0),
    fWidth_eta_1(0.1),
    fNBins_etaPhi_1(1440),
    fNBins_etaPhiPt_1(25920),
    fNBins_zEtaPhiPt_1(1036800),
    /* pT2 bins */
    fNBins_pt_2(18),
    fMin_pt_2(0.2),
    fMax_pt_2(2.0),
    fWidth_pt_2(0.1),
    /* phi2 bins */
    fNBins_phi_2(72),
    fMin_phi_2(0.0),
    fMax_phi_2(TMath::Pi()*2.0),
    fWidth_phi_2(TMath::Pi()*2.0/72.0),
    /* eta2 bins */
    fNBins_eta_2(20),
    fMin_eta_2(-1.0),
    fMax_eta_2(1.0),
    fWidth_eta_2(0.1),
    fNBins_etaPhi_2(1440),
    fNBins_etaPhiPt_2(25920),
    fNBins_zEtaPhiPt_2(1036800),
    /* pair bins */
    fNBins_deltaphi(72),
    fMin_deltaphi(0.0),
    fMax_deltaphi(TMath::Pi()*2.0),
    fWidth_deltaphi(TMath::Pi()*2.0/72.0),
    fNBins_deltaeta(20*2-1),
    fMin_deltaeta(-2.0),
    fMax_deltaeta(2.0),
    fWidth_deltaeta(4.0/(20*2-1)),

    /* accumulated values */
    fN1_1(0.0),
    fN1_2(0.0),
    fNnw1_1(0.0),
    fNnw1_2(0.0),
    fSum1Pt_1(0.0),
    fSum1Pt_2(0.0),
    fSum1Ptnw_1(0.0),
    fSum1Ptnw_2(0.0),
    /* histograms */
    fhN1_1_vsPt(nullptr),
    fhN1_1_vsEtaPhi(nullptr),
    fhSum1Pt_1_vsEtaPhi(nullptr),
    fhN1_1_vsZEtaPhiPt(nullptr),
    fhN1_2_vsPt(nullptr),
    fhN1_2_vsEtaPhi(nullptr),
    fhSum1Pt_2_vsEtaPhi(nullptr),
    fhN1_2_vsZEtaPhiPt(nullptr),
    fhN2_12_vsPtPt{nullptr,nullptr,nullptr,nullptr},
    fhN2_12_vsDEtaDPhi{nullptr,nullptr,nullptr,nullptr},
    fhSum2PtPt_12_vsDEtaDPhi{nullptr,nullptr,nullptr,nullptr},
    fhSum2DptDpt_12_vsDEtaDPhi{nullptr,nullptr,nullptr,nullptr},
    /* vs centrality profiles */
    fhN1_1_vsC(nullptr),
    fhSum1Pt_1_vsC(nullptr),
    fhN1nw_1_vsC(nullptr),
    fhSum1Ptnw_1_vsC(nullptr),
    fhN1_2_vsC(nullptr),
    fhSum1Pt_2_vsC(nullptr),
    fhN1nw_2_vsC(nullptr),
    fhSum1Ptnw_2_vsC(nullptr),
    fhN2_12_vsC{nullptr,nullptr,nullptr,nullptr},
    fhSum2PtPt_12_vsC{nullptr,nullptr,nullptr,nullptr},
    fhSum2DptDpt_12_vsC{nullptr,nullptr,nullptr,nullptr},
    fhN2nw_12_vsC{nullptr,nullptr,nullptr,nullptr},
    fhSum2PtPtnw_12_vsC{nullptr,nullptr,nullptr,nullptr},
    fhSum2DptDptnw_12_vsC{nullptr,nullptr,nullptr,nullptr},
    fhResonanceRoughMasses(nullptr),
    fhResonanceMasses(nullptr),
    fhDiscardedResonanceMasses(nullptr)
{
  for (Int_t ires = 0; ires < 16; ires++)
    fThresholdMult[ires] = 0x0;
}

/// \brief Default destructor
/// Deallocates the allocated memory
Ali2PCorrelations::~Ali2PCorrelations() {
  /* track 1 storage */
  delete [] fId_1;
  delete [] fCharge_1;
  delete [] fIxEta_1;
  delete [] fIxPhi_1;
  delete [] fIxPt_1;
  delete [] fFlags_1;
  delete [] fPt_1;
  delete [] fEta_1;
  delete [] fPhi_1;
  delete [] fCorrection_1;

  /* track 2 storage */
  delete [] fId_2;
  delete [] fCharge_2;
  delete [] fIxEta_2;
  delete [] fIxPhi_2;
  delete [] fIxPt_2;
  delete [] fFlags_2;
  delete [] fPt_2;
  delete [] fEta_2;
  delete [] fPhi_2;
  delete [] fCorrection_2;

  if (fCorrectionWeights_1 != NULL) delete [] fCorrectionWeights_1;
  if (fCorrectionWeights_2 != NULL) delete [] fCorrectionWeights_2;
  if (fEfficiencyCorrection_1 != NULL) delete [] fEfficiencyCorrection_1;
  if (fEfficiencyCorrection_2 != NULL) delete [] fEfficiencyCorrection_2;
  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
      delete fOutput;
  }
}


/// \brief Establishes the binning configuration
/// \param confstring string containing the binning configuration parameters
Bool_t Ali2PCorrelations::ConfigureBinning(const char *confstring) {
#define DPTDPTCORRBINCONFIGPAR 10

  Double_t min_pt = 0.0, max_pt = 0.0, width_pt = 0.0;
  Double_t min_eta = 0.0, max_eta = 0.0, width_eta = 0.0;
  Int_t    nBins_phi = 0;

  /* few sanity checks */
  TString str = confstring;
  TObjArray *array = str.Tokenize(";");
  for (Int_t item = 0; item < array->GetEntries(); item++) {
    const TString &stritem = ((TObjString*) array->At(item))->GetString();
    if (stritem.BeginsWith("phishift:")) {
      sscanf(stritem.Data(), "phishift:%lf", &fNBinsPhiShift);
    }
    else {
      TObjArray *a = stritem.Tokenize(",");
      if (a->GetEntries() != DPTDPTCORRBINCONFIGPAR) {
        delete a;
        return kFALSE;
      }
      sscanf(stritem, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d",
          &fMin_vertexZ, &fMax_vertexZ, &fWidth_vertexZ,
          &min_pt, &max_pt, &width_pt,
          &min_eta, &max_eta, &width_eta, &nBins_phi);
      delete a;
    }
  }
  delete array;

  fMin_pt_1 = fMin_pt_2 = min_pt;
  fMax_pt_1 = fMax_pt_2 = max_pt;
  fWidth_pt_1 = fWidth_pt_2 = width_pt;
  fMin_eta_1 = fMin_eta_2 = min_eta;
  fMax_eta_1 = fMax_eta_2 = max_eta;
  fWidth_eta_1 = fWidth_eta_2 = width_eta;
  fNBins_phi_1 = fNBins_phi_2 = nBins_phi;
  AliInfo("=====================================================");
  AliInfo(Form("Configured binning: %s", GetBinningConfigurationString().Data()));
  AliInfo("=====================================================");
  return kTRUE;
}

/// \brief Establishes the resonances rejection configuration
/// \param confstring string containing the resonances rejection configuration parameters
/// Basically one digit per supported resonance and the digit is the factor in one fourth
/// of the modulus around the resonance mass
void Ali2PCorrelations::ConfigureResonances(const char *confstring) {

  /* few sanity checks */
  TString str = confstring;
  if (!str.Contains("resonances:"))
    return;

  Int_t rescode;
  sscanf(str.Data(), "resonances:%d", &rescode);
  Int_t mult = 1;
  for (Int_t ires = 0; ires < fgkNoOfResonances; ires++) {
    fThresholdMult[ires] = Int_t(rescode / mult) % 10;
    mult *= 10;
  }
  AliInfo("=====================================================");
  AliInfo(Form("Configured resonance rejection cuts with string %s", GetResonancesConfigurationString().Data()));
  AliInfo("=====================================================");
}

/// \brief Build the configuration string
/// \return the configuration string corresponding to the current configuration
TString Ali2PCorrelations::GetBinningConfigurationString() const {
  if (fMin_pt_1 != fMin_pt_2 || fMax_pt_2 != fMax_pt_2 || fWidth_pt_1 != fWidth_pt_2 ||
      fMin_eta_1 != fMin_eta_2 || fMax_eta_1 != fMax_eta_2 || fWidth_eta_1 != fWidth_eta_2 ||
      fNBins_phi_1 != fNBins_phi_2) {
    return TString(Form("WrongAsymmetricBinning:phishift:%.1f;%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%d",
        fNBinsPhiShift,
        fMin_vertexZ, fMax_vertexZ, fWidth_vertexZ,
        fMin_pt_1, fMax_pt_1, fWidth_pt_1,
        fMin_eta_1, fMax_eta_1, fWidth_eta_1, fNBins_phi_1));
  }
  else {
    return TString(Form("Binning:phishift:%.1f;%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%d",
        fNBinsPhiShift,
        fMin_vertexZ, fMax_vertexZ, fWidth_vertexZ,
        fMin_pt_1, fMax_pt_1, fWidth_pt_1,
        fMin_eta_1, fMax_eta_1, fWidth_eta_1, fNBins_phi_1));
  }
}

/// \brief Build the resonances rejection configuration string
/// \return the configuration string corresponding to the current resonance rejection configuration
TString Ali2PCorrelations::GetResonancesConfigurationString() const {
  TString str = "resonances:";

  for (Int_t ires = 0; ires < fgkNoOfResonances; ires++) {
    str += TString::Format("%01d",fThresholdMult[fgkNoOfResonances - 1 - ires]);
  }
  return str;
}



/// \brief Initializes the member data structures
/// Allocates the needed memory an create the output histograms.
void Ali2PCorrelations::Initialize()
{
  Bool_t oldstatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fOutput = new TList();
  fOutput->SetName(GetName());
  fOutput->SetOwner();


  fNBins_vertexZ     = Int_t(0.5 + (fMax_vertexZ - fMin_vertexZ) / fWidth_vertexZ);
  fNBins_pt_1        = Int_t(0.5 + (fMax_pt_1 - fMin_pt_1 ) / fWidth_pt_1);
  fNBins_eta_1       = Int_t(0.5 + (fMax_eta_1 - fMin_eta_1) / fWidth_eta_1);
  fWidth_phi_1       = (fMax_phi_1  - fMin_phi_1) / fNBins_phi_1;
  fNBins_etaPhi_1    = fNBins_phi_1 * fNBins_eta_1;
  fNBins_etaPhiPt_1  = fNBins_etaPhi_1 * fNBins_pt_1;
  fNBins_zEtaPhiPt_1 = fNBins_vertexZ * fNBins_etaPhiPt_1;
  fNBins_pt_2        = Int_t(0.5 + (fMax_pt_2 - fMin_pt_2 ) / fWidth_pt_2);
  fNBins_eta_2       = Int_t(0.5 + (fMax_eta_2 - fMin_eta_2) / fWidth_eta_2);
  fWidth_phi_2       = (fMax_phi_2  - fMin_phi_2) / fNBins_phi_2;
  fNBins_etaPhi_2    = fNBins_phi_2 * fNBins_eta_2;
  fNBins_etaPhiPt_2  = fNBins_etaPhi_2 * fNBins_pt_2;
  fNBins_zEtaPhiPt_2 = fNBins_vertexZ * fNBins_etaPhiPt_2;
  fNBins_deltaeta    = fNBins_eta_1*2 -1;
  fMin_deltaeta      = fMin_eta_1 - fMax_eta_1;
  fMax_deltaeta      = fMax_eta_1 - fMin_eta_1;
  fWidth_deltaeta    = (fMax_deltaeta - fMin_deltaeta) / fNBins_deltaeta;
  fNBins_deltaphi    = fNBins_phi_1;
  fMin_deltaphi      = 0.0;
  fMax_deltaphi      = TMath::TwoPi();
  fWidth_deltaphi    = (fMax_deltaphi - fMin_deltaphi) / fNBins_deltaphi;

  /* incorporate configuration parameters */
  fOutput->Add(new TParameter<Bool_t>("DifferentialOutput",true,'f'));
  fOutput->Add(new TParameter<Int_t>("NoBinsVertexZ",fNBins_vertexZ,'f'));
  fOutput->Add(new TParameter<Int_t>("NoBinsPt",fNBins_pt_1,'f'));
  fOutput->Add(new TParameter<Int_t>("NoBinsEta",fNBins_eta_1,'f'));
  fOutput->Add(new TParameter<Int_t>("NoBinsPhi",fNBins_phi_1,'f'));
  fOutput->Add(new TParameter<Double_t>("MinVertexZ",fMin_vertexZ,'f'));
  fOutput->Add(new TParameter<Double_t>("MaxVertexZ",fMax_vertexZ,'f'));
  fOutput->Add(new TParameter<Double_t>("MinPt",fMin_pt_1,'f'));
  fOutput->Add(new TParameter<Double_t>("MaxPt",fMax_pt_1,'f'));
  fOutput->Add(new TParameter<Double_t>("MinEta",fMin_eta_1,'f'));
  fOutput->Add(new TParameter<Double_t>("MaxEta",fMax_eta_1,'f'));
  fOutput->Add(new TParameter<Double_t>("MinPhi",fMin_phi_1,'f'));
  fOutput->Add(new TParameter<Double_t>("MaxPhi",fMax_phi_1,'f'));

  /* incorporate the resonance rejection configuration parameter */
  Int_t rescode = 0;
  Int_t mult = 1;
  for (Int_t ires = 0; ires < fgkNoOfResonances; ires++) {
    rescode += fThresholdMult[ires] * mult;
    mult *= 10;
  }
  fOutput->Add(new TParameter<Int_t>("ResonancesCode",rescode,'f'));

  /* after the parameters dump the proper phi limits are set according to the phi shift */
  fMax_phi_1         = fMax_phi_1 - fWidth_phi_1 * fNBinsPhiShift;
  fMin_phi_1         = fMin_phi_1 - fWidth_phi_1 * fNBinsPhiShift;
  fMax_phi_2         = fMax_phi_2 - fWidth_phi_2 * fNBinsPhiShift;
  fMin_phi_2         = fMin_phi_2 - fWidth_phi_2 * fNBinsPhiShift;


  /* track 1 storage */
  fId_1 = new Int_t[fArraySize];
  fCharge_1 = new Int_t[fArraySize];
  fIxEta_1 = new Int_t[fArraySize];
  fIxPhi_1 = new Int_t[fArraySize];
  fIxPt_1 = new Int_t[fArraySize];
  fFlags_1 = new UInt_t[fArraySize];
  fPt_1 = new Float_t[fArraySize];
  fEta_1 = new Float_t[fArraySize];
  fPhi_1 = new Float_t[fArraySize];
  fCorrection_1 = new Float_t[fArraySize];

  /* track 2 storage */
  fId_2 = new Int_t[fArraySize];
  fCharge_2 = new Int_t[fArraySize];
  fIxEta_2 = new Int_t[fArraySize];
  fIxPhi_2 = new Int_t[fArraySize];
  fIxPt_2 = new Int_t[fArraySize];
  fFlags_2 = new UInt_t[fArraySize];
  fPt_2 = new Float_t[fArraySize];
  fEta_2 = new Float_t[fArraySize];
  fPhi_2 = new Float_t[fArraySize];
  fCorrection_2 = new Float_t[fArraySize];

  if (fSinglesOnly)  {
    fhN1_1_vsPt = new TH1F("n1_1_vsPt","#LT n_{1} #GT;p_{t,1} (GeV/c);#LT n_{1} #GT", fNBins_pt_1,  fMin_pt_1,  fMax_pt_1);
    fhN1_1_vsZEtaPhiPt = new TH3F("n1_1_vsZ_vsEtaPhi_vsPt","#LT n_{1} #GT;vtx_{z};#eta_{1}#times#varphi_{1};p_{t,1} (GeV/c)",
        fNBins_vertexZ,fMin_vertexZ,fMax_vertexZ, fNBins_etaPhi_1,0.0,Double_t(fNBins_etaPhi_1),fNBins_pt_1,fMin_pt_1,fMax_pt_1);

    fhN1_2_vsPt = new TH1F("n1_2_vsPt","#LT n_{1} #GT;p_{t,2} (GeV/c);#LT n_{1} #GT", fNBins_pt_2,  fMin_pt_2,  fMax_pt_2);
    fhN1_2_vsZEtaPhiPt = new TH3F("n1_2_vsZ_vsEtaPhi_vsPt","#LT n_{2} #GT;vtx_{z};#eta_{2}#times#varphi_{2};p_{t,2} (GeV/c)",
        fNBins_vertexZ,fMin_vertexZ,fMax_vertexZ, fNBins_etaPhi_2,0.0,Double_t(fNBins_etaPhi_2),fNBins_pt_2,fMin_pt_2,fMax_pt_2);

    fOutput->Add(fhN1_1_vsPt);
    fOutput->Add(fhN1_1_vsZEtaPhiPt);
    fOutput->Add(fhN1_2_vsPt);
    fOutput->Add(fhN1_2_vsZEtaPhiPt);
  }
  else {
    fhResonanceRoughMasses = new TH2F("ResonanceRoughMasses","Approx invariant mass;res id;Inv mass", fgkNoOfResonances, -0.5, fgkNoOfResonances - 0.5, 1000, 0.0, 5.0);
    fhResonanceMasses = new TH2F("ResonanceMasses","Invariant mass;res id;Inv mass", fgkNoOfResonances, -0.5, fgkNoOfResonances - 0.5, 1000, 0.0, 5.0);
    fhDiscardedResonanceMasses = new TH2F("DiscardedResonanceMasses","Discarded invariant mass;res id;Inv mass", fgkNoOfResonances, -0.5, fgkNoOfResonances - 0.5, 1000, 0.0, 5.0);

    /* histograms for track one */
    fhN1_1_vsEtaPhi = new TH2F("n1_1_vsEtaPhi","#LT n_{1} #GT;#eta_{1};#varphi_{1} (radian);#LT n_{1} #GT",
        fNBins_eta_1, fMin_eta_1, fMax_eta_1,  fNBins_phi_1, fMin_phi_1, fMax_phi_1);
    fhSum1Pt_1_vsEtaPhi = new TH2F("sumPt_1_vsEtaPhi","#LT #Sigma p_{t,1} #GT;#eta_{1};#varphi_{1} (radian);#LT #Sigma p_{t,1} #GT (GeV/c)",
        fNBins_eta_1,fMin_eta_1,fMax_eta_1,fNBins_phi_1,fMin_phi_1,fMax_phi_1);
    fhN1_1_vsC = new TProfile("n1_1_vsM","#LT n_{1} #GT (weighted);Centrality (%);#LT n_{1} #GT",100,0.0,100.0);
    fhSum1Pt_1_vsC = new TProfile("sumPt_1_vsM","#LT #Sigma p_{t,1} #GT (weighted);Centrality (%);#LT #Sigma p_{t,1} #GT (GeV/c)",100,0.0,100.0);
    fhN1nw_1_vsC = new TProfile("n1Nw_1_vsM","#LT n_{1} #GT;Centrality (%);#LT n_{1} #GT",100,0.0,100.0);
    fhSum1Ptnw_1_vsC = new TProfile("sumPtNw_1_vsM","#LT #Sigma p_{t,1} #GT;Centrality (%);#LT #Sigma p_{t,1} #GT (GeV/c)",100,0.0,100.0);

    /* histograms for track two */
    fhN1_2_vsEtaPhi = new TH2F("n1_2_vsEtaPhi","#LT n_{1} #GT;#eta_{2};#varphi_{2} (radian);#LT n_{1} #GT",
        fNBins_eta_2, fMin_eta_2, fMax_eta_2,  fNBins_phi_2, fMin_phi_2, fMax_phi_2);
    fhSum1Pt_2_vsEtaPhi = new TH2F("sumPt_2_vsEtaPhi","#LT #Sigma p_{t,2} #GT;#eta_{2};#varphi_{2} (radian);#LT #Sigma p_{t,2} #GT (GeV/c)",
        fNBins_eta_2,fMin_eta_2,fMax_eta_2,fNBins_phi_2,fMin_phi_2,fMax_phi_2);
    fhN1_2_vsC = new TProfile("n1_2_vsM","#LT n_{1} #GT (weighted);Centrality (%);#LT n_{1} #GT",100,0.0,100.0);
    fhSum1Pt_2_vsC = new TProfile("sumPt_2_vsM","#LT #Sigma p_{t,1} #GT (weighted);Centrality (%);#LT #Sigma p_{t,1} #GT (GeV/c)",100,0.0,100.0);
    fhN1nw_2_vsC = new TProfile("n1Nw_2_vsM","#LT n_{1} #GT;Centrality (%);#LT n_{1} #GT",100,0.0,100.0);
    fhSum1Ptnw_2_vsC = new TProfile("sumPtNw_2_vsM","#LT #Sigma p_{t,1} #GT;Centrality (%);#LT #Sigma p_{t,1} #GT (GeV/c)",100,0.0,100.0);

    /* histograms for track pairs */
    for (int i = 0; i<nTrackPairs; ++i) {
      const char *pname = TrackPairsNames[i];
      fhN2_12_vsDEtaDPhi[i] = new TH2F(TString::Format("n2_12_vsDEtaDPhi_%s",pname),TString::Format("#LT n_{2} #GT (%s);#Delta#eta;#Delta#varphi;#LT n_{2} #GT",pname),
          fNBins_deltaeta,fMin_deltaeta,fMax_deltaeta,fNBins_deltaphi,fMin_deltaphi,fMax_deltaphi);
      fhSum2PtPt_12_vsDEtaDPhi[i] = new TH2F(TString::Format("sumPtPt_12_vsDEtaDPhi_%s",pname),TString::Format("#LT #Sigma p_{t,1}p_{t,2} #GT (%s);#Delta#eta;#Delta#varphi;#LT #Sigma p_{t,1}p_{t,2} #GT (GeV^{2})",pname),
          fNBins_deltaeta,fMin_deltaeta,fMax_deltaeta,fNBins_deltaphi,fMin_deltaphi,fMax_deltaphi);
      fhSum2DptDpt_12_vsDEtaDPhi[i] = new TH2F(TString::Format("sumDptDpt_12_vsDEtaDPhi_%s",pname),TString::Format("#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (%s);#Delta#eta;#Delta#varphi;#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (GeV^{2})",pname),
          fNBins_deltaeta,fMin_deltaeta,fMax_deltaeta,fNBins_deltaphi,fMin_deltaphi,fMax_deltaphi);
      fhN2_12_vsPtPt[i] = new TH2F(TString::Format("n2_12_vsPtVsPt_%s",pname),TString::Format("#LT n_{2} #GT (%s);p_{t,1} (GeV/c);p_{t,2} (GeV/c);#LT n_{2} #GT",pname),
          fNBins_pt_1,fMin_pt_1,fMax_pt_1,fNBins_pt_2,fMin_pt_2,fMax_pt_2);
      fhN2_12_vsC[i] = new TProfile(TString::Format("n2_12_vsM_%s",pname),TString::Format("#LT n_{2} #GT (%s) (weighted);Centrality (%%);#LT n_{2} #GT",pname),100,0.0,100.0);
      fhSum2PtPt_12_vsC[i] = new TProfile(TString::Format("sumPtPt_12_vsM_%s",pname),TString::Format("#LT #Sigma p_{t,1}p_{t,2} #GT (%s) (weighted);Centrality (%%);#LT #Sigma p_{t,1}p_{t,2} #GT (GeV^{2})",pname),100,0.0,100.0);
      fhSum2DptDpt_12_vsC[i] = new TProfile(TString::Format("sumDptDpt_12_vsM_%s",pname),TString::Format("#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (%s) (weighted);Centrality (%%);#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (GeV^{2})",pname),100,0.0,100.0);
      fhN2nw_12_vsC[i] = new TProfile(TString::Format("n2Nw_12_vsM_%s",pname),TString::Format("#LT n_{2} #GT (%s);Centrality (%%);#LT n_{2} #GT",pname),100,0.0,100.0);
      fhSum2PtPtnw_12_vsC[i] = new TProfile(TString::Format("sumPtPtNw_12_vsM_%s",pname),TString::Format("#LT #Sigma p_{t,1}p_{t,2} #GT (%s);Centrality (%%);#LT #Sigma p_{t,1}p_{t,2} #GT (GeV^{2})",pname),100,0.0,100.0);
      fhSum2DptDptnw_12_vsC[i] = new TProfile(TString::Format("sumDptDptNw_12_vsM_%s",pname),TString::Format("#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (%s);Centrality (%%);#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (GeV^{2})",pname),100,0.0,100.0);
    }

    fOutput->Add(fhResonanceRoughMasses);
    fOutput->Add(fhResonanceMasses);
    fOutput->Add(fhDiscardedResonanceMasses);

    fOutput->Add(fhN1_1_vsEtaPhi);
    fOutput->Add(fhSum1Pt_1_vsEtaPhi);
    fOutput->Add(fhN1_1_vsC);
    fOutput->Add(fhSum1Pt_1_vsC);
    fOutput->Add(fhN1nw_1_vsC);
    fOutput->Add(fhSum1Ptnw_1_vsC);
    fOutput->Add(fhN1_2_vsEtaPhi);
    fOutput->Add(fhSum1Pt_2_vsEtaPhi);
    fOutput->Add(fhN1_2_vsC);
    fOutput->Add(fhSum1Pt_2_vsC);
    fOutput->Add(fhN1nw_2_vsC);
    fOutput->Add(fhSum1Ptnw_2_vsC);

    for (int i = 0; i<nTrackPairs; ++i) {
      fOutput->Add(fhN2_12_vsDEtaDPhi[i]);
      fOutput->Add(fhSum2PtPt_12_vsDEtaDPhi[i]);
      fOutput->Add(fhSum2DptDpt_12_vsDEtaDPhi[i]);
      fOutput->Add(fhN2_12_vsPtPt[i]);
      fOutput->Add(fhN2_12_vsC[i]);
      fOutput->Add(fhSum2PtPt_12_vsC[i]);
      fOutput->Add(fhSum2DptDpt_12_vsC[i]);
      fOutput->Add(fhN2nw_12_vsC[i]);
      fOutput->Add(fhSum2PtPtnw_12_vsC[i]);
      fOutput->Add(fhSum2DptDptnw_12_vsC[i]);
    }
  }

  TH1::AddDirectory(oldstatus);
}

/// \brief Stores the correction weights for both set of tracks
/// \param h3_1 the calibration weights for the first track
/// \param h3_2 the calibration weights for the second track
/// \return kTRUE if correctly done kFALSE otherwise
Bool_t Ali2PCorrelations::SetWeigths(const TH3F *h3_1, const TH3F *h3_2) {
  Bool_t done = kTRUE;

  if (fUseWeights){
    if (h3_1) {
      /* allocate memory for track 1 weights */
      fCorrectionWeights_1 = new Float_t[fNBins_vertexZ*fNBins_etaPhi_1*fNBins_pt_1];

      for (Int_t ixZ = 0; ixZ < fNBins_vertexZ; ixZ++) {
        Double_t zval = fMin_vertexZ + fWidth_vertexZ*(ixZ+0.5);
        for (Int_t ixEtaPhi=0; ixEtaPhi < fNBins_etaPhi_1; ixEtaPhi++) {
          for (Int_t ixPt=0; ixPt < fNBins_pt_1; ixPt++) {
            Double_t pTval = fMin_pt_1 + fWidth_pt_1*(ixPt+0.5);
            fCorrectionWeights_1[ixZ*fNBins_etaPhi_1*fNBins_pt_1+ixEtaPhi*fNBins_pt_1+ixPt] = h3_1->GetBinContent(h3_1->GetXaxis()->FindBin(zval),ixEtaPhi+1,h3_1->GetZaxis()->FindBin(pTval));
          }
        }
      }
    } // _weight_1
    else {
      AliFatal("The weights one histogram is a null pointer. ABORTING!!!");
      done = kFALSE;
    }

    if (h3_2) {
      /* allocate memory for track 1 weights */
      fCorrectionWeights_2 = new Float_t[fNBins_vertexZ*fNBins_etaPhi_2*fNBins_pt_2];

      for (Int_t ixZ = 0; ixZ < fNBins_vertexZ; ixZ++) {
        Double_t zval = fMin_vertexZ + fWidth_vertexZ*(ixZ+0.5);
        for (Int_t ixEtaPhi=0; ixEtaPhi < fNBins_etaPhi_2; ixEtaPhi++) {
          for (Int_t ixPt=0; ixPt < fNBins_pt_2; ixPt++) {
            Double_t pTval = fMin_pt_2 + fWidth_pt_2*(ixPt+0.5);
            fCorrectionWeights_2[ixZ*fNBins_etaPhi_2*fNBins_pt_2+ixEtaPhi*fNBins_pt_2+ixPt] = h3_2->GetBinContent(h3_2->GetXaxis()->FindBin(zval),ixEtaPhi+1,h3_2->GetZaxis()->FindBin(pTval));
          }
        }
      }
    } // _weight_2
    else {
      AliFatal("The weights two histogram is a null pointer. ABORTING!!!");
      done = kFALSE;
    }
  }
  else {
    AliError("Setting weights for a not use weights instance. Ignoring it");
    done = kFALSE;
  }
  return done;
}

/// \brief Stores the efficiency correction for both set of tracks
/// If the efficiency correction histograms are not present a efficiency correction value 1.0 is stored.
/// The correction value are always the inverse of the efficiency ones
///
/// Usually the efficiency is incorporated in the weights but this method provides
/// an additional way to incorporate pT dependent efficiency.
/// \param h1_1 histogram with efficiency values for track one
/// \param h1_2 histogram with efficiency values for track two
/// \return kTRUE if everything went ok otherwise kFALSE
Bool_t Ali2PCorrelations::SetEfficiencyCorrection(const TH1F *h1_1, const TH1F *h1_2) {
  Bool_t done = kTRUE;

  if (h1_1 != NULL || h1_2 != NULL){
    if (h1_1) {
      /* allocate memory for track 1 efficiency correction */
      fEfficiencyCorrection_1 = new Float_t[fNBins_pt_1];

      for (Int_t ixPt=0; ixPt < fNBins_pt_1; ixPt++) {
        Int_t bin = h1_1->GetXaxis()->FindBin(fMin_pt_1 + fWidth_pt_1/2.0 + ixPt*fWidth_pt_1);
        fEfficiencyCorrection_1[ixPt] = 1.0 / h1_1->GetBinContent(bin);
      }
    }
    else {
      AliFatal("The efficiency correction histogram one is a null pointer. ABORTING!!!");
      done = kFALSE;
    }

    if (h1_2) {
      /* allocate memory for track 1 efficiency correction */
      fEfficiencyCorrection_2 = new Float_t[fNBins_pt_2];

      for (Int_t ixPt=0; ixPt < fNBins_pt_2; ixPt++) {
        Int_t bin = h1_2->GetXaxis()->FindBin(fMin_pt_2 + fWidth_pt_2/2.0 + ixPt*fWidth_pt_2);
        fEfficiencyCorrection_2[ixPt] = 1.0 / h1_2->GetBinContent(bin);
      }
    }
    else {
      AliFatal("The efficiency correction histogram two is a null pointer. ABORTING!!!");
      done = kFALSE;
    }
  }
  else {
    AliInfo("Setting efficiency correction equal to one for both tracks");

    fEfficiencyCorrection_1 = new Float_t[fNBins_pt_1];
    for (Int_t ixPt=0; ixPt < fNBins_pt_1; ixPt++) {
      fEfficiencyCorrection_1[ixPt] = 1.0;
    }
    fEfficiencyCorrection_2 = new Float_t[fNBins_pt_2];
    for (Int_t ixPt=0; ixPt < fNBins_pt_2; ixPt++) {
      fEfficiencyCorrection_2[ixPt] = 1.0;
    }
    done = kTRUE;
  }
  AliInfo(Form("Configured efficiency correcton on %s",GetName()));
  printf("Track one:\n");
  for (Int_t bin = 0; bin < fNBins_pt_1; bin++) {
    printf("%f, ", fEfficiencyCorrection_1[bin]);
    if ((bin+1)%8 == 0) printf("\n");
  }
  printf("\n");
  printf("Track two:\n");
  for (Int_t bin = 0; bin < fNBins_pt_2; bin++) {
    printf("%f, ", fEfficiencyCorrection_2[bin]);
    if ((bin+1)%8 == 0) printf("\n");
  }
  printf("\n");
  return done;
}

/// \brief Stores the pairs efficiency correction for the four tracks combinations
/// The correction value is stored as the inverse of the efficiency ones
/// \param h11 histogram with efficiency values for the one-one pair
/// \param h12 histogram with efficiency values for the one-two pair
/// \param h22 histogram with efficiency values for the two-two pair
/// \param h21 histogram with efficiency values for the two-one pair
/// \return kTRUE if everything went ok otherwise kFALSE
Bool_t Ali2PCorrelations::SetPairEfficiencyCorrection(const THn *h11, const THn *h12, const THn *h22, const THn *h21) {
  AliInfo("");

  Bool_t done = kTRUE;

  if (h11 != NULL || h12 != NULL || h22 != NULL || h21 != NULL){
    if (h11) {
      if (h11->GetNdimensions() != 4) {
        AliFatal("Pair efficiency correction expected to have four dimensions. ABORTING!!!");
        done = kFALSE;
      }
      else {
        /* store the pair 11 efficiency histogram */
        fPairsEfficiency_PP = h11;
      }
    }
    else {
      AliFatal("The pair efficiency correction histogram one-one is a null pointer. ABORTING!!!");
      done = kFALSE;
    }

    if (h12) {
      if (h12->GetNdimensions() != 4) {
        AliFatal("Pair efficiency correction expected to have four dimensions. ABORTING!!!");
        done = kFALSE;
      }
      else {
        /* store the pair 12 efficiency histogram */
        fPairsEfficiency_PM = h12;
      }
    }
    else {
      AliFatal("The pair efficiency correction histogram one-two is a null pointer. ABORTING!!!");
      done = kFALSE;
    }

    if (h22) {
      if (h22->GetNdimensions() != 4) {
        AliFatal("Pair efficiency correction expected to have four dimensions. ABORTING!!!");
        done = kFALSE;
      }
      else {
        /* store the pair 22 efficiency histogram */
        fPairsEfficiency_MM = h22;
      }
    }
    else {
      AliFatal("The pair efficiency correction histogram two-two is a null pointer. ABORTING!!!");
      done = kFALSE;
    }

    if (h21) {
      if (h21->GetNdimensions() != 4) {
        AliFatal("Pair efficiency correction expected to have four dimensions. ABORTING!!!");
        done = kFALSE;
      }
      else {
        /* store the pair 11 efficiency histogram */
        fPairsEfficiency_MP = h21;
      }
    }
    else {
      AliFatal("The pair efficiency correction histogram one-one is a null pointer. ABORTING!!!");
      done = kFALSE;
    }
  }
  return done;
}

/// \brief Stores the pdf for plus and minus tracks depending on z vertex bin
///
/// Pdf are used to generate simulated tracks
/// \param pluspdf array with z vertex dependent pdf for positive tracks
/// \param minuspdf array with z vertex dependent pdf for negative tracks
/// \return kTRUE if everything went ok otherwise kFALSE
Bool_t Ali2PCorrelations::SetSimultationPdfs(const TObjArray *pluspdf, const TObjArray *minuspdf) {

  Bool_t done = kTRUE;

  if (fUseSimulation) {
    if (pluspdf != NULL) {
      if (pluspdf->GetEntriesFast() == fNBins_vertexZ) {
        fPositiveTrackPdfs = pluspdf;
      }
      else {
        AliFatal(Form("The positive tracks density histogram array number of pdfs %d differ from the number of zvertex bins %d. ABORTING!!!",
            pluspdf->GetEntriesFast(), fNBins_vertexZ));
        done = kFALSE;
      }
    }
    else {
      AliFatal("The positive tracks density histogram array is a null pointer. ABORTING!!!");
      done = kFALSE;
    }

    if (minuspdf != NULL) {
      if (minuspdf->GetEntriesFast() == fNBins_vertexZ) {
        fNegativeTrackPdfs = minuspdf;
      }
      else {
        AliFatal(Form("The negative tracks density histogram array number of pdfs %d differ from the number of zvertex bins %d. ABORTING!!!",
            minuspdf->GetEntriesFast(), fNBins_vertexZ));
        done = kFALSE;
      }
    }
    else {
      AliFatal("The negative tracks density histogram array is a null pointer. ABORTING!!!");
      done = kFALSE;
    }
  }
  else {
    AliError("Setting tracks density profiles for a not use simulation instance. Ignoring it");
    done = kFALSE;
  }
  return done;
}


/// \brief initializes the object instance to start a new event
/// \param centrality the new event centrality in percentage
/// \param vertexz the new event vertex \f$z\f$ coordinate
/// \return kTRUE if correctly initialized kFALSE otherwise
Bool_t Ali2PCorrelations::StartEvent(Float_t centrality, Float_t vertexz) {

  fVertexZ = vertexz;
  fCentrality = centrality;

  if (fVertexZ < fMin_vertexZ || fMax_vertexZ < fVertexZ) {
    AliInfo(Form("Event with z vertex: %.2f out of internal cuts", fVertexZ));
    return kFALSE;
  }

  fIxVertexZ = Int_t ((fVertexZ - fMin_vertexZ) / fWidth_vertexZ);
  if (fIxVertexZ < 0 || !(fIxVertexZ < fNBins_vertexZ)) {
    AliError(Form("Event z vertex bin %d out of bounds. Ignoring event.", fIxVertexZ));
    return kFALSE;
  }

  if (fUseSimulation) {
    fPositiveTrackCurrentPdf = (TH3F *) fPositiveTrackPdfs->At(fIxVertexZ);
    fNegativeTrackCurrentPdf = (TH3F *) fNegativeTrackPdfs->At(fIxVertexZ);
  }

  fNoOfTracks1 = 0;
  fNoOfTracks2 = 0;

  /* reset single counters */
  fN1_1 = fN1_2 = fSum1Pt_1 = fSum1Pt_2 = fNnw1_1 = fNnw1_2 = fSum1Ptnw_1 = fSum1Ptnw_2 = 0;

  return kTRUE;
}

/// \brief process a track and store its parameters if feasible
///
/// If simulation is ordered the actual track is discarded and a new one with the
/// same charge is produced out of the corresponding track pdf
/// \param trkId the external track Id
/// \param trk the passed track
/// \return kTRUE if the track is properly handled kFALSE otherwise
Bool_t Ali2PCorrelations::ProcessTrack(Int_t trkId, AliVTrack *trk) {

  if (fUseSimulation) {
    Double_t pT;
    Double_t eta;
    Double_t phi;

    if (trk->Charge() < 0)
      fNegativeTrackCurrentPdf->GetRandom3(eta,phi,pT);
    else
      fPositiveTrackCurrentPdf->GetRandom3(eta,phi,pT);

    return ProcessTrack(trkId, Int_t(trk->Charge()), pT, eta, phi);
  }
  else
    return ProcessTrack(trkId, Int_t(trk->Charge()), Float_t(trk->Pt()), Float_t(trk->Eta()), Float_t(trk->Phi()));
}


/// \brief process a true particle and store its parameters if feasible
///
/// If simulation is orderd the track is discarded and kFALSE is returned
/// \param trkId the external particle Id
/// \param part the passed particle
/// \return kTRUE if the particle is properly handled kFALSE otherwise
Bool_t Ali2PCorrelations::ProcessTrack(Int_t trkId, AliVParticle *part) {

  if (fUseSimulation) {
    return kFALSE;
  }
  else
    if (part->Charge() != 0) {
      return ProcessTrack(trkId, (part->Charge() > 0) ? 1 : -1, Float_t(part->Pt()), Float_t(part->Eta()), Float_t(part->Phi()));
    }
    else {
      return kFALSE;
    }
}


/// \brief process a track and store its parameters if feasible
/// \param trkId the external track Id
/// \param charge the track charge
/// \param pT the track \f$ p_T \f$
/// \param eta the track \f$ \eta \f$
/// \param phi the track \f$ \phi \f$
/// \return kTRUE if the track is properly handled kFALSE otherwise
/// For the time being positive tracks go to bank one and negative tracks go to bank two
Bool_t Ali2PCorrelations::ProcessTrack(Int_t trkId, Int_t charge, Float_t pT, Float_t eta, Float_t ophi) {

  if(charge == 0) return kFALSE;

  /* track 1 */
  if ((charge > 0) && pT > fMin_pt_1 && pT < fMax_pt_1)
  {
    if (!(fNoOfTracks1 < fArraySize)) {
      AliError(Form("Storage for track one full: %d", fArraySize));
      return kFALSE;
    }

    /* consider a potential phi origin shift */
    Float_t phi = ophi;
    if (!(phi < fMax_phi_1)) phi = phi - 2*TMath::Pi();
    Int_t ixPhi = Int_t ((phi - fMin_phi_1) / fWidth_phi_1);
    if (ixPhi < 0 || !(ixPhi < fNBins_phi_1)) {
      AliWarning("Track one phi out of bins");
      return kFALSE;
    }

    if (eta < fMin_eta_1 ||  fMax_eta_1 < eta) {
      AliWarning(Form("Wrongly passed track one with eta: %.2f", eta));
      return kFALSE;
    }

    Int_t ixEta = Int_t ((eta - fMin_eta_1) / fWidth_eta_1);
    if (ixEta < 0 || !(ixEta < fNBins_eta_1)) {
      AliWarning(Form("Track one eta bin %d out of bounds", ixEta));
      return kFALSE;
    }


    if (pT < fMin_pt_1 ||  fMax_pt_1 < pT) {
      AliWarning(Form("Wrongly passed track one with pT: %.2f", pT));
      return kFALSE;
    }

    Int_t ixPt = Int_t((pT - fMin_pt_1 ) / fWidth_pt_1);
    if (ixPt < 0  || !(ixPt < fNBins_pt_1)) {
      AliWarning(Form("Track pT bin %d out of bounds",ixPt));
      return kFALSE;
    }


    Int_t ixEtaPhi = ixEta*fNBins_phi_1+ixPhi;
    Int_t ixVertexP1 = fIxVertexZ*fNBins_etaPhiPt_1;
    Int_t ixZEtaPhiPt = ixVertexP1 + ixEtaPhi*fNBins_pt_1 + ixPt;
    if (ixZEtaPhiPt < 0 || !(ixZEtaPhiPt < fNBins_zEtaPhiPt_1)) {
      AliWarning(Form("Event zvertex and track eta phi and pt bin %d out of bounds", ixZEtaPhiPt));
      return kFALSE;
    }

    Float_t effcorr = fEfficiencyCorrection_1[ixPt];
    Float_t corr;
    if (fCorrectionWeights_1)
      corr = fCorrectionWeights_1[ixZEtaPhiPt];
    else
      corr = 1.0;

    /* the final correction incorporates also the efficiency correction */
    /* and affects to both, track densities and pT */
    corr *= effcorr;

    if (fSinglesOnly) {
      fhN1_1_vsPt->Fill(pT,corr);
      fhN1_1_vsZEtaPhiPt->Fill(fVertexZ,ixEtaPhi+0.5,pT,corr);
    }
    else {
      Float_t corrPt                = corr*pT;
      fId_1[fNoOfTracks1]           = trkId;
      fCharge_1[fNoOfTracks1]       = charge;
      fIxEta_1[fNoOfTracks1]        = ixEta;
      fIxPhi_1[fNoOfTracks1]        = ixPhi;
      fIxPt_1[fNoOfTracks1]         = ixPt;
      fFlags_1[fNoOfTracks1]        = 0x0;
      fPt_1[fNoOfTracks1]           = pT;
      fEta_1[fNoOfTracks1]          = eta;
      fPhi_1[fNoOfTracks1]          = ophi;
      fCorrection_1[fNoOfTracks1]   = corr;
      fN1_1                         += corr;
      fhN1_1_vsEtaPhi->Fill(eta,phi,corr);
      fSum1Pt_1                     += corrPt;
      fhSum1Pt_1_vsEtaPhi->Fill(eta,phi,corrPt);
      fNnw1_1                       += 1;
      fSum1Ptnw_1                   += pT;
      fNoOfTracks1++;
      if (!(fNoOfTracks1 < fArraySize)) {
        AliError(Form("Storage for track one full: %d", fArraySize));
        return kFALSE;
      }
    }
  }


  /* track 2 */
  if ((charge < 0) && pT > fMin_pt_2 && pT < fMax_pt_2)
  {
    if (!(fNoOfTracks2 < fArraySize)) {
      AliError(Form("Storage for track two full: %d", fArraySize));
      return kFALSE;
    }

    /* consider a potential phi origin shift */
    Float_t phi = ophi;
    if (!(phi < fMax_phi_2)) phi = phi - 2*TMath::Pi();
    Int_t ixPhi = Int_t ((phi - fMin_phi_2) / fWidth_phi_2);
    if (ixPhi <0 || !(ixPhi < fNBins_phi_2)) {
      AliWarning("Track two phi out of bins");
      return kFALSE;
    }

    if (eta < fMin_eta_2 ||  fMax_eta_2 < eta) {
      AliWarning(Form("Wrongly passed track two with eta: %.2f", eta));
      return kFALSE;
    }

    Int_t ixEta = Int_t((eta - fMin_eta_2) / fWidth_eta_2);
    if (ixEta < 0 || !(ixEta < fNBins_eta_2)) {
      AliWarning(Form("Track two eta bin %d out of bounds", ixEta));
      return kFALSE;
    }

    if (pT < fMin_pt_2 ||  fMax_pt_2 < pT) {
      AliWarning(Form("Wrongly passed track two with pT: %.2f", pT));
      return kFALSE;
    }

    Int_t ixPt = Int_t((pT - fMin_pt_2 ) / fWidth_pt_2 );
    if (ixPt < 0 || !(ixPt < fNBins_pt_2)) {
      AliWarning(Form("Track two pT bin %d out of bounds",ixPt));
      return kFALSE;
    }

    Int_t ixEtaPhi = ixEta*fNBins_phi_2+ixPhi;
    Int_t ixVertexP2 = fIxVertexZ*fNBins_etaPhiPt_2;
    Int_t ixZEtaPhiPt = ixVertexP2 + ixEtaPhi*fNBins_pt_2 + ixPt;
    if (ixZEtaPhiPt < 0 || !(ixZEtaPhiPt < fNBins_zEtaPhiPt_2)) {
      AliWarning(Form("Event zvertex and track eta phi and pt bin %d out of bounds", ixZEtaPhiPt));
      return kFALSE;
    }

    Float_t effcorr = fEfficiencyCorrection_2[ixPt];
    Float_t corr;
    if (fCorrectionWeights_2)
      corr = fCorrectionWeights_2[ixZEtaPhiPt];
    else
      corr = 1.0;

    /* the final correction incorporates also the efficiency correction */
    /* and affects to both, track densities and pT */
    corr *= effcorr;

    if (fSinglesOnly) {
      fhN1_2_vsPt->Fill(pT,corr);
      fhN1_2_vsZEtaPhiPt->Fill(fVertexZ,ixEtaPhi+0.5,pT,corr);
    }
    else {
      Float_t corrPt                = corr*pT;
      fId_2[fNoOfTracks2]           = trkId;
      fCharge_2[fNoOfTracks2]       = charge;
      fIxEta_2[fNoOfTracks2]        = ixEta;
      fIxPhi_2[fNoOfTracks2]        = ixPhi;
      fIxPt_2[fNoOfTracks2]         = ixPt;
      fFlags_2[fNoOfTracks2]        = 0x0;
      fPt_2[fNoOfTracks2]           = pT;
      fEta_2[fNoOfTracks2]          = eta;
      fPhi_2[fNoOfTracks2]          = ophi;
      fCorrection_2[fNoOfTracks2]   = corr;
      fN1_2                         += corr;
      fSum1Pt_2                     += corrPt;
      fNnw1_2                       += 1;
      fhN1_2_vsEtaPhi->Fill(eta,phi,corr);
      fhSum1Pt_2_vsEtaPhi->Fill(eta,phi,corrPt);
      fSum1Ptnw_2                   += pT;
      fNoOfTracks2++;
      if (!(fNoOfTracks2 < fArraySize)) {
        AliError(Form("Storage for track two full: %d", fArraySize));
        return kFALSE;
      }
    }
  }
  return kTRUE;
}

/// Process the event data when all tracks have been collected
///
/// Depending on the required track pair combination different processes
/// are invoked. Finally the proper profiles are filled.

void Ali2PCorrelations::ProcessEventData() {

  if (fSinglesOnly) {
    /* everything already done */
  }
  else {
    ProcessLikeSignPairs(1);
    ProcessLikeSignPairs(2);
    ProcessUnlikeSignPairs();

    /* finally fill the individual tracks profiles */
    fhN1_1_vsC->Fill(fCentrality, fN1_1);
    fhSum1Pt_1_vsC->Fill(fCentrality, fSum1Pt_1);
    fhN1nw_1_vsC->Fill(fCentrality, fNnw1_1);
    fhSum1Ptnw_1_vsC->Fill(fCentrality, fSum1Ptnw_1);
    fhN1_2_vsC->Fill(fCentrality, fN1_2);
    fhSum1Pt_2_vsC->Fill(fCentrality, fSum1Pt_2);
    fhN1nw_2_vsC->Fill(fCentrality, fNnw1_2);
    fhSum1Ptnw_2_vsC->Fill(fCentrality, fSum1Ptnw_2);
  }
}

/// \brief Process track combinations with the same charge
/// \param bank the tracks bank to use
void Ali2PCorrelations::ProcessLikeSignPairs(Int_t bank) {
  /* reset pair counters */
  double fN2_12 = 0;
  double fSum2PtPt_12 = 0;
  double fSum2DptDpt_12 = 0;
  double fNnw2_12 = 0;
  double fSum2PtPtnw_12 = 0;
  double fSum2DptDptnw_12 = 0;

  /* pair with same charge. The track list should be identical */
  if (bank == 1) {
    /* let's select the pair efficiency correction histogram */
    /* for the time being bank one is plus charge            */
    const THn *effcorr = fPairsEfficiency_PP;
    Int_t bins[4];

    /* we use only the bank one of tracks */
    for (Int_t ix1 = 0; ix1 < fNoOfTracks1; ix1++)
    {
      int ixEta_1    = fIxEta_1[ix1];
      int ixPhi_1    = fIxPhi_1[ix1];
      int ixPt_1     = fIxPt_1[ix1];
      float corr_1   = fCorrection_1[ix1];
      float pt_1     = fPt_1[ix1];
      float avgpt_1 = 0; /* TODO: load avg pt for eta1, phi1 */

      for (Int_t ix2 = ix1+1; ix2 < fNoOfTracks1; ix2++) {
        /* excluded self correlations */
        float corr       = corr_1 * fCorrection_1[ix2];
        float avgpt_2 = 0; /* TODO: load avg pt for eta2, phi2 */
        int ixDeltaEta_d   = ixEta_1-fIxEta_1[ix2]+fNBins_eta_1-1;
        int ixDeltaPhi_d   = ixPhi_1-fIxPhi_1[ix2]; if (ixDeltaPhi_d < 0) ixDeltaPhi_d += fNBins_phi_1;
        int ixDeltaEta_c   = fIxEta_1[ix2]-ixEta_1+fNBins_eta_1-1;
        int ixDeltaPhi_c   = fIxPhi_1[ix2]-ixPhi_1; if (ixDeltaPhi_c < 0) ixDeltaPhi_c += fNBins_phi_1;

        /* apply the pair correction if applicable */
        if (effcorr != NULL) {
          bins[0] = ixDeltaEta_d+1;
          bins[1] = ixDeltaPhi_d+1;
          bins[2] = ixPt_1+1;
          bins[3] = fIxPt_1[ix2]+1;
          corr = corr / effcorr->GetBinContent(bins);
        }

        int globalbin_d = fhN2_12_vsDEtaDPhi[kOO]->GetBin(ixDeltaEta_d+1,ixDeltaPhi_d+1);
        int globalbin_c = fhN2_12_vsDEtaDPhi[kOO]->GetBin(ixDeltaEta_c+1,ixDeltaPhi_c+1);

        fNnw2_12               += 2;
        fN2_12                 += 2*corr;
        fhN2_12_vsDEtaDPhi[kOO]->AddBinContent(globalbin_d,corr);
        fhN2_12_vsDEtaDPhi[kOO]->AddBinContent(globalbin_c,corr);
        float ptpt              = pt_1*fPt_1[ix2];
        fSum2PtPtnw_12         += 2*ptpt;
        fSum2PtPt_12           += 2*corr*ptpt;
        fhSum2PtPt_12_vsDEtaDPhi[kOO]->AddBinContent(globalbin_d,corr*ptpt);
        fhSum2PtPt_12_vsDEtaDPhi[kOO]->AddBinContent(globalbin_c,corr*ptpt);
        float dptdpt            = (pt_1-avgpt_1)*(fPt_1[ix2]-avgpt_2);
        fSum2DptDptnw_12       += 2*dptdpt;
        fSum2DptDpt_12         += 2*corr*dptdpt;
        fhSum2DptDpt_12_vsDEtaDPhi[kOO]->AddBinContent(globalbin_d,corr*dptdpt);
        fhSum2DptDpt_12_vsDEtaDPhi[kOO]->AddBinContent(globalbin_c,corr*dptdpt);
        fhN2_12_vsPtPt[kOO]->Fill(pt_1,fPt_1[ix2],corr);
        fhN2_12_vsPtPt[kOO]->Fill(fPt_1[ix2],pt_1,corr);
      } //ix2
    } //ix1
    fhN2_12_vsC[kOO]->Fill(fCentrality, fN2_12);
    fhSum2PtPt_12_vsC[kOO]->Fill(fCentrality, fSum2PtPt_12);
    fhSum2DptDpt_12_vsC[kOO]->Fill(fCentrality, fSum2DptDpt_12);
    fhN2nw_12_vsC[kOO]->Fill(fCentrality, fNnw2_12);
    fhSum2PtPtnw_12_vsC[kOO]->Fill(fCentrality, fSum2PtPtnw_12);
    fhSum2DptDptnw_12_vsC[kOO]->Fill(fCentrality, fSum2DptDptnw_12);

    /* update the number of entries in the differential histograms filled with AddBinContent */
    fhN2_12_vsDEtaDPhi[kOO]->SetEntries(fhN2_12_vsDEtaDPhi[kOO]->GetEntries() + fNnw2_12);
    fhSum2PtPt_12_vsDEtaDPhi[kOO]->SetEntries(fhSum2PtPt_12_vsDEtaDPhi[kOO]->GetEntries() + fNnw2_12);
    fhSum2DptDpt_12_vsDEtaDPhi[kOO]->SetEntries(fhSum2DptDpt_12_vsDEtaDPhi[kOO]->GetEntries() + fNnw2_12);
  }
  else if (bank == 2) {
    /* let's select the pair efficiency correction histogram */
    /* for the time being bank two is minus charge            */
    const THn *effcorr = fPairsEfficiency_MM;
    Int_t bins[4];

    /* we use only the bank two of tracks */
    for (Int_t ix1 = 0; ix1 < fNoOfTracks2; ix1++)
    {
      int ixEta_1    = fIxEta_2[ix1];
      int ixPhi_1    = fIxPhi_2[ix1];
      int ixPt_1     = fIxPt_2[ix1];
      float corr_1   = fCorrection_2[ix1];
      float pt_1     = fPt_2[ix1];
      float avgpt_1 = 0; /* TODO: load avg pt for eta1, phi1 */

      for (Int_t ix2 = ix1+1; ix2 < fNoOfTracks2; ix2++) {
        /* excluded self correlations */
        float corr      = corr_1 * fCorrection_2[ix2];
        float avgpt_2 = 0; /* TODO: load avg pt for eta2, phi2 */
        int ixDeltaEta_d   = ixEta_1-fIxEta_2[ix2]+fNBins_eta_1-1;
        int ixDeltaPhi_d   = ixPhi_1-fIxPhi_2[ix2]; if (ixDeltaPhi_d < 0) ixDeltaPhi_d += fNBins_phi_1;
        int ixDeltaEta_c   = fIxEta_2[ix2]-ixEta_1+fNBins_eta_1-1;
        int ixDeltaPhi_c   = fIxPhi_2[ix2]-ixPhi_1; if (ixDeltaPhi_c < 0) ixDeltaPhi_c += fNBins_phi_1;

        /* apply the pair correction if applicable */
        if (effcorr != NULL) {
          bins[0] = ixDeltaEta_d+1;
          bins[1] = ixDeltaPhi_d+1;
          bins[2] = ixPt_1+1;
          bins[3] = fIxPt_2[ix2]+1;
          corr = corr / effcorr->GetBinContent(bins);
        }

        int globalbin_d = fhN2_12_vsDEtaDPhi[kTT]->GetBin(ixDeltaEta_d+1,ixDeltaPhi_d+1);
        int globalbin_c = fhN2_12_vsDEtaDPhi[kTT]->GetBin(ixDeltaEta_c+1,ixDeltaPhi_c+1);

        fNnw2_12               += 2;
        fN2_12                 += 2*corr;
        fhN2_12_vsDEtaDPhi[kTT]->AddBinContent(globalbin_d,corr);
        fhN2_12_vsDEtaDPhi[kTT]->AddBinContent(globalbin_c,corr);
        float ptpt              = pt_1*fPt_2[ix2];
        fSum2PtPtnw_12         += 2*ptpt;
        fSum2PtPt_12           += 2*corr*ptpt;
        fhSum2PtPt_12_vsDEtaDPhi[kTT]->AddBinContent(globalbin_d,corr*ptpt);
        fhSum2PtPt_12_vsDEtaDPhi[kTT]->AddBinContent(globalbin_c,corr*ptpt);
        float dptdpt            = (pt_1-avgpt_1)*(fPt_2[ix2]-avgpt_2);
        fSum2DptDptnw_12       += 2*dptdpt;
        fSum2DptDpt_12         += 2*corr*dptdpt;
        fhSum2DptDpt_12_vsDEtaDPhi[kTT]->AddBinContent(globalbin_d,corr*dptdpt);
        fhSum2DptDpt_12_vsDEtaDPhi[kTT]->AddBinContent(globalbin_c,corr*dptdpt);
        fhN2_12_vsPtPt[kTT]->Fill(pt_1,fPt_2[ix2],corr);
        fhN2_12_vsPtPt[kTT]->Fill(fPt_2[ix2],pt_1,corr);
      } //ix2
    } //ix1
    fhN2_12_vsC[kTT]->Fill(fCentrality, fN2_12);
    fhSum2PtPt_12_vsC[kTT]->Fill(fCentrality, fSum2PtPt_12);
    fhSum2DptDpt_12_vsC[kTT]->Fill(fCentrality, fSum2DptDpt_12);
    fhN2nw_12_vsC[kTT]->Fill(fCentrality, fNnw2_12);
    fhSum2PtPtnw_12_vsC[kTT]->Fill(fCentrality, fSum2PtPtnw_12);
    fhSum2DptDptnw_12_vsC[kTT]->Fill(fCentrality, fSum2DptDptnw_12);

    /* update the number of entries in the differential histograms filled with AddBinContent */
    fhN2_12_vsDEtaDPhi[kTT]->SetEntries(fhN2_12_vsDEtaDPhi[kTT]->GetEntries() + fNnw2_12);
    fhSum2PtPt_12_vsDEtaDPhi[kTT]->SetEntries(fhSum2PtPt_12_vsDEtaDPhi[kTT]->GetEntries() + fNnw2_12);
    fhSum2DptDpt_12_vsDEtaDPhi[kTT]->SetEntries(fhSum2DptDpt_12_vsDEtaDPhi[kTT]->GetEntries() + fNnw2_12);
  }
}

/// \brief Flag the potential products of conversions and / or resonances
void Ali2PCorrelations::FlagConversionsAndResonances() {
  AliInfo("");
  /* pair with different charges. Both track list are different */
  for (Int_t ix1 = 0; ix1 < fNoOfTracks1; ix1++) {
    for (Int_t ix2 = 0; ix2 < fNoOfTracks2; ix2++) {

      for (Int_t ires = 0; ires < fgkNoOfResonances; ires++) {
        if (fThresholdMult[ires] != 0) {
          Float_t mass = checkIfResonance(ires, kTRUE, fPt_1[ix1], fEta_1[ix1], fPhi_1[ix1], fPt_2[ix2], fEta_2[ix2], fPhi_2[ix2]);

          if (0 < mass) {
            fFlags_1[ix1] |= (0x1 << ires);
            fFlags_2[ix2] |= (0x1 << ires);
            fhResonanceMasses->Fill(ires,TMath::Sqrt(mass));
          }
        }
      }
    } //ix2
  } //ix1
}

/// \brief Process track combinations with oposite charge
void Ali2PCorrelations::ProcessUnlikeSignPairs() {
  AliInfo("");
  /* flag the resonances / conversion candidates */
  FlagConversionsAndResonances();

  /* reset pair counters */
  double fN2_12 = 0;
  double fSum2PtPt_12 = 0;
  double fSum2DptDpt_12 = 0;
  double fNnw2_12 = 0;
  double fSum2PtPtnw_12 = 0;
  double fSum2DptDptnw_12 = 0;

  /* pair with different charges. Both track list are different */
  for (Int_t ix1 = 0; ix1 < fNoOfTracks1; ix1++)
  {
    /* let's select the pair efficiency correction histogram */
    /* we assume in principle that PM and MP corrections are */
    /* identical when the delta eta and delta phi bins are   */
    /* properly taken                                        */
    const THn *effcorr = fPairsEfficiency_PM;
    Int_t bins[4];

    int ixEta_1    = fIxEta_1[ix1];
    int ixPhi_1    = fIxPhi_1[ix1];
    int ixPt_1     = fIxPt_1[ix1];
    float corr_1   = fCorrection_1[ix1];
    float pt_1     = fPt_1[ix1];
    float avgpt_1 = 0; /* TODO: load avg pt for eta1, phi1 */

    for (Int_t ix2 = 0; ix2 < fNoOfTracks2; ix2++) {
      /* process the resonance suppression for this pair if needed */
      Bool_t processpair = kTRUE;
      for (Int_t ires = 0; ires < fgkNoOfResonances; ires++) {
        if (fThresholdMult[ires] != 0) {
          /* check if both tracks are flagged for the current resonance */
          if (((fFlags_1[ix1] & fFlags_2[ix2]) & UInt_t(0x1 << ires)) != UInt_t(0x1 << ires)) {
            /* no, they are not, continue */
            continue;
          }
          else {
            /* yes, check if applicable */
            Float_t mass = checkIfResonance(ires, kFALSE, fPt_1[ix1], fEta_1[ix1], fPhi_1[ix1], fPt_2[ix2], fEta_2[ix2], fPhi_2[ix2]);

            if (0 < mass) {
              fhDiscardedResonanceMasses->Fill(ires,TMath::Sqrt(mass));
              processpair = kFALSE;
              break;
            }
          }
        }
      }

      if (processpair) {
        Float_t corr      = corr_1 * fCorrection_2[ix2];
        float avgpt_2 = 0; /* TODO: load avg pt for eta2, phi2 */
        int ixDeltaEta_d   = ixEta_1-fIxEta_2[ix2]+fNBins_eta_1-1;
        int ixDeltaPhi_d   = ixPhi_1-fIxPhi_2[ix2]; if (ixDeltaPhi_d < 0) ixDeltaPhi_d += fNBins_phi_1;
        int ixDeltaEta_c   = fIxEta_2[ix2]-ixEta_1+fNBins_eta_1-1;
        int ixDeltaPhi_c   = fIxPhi_2[ix2]-ixPhi_1; if (ixDeltaPhi_c < 0) ixDeltaPhi_c += fNBins_phi_1;

        /* apply the pair correction if applicable */
        if (effcorr != NULL) {
          bins[0] = ixDeltaEta_d+1;
          bins[1] = ixDeltaPhi_d+1;
          bins[2] = ixPt_1+1;
          bins[3] = fIxPt_2[ix2]+1;
          corr = corr / effcorr->GetBinContent(bins);
        }

        int globalbin_d = fhN2_12_vsDEtaDPhi[kOT]->GetBin(ixDeltaEta_d+1,ixDeltaPhi_d+1);
        int globalbin_c = fhN2_12_vsDEtaDPhi[kTO]->GetBin(ixDeltaEta_c+1,ixDeltaPhi_c+1);

        fNnw2_12               += 1;
        fN2_12                 += corr;
        fhN2_12_vsDEtaDPhi[kOT]->AddBinContent(globalbin_d,corr);
        fhN2_12_vsDEtaDPhi[kTO]->AddBinContent(globalbin_c,corr);
        float ptpt              = pt_1*fPt_2[ix2];
        fSum2PtPtnw_12         += ptpt;
        fSum2PtPt_12           += corr*ptpt;
        fhSum2PtPt_12_vsDEtaDPhi[kOT]->AddBinContent(globalbin_d,corr*ptpt);
        fhSum2PtPt_12_vsDEtaDPhi[kTO]->AddBinContent(globalbin_c,corr*ptpt);
        float dptdpt            = (pt_1-avgpt_1)*(fPt_2[ix2]-avgpt_2);
        fSum2DptDptnw_12       += dptdpt;
        fSum2DptDpt_12         += corr*dptdpt;
        fhSum2DptDpt_12_vsDEtaDPhi[kOT]->AddBinContent(globalbin_d,corr*dptdpt);
        fhSum2DptDpt_12_vsDEtaDPhi[kTO]->AddBinContent(globalbin_c,corr*dptdpt);
        fhN2_12_vsPtPt[kOT]->Fill(pt_1,fPt_2[ix2],corr);
        fhN2_12_vsPtPt[kTO]->Fill(fPt_2[ix2],pt_1,corr);
      }
    } //ix2
  } //ix1
  fhN2_12_vsC[kOT]->Fill(fCentrality, fN2_12);
  fhN2_12_vsC[kTO]->Fill(fCentrality, fN2_12);
  fhSum2PtPt_12_vsC[kOT]->Fill(fCentrality, fSum2PtPt_12);
  fhSum2PtPt_12_vsC[kTO]->Fill(fCentrality, fSum2PtPt_12);
  fhSum2DptDpt_12_vsC[kOT]->Fill(fCentrality, fSum2DptDpt_12);
  fhSum2DptDpt_12_vsC[kTO]->Fill(fCentrality, fSum2DptDpt_12);
  fhN2nw_12_vsC[kOT]->Fill(fCentrality, fNnw2_12);
  fhN2nw_12_vsC[kTO]->Fill(fCentrality, fNnw2_12);
  fhSum2PtPtnw_12_vsC[kOT]->Fill(fCentrality, fSum2PtPtnw_12);
  fhSum2PtPtnw_12_vsC[kTO]->Fill(fCentrality, fSum2PtPtnw_12);
  fhSum2DptDptnw_12_vsC[kOT]->Fill(fCentrality, fSum2DptDptnw_12);
  fhSum2DptDptnw_12_vsC[kTO]->Fill(fCentrality, fSum2DptDptnw_12);

  /* update the number of entries in the differential histograms filled with AddBinContent */
  fhN2_12_vsDEtaDPhi[kOT]->SetEntries(fhN2_12_vsDEtaDPhi[kOT]->GetEntries() + fNnw2_12);
  fhN2_12_vsDEtaDPhi[kTO]->SetEntries(fhN2_12_vsDEtaDPhi[kTO]->GetEntries() + fNnw2_12);
  fhSum2PtPt_12_vsDEtaDPhi[kOT]->SetEntries(fhSum2PtPt_12_vsDEtaDPhi[kOT]->GetEntries() + fNnw2_12);
  fhSum2PtPt_12_vsDEtaDPhi[kTO]->SetEntries(fhSum2PtPt_12_vsDEtaDPhi[kTO]->GetEntries() + fNnw2_12);
  fhSum2DptDpt_12_vsDEtaDPhi[kOT]->SetEntries(fhSum2DptDpt_12_vsDEtaDPhi[kOT]->GetEntries() + fNnw2_12);
  fhSum2DptDpt_12_vsDEtaDPhi[kTO]->SetEntries(fhSum2DptDpt_12_vsDEtaDPhi[kTO]->GetEntries() + fNnw2_12);
}

/// \brief Fill the histograms once the whole process is finished
void  Ali2PCorrelations::FinalizeProcess()
{

  if (fSinglesOnly) {
    /* for the time being, the errors are trivial so, we do not use results file space */
    fhN1_1_vsPt->Sumw2(false);
    fhN1_1_vsZEtaPhiPt->Sumw2(false);
    fhN1_2_vsPt->Sumw2(false);
    fhN1_2_vsZEtaPhiPt->Sumw2(false);
  }
  else {
    /* for the time being, the errors are trivial or not valid so, we do not use results file space */
    fhN1_1_vsEtaPhi->Sumw2(false);
    fhSum1Pt_1_vsEtaPhi->Sumw2(false);
    fhN1_2_vsEtaPhi->Sumw2(false);
    fhSum1Pt_2_vsEtaPhi->Sumw2(false);
    for (int i = 0; i<nTrackPairs; ++i) {
      fhN2_12_vsDEtaDPhi[i]->Sumw2(false);
      fhSum2PtPt_12_vsDEtaDPhi[i]->Sumw2(false);
      fhSum2DptDpt_12_vsDEtaDPhi[i]->Sumw2(false);
      fhN2_12_vsPtPt[i]->Sumw2(false);
    }
  }
}


/// \cond CLASSIMP
ClassImp(Ali2PCorrelations);
/// \endcond
