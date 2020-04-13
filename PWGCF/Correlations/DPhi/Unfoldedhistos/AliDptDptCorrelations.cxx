/**************************************************************************
* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved.  *
*                                                                         *
* Permission to use, copy, modify and distribute this software and its    *
* documentation strictly for non-commercial purposes is hereby granted    *
* without fee, provided that the above copyright notice appears in all    *
* copies and that both the copyright notice and this permission notice    *
* appear in the supporting documentation. The authors make no claims      *
* about the suitability of this software for any purpose. It is           *
* provided "as is" without express or implied warranty.                   *
**************************************************************************/

/// \file AliDptDptCorrelations.cxx
/// \brief Implementation of the AliDptDptCorrelations class

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
#include "AliDptDptCorrelations.h"

Int_t AliDptDptCorrelations::fgkNoOfResonances = 4; ///< four resonances / conversions for the time being
Double_t AliDptDptCorrelations::fgkMass[16] = {/* photon */ 0.0, /* k0 */ 0.4976, /* lambda */ 1.115, /* rho */ 0.775};
Double_t AliDptDptCorrelations::fgkChildMass[2][16] = {
    {0.510e-3, 0.1396, 0.1396, 0.1396},
    {0.510e-3, 0.1396, 0.9383, 0.1396}
};
Double_t AliDptDptCorrelations::fgkMassThreshold[16] = {0.04,0.01,0.05,0.04};

/// Default constructor for object serialization
AliDptDptCorrelations::AliDptDptCorrelations() :
    TNamed(),
    fOutput(NULL),
    fVertexZ(0.0),
    fIxVertexZ(0),
    fCentrality(0.0),
    fHalfSymmetrize(kTRUE),
    fSinglesOnly(kTRUE),
    fUseWeights(kFALSE),
    fUseSimulation(kFALSE),
    fSameSign(kFALSE),
    fAllCombinations(kFALSE),
    fRequestedCharge_1(1),
    fRequestedCharge_2(-1),
    /* the arrays with tracks 1 and 2 information */
    fArraySize(5*1024),
    fNoOfTracks1(0),
    fId_1(NULL),
    fCharge_1(NULL),
    fIxEtaPhi_1(NULL),
    fIxPt_1(NULL),
    fFlags_1(NULL),
    fPt_1(NULL),
    fEta_1(NULL),
    fPhi_1(NULL),
    fCorrection_1(NULL),
    fNoOfTracks2(0),
    fId_2(NULL),
    fCharge_2(NULL),
    fIxEtaPhi_2(NULL),
    fIxPt_2(NULL),
    fFlags_2(NULL),
    fPt_2(NULL),
    fEta_2(NULL),
    fPhi_2(NULL),
    fCorrection_2(NULL),
    /* correction weights */
    fCorrectionWeights_1(NULL),
    fCorrectionWeights_2(NULL),
    /* efficiency correction */
    fEfficiencyCorrection_1(NULL),
    fEfficiencyCorrection_2(NULL),
    fPairsEfficiency_PP(NULL),
    fPairsEfficiency_PM(NULL),
    fPairsEfficiency_MM(NULL),
    fPairsEfficiency_MP(NULL),
    /* simulation pdfs */
    fPositiveTrackPdfs(NULL),
    fNegativeTrackPdfs(NULL),
    fPositiveTrackCurrentPdf(NULL),
    fNegativeTrackCurrentPdf(NULL),
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
    fNBins_etaPhi_12(2073600),
    /* accumulated values */
    fN1_1(0.0),
    fN1_2(0.0),
    fN2_12(0.0),
    fNnw1_1(0.0),
    fNnw1_2(0.0),
    fNnw2_12(0.0),
    fSum1Pt_1(0.0),
    fSum1Pt_2(0.0),
    fSum2PtPt_12(0.0),
    fSum1Ptnw_1(0.0),
    fSum1Ptnw_2(0.0),
    fSum2PtPtnw_12(0.0),
    fSum2NPt_12(0.0),
    fSum2PtN_12(0.0),
    fSum2NPtnw_12(0.0),
    fSum2PtNnw_12(0.0),
    /* storage arrays */
    fN1_1_vsPt(NULL),
    fN1_1_vsEtaPhi(NULL),
    fSum1Pt_1_vsEtaPhi(NULL),
    fN1_1_vsZEtaPhiPt(NULL),
    fN1_2_vsPt(NULL),
    fN1_2_vsEtaPhi(NULL),
    fSum1Pt_2_vsEtaPhi(NULL),
    fN1_2_vsZEtaPhiPt(NULL),
    fN2_12_vsPtPt(NULL),
    fN2_12_vsEtaPhi(NULL),
    fSum2PtPt_12_vsEtaPhi(NULL),
    fSum2PtN_12_vsEtaPhi(NULL),
    fSum2NPt_12_vsEtaPhi(NULL),
    /* histograms */
    fhN1_1_vsPt(NULL),
    fhN1_1_vsEtaPhi(NULL),
    fhSum1Pt_1_vsEtaPhi(NULL),
    fhN1_1_vsZEtaPhiPt(NULL),
    fhN1_2_vsPt(NULL),
    fhN1_2_vsEtaPhi(NULL),
    fhSum1Pt_2_vsEtaPhi(NULL),
    fhN1_2_vsZEtaPhiPt(NULL),
    fhN2_12_vsPtPt(NULL),
    fhN2_12_vsEtaPhi(NULL),
    fhSum2PtPt_12_vsEtaPhi(NULL),
    fhSum2PtN_12_vsEtaPhi(NULL),
    fhSum2NPt_12_vsEtaPhi(NULL),
    /* vs centrality profiles */
    fhN1_1_vsC(NULL),
    fhSum1Pt_1_vsC(NULL),
    fhN1nw_1_vsC(NULL),
    fhSum1Ptnw_1_vsC(NULL),
    fhN1_2_vsC(NULL),
    fhSum1Pt_2_vsC(NULL),
    fhN1nw_2_vsC(NULL),
    fhSum1Ptnw_2_vsC(NULL),
    fhN2_12_vsC(NULL),
    fhSum2PtPt_12_vsC(NULL),
    fhSum2PtN_12_vsC(NULL),
    fhSum2NPt_12_vsC(NULL),
    fhN2nw_12_vsC(NULL),
    fhSum2PtPtnw_12_vsC(NULL),
    fhSum2PtNnw_12_vsC(NULL),
    fhSum2NPtnw_12_vsC(NULL),
    fhResonanceRoughMasses(NULL),
    fhResonanceMasses(NULL),
    fhDiscardedResonanceMasses(NULL)
{
  for (Int_t ires = 0; ires < 16; ires++)
    fThresholdMult[ires] = 0x0;
}

/// Normal constructor
/// \param name the name for the object instance
AliDptDptCorrelations::AliDptDptCorrelations(const char *name) :
    TNamed(name,name),
    fOutput(NULL),
    fVertexZ(0.0),
    fIxVertexZ(0),
    fCentrality(0.0),
    fHalfSymmetrize(kTRUE),
    fSinglesOnly(kTRUE),
    fUseWeights(kFALSE),
    fUseSimulation(kFALSE),
    fSameSign(kFALSE),
    fAllCombinations(kFALSE),
    fRequestedCharge_1(1),
    fRequestedCharge_2(-1),
    /* the arrays with tracks 1 and 2 information */
    fArraySize(5*1024),
    fNoOfTracks1(0),
    fId_1(NULL),
    fCharge_1(NULL),
    fIxEtaPhi_1(NULL),
    fIxPt_1(NULL),
    fFlags_1(NULL),
    fPt_1(NULL),
    fEta_1(NULL),
    fPhi_1(NULL),
    fCorrection_1(NULL),
    fNoOfTracks2(0),
    fId_2(NULL),
    fCharge_2(NULL),
    fIxEtaPhi_2(NULL),
    fIxPt_2(NULL),
    fFlags_2(NULL),
    fPt_2(NULL),
    fEta_2(NULL),
    fPhi_2(NULL),
    fCorrection_2(NULL),
    /* correction weights */
    fCorrectionWeights_1(NULL),
    fCorrectionWeights_2(NULL),
    /* efficiency correction */
    fEfficiencyCorrection_1(NULL),
    fEfficiencyCorrection_2(NULL),
    fPairsEfficiency_PP(NULL),
    fPairsEfficiency_PM(NULL),
    fPairsEfficiency_MM(NULL),
    fPairsEfficiency_MP(NULL),
    /* simulation pdfs */
    fPositiveTrackPdfs(NULL),
    fNegativeTrackPdfs(NULL),
    fPositiveTrackCurrentPdf(NULL),
    fNegativeTrackCurrentPdf(NULL),
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
    fNBins_etaPhi_12(2073600),
    /* accumulated values */
    fN1_1(0.0),
    fN1_2(0.0),
    fN2_12(0.0),
    fNnw1_1(0.0),
    fNnw1_2(0.0),
    fNnw2_12(0.0),
    fSum1Pt_1(0.0),
    fSum1Pt_2(0.0),
    fSum2PtPt_12(0.0),
    fSum1Ptnw_1(0.0),
    fSum1Ptnw_2(0.0),
    fSum2PtPtnw_12(0.0),
    fSum2NPt_12(0.0),
    fSum2PtN_12(0.0),
    fSum2NPtnw_12(0.0),
    fSum2PtNnw_12(0.0),
    /* storage arrays */
    fN1_1_vsPt(NULL),
    fN1_1_vsEtaPhi(NULL),
    fSum1Pt_1_vsEtaPhi(NULL),
    fN1_1_vsZEtaPhiPt(NULL),
    fN1_2_vsPt(NULL),
    fN1_2_vsEtaPhi(NULL),
    fSum1Pt_2_vsEtaPhi(NULL),
    fN1_2_vsZEtaPhiPt(NULL),
    fN2_12_vsPtPt(NULL),
    fN2_12_vsEtaPhi(NULL),
    fSum2PtPt_12_vsEtaPhi(NULL),
    fSum2PtN_12_vsEtaPhi(NULL),
    fSum2NPt_12_vsEtaPhi(NULL),
    /* histograms */
    fhN1_1_vsPt(NULL),
    fhN1_1_vsEtaPhi(NULL),
    fhSum1Pt_1_vsEtaPhi(NULL),
    fhN1_1_vsZEtaPhiPt(NULL),
    fhN1_2_vsPt(NULL),
    fhN1_2_vsEtaPhi(NULL),
    fhSum1Pt_2_vsEtaPhi(NULL),
    fhN1_2_vsZEtaPhiPt(NULL),
    fhN2_12_vsPtPt(NULL),
    fhN2_12_vsEtaPhi(NULL),
    fhSum2PtPt_12_vsEtaPhi(NULL),
    fhSum2PtN_12_vsEtaPhi(NULL),
    fhSum2NPt_12_vsEtaPhi(NULL),
    /* vs centrality profiles */
    fhN1_1_vsC(NULL),
    fhSum1Pt_1_vsC(NULL),
    fhN1nw_1_vsC(NULL),
    fhSum1Ptnw_1_vsC(NULL),
    fhN1_2_vsC(NULL),
    fhSum1Pt_2_vsC(NULL),
    fhN1nw_2_vsC(NULL),
    fhSum1Ptnw_2_vsC(NULL),
    fhN2_12_vsC(NULL),
    fhSum2PtPt_12_vsC(NULL),
    fhSum2PtN_12_vsC(NULL),
    fhSum2NPt_12_vsC(NULL),
    fhN2nw_12_vsC(NULL),
    fhSum2PtPtnw_12_vsC(NULL),
    fhSum2PtNnw_12_vsC(NULL),
    fhSum2NPtnw_12_vsC(NULL),
    fhResonanceRoughMasses(NULL),
    fhResonanceMasses(NULL),
    fhDiscardedResonanceMasses(NULL)
{
  for (Int_t ires = 0; ires < 16; ires++)
    fThresholdMult[ires] = 0x0;
}

/// \brief Default destructor
/// Deallocates the allocated memory
AliDptDptCorrelations::~AliDptDptCorrelations() {
  if (fN1_1_vsPt != NULL) delete [] fN1_1_vsPt;
  if (fN1_1_vsEtaPhi != NULL) delete [] fN1_1_vsEtaPhi;
  if (fSum1Pt_1_vsEtaPhi != NULL) delete [] fSum1Pt_1_vsEtaPhi;
  if (fN1_1_vsZEtaPhiPt != NULL) delete [] fN1_1_vsZEtaPhiPt;
  if (fN1_2_vsPt != NULL) delete [] fN1_2_vsPt;
  if (fN1_2_vsEtaPhi != NULL) delete [] fN1_2_vsEtaPhi;
  if (fSum1Pt_2_vsEtaPhi != NULL) delete [] fSum1Pt_2_vsEtaPhi;
  if (fN1_2_vsZEtaPhiPt != NULL) delete [] fN1_2_vsZEtaPhiPt;
  if (fN2_12_vsPtPt != NULL) delete [] fN2_12_vsPtPt;
  if (fN2_12_vsEtaPhi != NULL) delete [] fN2_12_vsEtaPhi;
  if (fSum2PtPt_12_vsEtaPhi != NULL) delete [] fSum2PtPt_12_vsEtaPhi;
  if (fSum2PtN_12_vsEtaPhi != NULL) delete [] fSum2PtN_12_vsEtaPhi;
  if (fSum2NPt_12_vsEtaPhi != NULL) delete [] fSum2NPt_12_vsEtaPhi;
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
Bool_t AliDptDptCorrelations::ConfigureBinning(const char *confstring) {
#define DPTDPTCORRBINCONFIGPAR 10

  Double_t min_pt = 0.0, max_pt = 0.0, width_pt = 0.0;
  Double_t min_eta = 0.0, max_eta = 0.0, width_eta = 0.0;
  Int_t    nBins_phi = 0;

  /* few sanity checks */
  TString str = confstring;
  if (!str.Contains("halfsymm") || !str.Contains("phishift"))
    return kFALSE;

  TObjArray *array = str.Tokenize(";");
  for (Int_t item = 0; item < array->GetEntries(); item++) {
    const TString &stritem = ((TObjString*) array->At(item))->GetString();
    if (stritem.BeginsWith("halfsymm:")) {
      fHalfSymmetrize = stritem.Contains("yes");
    }
    else if (stritem.BeginsWith("phishift:")) {
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
void AliDptDptCorrelations::ConfigureResonances(const char *confstring) {

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
TString AliDptDptCorrelations::GetBinningConfigurationString() const {
  if (fMin_pt_1 != fMin_pt_2 || fMax_pt_2 != fMax_pt_2 || fWidth_pt_1 != fWidth_pt_2 ||
      fMin_eta_1 != fMin_eta_2 || fMax_eta_1 != fMax_eta_2 || fWidth_eta_1 != fWidth_eta_2 ||
      fNBins_phi_1 != fNBins_phi_2) {
    return TString(Form("WrongAsymmetricBinning:halfsymm:%s;phishift:%.1f;%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%d",
        (fHalfSymmetrize ? "yes" : "not"),
        fNBinsPhiShift,
        fMin_vertexZ, fMax_vertexZ, fWidth_vertexZ,
        fMin_pt_1, fMax_pt_1, fWidth_pt_1,
        fMin_eta_1, fMax_eta_1, fWidth_eta_1, fNBins_phi_1));
  }
  else {
    return TString(Form("Binning:halfsymm:%s;phishift:%.1f;%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%d",
        (fHalfSymmetrize ? "yes" : "not"),
        fNBinsPhiShift,
        fMin_vertexZ, fMax_vertexZ, fWidth_vertexZ,
        fMin_pt_1, fMax_pt_1, fWidth_pt_1,
        fMin_eta_1, fMax_eta_1, fWidth_eta_1, fNBins_phi_1));
  }
}

/// \brief Build the resonances rejection configuration string
/// \return the configuration string corresponding to the current resonance rejection configuration
TString AliDptDptCorrelations::GetResonancesConfigurationString() const {
  TString str = "resonances:";

  for (Int_t ires = 0; ires < fgkNoOfResonances; ires++) {
    str += TString::Format("%01d",fThresholdMult[fgkNoOfResonances - 1 - ires]);
  }
  return str;
}



/// \brief Initializes the member data structures
/// Allocates the needed memory an create the output histograms.
void AliDptDptCorrelations::Initialize()
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
  fNBins_etaPhi_12   = fNBins_etaPhi_1 * fNBins_etaPhi_2;

  /* incorporate configuration parameters */
  fOutput->Add(new TParameter<Bool_t>("HalfSymmetricResults",fHalfSymmetrize,'f'));
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
  fIxEtaPhi_1 = new Int_t[fArraySize];
  fIxPt_1 = new Int_t[fArraySize];
  fFlags_1 = new UInt_t[fArraySize];
  fPt_1 = new Float_t[fArraySize];
  fEta_1 = new Float_t[fArraySize];
  fPhi_1 = new Float_t[fArraySize];
  fCorrection_1 = new Float_t[fArraySize];
  fN1_1_vsPt = new Double_t[fNBins_pt_1]; for (Int_t i = 0; i < fNBins_pt_1; i++) fN1_1_vsPt[i] = 0.0;
  fN1_1_vsEtaPhi = new Double_t[fNBins_etaPhi_1]; for (Int_t i = 0; i < fNBins_etaPhi_1; i++) fN1_1_vsEtaPhi[i] = 0.0;
  fSum1Pt_1_vsEtaPhi = new Double_t[fNBins_etaPhi_1]; for (Int_t i = 0; i < fNBins_etaPhi_1; i++) fSum1Pt_1_vsEtaPhi[i] = 0.0;
  fN1_1_vsZEtaPhiPt = new Float_t[fNBins_zEtaPhiPt_1]; for (Int_t i = 0; i < fNBins_zEtaPhiPt_1; i++) fN1_1_vsZEtaPhiPt[i] = 0.0;

  /* track 2 storage */
  fId_2 = new Int_t[fArraySize];
  fCharge_2 = new Int_t[fArraySize];
  fIxEtaPhi_2 = new Int_t[fArraySize];
  fIxPt_2 = new Int_t[fArraySize];
  fFlags_2 = new UInt_t[fArraySize];
  fPt_2 = new Float_t[fArraySize];
  fEta_2 = new Float_t[fArraySize];
  fPhi_2 = new Float_t[fArraySize];
  fCorrection_2 = new Float_t[fArraySize];
  fN1_2_vsPt = new Double_t[fNBins_pt_2]; for (Int_t i = 0; i < fNBins_pt_2; i++) fN1_2_vsPt[i] = 0.0;
  fN1_2_vsEtaPhi = new Double_t[fNBins_etaPhi_2]; for (Int_t i = 0; i < fNBins_etaPhi_2; i++) fN1_2_vsEtaPhi[i] = 0.0;
  fSum1Pt_2_vsEtaPhi = new Double_t[fNBins_etaPhi_2]; for (Int_t i = 0; i < fNBins_etaPhi_2; i++) fSum1Pt_2_vsEtaPhi[i] = 0.0;
  fN1_2_vsZEtaPhiPt = new Float_t[fNBins_zEtaPhiPt_2]; for (Int_t i = 0; i < fNBins_zEtaPhiPt_2; i++) fN1_2_vsZEtaPhiPt[i] = 0.0;

  fN2_12_vsPtPt = new Double_t[fNBins_pt_1*fNBins_pt_2]; for (Int_t i = 0; i < fNBins_pt_1*fNBins_pt_2; i++) fN2_12_vsPtPt[i] = 0.0;
  fN2_12_vsEtaPhi = new Float_t[fNBins_etaPhi_12]; for (Int_t i = 0; i < fNBins_etaPhi_12; i++) fN2_12_vsEtaPhi[i] = 0.0;
  fSum2PtPt_12_vsEtaPhi = new Float_t[fNBins_etaPhi_12]; for (Int_t i = 0; i < fNBins_etaPhi_12; i++) fSum2PtPt_12_vsEtaPhi[i] = 0.0;
  fSum2PtN_12_vsEtaPhi = new Float_t[fNBins_etaPhi_12]; for (Int_t i = 0; i < fNBins_etaPhi_12; i++) fSum2PtN_12_vsEtaPhi[i] = 0.0;
  fSum2NPt_12_vsEtaPhi = new Float_t[fNBins_etaPhi_12]; for (Int_t i = 0; i < fNBins_etaPhi_12; i++) fSum2NPt_12_vsEtaPhi[i] = 0.0;

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

    fhN1_1_vsEtaPhi = new TH2F("n1_1_vsEtaPhi","#LT n_{1} #GT;#eta_{1};#varphi_{1} (radian);#LT n_{1} #GT",
        fNBins_eta_1, fMin_eta_1, fMax_eta_1,  fNBins_phi_1, fMin_phi_1, fMax_phi_1);
    fhSum1Pt_1_vsEtaPhi = new TH2F("sumPt_1_vsEtaPhi","#LT #Sigma p_{t,1} #GT;#eta_{1};#varphi_{1} (radian);#LT #Sigma p_{t,1} #GT (GeV/c)",
        fNBins_eta_1,fMin_eta_1,fMax_eta_1,fNBins_phi_1,fMin_phi_1,fMax_phi_1);
    fhN1_1_vsC = new TProfile("n1_1_vsM","#LT n_{1} #GT (weighted);Centrality (%);#LT n_{1} #GT",100,0.0,100.0);
    fhSum1Pt_1_vsC = new TProfile("sumPt_1_vsM","#LT #Sigma p_{t,1} #GT (weighted);Centrality (%);#LT #Sigma p_{t,1} #GT (GeV/c)",100,0.0,100.0);
    fhN1nw_1_vsC = new TProfile("n1Nw_1_vsM","#LT n_{1} #GT;Centrality (%);#LT n_{1} #GT",100,0.0,100.0);
    fhSum1Ptnw_1_vsC = new TProfile("sumPtNw_1_vsM","#LT #Sigma p_{t,1} #GT;Centrality (%);#LT #Sigma p_{t,1} #GT (GeV/c)",100,0.0,100.0);

    fhN1_2_vsEtaPhi = new TH2F("n1_2_vsEtaPhi","#LT n_{1} #GT;#eta_{2};#varphi_{2} (radian);#LT n_{1} #GT",
        fNBins_eta_2, fMin_eta_2, fMax_eta_2,  fNBins_phi_2, fMin_phi_2, fMax_phi_2);
    fhSum1Pt_2_vsEtaPhi = new TH2F("sumPt_2_vsEtaPhi","#LT #Sigma p_{t,2} #GT;#eta_{2};#varphi_{2} (radian);#LT #Sigma p_{t,2} #GT (GeV/c)",
        fNBins_eta_2,fMin_eta_2,fMax_eta_2,fNBins_phi_2,fMin_phi_2,fMax_phi_2);
    fhN1_2_vsC = new TProfile("n1_2_vsM","#LT n_{1} #GT (weighted);Centrality (%);#LT n_{1} #GT",100,0.0,100.0);
    fhSum1Pt_2_vsC = new TProfile("sumPt_2_vsM","#LT #Sigma p_{t,1} #GT (weighted);Centrality (%);#LT #Sigma p_{t,1} #GT (GeV/c)",100,0.0,100.0);
    fhN1nw_2_vsC = new TProfile("n1Nw_2_vsM","#LT n_{1} #GT;Centrality (%);#LT n_{1} #GT",100,0.0,100.0);
    fhSum1Ptnw_2_vsC = new TProfile("sumPtNw_2_vsM","#LT #Sigma p_{t,1} #GT;Centrality (%);#LT #Sigma p_{t,1} #GT (GeV/c)",100,0.0,100.0);

    fhN2_12_vsEtaPhi = new TH1F("n2_12_vsEtaPhi","#LT n_{2} #GT;#eta_{1}#times#varphi_{1}#times#eta_{2}#times#varphi_{2};#LT n_{2} #GT",
        fNBins_etaPhi_12,0.,Double_t(fNBins_etaPhi_12));
    fhSum2PtPt_12_vsEtaPhi = new TH1F("sumPtPt_12_vsEtaPhi","#LT #Sigma p_{t,1}p_{t,2} #GT;#eta_{1}#times#varphi_{1}#times#eta_{2}#times#varphi_{2};#LT #Sigma p_{t,1}p_{t,2} #GT (GeV)^{2}",
        fNBins_etaPhi_12,0.,Double_t(fNBins_etaPhi_12));
    fhSum2PtN_12_vsEtaPhi = new TH1F("sumPtN_12_vsEtaPhi","#LT #Sigma p_{t,1}N #GT;#eta_{1}#times#varphi_{1}#times#eta_{2}#times#varphi_{2};#LT #Sigma p_{t,1}N #GT (GeV)",
        fNBins_etaPhi_12,0.,Double_t(fNBins_etaPhi_12));
    fhSum2NPt_12_vsEtaPhi = new TH1F("sumNPt_12_vsEtaPhi","#LT N#Sigma p_{t,2} #GT;#eta_{1}#times#varphi_{1}#times#eta_{2}#times#varphi_{2};#LT N#Sigma p_{t,2} #GT (GeV)",
        fNBins_etaPhi_12,0.,Double_t(fNBins_etaPhi_12));
    fhN2_12_vsPtPt = new TH2F("n2_12_vsPtVsPt","#LT n_{2} #GT;p_{t,1} (GeV/c);p_{t,2} (GeV/c);#LT n_{2} #GT",
        fNBins_pt_1,fMin_pt_1,fMax_pt_1,fNBins_pt_2,fMin_pt_2,fMax_pt_2);
    fhN2_12_vsC = new TProfile("n2_12_vsM","#LT n_{2} #GT (weighted);Centrality (%);#LT n_{2} #GT",100,0.0,100.0);
    fhSum2PtPt_12_vsC = new TProfile("sumPtPt_12_vsM","#LT #Sigma p_{t,1}p_{t,2} #GT (weighted);Centrality (%);#LT #Sigma p_{t,1}p_{t,2} #GT (GeV)^{2}",100,0.0,100.0);
    fhSum2PtN_12_vsC = new TProfile("sumPtN_12_vsM","#LT #Sigma p_{t,1}N #GT (weighted);Centrality (%);#LT #Sigma p_{t,1}N #GT (GeV)",100,0.0,100.0);
    fhSum2NPt_12_vsC = new TProfile("sumNPt_12_vsM","#LT N#Sigma p_{t,2} #GT (weighted);Centrality (%);#LT N#Sigma p_{t,2} #GT (GeV)",100,0.0,100.0);
    fhN2nw_12_vsC = new TProfile("n2Nw_12_vsM","#LT n_{2} #GT;Centrality (%);#LT n_{2} #GT",100,0.0,100.0);
    fhSum2PtPtnw_12_vsC = new TProfile("sumPtPtNw_12_vsM","#LT #Sigma p_{t,1}p_{t,2} #GT;Centrality (%);#LT #Sigma p_{t,1}p_{t,2} #GT (GeV)^{2}",100,0.0,100.0);
    fhSum2PtNnw_12_vsC = new TProfile("sumPtNNw_12_vsM","#LT #Sigma p_{t,1}N #GT;Centrality (%);#LT #Sigma p_{t,1}N #GT (GeV)",100,0.0,100.0);
    fhSum2NPtnw_12_vsC = new TProfile("sumNPtNw_12_vsM","#LT N#Sigma p_{t,2} #GT;Centrality (%);#LT N#Sigma p_{t,2} #GT (GeV)",100,0.0,100.0);

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

    fOutput->Add(fhN2_12_vsEtaPhi);
    fOutput->Add(fhSum2PtPt_12_vsEtaPhi);
    fOutput->Add(fhSum2PtN_12_vsEtaPhi);
    fOutput->Add(fhSum2NPt_12_vsEtaPhi);
    fOutput->Add(fhN2_12_vsPtPt);
    fOutput->Add(fhN2_12_vsC);
    fOutput->Add(fhSum2PtPt_12_vsC);
    fOutput->Add(fhSum2PtN_12_vsC);
    fOutput->Add(fhSum2NPt_12_vsC);
    fOutput->Add(fhN2nw_12_vsC);
    fOutput->Add(fhSum2PtPtnw_12_vsC);
    fOutput->Add(fhSum2PtNnw_12_vsC);
    fOutput->Add(fhSum2NPtnw_12_vsC);
  }

  TH1::AddDirectory(oldstatus);
}

/// \brief Stores the correction weights for both set of tracks
/// \param h3_1 the calibration weights for the first track
/// \param h3_2 the calibration weights for the second track
/// \return kTRUE if correctly done kFALSE otherwise
Bool_t AliDptDptCorrelations::SetWeigths(const TH3F *h3_1, const TH3F *h3_2) {
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
Bool_t AliDptDptCorrelations::SetEfficiencyCorrection(const TH1F *h1_1, const TH1F *h1_2) {
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
Bool_t AliDptDptCorrelations::SetPairEfficiencyCorrection(const THn *h11, const THn *h12, const THn *h22, const THn *h21) {
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
Bool_t AliDptDptCorrelations::SetSimultationPdfs(const TObjArray *pluspdf, const TObjArray *minuspdf) {

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
Bool_t AliDptDptCorrelations::StartEvent(Float_t centrality, Float_t vertexz) {

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
Bool_t AliDptDptCorrelations::ProcessTrack(Int_t trkId, AliVTrack *trk) {

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
Bool_t AliDptDptCorrelations::ProcessTrack(Int_t trkId, AliVParticle *part) {

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
Bool_t AliDptDptCorrelations::ProcessTrack(Int_t trkId, Int_t charge, Float_t pT, Float_t eta, Float_t ophi) {

  if(charge == 0) return kFALSE;

  /* track 1 */
  if (fRequestedCharge_1 == charge && pT > fMin_pt_1 && pT < fMax_pt_1)
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
      fN1_1_vsPt[ixPt] += corr;
      fN1_1_vsZEtaPhiPt[ixZEtaPhiPt] += corr;
    }
    else {
      Float_t corrPt                = corr*pT;
      fId_1[fNoOfTracks1]           = trkId;
      fCharge_1[fNoOfTracks1]       = charge;
      fIxEtaPhi_1[fNoOfTracks1]     = ixEtaPhi;
      fIxPt_1[fNoOfTracks1]         = ixPt;
      fFlags_1[fNoOfTracks1]        = 0x0;
      fPt_1[fNoOfTracks1]           = pT;
      fEta_1[fNoOfTracks1]          = eta;
      fPhi_1[fNoOfTracks1]          = ophi;
      fCorrection_1[fNoOfTracks1]   = corr;
      fN1_1                         += corr;
      fN1_1_vsEtaPhi[ixEtaPhi]      += corr;
      fSum1Pt_1                     += corrPt;
      fSum1Pt_1_vsEtaPhi[ixEtaPhi]  += corrPt;
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
  if (fRequestedCharge_2 == charge && pT > fMin_pt_2 && pT < fMax_pt_2)
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
      fN1_2_vsPt[ixPt] += corr;
      fN1_2_vsZEtaPhiPt[ixZEtaPhiPt] += corr;
    }
    else {
      Float_t corrPt                = corr*pT;
      fId_2[fNoOfTracks2]           = trkId;
      fCharge_2[fNoOfTracks2]       = charge;
      fIxEtaPhi_2[fNoOfTracks2]     = ixEtaPhi;
      fIxPt_2[fNoOfTracks2]         = ixPt;
      fFlags_2[fNoOfTracks2]        = 0x0;
      fPt_2[fNoOfTracks2]           = pT;
      fEta_2[fNoOfTracks2]          = eta;
      fPhi_2[fNoOfTracks2]          = ophi;
      fCorrection_2[fNoOfTracks2]   = corr;
      fN1_2                         += corr;
      fSum1Pt_2                     += corrPt;
      fNnw1_2                       += 1;
      fN1_2_vsEtaPhi[ixEtaPhi]      += corr;
      fSum1Pt_2_vsEtaPhi[ixEtaPhi]  += corrPt;
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

void AliDptDptCorrelations::ProcessEventData() {

  if (fSinglesOnly) {
    /* everything already done */
  }
  else {
    /* reset pair counters */
    fN2_12   = fSum2PtPt_12   = fSum2NPt_12    = fSum2PtN_12    = 0;
    fNnw2_12 = fSum2PtPtnw_12 = fSum2NPtnw_12  = fSum2PtNnw_12  = 0;

    if (fSameSign) {
      if (fHalfSymmetrize) {
        ProcessLikeSignPairs(1);
      }
      else {
        ProcessNotHalfSymmLikeSignPairs(1);
      }
    }
    else {
      if (fAllCombinations) {
        if (fHalfSymmetrize) {
          ProcessLikeSignPairs(1);
          ProcessLikeSignPairs(2);
          ProcessUnlikeSignPairs();
        }
        else {
          ProcessNotHalfSymmLikeSignPairs(1);
          ProcessNotHalfSymmLikeSignPairs(2);
          ProcessNotHalfSymmUnlikeSignPairs();
        }
      }
      else {
        if (fHalfSymmetrize) {
          ProcessUnlikeSignPairs();
        }
        else {
          ProcessNotHalfSymmUnlikeSignPairs();
        }
      }
    }
    /* now fill the profiles */
    fhN1_1_vsC->Fill(fCentrality, fN1_1);
    fhSum1Pt_1_vsC->Fill(fCentrality, fSum1Pt_1);
    fhN1nw_1_vsC->Fill(fCentrality, fNnw1_1);
    fhSum1Ptnw_1_vsC->Fill(fCentrality, fSum1Ptnw_1);
    fhN1_2_vsC->Fill(fCentrality, fN1_2);
    fhSum1Pt_2_vsC->Fill(fCentrality, fSum1Pt_2);
    fhN1nw_2_vsC->Fill(fCentrality, fNnw1_2);
    fhSum1Ptnw_2_vsC->Fill(fCentrality, fSum1Ptnw_2);
    fhN2_12_vsC->Fill(fCentrality, fN2_12);
    fhSum2PtPt_12_vsC->Fill(fCentrality, fSum2PtPt_12);
    fhSum2PtN_12_vsC->Fill(fCentrality, fSum2PtN_12);
    fhSum2NPt_12_vsC->Fill(fCentrality, fSum2NPt_12);
    fhN2nw_12_vsC->Fill(fCentrality, fNnw2_12);
    fhSum2PtPtnw_12_vsC->Fill(fCentrality, fSum2PtPtnw_12);
    fhSum2PtNnw_12_vsC->Fill(fCentrality, fSum2PtNnw_12);
    fhSum2NPtnw_12_vsC->Fill(fCentrality, fSum2NPtnw_12);
  }
}

/// Helper function to compute an index for a two dimensional
/// symmetric matrix linearized in a single array to reduce
/// memory usage
/// \param ix the index on the x axis (zero based)
/// \param nxbins the number of bins on the x axis
/// \param iy the index on the y axis (zero based)
/// \param nybins the number of bins on the y axis
/// \return the linearized index

inline Int_t getLinearSymmIndex(Int_t ix, Int_t nxbins, Int_t iy, Int_t nybins) {

  if (iy < ix)
      return iy*nxbins - Int_t((iy-1)*iy / 2) + (ix-iy);
  else
      return ix*nybins - Int_t((ix-1)*ix / 2) + (iy-ix);
}

/// \brief Process track combinations with the same charge
/// \param bank the tracks bank to use
void AliDptDptCorrelations::ProcessLikeSignPairs(Int_t bank) {
  /* pair with same charge. The track list should be identical */
  if (bank == 1) {
    /* let's select the pair efficiency correction histogram */
    const THn *effcorr = NULL;
    if (fRequestedCharge_1 > 0)
      effcorr = fPairsEfficiency_PP;
    else
      effcorr = fPairsEfficiency_MM;
    Int_t bins[4];

    /* we use only the bank one of tracks */
    for (Int_t ix1 = 0; ix1 < fNoOfTracks1; ix1++)
    {
      Int_t ixEtaPhi_1 = fIxEtaPhi_1[ix1];
      Int_t ixPt_1     = fIxPt_1[ix1];
      Float_t corr_1   = fCorrection_1[ix1];
      Float_t pt_1     = fPt_1[ix1];

      for (Int_t ix2 = ix1+1; ix2 < fNoOfTracks1; ix2++) {
        /* excluded self correlations */
        Float_t corr      = corr_1 * fCorrection_1[ix2];

        /* apply the pair correction if applicable */
        if (effcorr != NULL) {
          Int_t ieta1 = Int_t(ixEtaPhi_1/fNBins_phi_1);
          Int_t iphi1 = ixEtaPhi_1 % fNBins_phi_1;
          Int_t ieta2 = Int_t(fIxEtaPhi_1[ix2]/fNBins_phi_1);
          Int_t iphi2 = fIxEtaPhi_1[ix2] % fNBins_phi_1;
          Int_t deltaetabin = ieta1-ieta2+fNBins_eta_1;
          Int_t ixdeltaphi = iphi1-iphi2; if (ixdeltaphi < 0) ixdeltaphi += fNBins_phi_1;
          bins[0] = deltaetabin;
          bins[1] = ixdeltaphi+1;
          bins[2] = ixPt_1+1;
          bins[3] = fIxPt_1[ix2]+1;
          corr = corr / effcorr->GetBinContent(bins);
        }

        Int_t ij        = getLinearSymmIndex(ixEtaPhi_1, fNBins_etaPhi_1, fIxEtaPhi_1[ix2], fNBins_etaPhi_2);

        fN2_12                                              += 2*corr;
        fN2_12_vsEtaPhi[ij]                                 += corr;
        Float_t ptpt                                         = pt_1*fPt_1[ix2];
        fSum2PtPt_12                                        += 2*corr*ptpt;
        fSum2PtN_12                                         += corr*pt_1 + corr*fPt_1[ix2];
        fSum2NPt_12                                         += corr*fPt_1[ix2] + corr*pt_1;
        fSum2PtPt_12_vsEtaPhi[ij]                           += corr*ptpt;
        fSum2PtN_12_vsEtaPhi[ij]                            += corr*pt_1;
        fSum2NPt_12_vsEtaPhi[ij]                            += corr*fPt_1[ix2];
        fN2_12_vsPtPt[ixPt_1*fNBins_pt_1 + fIxPt_1[ix2]]    += corr;
        fN2_12_vsPtPt[fIxPt_1[ix2]*fNBins_pt_1 + ixPt_1]    += corr;

        fNnw2_12                    += 2;
        fSum2PtPtnw_12              += 2*ptpt;
        fSum2PtNnw_12               += pt_1;
        fSum2PtNnw_12               += fPt_1[ix2];
        fSum2NPtnw_12               += fPt_1[ix2];
        fSum2NPtnw_12               += pt_1;
      } //ix2
    } //ix1
  }
  else if (bank == 2) {
    /* let's select the pair efficiency correction histogram */
    const THn *effcorr = NULL;
    if (fRequestedCharge_2 > 0)
      effcorr = fPairsEfficiency_PP;
    else
      effcorr = fPairsEfficiency_MM;
    Int_t bins[4];

    for (Int_t ix1 = 0; ix1 < fNoOfTracks2; ix1++)
    {
      Int_t ixEtaPhi_1 = fIxEtaPhi_2[ix1];
      Int_t ixPt_1     = fIxPt_2[ix1];
      Float_t corr_1   = fCorrection_2[ix1];
      Float_t pt_1     = fPt_2[ix1];

      for (Int_t ix2 = ix1+1; ix2 < fNoOfTracks2; ix2++) {
        /* excluded self correlations */
        Float_t corr      = corr_1 * fCorrection_2[ix2];

        /* apply the pair correction if applicable */
        if (effcorr != NULL) {
          Int_t ieta1 = Int_t(ixEtaPhi_1/fNBins_phi_1);
          Int_t iphi1 = ixEtaPhi_1 % fNBins_phi_1;
          Int_t ieta2 = Int_t(fIxEtaPhi_2[ix2]/fNBins_phi_1);
          Int_t iphi2 = fIxEtaPhi_2[ix2] % fNBins_phi_1;
          Int_t deltaetabin = ieta1-ieta2+fNBins_eta_1;
          Int_t ixdeltaphi = iphi1-iphi2; if (ixdeltaphi < 0) ixdeltaphi += fNBins_phi_1;
          bins[0] = deltaetabin;
          bins[1] = ixdeltaphi+1;
          bins[2] = ixPt_1+1;
          bins[3] = fIxPt_2[ix2]+1;
          corr = corr / effcorr->GetBinContent(bins);
        }

        Int_t ij        = getLinearSymmIndex(ixEtaPhi_1, fNBins_etaPhi_1, fIxEtaPhi_2[ix2], fNBins_etaPhi_2);

        fN2_12                                              += 2*corr;
        fN2_12_vsEtaPhi[ij]                                 += corr;
        Float_t ptpt                                         = pt_1*fPt_2[ix2];
        fSum2PtPt_12                                        += 2*corr*ptpt;
        fSum2PtN_12                                         += corr*pt_1 + corr*fPt_2[ix2];
        fSum2NPt_12                                         += corr*fPt_2[ix2] + corr*pt_1;
        fSum2PtPt_12_vsEtaPhi[ij]                           += corr*ptpt;
        fSum2PtN_12_vsEtaPhi[ij]                            += corr*pt_1;
        fSum2NPt_12_vsEtaPhi[ij]                            += corr*fPt_2[ix2];
        fN2_12_vsPtPt[ixPt_1*fNBins_pt_2 + fIxPt_2[ix2]]    += corr;
        fN2_12_vsPtPt[fIxPt_2[ix2]*fNBins_pt_2 + ixPt_1]    += corr;

        fNnw2_12                    += 2;
        fSum2PtPtnw_12              += 2*ptpt;
        fSum2PtNnw_12               += pt_1;
        fSum2PtNnw_12               += fPt_2[ix2];
        fSum2NPtnw_12               += fPt_2[ix2];
        fSum2NPtnw_12               += pt_1;
      } //ix2
    } //ix1
  }
}

/// \brief Process track combinations with the same charge
/// The half symmetrizing process for memory reduction is not used
/// \param bank the tracks bank to use
void AliDptDptCorrelations::ProcessNotHalfSymmLikeSignPairs(Int_t bank) {
  /* pair with same charge. The track list should be identical */
  if (bank == 1) {
    /* let's select the pair efficiency correction histogram */
    const THn *effcorr = NULL;
    if (fRequestedCharge_1 > 0)
      effcorr = fPairsEfficiency_PP;
    else
      effcorr = fPairsEfficiency_MM;
    Int_t bins[4];

    /* we use only the bank one of tracks */
    for (Int_t ix1 = 0; ix1 < fNoOfTracks1; ix1++)
    {
      Int_t ixEtaPhi_1 = fIxEtaPhi_1[ix1];
      Int_t ixPt_1     = fIxPt_1[ix1];
      Float_t corr_1   = fCorrection_1[ix1];
      Float_t pt_1     = fPt_1[ix1];

      for (Int_t ix2 = ix1+1; ix2 < fNoOfTracks1; ix2++) {
        /* excluded self correlations */
        Float_t corr      = corr_1 * fCorrection_1[ix2];

        /* apply the pair correction if applicable */
        if (effcorr != NULL) {
          Int_t ieta1 = Int_t(ixEtaPhi_1/fNBins_phi_1);
          Int_t iphi1 = ixEtaPhi_1 % fNBins_phi_1;
          Int_t ieta2 = Int_t(fIxEtaPhi_1[ix2]/fNBins_phi_1);
          Int_t iphi2 = fIxEtaPhi_1[ix2] % fNBins_phi_1;
          Int_t deltaetabin = ieta1-ieta2+fNBins_eta_1;
          Int_t ixdeltaphi = iphi1-iphi2; if (ixdeltaphi < 0) ixdeltaphi += fNBins_phi_1;
          bins[0] = deltaetabin;
          bins[1] = ixdeltaphi+1;
          bins[2] = ixPt_1+1;
          bins[3] = fIxPt_1[ix2]+1;
          corr = corr / effcorr->GetBinContent(bins);
        }

        Int_t ij                                             = ixEtaPhi_1*fNBins_etaPhi_1 + fIxEtaPhi_1[ix2];

        fN2_12                                              += corr;
        fN2_12_vsEtaPhi[ij]                                 += corr;
        Float_t ptpt                                         = pt_1*fPt_1[ix2];
        fSum2PtPt_12                                        += corr*ptpt;
        fSum2PtN_12                                         += corr*pt_1;
        fSum2NPt_12                                         += corr*fPt_1[ix2];
        fSum2PtPt_12_vsEtaPhi[ij]                           += corr*ptpt;
        fSum2PtN_12_vsEtaPhi[ij]                            += corr*pt_1;
        fSum2NPt_12_vsEtaPhi[ij]                            += corr*fPt_1[ix2];
        fN2_12_vsPtPt[ixPt_1*fNBins_pt_1 + fIxPt_1[ix2]]    += corr;

        fNnw2_12                    += 1;
        fSum2PtPtnw_12              += ptpt;
        fSum2PtNnw_12               += pt_1;
        fSum2NPtnw_12               += fPt_1[ix2];
      } //ix2
    } //ix1
  }
  else if (bank == 2) {
    /* let's select the pair efficiency correction histogram */
    const THn *effcorr = NULL;
    if (fRequestedCharge_2 > 0)
      effcorr = fPairsEfficiency_PP;
    else
      effcorr = fPairsEfficiency_MM;
    Int_t bins[4];

    for (Int_t ix1 = 0; ix1 < fNoOfTracks2; ix1++)
    {
      Int_t ixEtaPhi_1 = fIxEtaPhi_2[ix1];
      Int_t ixPt_1     = fIxPt_2[ix1];
      Float_t corr_1   = fCorrection_2[ix1];
      Float_t pt_1     = fPt_2[ix1];

      for (Int_t ix2 = ix1+1; ix2 < fNoOfTracks2; ix2++) {
        /* excluded self correlations */
        Float_t corr      = corr_1 * fCorrection_2[ix2];

        /* apply the pair correction if applicable */
        if (effcorr != NULL) {
          Int_t ieta1 = Int_t(ixEtaPhi_1/fNBins_phi_1);
          Int_t iphi1 = ixEtaPhi_1 % fNBins_phi_1;
          Int_t ieta2 = Int_t(fIxEtaPhi_2[ix2]/fNBins_phi_1);
          Int_t iphi2 = fIxEtaPhi_2[ix2] % fNBins_phi_1;
          Int_t deltaetabin = ieta1-ieta2+fNBins_eta_1;
          Int_t ixdeltaphi = iphi1-iphi2; if (ixdeltaphi < 0) ixdeltaphi += fNBins_phi_1;
          bins[0] = deltaetabin;
          bins[1] = ixdeltaphi+1;
          bins[2] = ixPt_1+1;
          bins[3] = fIxPt_2[ix2]+1;
          corr = corr / effcorr->GetBinContent(bins);
        }

        Int_t ij                                             = ixEtaPhi_1*fNBins_etaPhi_2 + fIxEtaPhi_2[ix2];

        fN2_12                                              += corr;
        fN2_12_vsEtaPhi[ij]                                 += corr;
        Float_t ptpt                                         = pt_1*fPt_2[ix2];
        fSum2PtPt_12                                        += corr*ptpt;
        fSum2PtN_12                                         += corr*pt_1;
        fSum2NPt_12                                         += corr*fPt_2[ix2];
        fSum2PtPt_12_vsEtaPhi[ij]                           += corr*ptpt;
        fSum2PtN_12_vsEtaPhi[ij]                            += corr*pt_1;
        fSum2NPt_12_vsEtaPhi[ij]                            += corr*fPt_2[ix2];
        fN2_12_vsPtPt[ixPt_1*fNBins_pt_2 + fIxPt_2[ix2]]    += corr;

        fNnw2_12                    += 1;
        fSum2PtPtnw_12              += ptpt;
        fSum2PtNnw_12               += pt_1;
        fSum2NPtnw_12               += fPt_2[ix2];
      } //ix2
    } //ix1
  }
}

/// \brief Flag the potential products of conversions and / or resonances
void AliDptDptCorrelations::FlagConversionsAndResonances() {
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
void AliDptDptCorrelations::ProcessUnlikeSignPairs() {
  AliInfo("");
  /* flag the resonances / conversion candidates */
  FlagConversionsAndResonances();

  /* pair with different charges. Both track list are different */
  for (Int_t ix1 = 0; ix1 < fNoOfTracks1; ix1++)
  {
    /* let's select the pair efficiency correction histogram */
    const THn *effcorr = NULL;
    if (fRequestedCharge_1 > 0)
      effcorr = fPairsEfficiency_PM;
    else
      effcorr = fPairsEfficiency_MP;
    Int_t bins[4];

    Int_t ixEtaPhi_1 = fIxEtaPhi_1[ix1];
    Int_t ixPt_1     = fIxPt_1[ix1];
    Float_t corr_1   = fCorrection_1[ix1];
    Float_t pt_1     = fPt_1[ix1];

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

        /* apply the pair correction if applicable */
        if (effcorr != NULL) {
          Int_t ieta1 = Int_t(ixEtaPhi_1/fNBins_phi_1);
          Int_t iphi1 = ixEtaPhi_1 % fNBins_phi_1;
          Int_t ieta2 = Int_t(fIxEtaPhi_2[ix2]/fNBins_phi_1);
          Int_t iphi2 = fIxEtaPhi_2[ix2] % fNBins_phi_1;
          Int_t deltaetabin = ieta1-ieta2+fNBins_eta_1;
          Int_t ixdeltaphi = iphi1-iphi2; if (ixdeltaphi < 0) ixdeltaphi += fNBins_phi_1;
          bins[0] = deltaetabin;
          bins[1] = ixdeltaphi+1;
          bins[2] = ixPt_1+1;
          bins[3] = fIxPt_2[ix2]+1;
          corr = corr / effcorr->GetBinContent(bins);
        }

        Int_t ij        = getLinearSymmIndex(ixEtaPhi_1, fNBins_etaPhi_1, fIxEtaPhi_2[ix2], fNBins_etaPhi_2);

        fN2_12                                              += 2*corr;
        fN2_12_vsEtaPhi[ij]                                 += corr;
        Float_t ptpt                                         = pt_1*fPt_2[ix2];
        fSum2PtPt_12                                        += 2*corr*ptpt;
        fSum2PtN_12                                         += corr*pt_1 + corr*fPt_2[ix2];
        fSum2NPt_12                                         += corr*fPt_2[ix2] + corr*pt_1 ;
        fSum2PtPt_12_vsEtaPhi[ij]                           += corr*ptpt;
        fSum2PtN_12_vsEtaPhi[ij]                            += corr*pt_1;
        fSum2NPt_12_vsEtaPhi[ij]                            += corr*fPt_2[ix2];
        fN2_12_vsPtPt[ixPt_1*fNBins_pt_2 + fIxPt_2[ix2]]    += corr;
        fN2_12_vsPtPt[fIxPt_2[ix2]*fNBins_pt_1 + ixPt_1]    += corr;

        fNnw2_12                    += 2;
        fSum2PtPtnw_12              += 2*ptpt;
        fSum2PtNnw_12               += pt_1;
        fSum2PtNnw_12               += fPt_2[ix2];
        fSum2NPtnw_12               += fPt_2[ix2];
        fSum2NPtnw_12               += pt_1;
      }
    } //ix2
  } //ix1
}

/// \brief Process track combinations with oposite charge
/// The half symmetrizing process for memory reduction is not used
void AliDptDptCorrelations::ProcessNotHalfSymmUnlikeSignPairs() {
  AliInfo("");
  /* flag the resonances / conversion candidates */
  FlagConversionsAndResonances();

  /* pair with different charges. Both track list are different */
  for (Int_t ix1 = 0; ix1 < fNoOfTracks1; ix1++)
  {
    /* let's select the pair efficiency correction histogram */
    const THn *effcorr = NULL;
    if (fRequestedCharge_1 > 0)
      effcorr = fPairsEfficiency_PM;
    else
      effcorr = fPairsEfficiency_MP;
    Int_t bins[4];

    Int_t ixEtaPhi_1 = fIxEtaPhi_1[ix1];
    Int_t ixPt_1     = fIxPt_1[ix1];
    Float_t corr_1   = fCorrection_1[ix1];
    Float_t pt_1     = fPt_1[ix1];

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

        /* apply the pair correction if applicable */
        if (effcorr != NULL) {
          Int_t ieta1 = Int_t(ixEtaPhi_1/fNBins_phi_1);
          Int_t iphi1 = ixEtaPhi_1 % fNBins_phi_1;
          Int_t ieta2 = Int_t(fIxEtaPhi_2[ix2]/fNBins_phi_1);
          Int_t iphi2 = fIxEtaPhi_2[ix2] % fNBins_phi_1;
          Int_t deltaetabin = ieta1-ieta2+fNBins_eta_1;
          Int_t ixdeltaphi = iphi1-iphi2; if (ixdeltaphi < 0) ixdeltaphi += fNBins_phi_1;
          bins[0] = deltaetabin;
          bins[1] = ixdeltaphi+1;
          bins[2] = ixPt_1+1;
          bins[3] = fIxPt_2[ix2]+1;
          corr = corr / effcorr->GetBinContent(bins);
        }

        Int_t ij                                             = ixEtaPhi_1*fNBins_etaPhi_1 + fIxEtaPhi_2[ix2];

        fN2_12                                              += corr;
        fN2_12_vsEtaPhi[ij]                                 += corr;
        Float_t ptpt                                         = pt_1*fPt_2[ix2];
        fSum2PtPt_12                                        += corr*ptpt;
        fSum2PtN_12                                         += corr*pt_1;
        fSum2NPt_12                                         += corr*fPt_2[ix2];
        fSum2PtPt_12_vsEtaPhi[ij]                           += corr*ptpt;
        fSum2PtN_12_vsEtaPhi[ij]                            += corr*pt_1;
        fSum2NPt_12_vsEtaPhi[ij]                            += corr*fPt_2[ix2];
        fN2_12_vsPtPt[ixPt_1*fNBins_pt_2 + fIxPt_2[ix2]]    += corr;

        fNnw2_12                    += 1;
        fSum2PtPtnw_12              += ptpt;
        fSum2PtNnw_12               += pt_1;
        fSum2NPtnw_12               += fPt_2[ix2];
      }
    } //ix2
  } //ix1
}

/// \brief Fill the histograms once the whole process is finished
void  AliDptDptCorrelations::FinalizeProcess()
{

  if (fSinglesOnly) {
    fillHistoWithArray(fhN1_1_vsPt,              fN1_1_vsPt,            fNBins_pt_1);
    fillHistoWithArray(fhN1_1_vsZEtaPhiPt,       fN1_1_vsZEtaPhiPt,     fNBins_vertexZ, fNBins_etaPhi_1, fNBins_pt_1);
    fillHistoWithArray(fhN1_2_vsPt,              fN1_2_vsPt,            fNBins_pt_2);
    fillHistoWithArray(fhN1_2_vsZEtaPhiPt,       fN1_2_vsZEtaPhiPt,     fNBins_vertexZ, fNBins_etaPhi_2, fNBins_pt_2);

    /* for the time being, the errors are trivial so, we do not use results file space */
    fhN1_1_vsPt->Sumw2(false);
    fhN1_1_vsZEtaPhiPt->Sumw2(false);
    fhN1_2_vsPt->Sumw2(false);
    fhN1_2_vsZEtaPhiPt->Sumw2(false);
  }
  else {
    fillHistoWithArray(fhN1_1_vsEtaPhi,          fN1_1_vsEtaPhi,        fNBins_eta_1,   fNBins_phi_1);
    fillHistoWithArray(fhSum1Pt_1_vsEtaPhi,      fSum1Pt_1_vsEtaPhi,    fNBins_eta_1,   fNBins_phi_1);
    fillHistoWithArray(fhN1_2_vsEtaPhi,          fN1_2_vsEtaPhi,        fNBins_eta_2,   fNBins_phi_2);
    fillHistoWithArray(fhSum1Pt_2_vsEtaPhi,      fSum1Pt_2_vsEtaPhi,    fNBins_eta_2,   fNBins_phi_2);

    fillHistoWithArray(fhN2_12_vsEtaPhi,         fN2_12_vsEtaPhi,       fNBins_etaPhi_12);
    fillHistoWithArray(fhSum2PtPt_12_vsEtaPhi,   fSum2PtPt_12_vsEtaPhi, fNBins_etaPhi_12);
    fillHistoWithArray(fhSum2PtN_12_vsEtaPhi,    fSum2PtN_12_vsEtaPhi,  fNBins_etaPhi_12);
    fillHistoWithArray(fhSum2NPt_12_vsEtaPhi,    fSum2NPt_12_vsEtaPhi,  fNBins_etaPhi_12);
    fillHistoWithArray(fhN2_12_vsPtPt,           fN2_12_vsPtPt,         fNBins_pt_1,    fNBins_pt_2);

    /* for the time being, the errors are trivial so, we do not use results file space */
    fhN1_1_vsEtaPhi->Sumw2(false);
    fhSum1Pt_1_vsEtaPhi->Sumw2(false);
    fhN1_2_vsEtaPhi->Sumw2(false);
    fhSum1Pt_2_vsEtaPhi->Sumw2(false);
    fhN2_12_vsEtaPhi->Sumw2(false);
    fhSum2PtPt_12_vsEtaPhi->Sumw2(false);
    fhSum2PtN_12_vsEtaPhi->Sumw2(false);
    fhSum2NPt_12_vsEtaPhi->Sumw2(false);
    fhN2_12_vsPtPt->Sumw2(false);
  }
}

//Tools
//===================================================================================================
/// \brief Fills one dimension histogram from an array of float
/// \param h one dimension histogram to fill
/// \param array the source array
/// \param size the number of bins in the histogram
void  AliDptDptCorrelations::fillHistoWithArray(TH1 * h, float * array, int size)
{
  int i, i1;
  float v1, ev1, v2, ev2, sum, esum;
  for (i=0, i1=1; i<size; ++i,++i1) {
    v1  = array[i]; ev1 = sqrt(v1);
    v2  = h->GetBinContent(i1);
    ev2 = h->GetBinError(i1);
    sum = v1 + v2;
    esum = sqrt(ev1*ev1+ev2*ev2);
    h->SetBinContent(i1,sum);
    h->SetBinError(i1,esum);
  }
}

/// \brief Fills two dimensions histogram from an array of float
/// \param h two dimensions histogram to fill
/// \param array the source array
/// \param size1 the number of bins in the histogram x dimension
/// \param size2 the number of bins in the histogram y dimension
void  AliDptDptCorrelations::fillHistoWithArray(TH2 * h, float * array, int size1, int size2)
{
  int i, i1;
  int j, j1;
  float v1, ev1, v2, ev2, sum, esum;
  for (i=0, i1=1; i<size1; ++i,++i1) {
    for (j=0, j1=1; j<size2; ++j,++j1) {
      v1  = array[i*size2+j]; ev1 = sqrt(v1);
      v2  = h->GetBinContent(i1,j1);
      ev2 = h->GetBinError(i1,j1);
      sum = v1 + v2;
      esum = sqrt(ev1*ev1+ev2*ev2);
      h->SetBinContent(i1,j1,sum);
      h->SetBinError(i1,j1,esum);
    }
  }
}

/// \brief Fills three dimensions histogram from an array of float
/// \param h three dimensions histogram to fill
/// \param array the source array
/// \param size1 the number of bins in the histogram x dimension
/// \param size2 the number of bins in the histogram y dimension
/// \param size3 the number of bins in the histogram z dimension
void  AliDptDptCorrelations::fillHistoWithArray(TH3 * h, float * array, int size1, int size2, int size3)
{
  int i, i1;
  int j, j1;
  int k, k1;
  float v1, ev1, v2, ev2, sum, esum;
  int size23 = size2*size3;
  for (i=0, i1=1; i<size1; ++i,++i1) {
    for (j=0, j1=1; j<size2; ++j,++j1) {
      for (k=0, k1=1; k<size3; ++k,++k1) {
        v1  = array[i*size23+j*size3+k]; ev1 = sqrt(v1);
        v2  = h->GetBinContent(i1,j1,k1);
        ev2 = h->GetBinError(i1,j1,k1);
        sum = v1 + v2;
        esum = sqrt(ev1*ev1+ev2*ev2);
        h->SetBinContent(i1,j1,k1,sum);
        h->SetBinError(i1,j1,k1,esum);
      }
    }
  }
}

/// \brief Fills one dimension histogram from an array of double
/// \param h one dimension histogram to fill
/// \param array the source array
/// \param size the number of bins in the histogram x dimension
void  AliDptDptCorrelations::fillHistoWithArray(TH1 * h, double * array, int size)
{
  int i, i1;
  double v1, ev1, v2, ev2, sum, esum;
  for (i=0, i1=1; i<size; ++i,++i1) {
    v1  = array[i]; ev1 = sqrt(v1);
    v2  = h->GetBinContent(i1);
    ev2 = h->GetBinError(i1);
    sum = v1 + v2;
    esum = sqrt(ev1*ev1+ev2*ev2);
    h->SetBinContent(i1,sum);
    h->SetBinError(i1,esum);
  }
}

/// \brief Fills two dimensions histogram from an array of double
/// \param h two dimensions histogram to fill
/// \param array the source array
/// \param size1 the number of bins in the histogram x dimension
/// \param size2 the number of bins in the histogram y dimension
void  AliDptDptCorrelations::fillHistoWithArray(TH2 * h, double * array, int size1, int size2)
{
  int i, i1;
  int j, j1;
  double v1, ev1, v2, ev2, sum, esum;
  for (i=0, i1=1; i<size1; ++i,++i1) {
    for (j=0, j1=1; j<size2; ++j,++j1) {
      v1  = array[i*size2+j]; ev1 = sqrt(v1);
      v2  = h->GetBinContent(i1,j1);
      ev2 = h->GetBinError(i1,j1);
      sum = v1 + v2;
      esum = sqrt(ev1*ev1+ev2*ev2);
      h->SetBinContent(i1,j1,sum);
      h->SetBinError(i1,j1,esum);
    }
  }
}

/// \brief Fills three dimensions histogram from an array of double
/// \param h three dimensions histogram to fill
/// \param array the source array
/// \param size1 the number of bins in the histogram x dimension
/// \param size2 the number of bins in the histogram y dimension
/// \param size3 the number of bins in the histogram z dimension
void  AliDptDptCorrelations::fillHistoWithArray(TH3 * h, double * array, int size1, int size2, int size3)
{
  int i, i1;
  int j, j1;
  int k, k1;
  double v1, ev1, v2, ev2, sum, esum;
  int size23 = size2*size3;
  for (i=0, i1=1; i<size1; ++i,++i1) {
    for (j=0, j1=1; j<size2; ++j,++j1) {
      for (k=0, k1=1; k<size3; ++k,++k1) {
        v1  = array[i*size23+j*size3+k]; ev1 = sqrt(v1);
        v2  = h->GetBinContent(i1,j1,k1);
        ev2 = h->GetBinError(i1,j1,k1);
        sum = v1 + v2;
        esum = sqrt(ev1*ev1+ev2*ev2);
        h->SetBinContent(i1,j1,k1,sum);
        h->SetBinError(i1,j1,k1,esum);
      }
    }
  }
}



/// \cond CLASSIMP
ClassImp(AliDptDptCorrelations);
/// \endcond
