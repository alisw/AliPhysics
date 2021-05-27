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

/// \file AliTwoParticleCorrelationsBase.cxx
/// \brief Implementation of the AliTwoParticleCorrelationsBase base class

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
#include "AliTwoParticleCorrelationsBase.h"

Int_t AliTwoParticleCorrelationsBase::fgkNoOfResonances = 4; ///< four resonances / conversions for the time being
Double_t AliTwoParticleCorrelationsBase::fgkMass[16] = {/* photon */ 0.0, /* k0 */ 0.4976, /* lambda */ 1.115, /* rho */ 0.775};
Double_t AliTwoParticleCorrelationsBase::fgkChildMass[2][16] = {
    {0.510e-3, 0.1396, 0.1396, 0.1396},
    {0.510e-3, 0.1396, 0.9383, 0.1396}
};
Double_t AliTwoParticleCorrelationsBase::fgkMassThreshold[16] = {0.04,0.01,0.05,0.04};

/// Default constructor for object serialization
AliTwoParticleCorrelationsBase::AliTwoParticleCorrelationsBase() :
    TNamed(),
    fOutput(NULL),
    fVertexZ(0.0),
    fIxVertexZ(0),
    fCentrality(0.0),
    fSinglesOnly(kTRUE),
    fUseWeights(kFALSE),
    fUseSimulation(kFALSE),
    /* the arrays with tracks 1 and 2 information */
    fArraySize(5*1024),
    fNoOfTracks1(0),
    fFlags_1(nullptr),
    fPt_1(nullptr),
    fEta_1(nullptr),
    fPhi_1(nullptr),
    fNoOfTracks2(0),
    fFlags_2(nullptr),
    fPt_2(nullptr),
    fEta_2(nullptr),
    fPhi_2(nullptr),
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
    fhN1_1_vsPt(NULL),
    fhN1_1_vsEtaPhi(NULL),
    fhSum1Pt_1_vsEtaPhi(NULL),
    fhN1_1_vsZEtaPhiPt(NULL),
    fhSum1Pt_1_vsZEtaPhiPt(nullptr),
    fhN1_2_vsPt(NULL),
    fhN1_2_vsEtaPhi(NULL),
    fhSum1Pt_2_vsEtaPhi(NULL),
    fhN1_2_vsZEtaPhiPt(NULL),
    fhSum1Pt_2_vsZEtaPhiPt(nullptr),
    /* vs centrality profiles */
    fhN1_1_vsC(NULL),
    fhSum1Pt_1_vsC(NULL),
    fhN1nw_1_vsC(NULL),
    fhSum1Ptnw_1_vsC(NULL),
    fhN1_2_vsC(NULL),
    fhSum1Pt_2_vsC(NULL),
    fhN1nw_2_vsC(NULL),
    fhSum1Ptnw_2_vsC(NULL),
    fhResonanceRoughMasses(NULL),
    fhResonanceMasses(NULL),
    fhDiscardedResonanceMasses(NULL)
{
  for (Int_t ires = 0; ires < 16; ires++)
    fThresholdMult[ires] = 0x0;
}

/// Normal constructor
/// \param name the name for the object instance
AliTwoParticleCorrelationsBase::AliTwoParticleCorrelationsBase(const char *name) :
    TNamed(name,name),
    fOutput(NULL),
    fVertexZ(0.0),
    fIxVertexZ(0),
    fCentrality(0.0),
    fSinglesOnly(kTRUE),
    fUseWeights(kFALSE),
    fUseSimulation(kFALSE),
    /* the arrays with tracks 1 and 2 information */
    fArraySize(5*1024),
    fNoOfTracks1(0),
    fFlags_1(nullptr),
    fPt_1(nullptr),
    fEta_1(nullptr),
    fPhi_1(nullptr),
    fNoOfTracks2(0),
    fFlags_2(nullptr),
    fPt_2(nullptr),
    fEta_2(nullptr),
    fPhi_2(nullptr),
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
    /* accumulated values */
    fN1_1(0.0),
    fN1_2(0.0),
    fNnw1_1(0.0),
    fNnw1_2(0.0),
    fSum1Pt_1(0.0),
    fSum1Pt_2(0.0),
    fSum1Ptnw_1(0.0),
    fSum1Ptnw_2(0.0),
    /* storage arrays */
    /* histograms */
    fhN1_1_vsPt(NULL),
    fhN1_1_vsEtaPhi(NULL),
    fhSum1Pt_1_vsEtaPhi(NULL),
    fhN1_1_vsZEtaPhiPt(NULL),
    fhSum1Pt_1_vsZEtaPhiPt(nullptr),
    fhN1_2_vsPt(NULL),
    fhN1_2_vsEtaPhi(NULL),
    fhSum1Pt_2_vsEtaPhi(NULL),
    fhN1_2_vsZEtaPhiPt(NULL),
    fhSum1Pt_2_vsZEtaPhiPt(nullptr),
    /* vs centrality profiles */
    fhN1_1_vsC(NULL),
    fhSum1Pt_1_vsC(NULL),
    fhN1nw_1_vsC(NULL),
    fhSum1Ptnw_1_vsC(NULL),
    fhN1_2_vsC(NULL),
    fhSum1Pt_2_vsC(NULL),
    fhN1nw_2_vsC(NULL),
    fhSum1Ptnw_2_vsC(NULL),
    fhResonanceRoughMasses(NULL),
    fhResonanceMasses(NULL),
    fhDiscardedResonanceMasses(NULL)
{
  for (Int_t ires = 0; ires < 16; ires++)
    fThresholdMult[ires] = 0x0;
}

/// \brief Default destructor
/// Deallocates the allocated memory
AliTwoParticleCorrelationsBase::~AliTwoParticleCorrelationsBase() {
  /* track 1 storage */
  delete [] fFlags_1;
  delete [] fPt_1;
  delete [] fEta_1;
  delete [] fPhi_1;

  /* track 2 storage */
  delete [] fFlags_2;
  delete [] fPt_2;
  delete [] fEta_2;
  delete [] fPhi_2;

  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
      delete fOutput;
  }
}


/// \brief Establishes the binning configuration
/// \param confstring string containing the binning configuration parameters
Bool_t AliTwoParticleCorrelationsBase::ConfigureBinning(const char *confstring) {
#define DPTDPTCORRBINCONFIGPAR 10

  Double_t min_pt = 0.0, max_pt = 0.0, width_pt = 0.0;
  Double_t min_eta = 0.0, max_eta = 0.0, width_eta = 0.0;
  Int_t    nBins_phi = 0;

  TString stritem = TString(confstring);

  TObjArray *a = stritem.Tokenize(",");
  if (a->GetEntries() != DPTDPTCORRBINCONFIGPAR) {
    delete a;
    return false;
  }
  sscanf(stritem, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d",
      &fMin_vertexZ, &fMax_vertexZ, &fWidth_vertexZ,
      &min_pt, &max_pt, &width_pt,
      &min_eta, &max_eta, &width_eta, &nBins_phi);
  delete a;

  fMin_pt_1 = fMin_pt_2 = min_pt;
  fMax_pt_1 = fMax_pt_2 = max_pt;
  fWidth_pt_1 = fWidth_pt_2 = width_pt;
  fMin_eta_1 = fMin_eta_2 = min_eta;
  fMax_eta_1 = fMax_eta_2 = max_eta;
  fWidth_eta_1 = fWidth_eta_2 = width_eta;
  fNBins_phi_1 = fNBins_phi_2 = nBins_phi;
  return true;
}

/// \brief Establishes the resonances rejection configuration
/// \param confstring string containing the resonances rejection configuration parameters
/// Basically one digit per supported resonance and the digit is the factor in one fourth
/// of the modulus around the resonance mass
void AliTwoParticleCorrelationsBase::ConfigureResonances(const char *confstring) {

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
TString AliTwoParticleCorrelationsBase::GetBinningConfigurationString() const {
  if (fMin_pt_1 != fMin_pt_2 || fMax_pt_2 != fMax_pt_2 || fWidth_pt_1 != fWidth_pt_2 ||
      fMin_eta_1 != fMin_eta_2 || fMax_eta_1 != fMax_eta_2 || fWidth_eta_1 != fWidth_eta_2 ||
      fNBins_phi_1 != fNBins_phi_2) {
    return TString::Format("WrongAsymmetricBinning:%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%d",
        fMin_vertexZ, fMax_vertexZ, fWidth_vertexZ,
        fMin_pt_1, fMax_pt_1, fWidth_pt_1,
        fMin_eta_1, fMax_eta_1, fWidth_eta_1, fNBins_phi_1);
  }
  else {
    return TString::Format("%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%d",
        fMin_vertexZ, fMax_vertexZ, fWidth_vertexZ,
        fMin_pt_1, fMax_pt_1, fWidth_pt_1,
        fMin_eta_1, fMax_eta_1, fWidth_eta_1, fNBins_phi_1);
  }
}

/// \brief Build the resonances rejection configuration string
/// \return the configuration string corresponding to the current resonance rejection configuration
TString AliTwoParticleCorrelationsBase::GetResonancesConfigurationString() const {
  TString str = "resonances:";

  for (Int_t ires = 0; ires < fgkNoOfResonances; ires++) {
    str += TString::Format("%01d",fThresholdMult[fgkNoOfResonances - 1 - ires]);
  }
  return str;
}



/// \brief Initializes the member data structures
/// Allocates the needed memory an create the output histograms.
void AliTwoParticleCorrelationsBase::Initialize()
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

  /* incorporate configuration parameters */
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
  fFlags_1 = new UInt_t[fArraySize];
  fPt_1 = new Float_t[fArraySize];
  fEta_1 = new Float_t[fArraySize];
  fPhi_1 = new Float_t[fArraySize];

  /* track 2 storage */
  fFlags_2 = new UInt_t[fArraySize];
  fPt_2 = new Float_t[fArraySize];
  fEta_2 = new Float_t[fArraySize];
  fPhi_2 = new Float_t[fArraySize];

  if (fSinglesOnly)  {
    fhN1_1_vsPt = new TH1F("n1_1_vsPt","#LT n_{1} #GT;p_{t,1} (GeV/c);#LT n_{1} #GT", fNBins_pt_1,  fMin_pt_1,  fMax_pt_1);
    fhN1_1_vsZEtaPhiPt = new TH3F("n1_1_vsZ_vsEtaPhi_vsPt","#LT n_{1} #GT;vtx_{z};#eta_{1}#times#varphi_{1};p_{t,1} (GeV/c)",
        fNBins_vertexZ,fMin_vertexZ,fMax_vertexZ, fNBins_etaPhi_1,0.0,Double_t(fNBins_etaPhi_1),fNBins_pt_1,fMin_pt_1,fMax_pt_1);
    fhSum1Pt_1_vsZEtaPhiPt = new TH3F("sumPt1_1_vsZ_vsEtaPhi_vsPt","#LT #Sigma p_{t,1} #GT;vtx_{z};#eta_{1}#times#varphi_{1};p_{t,1} (GeV/c)",
        fNBins_vertexZ,fMin_vertexZ,fMax_vertexZ, fNBins_etaPhi_1,0.0,Double_t(fNBins_etaPhi_1),fNBins_pt_1,fMin_pt_1,fMax_pt_1);
    fhN1_1_vsEtaPhi = new TH2F("n1_1_vsEtaPhi","#LT n_{1} #GT;#eta_{1};#varphi_{1} (radian);#LT n_{1} #GT",
        fNBins_eta_1, fMin_eta_1, fMax_eta_1,  fNBins_phi_1, fMin_phi_1, fMax_phi_1);
    fhSum1Pt_1_vsEtaPhi = new TH2F("sumPt_1_vsEtaPhi","#LT #Sigma p_{t,1} #GT;#eta_{1};#varphi_{1} (radian);#LT #Sigma p_{t,1} #GT (GeV/c)",
        fNBins_eta_1,fMin_eta_1,fMax_eta_1,fNBins_phi_1,fMin_phi_1,fMax_phi_1);

    fhN1_2_vsPt = new TH1F("n1_2_vsPt","#LT n_{1} #GT;p_{t,2} (GeV/c);#LT n_{1} #GT", fNBins_pt_2,  fMin_pt_2,  fMax_pt_2);
    fhN1_2_vsZEtaPhiPt = new TH3F("n1_2_vsZ_vsEtaPhi_vsPt","#LT n_{2} #GT;vtx_{z};#eta_{2}#times#varphi_{2};p_{t,2} (GeV/c)",
        fNBins_vertexZ,fMin_vertexZ,fMax_vertexZ, fNBins_etaPhi_2,0.0,Double_t(fNBins_etaPhi_2),fNBins_pt_2,fMin_pt_2,fMax_pt_2);
    fhSum1Pt_2_vsZEtaPhiPt = new TH3F("sumPt1_2_vsZ_vsEtaPhi_vsPt","#LT #Sigma p_{t,2} #GT;vtx_{z};#eta_{2}#times#varphi_{2};p_{t,2} (GeV/c)",
        fNBins_vertexZ,fMin_vertexZ,fMax_vertexZ, fNBins_etaPhi_1,0.0,Double_t(fNBins_etaPhi_1),fNBins_pt_1,fMin_pt_1,fMax_pt_1);
    fhN1_2_vsEtaPhi = new TH2F("n1_2_vsEtaPhi","#LT n_{1} #GT;#eta_{2};#varphi_{2} (radian);#LT n_{1} #GT",
        fNBins_eta_2, fMin_eta_2, fMax_eta_2,  fNBins_phi_2, fMin_phi_2, fMax_phi_2);
    fhSum1Pt_2_vsEtaPhi = new TH2F("sumPt_2_vsEtaPhi","#LT #Sigma p_{t,2} #GT;#eta_{2};#varphi_{2} (radian);#LT #Sigma p_{t,2} #GT (GeV/c)",
        fNBins_eta_2,fMin_eta_2,fMax_eta_2,fNBins_phi_2,fMin_phi_2,fMax_phi_2);

    fOutput->Add(fhN1_1_vsPt);
    fOutput->Add(fhN1_1_vsZEtaPhiPt);
    fOutput->Add(fhSum1Pt_1_vsZEtaPhiPt);
    fOutput->Add(fhN1_2_vsPt);
    fOutput->Add(fhN1_2_vsZEtaPhiPt);
    fOutput->Add(fhSum1Pt_2_vsZEtaPhiPt);
    fOutput->Add(fhN1_1_vsEtaPhi);
    fOutput->Add(fhSum1Pt_1_vsEtaPhi);
    fOutput->Add(fhN1_2_vsEtaPhi);
    fOutput->Add(fhSum1Pt_2_vsEtaPhi);
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
  }
  TH1::AddDirectory(oldstatus);
}

/// \brief Stores the correction weights for both set of tracks
/// \param h3_1 the calibration weights for the first track
/// \param h3_2 the calibration weights for the second track
/// \return kTRUE if correctly done kFALSE otherwise
Bool_t AliTwoParticleCorrelationsBase::SetWeigths(const TH3F *h3_1, const TH3F *h3_2) {
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
Bool_t AliTwoParticleCorrelationsBase::SetEfficiencyCorrection(const TH1F *h1_1, const TH1F *h1_2) {
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
Bool_t AliTwoParticleCorrelationsBase::SetPairEfficiencyCorrection(const THn *h11, const THn *h12, const THn *h22, const THn *h21) {
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
Bool_t AliTwoParticleCorrelationsBase::SetSimultationPdfs(const TObjArray *pluspdf, const TObjArray *minuspdf) {

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
Bool_t AliTwoParticleCorrelationsBase::StartEvent(Float_t centrality, Float_t vertexz) {

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
Bool_t AliTwoParticleCorrelationsBase::ProcessTrack(Int_t trkId, AliVTrack *trk) {

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
Bool_t AliTwoParticleCorrelationsBase::ProcessTrack(Int_t trkId, AliVParticle *part) {

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


/// \brief Flag the potential products of conversions and / or resonances
void AliTwoParticleCorrelationsBase::FlagConversionsAndResonances() {
  AliInfo("");
  /* pair with different charges. Both track list are different */
  for (Int_t ix1 = 0; ix1 < fNoOfTracks1; ix1++) {
    for (Int_t ix2 = 0; ix2 < fNoOfTracks2; ix2++) {

      for (Int_t ires = 0; ires < fgkNoOfResonances; ires++) {
        if (fThresholdMult[ires] != 0) {
          Float_t mass = checkIfResonance(ires, kTRUE, ix1, ix2);

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

/// \cond CLASSIMP
ClassImp(AliTwoParticleCorrelationsBase);
/// \endcond