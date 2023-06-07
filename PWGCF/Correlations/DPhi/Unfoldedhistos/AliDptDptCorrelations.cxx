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
#include "AliDptDptCorrelations.h"


/// Default constructor for object serialization
AliDptDptCorrelations::AliDptDptCorrelations()
  : AliTwoParticleCorrelationsBase(),
    fHalfSymmetrize(kTRUE),
    fSameSign(kFALSE),
    fAllCombinations(kFALSE),
    fRequestedCharge_1(1),
    fRequestedCharge_2(-1),
    /* the arrays with tracks 1 and 2 information */
    fIxEtaPhi(nullptr),
    fIxPt(nullptr),
    fCorrection(nullptr),
    fNBins_etaPhi_12(2073600),
    /* accumulated values */
    fN2_12(0.0),
    fNnw2_12(0.0),
    fSum2PtPt_12(0.0),
    fSum2PtPtnw_12(0.0),
    fSum2NPt_12(0.0),
    fSum2PtN_12(0.0),
    fSum2NPtnw_12(0.0),
    fSum2PtNnw_12(0.0),
    /* storage arrays */
    fN1_vsPt{},
    fN1_vsEtaPhi{},
    fSum1Pt_vsEtaPhi{},
    fN1_vsZEtaPhiPt{},
    fN2_12_vsPtPt(NULL),
    fN2_12_vsEtaPhi(NULL),
    fSum2PtPt_12_vsEtaPhi(NULL),
    fSum2PtN_12_vsEtaPhi(NULL),
    fSum2NPt_12_vsEtaPhi(NULL),
    /* histograms */
    fhN2_12_vsPtPt(NULL),
    fhN2_12_vsEtaPhi(NULL),
    fhSum2PtPt_12_vsEtaPhi(NULL),
    fhSum2PtN_12_vsEtaPhi(NULL),
    fhSum2NPt_12_vsEtaPhi(NULL),
    /* vs centrality profiles */
    fhN2_12_vsC(NULL),
    fhSum2PtPt_12_vsC(NULL),
    fhSum2PtN_12_vsC(NULL),
    fhSum2NPt_12_vsC(NULL),
    fhN2nw_12_vsC(NULL),
    fhSum2PtPtnw_12_vsC(NULL),
    fhSum2PtNnw_12_vsC(NULL),
    fhSum2NPtnw_12_vsC(NULL)
{
}

/// Normal constructor
/// \param name the name for the object instance
AliDptDptCorrelations::AliDptDptCorrelations(const char* name)
  : AliTwoParticleCorrelationsBase(name),
    fHalfSymmetrize(kTRUE),
    fSameSign(kFALSE),
    fAllCombinations(kFALSE),
    fRequestedCharge_1(1),
    fRequestedCharge_2(-1),
    /* the arrays with tracks 1 and 2 information */
    fIxEtaPhi(nullptr),
    fIxPt(nullptr),
    fCorrection(nullptr),
    fNBins_etaPhi_12(2073600),
    /* accumulated values */
    fN2_12(0.0),
    fNnw2_12(0.0),
    fSum2PtPt_12(0.0),
    fSum2PtPtnw_12(0.0),
    fSum2NPt_12(0.0),
    fSum2PtN_12(0.0),
    fSum2NPtnw_12(0.0),
    fSum2PtNnw_12(0.0),
    /* storage arrays */
    fN1_vsPt{},
    fN1_vsEtaPhi{},
    fSum1Pt_vsEtaPhi{},
    fN1_vsZEtaPhiPt{},
    fN2_12_vsPtPt(NULL),
    fN2_12_vsEtaPhi(NULL),
    fSum2PtPt_12_vsEtaPhi(NULL),
    fSum2PtN_12_vsEtaPhi(NULL),
    fSum2NPt_12_vsEtaPhi(NULL),
    /* histograms */
    fhN2_12_vsPtPt(NULL),
    fhN2_12_vsEtaPhi(NULL),
    fhSum2PtPt_12_vsEtaPhi(NULL),
    fhSum2PtN_12_vsEtaPhi(NULL),
    fhSum2NPt_12_vsEtaPhi(NULL),
    /* vs centrality profiles */
    fhN2_12_vsC(NULL),
    fhSum2PtPt_12_vsC(NULL),
    fhSum2PtN_12_vsC(NULL),
    fhSum2NPt_12_vsC(NULL),
    fhN2nw_12_vsC(NULL),
    fhSum2PtPtnw_12_vsC(NULL),
    fhSum2PtNnw_12_vsC(NULL),
    fhSum2NPtnw_12_vsC(NULL)
{
}

/// \brief Default destructor
/// Deallocates the allocated memory
AliDptDptCorrelations::~AliDptDptCorrelations() {
  for (int i = 0; i < 2; ++i) {
    /* for both h+ and h- tracks */
    if (fN1_vsPt[i] != nullptr)
      delete[] fN1_vsPt[i];
    if (fN1_vsEtaPhi[i] != nullptr)
      delete[] fN1_vsEtaPhi[i];
    if (fSum1Pt_vsEtaPhi[i] != nullptr)
      delete[] fSum1Pt_vsEtaPhi[i];
    if (fN1_vsZEtaPhiPt[i] != nullptr)
      delete[] fN1_vsZEtaPhiPt[i];
  }
  if (fN2_12_vsPtPt != NULL) delete [] fN2_12_vsPtPt;
  if (fN2_12_vsEtaPhi != NULL) delete [] fN2_12_vsEtaPhi;
  if (fSum2PtPt_12_vsEtaPhi != NULL) delete [] fSum2PtPt_12_vsEtaPhi;
  if (fSum2PtN_12_vsEtaPhi != NULL) delete [] fSum2PtN_12_vsEtaPhi;
  if (fSum2NPt_12_vsEtaPhi != NULL) delete [] fSum2NPt_12_vsEtaPhi;

  /* track storage */
  delete[] fIxEtaPhi;
  delete[] fIxPt;
  delete[] fCorrection;
}


/// \brief Establishes the binning configuration
/// \param confstring string containing the binning configuration parameters
Bool_t AliDptDptCorrelations::ConfigureBinning(const char *confstring) {
  /* few sanity checks */
  TString str = confstring;
  if (!str.Contains("halfsymm") || !str.Contains("phishift"))
    return false;

  TObjArray *array = str.Tokenize(";");
  bool parentres = false;
  for (Int_t item = 0; item < array->GetEntries(); item++) {
    const TString &stritem = ((TObjString*) array->At(item))->GetString();
    if (stritem.BeginsWith("halfsymm:")) {
      fHalfSymmetrize = stritem.Contains("yes");
    }
    else if (stritem.BeginsWith("phishift:")) {
      sscanf(stritem.Data(), "phishift:%lf", &fNBinsPhiShift);
    }
    else {
      parentres = AliTwoParticleCorrelationsBase::ConfigureBinning(stritem.Data());
    }
  }
  delete array;
  if (!parentres) return false;

  AliInfo("=====================================================");
  AliInfo(Form("Configured binning: %s", GetBinningConfigurationString().Data()));
  AliInfo("=====================================================");
  return true;
}

/// \brief Build the configuration string
/// \return the configuration string corresponding to the current configuration
TString AliDptDptCorrelations::GetBinningConfigurationString() const {
  TString parentcfg = AliTwoParticleCorrelationsBase::GetBinningConfigurationString();
  return TString::Format("Binning:halfsymm:%s;phishift:%.1f;%s",
      (fHalfSymmetrize ? "yes" : "not"),fNBinsPhiShift,parentcfg.Data());
}

/// \brief Initializes the member data structures
/// Allocates the needed memory an create the output histograms.
void AliDptDptCorrelations::Initialize()
{
  Bool_t oldstatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  AliTwoParticleCorrelationsBase::Initialize();

  /* incorporate configuration parameters */
  fOutput->Add(new TParameter<Bool_t>("HalfSymmetricResults",fHalfSymmetrize,'f'));

  /* track storage */
  fIxEtaPhi = new Int_t[fArraySize];
  fIxPt = new Int_t[fArraySize];
  fCorrection = new Float_t[fArraySize];
  for (int i = 0; i < 2; ++i) {
    /* for both h+ and h- tracks */
    fN1_vsPt.push_back(new Double_t[fNBins_pt]);
    for (Int_t j = 0; j < fNBins_pt; j++)
      fN1_vsPt[i][j] = 0.0;
    fN1_vsEtaPhi.push_back(new Double_t[fNBins_etaPhi]);
    for (Int_t j = 0; j < fNBins_etaPhi; j++)
      fN1_vsEtaPhi[i][j] = 0.0;
    fSum1Pt_vsEtaPhi.push_back(new Double_t[fNBins_etaPhi]);
    for (Int_t j = 0; j < fNBins_etaPhi; j++)
      fSum1Pt_vsEtaPhi[i][j] = 0.0;
    fN1_vsZEtaPhiPt.push_back(new Float_t[fNBins_zEtaPhiPt]);
    for (Int_t j = 0; j < fNBins_zEtaPhiPt; j++)
      fN1_vsZEtaPhiPt[i][j] = 0.0;
  }

  fN2_12_vsPtPt = new Double_t[fNBins_pt * fNBins_pt];
  for (Int_t i = 0; i < fNBins_pt * fNBins_pt; i++)
    fN2_12_vsPtPt[i] = 0.0;
  fN2_12_vsEtaPhi = new Float_t[fNBins_etaPhi_12]; for (Int_t i = 0; i < fNBins_etaPhi_12; i++) fN2_12_vsEtaPhi[i] = 0.0;
  fSum2PtPt_12_vsEtaPhi = new Float_t[fNBins_etaPhi_12]; for (Int_t i = 0; i < fNBins_etaPhi_12; i++) fSum2PtPt_12_vsEtaPhi[i] = 0.0;
  fSum2PtN_12_vsEtaPhi = new Float_t[fNBins_etaPhi_12]; for (Int_t i = 0; i < fNBins_etaPhi_12; i++) fSum2PtN_12_vsEtaPhi[i] = 0.0;
  fSum2NPt_12_vsEtaPhi = new Float_t[fNBins_etaPhi_12]; for (Int_t i = 0; i < fNBins_etaPhi_12; i++) fSum2NPt_12_vsEtaPhi[i] = 0.0;

  if (!fSinglesOnly)  {
    fhN2_12_vsEtaPhi = new TH1F("n2_12_vsEtaPhi","#LT n_{2} #GT;#eta_{1}#times#varphi_{1}#times#eta_{2}#times#varphi_{2};#LT n_{2} #GT",
        fNBins_etaPhi_12,0.,Double_t(fNBins_etaPhi_12));
    fhSum2PtPt_12_vsEtaPhi = new TH1F("sumPtPt_12_vsEtaPhi","#LT #Sigma p_{t,1}p_{t,2} #GT;#eta_{1}#times#varphi_{1}#times#eta_{2}#times#varphi_{2};#LT #Sigma p_{t,1}p_{t,2} #GT (GeV)^{2}",
        fNBins_etaPhi_12,0.,Double_t(fNBins_etaPhi_12));
    fhSum2PtN_12_vsEtaPhi = new TH1F("sumPtN_12_vsEtaPhi","#LT #Sigma p_{t,1}N #GT;#eta_{1}#times#varphi_{1}#times#eta_{2}#times#varphi_{2};#LT #Sigma p_{t,1}N #GT (GeV)",
        fNBins_etaPhi_12,0.,Double_t(fNBins_etaPhi_12));
    fhSum2NPt_12_vsEtaPhi = new TH1F("sumNPt_12_vsEtaPhi","#LT N#Sigma p_{t,2} #GT;#eta_{1}#times#varphi_{1}#times#eta_{2}#times#varphi_{2};#LT N#Sigma p_{t,2} #GT (GeV)",
        fNBins_etaPhi_12,0.,Double_t(fNBins_etaPhi_12));
    fhN2_12_vsPtPt = new TH2F("n2_12_vsPtVsPt", "#LT n_{2} #GT;p_{t,1} (GeV/c);p_{t,2} (GeV/c);#LT n_{2} #GT",
                              fNBins_pt, fMin_pt, fMax_pt, fNBins_pt, fMin_pt, fMax_pt);
    fhN2_12_vsC = new TProfile("n2_12_vsM","#LT n_{2} #GT (weighted);Centrality (%);#LT n_{2} #GT",100,0.0,100.0);
    fhSum2PtPt_12_vsC = new TProfile("sumPtPt_12_vsM","#LT #Sigma p_{t,1}p_{t,2} #GT (weighted);Centrality (%);#LT #Sigma p_{t,1}p_{t,2} #GT (GeV)^{2}",100,0.0,100.0);
    fhSum2PtN_12_vsC = new TProfile("sumPtN_12_vsM","#LT #Sigma p_{t,1}N #GT (weighted);Centrality (%);#LT #Sigma p_{t,1}N #GT (GeV)",100,0.0,100.0);
    fhSum2NPt_12_vsC = new TProfile("sumNPt_12_vsM","#LT N#Sigma p_{t,2} #GT (weighted);Centrality (%);#LT N#Sigma p_{t,2} #GT (GeV)",100,0.0,100.0);
    fhN2nw_12_vsC = new TProfile("n2Nw_12_vsM","#LT n_{2} #GT;Centrality (%);#LT n_{2} #GT",100,0.0,100.0);
    fhSum2PtPtnw_12_vsC = new TProfile("sumPtPtNw_12_vsM","#LT #Sigma p_{t,1}p_{t,2} #GT;Centrality (%);#LT #Sigma p_{t,1}p_{t,2} #GT (GeV)^{2}",100,0.0,100.0);
    fhSum2PtNnw_12_vsC = new TProfile("sumPtNNw_12_vsM","#LT #Sigma p_{t,1}N #GT;Centrality (%);#LT #Sigma p_{t,1}N #GT (GeV)",100,0.0,100.0);
    fhSum2NPtnw_12_vsC = new TProfile("sumNPtNw_12_vsM","#LT N#Sigma p_{t,2} #GT;Centrality (%);#LT N#Sigma p_{t,2} #GT (GeV)",100,0.0,100.0);

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

/// \brief initializes the object instance to start a new event
/// \param centrality the new event centrality in percentage
/// \param vertexz the new event vertex \f$z\f$ coordinate
/// \return kTRUE if correctly initialized kFALSE otherwise
Bool_t AliDptDptCorrelations::StartEvent(Float_t centrality, Float_t vertexz) {

  return AliTwoParticleCorrelationsBase::StartEvent(centrality, vertexz);
}

/// \brief process a track and store its parameters if feasible
/// \param trkId the external track Id
/// \param charge the track charge
/// \param pT the track \f$ p_T \f$
/// \param eta the track \f$ \eta \f$
/// \param phi the track \f$ \phi \f$
/// \return kTRUE if the track is properly handled kFALSE otherwise
bool AliDptDptCorrelations::ProcessTrack(int pid, float pT, float eta, float ophi)
{

  if (pid == 0) {
    /* we got a h+ track */
    if (!(fRequestedCharge_1 > 0 || fRequestedCharge_2 > 0)) {
      /* but not positive charge analysis was required */
      return false;
    }
  } else if (pid == 1) {
    /* we got a h- track */
    if (!(fRequestedCharge_1 < 0 || fRequestedCharge_2 < 0)) {
      /* but not negative charge analysis was required */
      return false;
    }
  } else {
    /* only prepared for hadron +/- correlations */
    return false;
  }

  if (pT > fMin_pt && pT < fMax_pt) {
    if (!(fNoOfTracks < fArraySize)) {
      AliError(Form("Storage for tracks full: %d", fArraySize));
      return kFALSE;
    }

    /* consider a potential phi origin shift */
    Float_t phi = ophi;
    if (!(phi < fMax_phi))
      phi = phi - 2 * TMath::Pi();
    Int_t ixPhi = Int_t((phi - fMin_phi) / fWidth_phi);
    if (ixPhi < 0 || !(ixPhi < fNBins_phi)) {
      AliWarning("Track one phi out of bins");
      return kFALSE;
    }

    if (eta < fMin_eta || fMax_eta < eta) {
      AliWarning(Form("Wrongly passed track one with eta: %.2f", eta));
      return kFALSE;
    }

    Int_t ixEta = Int_t((eta - fMin_eta) / fWidth_eta);
    if (ixEta < 0 || !(ixEta < fNBins_eta)) {
      AliWarning(Form("Track one eta bin %d out of bounds", ixEta));
      return kFALSE;
    }

    if (pT < fMin_pt || fMax_pt < pT) {
      AliWarning(Form("Wrongly passed track one with pT: %.2f", pT));
      return kFALSE;
    }

    Int_t ixPt = Int_t((pT - fMin_pt) / fWidth_pt);
    if (ixPt < 0 || !(ixPt < fNBins_pt)) {
      AliWarning(Form("Track pT bin %d out of bounds", ixPt));
      return kFALSE;
    }

    Int_t ixEtaPhi = ixEta * fNBins_phi + ixPhi;
    Int_t ixVertexP1 = fIxVertexZ * fNBins_etaPhiPt;
    Int_t ixZEtaPhiPt = ixVertexP1 + ixEtaPhi * fNBins_pt + ixPt;
    if (ixZEtaPhiPt < 0 || !(ixZEtaPhiPt < fNBins_zEtaPhiPt)) {
      AliWarning(Form("Event zvertex and track eta phi and pt bin %d out of bounds", ixZEtaPhiPt));
      return kFALSE;
    }

    Float_t effcorr = fEfficiencyCorrection[pid][ixPt];
    Float_t corr;
    if (fCorrectionWeights[pid] != nullptr)
      corr = fCorrectionWeights[pid][ixZEtaPhiPt];
    else
      corr = 1.0;

    /* the final correction incorporates also the efficiency correction */
    /* and affects to both, track densities and pT */
    corr *= effcorr;

    if (fSinglesOnly) {
      fN1_vsPt[pid][ixPt] += corr;
      fN1_vsZEtaPhiPt[pid][ixZEtaPhiPt] += corr;
    } else {
      float corrPt = corr * pT;
      fPID[fNoOfTracks] = pid;
      fIxEtaPhi[fNoOfTracks] = ixEtaPhi;
      fIxPt[fNoOfTracks] = ixPt;
      fFlags[fNoOfTracks] = 0x0;
      fPt[fNoOfTracks] = pT;
      fEta[fNoOfTracks] = eta;
      fPhi[fNoOfTracks] = ophi;
      fCorrection[fNoOfTracks] = corr;
      fN1[pid] += corr;
      fN1_vsEtaPhi[pid][ixEtaPhi] += corr;
      fSum1Pt[pid] += corrPt;
      fSum1Pt_vsEtaPhi[pid][ixEtaPhi] += corrPt;
      fNnw1[pid] += 1;
      fSum1Ptnw[pid] += pT;
      fNoOfTracks++;
      if (!(fNoOfTracks < fArraySize)) {
        AliError(Form("Storage for tracks full: %d", fArraySize));
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
        ProcessPairs<kLS, kHALF>();
      }
      else {
        ProcessPairs<kLS, kFULL>();
      }
    }
    else {
      if (fAllCombinations) {
        if (fHalfSymmetrize) {
          ProcessPairs<kALL, kHALF>();
        }
        else {
          ProcessPairs<kALL, kFULL>();
        }
      }
      else {
        if (fHalfSymmetrize) {
          ProcessPairs<kUS, kHALF>();
        }
        else {
          ProcessPairs<kUS, kFULL>();
        }
      }
    }
    /* now fill the profiles */
    for (int pid = 0; pid < 2; ++pid) {
      fhN1_vsC[pid]->Fill(fCentrality, fN1[pid]);
      fhSum1Pt_vsC[pid]->Fill(fCentrality, fSum1Pt[pid]);
      fhN1nw_vsC[pid]->Fill(fCentrality, fNnw1[pid]);
      fhSum1Ptnw_vsC[pid]->Fill(fCentrality, fSum1Ptnw[pid]);
    }
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
    return iy * nxbins - Int_t((iy - 1) * iy / 2) + (ix - iy);
  else
    return ix * nybins - Int_t((ix - 1) * ix / 2) + (iy - ix);
}

/// \brief Process track combinations with the same charge
template <AliDptDptCorrelations::kLSUS kind, AliDptDptCorrelations::kSYMM half>
void AliDptDptCorrelations::ProcessPairs()
{
  AliInfo("");

  if constexpr (kind != kLS) {
    if (fPostRejectResonances) {
      FlagConversionsAndResonances();
    }
  }

  for (int ix1 = 0; ix1 < fNoOfTracks; ++ix1) {
    if constexpr (kind == kLS) {
      /* for like sign all track should have the same charge/pid */
      if (fPID[0] != fPID[ix1]) {
        AliFatal("Something went really bad in here. ABORTING!!");
      }
    }

    for (Int_t ix2 = ix1 + 1; ix2 < fNoOfTracks; ix2++) {
      if constexpr (kind == kUS) {
        if (fPID[ix1] == fPID[ix2]) {
          /* same charge skeep it */
          continue;
        }
      }

      bool processpair = true;
      if constexpr (kind != kLS) {
        if (fPostRejectResonances) {
          /* process post particle selection resonance rejection */
          if (fPID[ix1] != fPID[ix2]) {
            /* process the resonance suppression for this pair if needed */
            for (Int_t ires = 0; ires < fgkNoOfResonances; ires++) {
              if (fThresholdMult[ires] != 0) {
                /* check if both tracks are flagged for the current resonance */
                if (((fFlags[ix1] & fFlags[ix2]) & UInt_t(0x1 << ires)) != UInt_t(0x1 << ires)) {
                  /* no, they are not, continue */
                  continue;
                } else {
                  /* yes, check if applicable */
                  Float_t mass = checkIfResonance(ires,
                                                  kFALSE,
                                                  fPt[ix1],
                                                  fEta[ix1],
                                                  fPhi[ix1],
                                                  fPt[ix2],
                                                  fEta[ix2],
                                                  fPhi[ix2]);

                  if (0 < mass) {
                    fhDiscardedResonanceMasses->Fill(ires, TMath::Sqrt(mass));
                    processpair = false;
                    break;
                  }
                }
              }
            }
          }
        }
      }

      if (processpair) {
        float corr = fCorrection[ix1] * fCorrection[ix2];

        /* apply the pair correction if applicable */
        const THn* effcorr = fPairsEfficiency[fPID[ix1]][fPID[ix2]];
        int bins[4];
        if (effcorr != NULL) {
          Int_t ieta1 = Int_t(fIxEtaPhi[ix1] / fNBins_phi);
          Int_t iphi1 = fIxEtaPhi[ix1] % fNBins_phi;
          Int_t ieta2 = Int_t(fIxEtaPhi[ix2] / fNBins_phi);
          Int_t iphi2 = fIxEtaPhi[ix2] % fNBins_phi;
          Int_t deltaetabin = ieta1 - ieta2 + fNBins_eta;
          Int_t ixdeltaphi = iphi1 - iphi2;
          if (ixdeltaphi < 0)
            ixdeltaphi += fNBins_phi;
          bins[0] = deltaetabin;
          bins[1] = ixdeltaphi + 1;
          bins[2] = fIxPt[ix1] + 1;
          bins[3] = fIxPt[ix2] + 1;
          corr = corr / effcorr->GetBinContent(bins);
        }

        if constexpr (half == kHALF) {
          int ij = getLinearSymmIndex(fIxEtaPhi[ix1], fNBins_etaPhi, fIxEtaPhi[ix2], fNBins_etaPhi);

          fN2_12 += 2 * corr;
          fN2_12_vsEtaPhi[ij] += corr;
          float ptpt = fPt[ix1] * fPt[ix2];
          fSum2PtPt_12 += 2 * corr * ptpt;
          fSum2PtN_12 += corr * fPt[ix1] + corr * fPt[ix2];
          fSum2NPt_12 += corr * fPt[ix2] + corr * fPt[ix1];
          fSum2PtPt_12_vsEtaPhi[ij] += corr * ptpt;
          fSum2PtN_12_vsEtaPhi[ij] += corr * fPt[ix1];
          fSum2NPt_12_vsEtaPhi[ij] += corr * fPt[ix2];
          fN2_12_vsPtPt[fIxPt[ix1] * fNBins_pt + fIxPt[ix2]] += corr;
          fN2_12_vsPtPt[fIxPt[ix2] * fNBins_pt + fIxPt[ix1]] += corr;

          fNnw2_12 += 2;
          fSum2PtPtnw_12 += 2 * ptpt;
          fSum2PtNnw_12 += fPt[ix1];
          fSum2PtNnw_12 += fPt[ix2];
          fSum2NPtnw_12 += fPt[ix2];
          fSum2NPtnw_12 += fPt[ix1];
        } else {
          int ij = fIxEtaPhi[ix1] * fNBins_etaPhi + fIxEtaPhi[ix2];

          fN2_12 += corr;
          fN2_12_vsEtaPhi[ij] += corr;
          Float_t ptpt = fPt[ix1] * fPt[ix2];
          fSum2PtPt_12 += corr * ptpt;
          fSum2PtN_12 += corr * fPt[ix1];
          fSum2NPt_12 += corr * fPt[ix2];
          fSum2PtPt_12_vsEtaPhi[ij] += corr * ptpt;
          fSum2PtN_12_vsEtaPhi[ij] += corr * fPt[ix1];
          fSum2NPt_12_vsEtaPhi[ij] += corr * fPt[ix2];
          fN2_12_vsPtPt[fIxPt[ix1] * fNBins_pt + fIxPt[ix2]] += corr;

          fNnw2_12 += 1;
          fSum2PtPtnw_12 += ptpt;
          fSum2PtNnw_12 += fPt[ix1];
          fSum2NPtnw_12 += fPt[ix2];
        }
      }
    } // ix2
  }   // ix1
}

/// \brief Fill the histograms once the whole process is finished
void  AliDptDptCorrelations::FinalizeProcess()
{

  if (fSinglesOnly) {
    for (int pid = 0; pid < 2; ++pid) {
      fillHistoWithArray(fhN1_vsPt[pid], fN1_vsPt[pid], fNBins_pt);
      fillHistoWithArray(fhN1_vsZEtaPhiPt[pid], fN1_vsZEtaPhiPt[pid], fNBins_vertexZ, fNBins_etaPhi, fNBins_pt);

      /* for the time being, the errors are trivial so, we do not use results file space */
      fhN1_vsPt[pid]->Sumw2(false);
      fhN1_vsZEtaPhiPt[pid]->Sumw2(false);
    }
  } else {
    for (int pid = 0; pid < 2; ++pid) {
      fillHistoWithArray(fhN1_vsEtaPhi[pid], fN1_vsEtaPhi[pid], fNBins_eta, fNBins_phi);
      fillHistoWithArray(fhSum1Pt_vsEtaPhi[pid], fSum1Pt_vsEtaPhi[pid], fNBins_eta, fNBins_phi);
      fhN1_vsEtaPhi[pid]->Sumw2(false);
      fhSum1Pt_vsEtaPhi[pid]->Sumw2(false);
    }

    fillHistoWithArray(fhN2_12_vsEtaPhi, fN2_12_vsEtaPhi, fNBins_etaPhi_12);
    fillHistoWithArray(fhSum2PtPt_12_vsEtaPhi, fSum2PtPt_12_vsEtaPhi, fNBins_etaPhi_12);
    fillHistoWithArray(fhSum2PtN_12_vsEtaPhi, fSum2PtN_12_vsEtaPhi, fNBins_etaPhi_12);
    fillHistoWithArray(fhSum2NPt_12_vsEtaPhi, fSum2NPt_12_vsEtaPhi, fNBins_etaPhi_12);
    fillHistoWithArray(fhN2_12_vsPtPt, fN2_12_vsPtPt, fNBins_pt, fNBins_pt);

    /* for the time being, the errors are trivial so, we do not use results file space */
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
