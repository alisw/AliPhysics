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
AliDptDptCorrelations::AliDptDptCorrelations() :
    AliTwoParticleCorrelationsBase(),
    fHalfSymmetrize(kTRUE),
    fSameSign(kFALSE),
    fAllCombinations(kFALSE),
    fRequestedCharge_1(1),
    fRequestedCharge_2(-1),
    /* the arrays with tracks 1 and 2 information */
    fId_1(NULL),
    fCharge_1(NULL),
    fIxEtaPhi_1(NULL),
    fIxPt_1(NULL),
    fCorrection_1(NULL),
    fId_2(NULL),
    fCharge_2(NULL),
    fIxEtaPhi_2(NULL),
    fIxPt_2(NULL),
    fCorrection_2(NULL),
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
AliDptDptCorrelations::AliDptDptCorrelations(const char *name) :
    AliTwoParticleCorrelationsBase(name),
    fHalfSymmetrize(kTRUE),
    fSameSign(kFALSE),
    fAllCombinations(kFALSE),
    fRequestedCharge_1(1),
    fRequestedCharge_2(-1),
    /* the arrays with tracks 1 and 2 information */
    fId_1(NULL),
    fCharge_1(NULL),
    fIxEtaPhi_1(NULL),
    fIxPt_1(NULL),
    fCorrection_1(NULL),
    fId_2(NULL),
    fCharge_2(NULL),
    fIxEtaPhi_2(NULL),
    fIxPt_2(NULL),
    fCorrection_2(NULL),
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

  /* track 1 storage */
  delete [] fId_1;
  delete [] fCharge_1;
  delete [] fIxEtaPhi_1;
  delete [] fIxPt_1;
  delete [] fCorrection_1;

  /* track 2 storage */
  delete [] fId_2;
  delete [] fCharge_2;
  delete [] fIxEtaPhi_2;
  delete [] fIxPt_2;
  delete [] fCorrection_2;
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

  /* track 1 storage */
  fId_1 = new Int_t[fArraySize];
  fCharge_1 = new Int_t[fArraySize];
  fIxEtaPhi_1 = new Int_t[fArraySize];
  fIxPt_1 = new Int_t[fArraySize];
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

  if (!fSinglesOnly)  {
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
            Float_t mass = checkIfResonance(ires, kFALSE, ix1, ix2);

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
            Float_t mass = checkIfResonance(ires, kFALSE, ix1, ix2);

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
