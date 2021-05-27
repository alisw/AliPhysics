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

/// Default constructor for object serialization
Ali2PCorrelations::Ali2PCorrelations() :
    AliTwoParticleCorrelationsBase(),
    fhPositivePtAverage(nullptr),
    fhNegativePtAverage(nullptr),
    fId_1(nullptr),
    fCharge_1(nullptr),
    fIxEta_1(nullptr),
    fIxPhi_1(nullptr),
    fIxPt_1(nullptr),
    fAvgPt_1(nullptr),
    fCorrection_1(nullptr),
    fId_2(nullptr),
    fCharge_2(nullptr),
    fIxEta_2(nullptr),
    fIxPhi_2(nullptr),
    fIxPt_2(nullptr),
    fAvgPt_2(nullptr),
    fCorrection_2(nullptr),
    /* pair bins */
    fNBins_deltaphi(72),
    fMin_deltaphi(0.0),
    fMax_deltaphi(TMath::Pi()*2.0),
    fWidth_deltaphi(TMath::Pi()*2.0/72.0),
    fNBins_deltaeta(20*2-1),
    fMin_deltaeta(-2.0),
    fMax_deltaeta(2.0),
    fWidth_deltaeta(4.0/(20*2-1)),
    /* histograms */
    fhN2_12_vsPtPt{nullptr,nullptr,nullptr,nullptr},
    fhN2_12_vsDEtaDPhi{nullptr,nullptr,nullptr,nullptr},
    fhSum2PtPt_12_vsDEtaDPhi{nullptr,nullptr,nullptr,nullptr},
    fhSum2DptDpt_12_vsDEtaDPhi{nullptr,nullptr,nullptr,nullptr},
    /* vs centrality profiles */
    fhN2_12_vsC{nullptr,nullptr,nullptr,nullptr},
    fhSum2PtPt_12_vsC{nullptr,nullptr,nullptr,nullptr},
    fhSum2DptDpt_12_vsC{nullptr,nullptr,nullptr,nullptr},
    fhN2nw_12_vsC{nullptr,nullptr,nullptr,nullptr},
    fhSum2PtPtnw_12_vsC{nullptr,nullptr,nullptr,nullptr},
    fhSum2DptDptnw_12_vsC{nullptr,nullptr,nullptr,nullptr}
{
}

/// Normal constructor
/// \param name the name for the object instance
Ali2PCorrelations::Ali2PCorrelations(const char *name) :
    AliTwoParticleCorrelationsBase(name),
    fhPositivePtAverage(nullptr),
    fhNegativePtAverage(nullptr),
    fId_1(nullptr),
    fCharge_1(nullptr),
    fIxEta_1(nullptr),
    fIxPhi_1(nullptr),
    fIxPt_1(nullptr),
    fAvgPt_1(nullptr),
    fCorrection_1(nullptr),
    fId_2(nullptr),
    fCharge_2(nullptr),
    fIxEta_2(nullptr),
    fIxPhi_2(nullptr),
    fIxPt_2(nullptr),
    fAvgPt_2(nullptr),
    fCorrection_2(nullptr),
    /* pair bins */
    fNBins_deltaphi(72),
    fMin_deltaphi(0.0),
    fMax_deltaphi(TMath::Pi()*2.0),
    fWidth_deltaphi(TMath::Pi()*2.0/72.0),
    fNBins_deltaeta(20*2-1),
    fMin_deltaeta(-2.0),
    fMax_deltaeta(2.0),
    fWidth_deltaeta(4.0/(20*2-1)),
    /* histograms */
    fhN2_12_vsPtPt{nullptr,nullptr,nullptr,nullptr},
    fhN2_12_vsDEtaDPhi{nullptr,nullptr,nullptr,nullptr},
    fhSum2PtPt_12_vsDEtaDPhi{nullptr,nullptr,nullptr,nullptr},
    fhSum2DptDpt_12_vsDEtaDPhi{nullptr,nullptr,nullptr,nullptr},
    fhN2_12_vsC{nullptr,nullptr,nullptr,nullptr},
    fhSum2PtPt_12_vsC{nullptr,nullptr,nullptr,nullptr},
    fhSum2DptDpt_12_vsC{nullptr,nullptr,nullptr,nullptr},
    fhN2nw_12_vsC{nullptr,nullptr,nullptr,nullptr},
    fhSum2PtPtnw_12_vsC{nullptr,nullptr,nullptr,nullptr},
    fhSum2DptDptnw_12_vsC{nullptr,nullptr,nullptr,nullptr}
{
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
  delete [] fAvgPt_1;
  delete [] fCorrection_1;

  /* track 2 storage */
  delete [] fId_2;
  delete [] fCharge_2;
  delete [] fIxEta_2;
  delete [] fIxPhi_2;
  delete [] fIxPt_2;
  delete [] fAvgPt_2;
  delete [] fCorrection_2;
}


/// \brief Establishes the binning configuration
/// \param confstring string containing the binning configuration parameters
Bool_t Ali2PCorrelations::ConfigureBinning(const char *confstring) {
  /* few sanity checks */
  TString str = confstring;
  if (!str.Contains("phishift"))
    return false;

  TObjArray *array = str.Tokenize(";");
  bool parentres = false;
  for (Int_t item = 0; item < array->GetEntries(); item++) {
    const TString &stritem = ((TObjString*) array->At(item))->GetString();
    if (stritem.BeginsWith("phishift:")) {
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
TString Ali2PCorrelations::GetBinningConfigurationString() const {
  TString parentcfg = AliTwoParticleCorrelationsBase::GetBinningConfigurationString();
  return TString::Format("Binning:phishift:%.1f;%s",fNBinsPhiShift,parentcfg.Data());
}

/// \brief Stores the average pt distribution in eta phi
/// \param h2_1 the avg pT for positive tracks
/// \param h2_2 the avg pT for negative tracks
bool Ali2PCorrelations::SetPtAvg(const TH2 *h2_1, const TH2 *h2_2) {
  if (h2_1 != nullptr && h2_2 != nullptr) {
    AliInfo("Stored transvers momentum average distributions");
    fhPositivePtAverage = h2_1;
    fhNegativePtAverage = h2_2;
    return true;
  }
  else {
    AliError("Transvers momentum average distributions empty!!! P2 will not work!!!");
    return false;
  }
}

/// \brief Initializes the member data structures
/// Allocates the needed memory an create the output histograms.
void Ali2PCorrelations::Initialize()
{
  Bool_t oldstatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  AliTwoParticleCorrelationsBase::Initialize();

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

  /* track 1 storage */
  fId_1 = new Int_t[fArraySize];
  fCharge_1 = new Int_t[fArraySize];
  fIxEta_1 = new Int_t[fArraySize];
  fIxPhi_1 = new Int_t[fArraySize];
  fIxPt_1 = new Int_t[fArraySize];
  fAvgPt_1 = new float[fArraySize];
  fCorrection_1 = new Float_t[fArraySize];

  /* track 2 storage */
  fId_2 = new Int_t[fArraySize];
  fCharge_2 = new Int_t[fArraySize];
  fIxEta_2 = new Int_t[fArraySize];
  fIxPhi_2 = new Int_t[fArraySize];
  fIxPt_2 = new Int_t[fArraySize];
  fAvgPt_2 = new float[fArraySize];
  fCorrection_2 = new Float_t[fArraySize];

  if (!fSinglesOnly)  {
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

/// \brief initializes the object instance to start a new event
/// \param centrality the new event centrality in percentage
/// \param vertexz the new event vertex \f$z\f$ coordinate
/// \return kTRUE if correctly initialized kFALSE otherwise
Bool_t Ali2PCorrelations::StartEvent(Float_t centrality, Float_t vertexz) {

  return AliTwoParticleCorrelationsBase::StartEvent(centrality, vertexz);
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
      fhSum1Pt_1_vsZEtaPhiPt->Fill(fVertexZ,ixEtaPhi+0.5,pT,corr*pT);
      fhN1_1_vsEtaPhi->Fill(eta,phi,corr);
      fhSum1Pt_1_vsEtaPhi->Fill(eta,phi,corr*pT);
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
      fAvgPt_1[fNoOfTracks1]        = (fhPositivePtAverage != nullptr) ? fhPositivePtAverage->GetBinContent(ixEta+1,ixPhi+1) : 0.0;
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
      fhSum1Pt_2_vsZEtaPhiPt->Fill(fVertexZ,ixEtaPhi+0.5,pT,corr*pT);
      fhN1_2_vsEtaPhi->Fill(eta,phi,corr);
      fhSum1Pt_2_vsEtaPhi->Fill(eta,phi,corr*pT);
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
      fAvgPt_2[fNoOfTracks2]        = (fhNegativePtAverage != nullptr) ? fhNegativePtAverage->GetBinContent(ixEta+1,ixPhi+1) : 0.0;
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
      float avgpt_1  = fAvgPt_1[ix1];

      for (Int_t ix2 = ix1+1; ix2 < fNoOfTracks1; ix2++) {
        /* excluded self correlations */
        float corr       = corr_1 * fCorrection_1[ix2];
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
        float dptdpt            = (pt_1-avgpt_1)*(fPt_1[ix2]-fAvgPt_1[ix2]);
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
      float avgpt_1  = fAvgPt_2[ix1];

      for (Int_t ix2 = ix1+1; ix2 < fNoOfTracks2; ix2++) {
        /* excluded self correlations */
        float corr      = corr_1 * fCorrection_2[ix2];
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
        float dptdpt            = (pt_1-avgpt_1)*(fPt_2[ix2]-fAvgPt_2[ix2]);
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
    float avgpt_1  = fAvgPt_1[ix1];

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
        float dptdpt            = (pt_1-avgpt_1)*(fPt_2[ix2]-fAvgPt_2[ix2]);
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
