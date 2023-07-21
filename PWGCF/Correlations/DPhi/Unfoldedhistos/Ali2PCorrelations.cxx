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

/// Default constructor for object serialization
Ali2PCorrelations::Ali2PCorrelations()
  : AliTwoParticleCorrelationsBase(),
    fhPtAverage{},
    fIxEta(nullptr),
    fIxPhi(nullptr),
    fIxPt(nullptr),
    fAvgPt(nullptr),
    fCorrection(nullptr),
    /* pair bins */
    fNBins_deltaphi(72),
    fMin_deltaphi(0.0),
    fMax_deltaphi(TMath::Pi() * 2.0),
    fWidth_deltaphi(TMath::Pi() * 2.0 / 72.0),
    fNBins_deltaeta(20 * 2 - 1),
    fMin_deltaeta(-2.0),
    fMax_deltaeta(2.0),
    fWidth_deltaeta(4.0 / (20 * 2 - 1)),
    /* histograms */
    fhN2_12_vsPtPt{},
    fhN2_12_vsDEtaDPhi{},
    fhN2_12_vsDEtaDPhi_na{},
    fhSum2PtPt_12_vsDEtaDPhi{},
    fhSum2DptDpt_12_vsDEtaDPhi{},
    fhSum2DptDpt_12_vsDEtaDPhi_na{},
    /* vs centrality profiles */
    fhN2_12_vsC{},
    fhSum2PtPt_12_vsC{},
    fhSum2DptDpt_12_vsC{},
    fhN2nw_12_vsC{},
    fhSum2PtPtnw_12_vsC{},
    fhSum2DptDptnw_12_vsC{}
{
}

/// Normal constructor
/// \param name the name for the object instance
Ali2PCorrelations::Ali2PCorrelations(const char* name)
  : AliTwoParticleCorrelationsBase(name),
    fhPtAverage{},
    fIxEta(nullptr),
    fIxPhi(nullptr),
    fIxPt(nullptr),
    fAvgPt(nullptr),
    fCorrection(nullptr),
    /* pair bins */
    fNBins_deltaphi(72),
    fMin_deltaphi(0.0),
    fMax_deltaphi(TMath::Pi() * 2.0),
    fWidth_deltaphi(TMath::Pi() * 2.0 / 72.0),
    fNBins_deltaeta(20 * 2 - 1),
    fMin_deltaeta(-2.0),
    fMax_deltaeta(2.0),
    fWidth_deltaeta(4.0 / (20 * 2 - 1)),
    /* histograms */
    fhN2_12_vsPtPt{},
    fhN2_12_vsDEtaDPhi{},
    fhN2_12_vsDEtaDPhi_na{},
    fhSum2PtPt_12_vsDEtaDPhi{},
    fhSum2DptDpt_12_vsDEtaDPhi{},
    fhSum2DptDpt_12_vsDEtaDPhi_na{},
    /* vs centrality profiles */
    fhN2_12_vsC{},
    fhSum2PtPt_12_vsC{},
    fhSum2DptDpt_12_vsC{},
    fhN2nw_12_vsC{},
    fhSum2PtPtnw_12_vsC{},
    fhSum2DptDptnw_12_vsC{}
{
}

/// \brief Default destructor
/// Deallocates the allocated memory
Ali2PCorrelations::~Ali2PCorrelations() {
  /* track storage */
  delete[] fIxEta;
  delete[] fIxPhi;
  delete[] fIxPt;
  delete[] fAvgPt;
  delete[] fCorrection;
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
/// \param h2 the avg pT for the different species
bool Ali2PCorrelations::SetPtAvg(std::vector<const TH2*> h2)
{
  if (h2.size() > 0) {
    AliInfo("Stored transverse momentum average distributions");
    for (uint i = 0; i < h2.size(); ++i) {
      fhPtAverage.push_back(h2[i]);
    }
    return true;
  } else {
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

  fNBins_deltaeta = fNBins_eta * 2 - 1;
  fMin_deltaeta = fMin_eta - fMax_eta;
  fMax_deltaeta = fMax_eta - fMin_eta;
  fWidth_deltaeta    = (fMax_deltaeta - fMin_deltaeta) / fNBins_deltaeta;
  fNBins_deltaphi = fNBins_phi;
  fWidth_deltaphi    = TMath::TwoPi() / fNBins_deltaphi;
  fMin_deltaphi      = 0.0 - fWidth_deltaphi / 2.0;
  fMax_deltaphi      = TMath::TwoPi() - fWidth_deltaphi / 2.0;

  /* incorporate configuration parameters */
  fOutput->Add(new TParameter<Bool_t>("DifferentialOutput",true,'f'));

  /* track storage */
  fIxEta = new Int_t[fArraySize];
  fIxPhi = new Int_t[fArraySize];
  fIxPt = new Int_t[fArraySize];
  fAvgPt = new float[fArraySize];
  fCorrection = new Float_t[fArraySize];

  if (!fSinglesOnly) {
    /* initialize the accumlators */
    fN2_12 = std::vector<std::vector<double>>(fSpeciesNames.size(), std::vector<double>(fSpeciesNames.size(), 0.0));
    fSum2PtPt_12 = std::vector<std::vector<double>>(fSpeciesNames.size(), std::vector<double>(fSpeciesNames.size(), 0.0));
    fSum2DptDpt_12 = std::vector<std::vector<double>>(fSpeciesNames.size(), std::vector<double>(fSpeciesNames.size(), 0.0));
    fNnw2_12 = std::vector<std::vector<double>>(fSpeciesNames.size(), std::vector<double>(fSpeciesNames.size(), 0.0));
    fSum2PtPtnw_12 = std::vector<std::vector<double>>(fSpeciesNames.size(), std::vector<double>(fSpeciesNames.size(), 0.0));
    fSum2DptDptnw_12 = std::vector<std::vector<double>>(fSpeciesNames.size(), std::vector<double>(fSpeciesNames.size(), 0.0));

    /* histograms for track pairs */
    for (uint pidx = 0; pidx < fSpeciesNames.size(); ++pidx) {
      fhN2_12_vsDEtaDPhi.resize(fSpeciesNames.size());
      fhN2_12_vsDEtaDPhi_na.resize(fSpeciesNames.size());
      fhSum2PtPt_12_vsDEtaDPhi.resize(fSpeciesNames.size());
      fhSum2DptDpt_12_vsDEtaDPhi.resize(fSpeciesNames.size());
      fhSum2DptDpt_12_vsDEtaDPhi_na.resize(fSpeciesNames.size());
      fhN2_12_vsPtPt.resize(fSpeciesNames.size());
      fhN2_12_vsC.resize(fSpeciesNames.size());
      fhSum2PtPt_12_vsC.resize(fSpeciesNames.size());
      fhSum2DptDpt_12_vsC.resize(fSpeciesNames.size());
      fhN2nw_12_vsC.resize(fSpeciesNames.size());
      fhSum2PtPtnw_12_vsC.resize(fSpeciesNames.size());
      fhSum2DptDptnw_12_vsC.resize(fSpeciesNames.size());
      for (uint pidy = 0; pidy < fSpeciesNames.size(); ++pidy) {
        char pname[128];
        sprintf(pname, "%s%s", fSpeciesNames[pidx].c_str(), fSpeciesNames[pidy].c_str());
        fhN2_12_vsDEtaDPhi[pidx].push_back(new TH2F(TString::Format("n2_12_vsDEtaDPhi_%s", pname), TString::Format("#LT n_{2} #GT (%s);#Delta#eta;#Delta#varphi;#LT n_{2} #GT", pname),
                                                    fNBins_deltaeta, fMin_deltaeta, fMax_deltaeta, fNBins_deltaphi, fMin_deltaphi, fMax_deltaphi));
        fhN2_12_vsDEtaDPhi_na[pidx].push_back(new TH2F(TString::Format("n2_12_vsDEtaDPhiNa_%s", pname), TString::Format("#LT n_{2} #GT na (%s);#Delta#eta;#Delta#varphi;#LT n_{2} #GT", pname),
                                                       fNBins_deltaeta, fMin_deltaeta, fMax_deltaeta, fNBins_deltaphi, fMin_deltaphi, fMax_deltaphi));
        fhSum2PtPt_12_vsDEtaDPhi[pidx].push_back(new TH2F(TString::Format("sumPtPt_12_vsDEtaDPhi_%s", pname), TString::Format("#LT #Sigma p_{t,1}p_{t,2} #GT (%s);#Delta#eta;#Delta#varphi;#LT #Sigma p_{t,1}p_{t,2} #GT (GeV^{2})", pname),
                                                          fNBins_deltaeta, fMin_deltaeta, fMax_deltaeta, fNBins_deltaphi, fMin_deltaphi, fMax_deltaphi));
        fhSum2DptDpt_12_vsDEtaDPhi[pidx].push_back(new TH2F(TString::Format("sumDptDpt_12_vsDEtaDPhi_%s", pname), TString::Format("#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (%s);#Delta#eta;#Delta#varphi;#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (GeV^{2})", pname),
                                                            fNBins_deltaeta, fMin_deltaeta, fMax_deltaeta, fNBins_deltaphi, fMin_deltaphi, fMax_deltaphi));
        fhSum2DptDpt_12_vsDEtaDPhi_na[pidx].push_back(new TH2F(TString::Format("sumDptDpt_12_vsDEtaDPhiNa_%s", pname), TString::Format("#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT na (%s);#Delta#eta;#Delta#varphi;#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (GeV^{2})", pname),
                                                               fNBins_deltaeta, fMin_deltaeta, fMax_deltaeta, fNBins_deltaphi, fMin_deltaphi, fMax_deltaphi));
        fhN2_12_vsPtPt[pidx].push_back(new TH2F(TString::Format("n2_12_vsPtVsPt_%s", pname), TString::Format("#LT n_{2} #GT (%s);p_{t,1} (GeV/c);p_{t,2} (GeV/c);#LT n_{2} #GT", pname),
                                                fNBins_pt, fMin_pt, fMax_pt, fNBins_pt, fMin_pt, fMax_pt));
        fhN2_12_vsC[pidx].push_back(new TProfile(TString::Format("n2_12_vsM_%s", pname), TString::Format("#LT n_{2} #GT (%s) (weighted);Centrality (%%);#LT n_{2} #GT", pname), 100, 0.0, 100.0));
        fhSum2PtPt_12_vsC[pidx].push_back(new TProfile(TString::Format("sumPtPt_12_vsM_%s", pname), TString::Format("#LT #Sigma p_{t,1}p_{t,2} #GT (%s) (weighted);Centrality (%%);#LT #Sigma p_{t,1}p_{t,2} #GT (GeV^{2})", pname), 100, 0.0, 100.0));
        fhSum2DptDpt_12_vsC[pidx].push_back(new TProfile(TString::Format("sumDptDpt_12_vsM_%s", pname), TString::Format("#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (%s) (weighted);Centrality (%%);#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (GeV^{2})", pname), 100, 0.0, 100.0));
        fhN2nw_12_vsC[pidx].push_back(new TProfile(TString::Format("n2Nw_12_vsM_%s", pname), TString::Format("#LT n_{2} #GT (%s);Centrality (%%);#LT n_{2} #GT", pname), 100, 0.0, 100.0));
        fhSum2PtPtnw_12_vsC[pidx].push_back(new TProfile(TString::Format("sumPtPtNw_12_vsM_%s", pname), TString::Format("#LT #Sigma p_{t,1}p_{t,2} #GT (%s);Centrality (%%);#LT #Sigma p_{t,1}p_{t,2} #GT (GeV^{2})", pname), 100, 0.0, 100.0));
        fhSum2DptDptnw_12_vsC[pidx].push_back(new TProfile(TString::Format("sumDptDptNw_12_vsM_%s", pname), TString::Format("#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (%s);Centrality (%%);#LT #Sigma (p_{t,1} - #LT p_{t,1} #GT)(p_{t,2} - #LT p_{t,2} #GT) #GT (GeV^{2})", pname), 100, 0.0, 100.0));

        fOutput->Add(fhN2_12_vsDEtaDPhi[pidx][pidy]);
        fOutput->Add(fhN2_12_vsDEtaDPhi_na[pidx][pidy]);
        fOutput->Add(fhSum2PtPt_12_vsDEtaDPhi[pidx][pidy]);
        fOutput->Add(fhSum2DptDpt_12_vsDEtaDPhi[pidx][pidy]);
        fOutput->Add(fhSum2DptDpt_12_vsDEtaDPhi_na[pidx][pidy]);
        fOutput->Add(fhN2_12_vsPtPt[pidx][pidy]);
        fOutput->Add(fhN2_12_vsC[pidx][pidy]);
        fOutput->Add(fhSum2PtPt_12_vsC[pidx][pidy]);
        fOutput->Add(fhSum2DptDpt_12_vsC[pidx][pidy]);
        fOutput->Add(fhN2nw_12_vsC[pidx][pidy]);
        fOutput->Add(fhSum2PtPtnw_12_vsC[pidx][pidy]);
        fOutput->Add(fhSum2DptDptnw_12_vsC[pidx][pidy]);
      }
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
/// \param pid the track PID
/// \param pT the track \f$ p_T \f$
/// \param eta the track \f$ \eta \f$
/// \param phi the track \f$ \phi \f$
/// \return kTRUE if the track is properly handled kFALSE otherwise
/// For the time being positive tracks go to bank one and negative tracks go to bank two
bool Ali2PCorrelations::ProcessTrack(int pid, float pT, float eta, float ophi)
{
  if (!(fNoOfTracks < fArraySize)) {
    AliError(Form("Storage for track one full: %d", fArraySize));
    return kFALSE;
  }

  /* consider a potential phi origin shift */
  float phi = ophi;
  if (!(phi < fMax_phi))
    phi = phi - 2 * TMath::Pi();
  if (phi < fMin_phi) {
    return kFALSE;
  }
  float ixPhi = int((phi - fMin_phi) / fWidth_phi);
  if (ixPhi < 0 || !(ixPhi < fNBins_phi)) {
    AliWarning("Track one phi out of bins");
    return kFALSE;
  }

  if (eta < fMin_eta || fMax_eta < eta) {
    AliWarning(Form("Wrongly passed track one with eta: %.2f", eta));
    return kFALSE;
  }

  int ixEta = int((eta - fMin_eta) / fWidth_eta);
  if (ixEta < 0 || !(ixEta < fNBins_eta)) {
    AliWarning(Form("Track one eta bin %d out of bounds", ixEta));
    return kFALSE;
  }

  if (pT < fMin_pt || fMax_pt < pT) {
    AliWarning(Form("Wrongly passed track one with pT: %.2f", pT));
    return kFALSE;
  }

  int ixPt = int((pT - fMin_pt) / fWidth_pt);
  if (ixPt < 0 || !(ixPt < fNBins_pt)) {
    AliWarning(Form("Track pT bin %d out of bounds",ixPt));
    return kFALSE;
  }

  int ixEtaPhi = ixEta * fNBins_phi + ixPhi;
  int ixVertexP1 = fIxVertexZ * fNBins_etaPhiPt;
  int ixZEtaPhiPt = ixVertexP1 + ixEtaPhi * fNBins_pt + ixPt;
  if (ixZEtaPhiPt < 0 || !(ixZEtaPhiPt < fNBins_zEtaPhiPt)) {
    AliWarning(Form("Event zvertex and track eta phi and pt bin %d out of bounds", ixZEtaPhiPt));
    return kFALSE;
  }

  float effcorr = fEfficiencyCorrection[pid][ixPt];
  float corr;
  /* for the time being the weights are only applied at the hadron level */
  if (fCorrectionWeights[pid % 2] != nullptr)
    corr = fCorrectionWeights[pid % 2][ixZEtaPhiPt];
  else
    corr = 1.0;

  /* the final correction incorporates also the efficiency correction */
  /* and affects to both, track densities and pT */
  corr *= effcorr;

  if (fSinglesOnly) {
    auto fillsingles =[&](int ix) {
      fhN1_vsPt[ix]->Fill(pT, corr);
      fhN1_vsZEtaPhiPt[ix]->Fill(fVertexZ, ixEtaPhi + 0.5, pT, corr);
      fhSum1Pt_vsZEtaPhiPt[ix]->Fill(fVertexZ, ixEtaPhi + 0.5, pT, corr * pT);
      fhN1_vsEtaPhi[ix]->Fill(eta, phi, corr);
      fhSum1Pt_vsEtaPhi[ix]->Fill(eta, phi, corr * pT);
    };
    fillsingles(pid);
    fillsingles(pid % 2); /* always the whole charged hadrons as well */
  }
  else {
    float corrPt = corr * pT;
    auto fillsingles =[&](int ix) {
      fN1[ix] += corr;
      fhN1_vsEtaPhi[ix]->Fill(eta, phi, corr);
      fSum1Pt[ix] += corrPt;
      fhSum1Pt_vsEtaPhi[ix]->Fill(eta, phi, corrPt);
      fNnw1[ix] += 1;
      fSum1Ptnw[ix] += pT;
    };
    fPID[fNoOfTracks] = pid;
    fIxEta[fNoOfTracks] = ixEta;
    fIxPhi[fNoOfTracks] = ixPhi;
    fIxPt[fNoOfTracks] = ixPt;
    fFlags[fNoOfTracks] = 0x0;
    fPt[fNoOfTracks] = pT;
    fEta[fNoOfTracks] = eta;
    fPhi[fNoOfTracks] = phi;
    fAvgPt[fNoOfTracks] = (fhPtAverage[pid] != nullptr) ? fhPtAverage[pid]->GetBinContent(ixEta + 1, ixPhi + 1) : 0.0;
    fCorrection[fNoOfTracks] = corr;
    fillsingles(pid);
    fillsingles(pid % 2); /* always the whole charged hadrons as well */
    fNoOfTracks++;
    if (!(fNoOfTracks < fArraySize)) {
      AliError(Form("Storage for track one full: %d", fArraySize));
      return kFALSE;
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
    if (fSpeciesNames.size() > 2) {
      /* identified, incorporate charged hadrons */
      ProcessPairs<false>();
    } else {
      /* not identified, just charged hadrons */
      ProcessPairs<true>();
    }

    /* finally fill the individual tracks profiles */
    for (uint pid = 0; pid < fSpeciesNames.size(); ++pid) {
      fhN1_vsC[pid]->Fill(fCentrality, fN1[pid]);
      fhSum1Pt_vsC[pid]->Fill(fCentrality, fSum1Pt[pid]);
      fhN1nw_vsC[pid]->Fill(fCentrality, fNnw1[pid]);
      fhSum1Ptnw_vsC[pid]->Fill(fCentrality, fSum1Ptnw[pid]);
    }
  }
}

/// \brief Process track combinations with the same charge
/// \param bank the tracks bank to use
template<bool chargedhadrons>
void Ali2PCorrelations::ProcessPairs()
{
  /* flag resonances candidates if needed */
  if (fPostRejectResonances) {
    FlagConversionsAndResonances();
  }

  /* reset pair counters */
  for (unsigned int icx = 0; icx < fSpeciesNames.size(); ++icx) {
    for (unsigned int icy = 0; icy < fSpeciesNames.size(); ++icy) {
      fN2_12[icx][icy] = 0;
      fSum2PtPt_12[icx][icy] = 0;
      fSum2DptDpt_12[icx][icy] = 0;
      fNnw2_12[icx][icy] = 0;
      fSum2PtPtnw_12[icx][icy] = 0;
      fSum2DptDptnw_12[icx][icy] = 0;
    }
  }

  for (int ix1 = 0; ix1 < fNoOfTracks; ++ix1) {
    for (int ix2 = ix1 + 1; ix2 < fNoOfTracks; ++ix2) {
      /* self correlations already excluded */

      bool processpair = kTRUE;
      /* TODO: probably a table with the crossed PID which could give a certain resonance will help */
      if (fPostRejectResonances) {
        /* process post particle selection resonance rejection */
        if (fPID[ix1] != fPID[ix2]) {
          /* for different species/charges process the resonance suppression for this pair if needed */
          for (int ires = 0; ires < fgkNoOfResonances; ires++) {
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
                  processpair = kFALSE;
                  break;
                }
              }
            }
          }
        }
      }
      if (processpair) {
        float corr = fCorrection[ix1] * fCorrection[ix2];
        int ixDeltaEta_d = fIxEta[ix1] - fIxEta[ix2] + fNBins_eta - 1;
        int ixDeltaPhi_d = fIxPhi[ix1] - fIxPhi[ix2];
        if (ixDeltaPhi_d < 0)
          ixDeltaPhi_d += fNBins_phi;
        int ixDeltaEta_c = fIxEta[ix2] - fIxEta[ix1] + fNBins_eta - 1;
        int ixDeltaPhi_c = fIxPhi[ix2] - fIxPhi[ix1];
        if (ixDeltaPhi_c < 0)
          ixDeltaPhi_c += fNBins_phi;
        float deltaEta_d = fEta[ix1] - fEta[ix2];
        float deltaEta_c = -deltaEta_d;
        float deltaPhi_d = float(TVector2::Phi_0_2pi(fPhi[ix1] - fPhi[ix2]));
        if (not(deltaPhi_d < fMax_deltaphi))
          deltaPhi_d -= TMath::TwoPi();
        float deltaPhi_c = float(TVector2::Phi_0_2pi(fPhi[ix2] - fPhi[ix1]));
        if (not(deltaPhi_c < fMax_deltaphi))
          deltaPhi_c -= TMath::TwoPi();

        /* apply the pair correction if applicable */
        const THn* effcorr = fPairsEfficiency[fPID[ix1]][fPID[ix2]];
        int bins[4];
        if (effcorr != NULL) {
          bins[0] = ixDeltaEta_d + 1;
          bins[1] = ixDeltaPhi_d + 1;
          bins[2] = fIxPt[ix1] + 1;
          bins[3] = fIxPt[ix2] + 1;
          corr = corr / effcorr->GetBinContent(bins);
        }
        int globalbin_d = fhN2_12_vsDEtaDPhi[fPID[ix1]][fPID[ix2]]->GetBin(ixDeltaEta_d + 1, ixDeltaPhi_d + 1);
        int globalbin_c = fhN2_12_vsDEtaDPhi[fPID[ix2]][fPID[ix1]]->GetBin(ixDeltaEta_c + 1, ixDeltaPhi_c + 1);
        int globalbin_na_d = fhN2_12_vsDEtaDPhi_na[fPID[ix1]][fPID[ix2]]->FindBin(deltaEta_d, deltaPhi_d);
        int globalbin_na_c = fhN2_12_vsDEtaDPhi_na[fPID[ix2]][fPID[ix1]]->FindBin(deltaEta_c, deltaPhi_c);

        auto fillpairs = [&](int pidx1,int pidx2) {
          fNnw2_12[pidx1][pidx2] += 1;
          fNnw2_12[pidx2][pidx1] += 1;
          fN2_12[pidx1][pidx2] += corr;
          fN2_12[pidx2][pidx1] += corr;
          fhN2_12_vsDEtaDPhi[pidx1][pidx2]->AddBinContent(globalbin_d, corr);
          fhN2_12_vsDEtaDPhi[pidx2][pidx1]->AddBinContent(globalbin_c, corr);
          fhN2_12_vsDEtaDPhi_na[pidx1][pidx2]->AddBinContent(globalbin_na_d, corr);
          fhN2_12_vsDEtaDPhi_na[pidx2][pidx1]->AddBinContent(globalbin_na_c, corr);
          float ptpt = fPt[ix1] * fPt[ix2];
          fSum2PtPtnw_12[pidx1][pidx2] += ptpt;
          fSum2PtPtnw_12[pidx2][pidx1] += ptpt;
          fSum2PtPt_12[pidx1][pidx2] += corr * ptpt;
          fSum2PtPt_12[pidx2][pidx1] += corr * ptpt;
          fhSum2PtPt_12_vsDEtaDPhi[pidx1][pidx2]->AddBinContent(globalbin_d, corr * ptpt);
          fhSum2PtPt_12_vsDEtaDPhi[pidx2][pidx1]->AddBinContent(globalbin_c, corr * ptpt);
          float dptdptnw = (fPt[ix1] - fAvgPt[ix1]) * (fPt[ix2] - fAvgPt[ix2]);
          float dptdpt = corr * (fPt[ix1] - fAvgPt[ix1]) * (fPt[ix2] - fAvgPt[ix2]);
          fSum2DptDptnw_12[pidx1][pidx2] += dptdptnw;
          fSum2DptDptnw_12[pidx2][pidx1] += dptdptnw;
          fSum2DptDpt_12[pidx1][pidx2] += dptdpt;
          fSum2DptDpt_12[pidx2][pidx1] += dptdpt;
          fhSum2DptDpt_12_vsDEtaDPhi[pidx1][pidx2]->AddBinContent(globalbin_d, dptdpt);
          fhSum2DptDpt_12_vsDEtaDPhi[pidx2][pidx1]->AddBinContent(globalbin_c, dptdpt);
          fhSum2DptDpt_12_vsDEtaDPhi_na[pidx1][pidx2]->AddBinContent(globalbin_na_d, dptdpt);
          fhSum2DptDpt_12_vsDEtaDPhi_na[pidx2][pidx1]->AddBinContent(globalbin_na_c, dptdpt);
          fhN2_12_vsPtPt[pidx1][pidx2]->Fill(fPt[ix1], fPt[ix2], corr);
          fhN2_12_vsPtPt[pidx2][pidx1]->Fill(fPt[ix2], fPt[ix1], corr);
        };
        fillpairs(fPID[ix1],fPID[ix2]);
        if constexpr (!chargedhadrons) {
          fillpairs(fPID[ix1] % 2,fPID[ix2] % 2); /* always fill the charged hadron correlations if identified */
        }
      }
    }
  }
  for (uint pidx = 0; pidx < fSpeciesNames.size(); ++pidx) {
    for (uint pidy = 0; pidy < fSpeciesNames.size(); ++pidy) {
      fhN2_12_vsC[pidx][pidy]->Fill(fCentrality, fN2_12[pidx][pidy]);
      fhSum2PtPt_12_vsC[pidx][pidy]->Fill(fCentrality, fSum2PtPt_12[pidx][pidy]);
      fhSum2DptDpt_12_vsC[pidx][pidy]->Fill(fCentrality, fSum2DptDpt_12[pidx][pidy]);
      fhN2nw_12_vsC[pidx][pidy]->Fill(fCentrality, fNnw2_12[pidx][pidy]);
      fhSum2PtPtnw_12_vsC[pidx][pidy]->Fill(fCentrality, fSum2PtPtnw_12[pidx][pidy]);
      fhSum2DptDptnw_12_vsC[pidx][pidy]->Fill(fCentrality, fSum2DptDptnw_12[pidx][pidy]);

      /* update the number of entries in the differential histograms filled with AddBinContent */
      fhN2_12_vsDEtaDPhi[pidx][pidy]->SetEntries(fhN2_12_vsDEtaDPhi[pidx][pidy]->GetEntries() + fNnw2_12[pidx][pidy]);
      fhN2_12_vsDEtaDPhi_na[pidx][pidy]->SetEntries(fhN2_12_vsDEtaDPhi_na[pidx][pidy]->GetEntries() + fNnw2_12[pidx][pidy]);
      fhSum2PtPt_12_vsDEtaDPhi[pidx][pidy]->SetEntries(fhSum2PtPt_12_vsDEtaDPhi[pidx][pidy]->GetEntries() + fNnw2_12[pidx][pidy]);
      fhSum2DptDpt_12_vsDEtaDPhi[pidx][pidy]->SetEntries(fhSum2DptDpt_12_vsDEtaDPhi[pidx][pidy]->GetEntries() + fNnw2_12[pidx][pidy]);
      fhSum2DptDpt_12_vsDEtaDPhi_na[pidx][pidy]->SetEntries(fhSum2DptDpt_12_vsDEtaDPhi_na[pidx][pidy]->GetEntries() + fNnw2_12[pidx][pidy]);
    }
  }
}

/// \brief Fill the histograms once the whole process is finished
void  Ali2PCorrelations::FinalizeProcess()
{
  if (fSinglesOnly) {
    for (uint pid = 0; pid < fSpeciesNames.size(); ++pid) {
      /* for the time being, the errors are trivial so, we do not use results file space */
      fhN1_vsPt[pid]->Sumw2(false);
      fhN1_vsZEtaPhiPt[pid]->Sumw2(false);
    }
  } else {
    for (uint pidx = 0; pidx < fSpeciesNames.size(); ++pidx) {
      /* for the time being, the errors are trivial or not valid so, we do not use results file space */
      fhN1_vsEtaPhi[pidx]->Sumw2(false);
      fhSum1Pt_vsEtaPhi[pidx]->Sumw2(false);
      for (uint pidy = 0; pidy < fSpeciesNames.size(); ++pidy) {
        fhN2_12_vsDEtaDPhi[pidx][pidy]->Sumw2(false);
        fhN2_12_vsDEtaDPhi_na[pidx][pidy]->Sumw2(false);
        fhSum2PtPt_12_vsDEtaDPhi[pidx][pidy]->Sumw2(false);
        fhSum2DptDpt_12_vsDEtaDPhi[pidx][pidy]->Sumw2(false);
        fhSum2DptDpt_12_vsDEtaDPhi_na[pidx][pidy]->Sumw2(false);
        fhN2_12_vsPtPt[pidx][pidy]->Sumw2(false);
      }
    }
  }
}

/// \cond CLASSIMP
ClassImp(Ali2PCorrelations);
/// \endcond
