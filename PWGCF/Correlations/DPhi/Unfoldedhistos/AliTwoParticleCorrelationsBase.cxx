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
AliTwoParticleCorrelationsBase::AliTwoParticleCorrelationsBase()
  : TNamed(),
    fOutput(NULL),
    fVertexZ(0.0),
    fIxVertexZ(0),
    fCentrality(0.0),
    fSinglesOnly(kTRUE),
    fUseWeights(kFALSE),
    fUseSimulation(kFALSE),
    /* the arrays with tracks 1 and 2 information */
    fArraySize(10 * 1024),
    fNoOfTracks(0),
    fFlags(nullptr),
    fPt(nullptr),
    fEta(nullptr),
    fPhi(nullptr),
    fPID(nullptr),
    /* correction weights */
    fCorrectionWeights{},
    /* efficiency correction */
    fEfficiencyCorrection{},
    fPairsEfficiency{},
    /* simulation pdfs */
    fTrackPdfs{},
    fTrackCurrentPdf{},
    fSimEventsPerEvent(1),
    /* vertex bins */
    fNBins_vertexZ(40),
    fMin_vertexZ(-10.0),
    fMax_vertexZ(10.0),
    fWidth_vertexZ(0.5),
    /* phi origin shift */
    fNBinsPhiShift(0.0),
    /* pT bins */
    fNBins_pt(18),
    fMin_pt(0.2),
    fMax_pt(2.0),
    fWidth_pt(0.1),
    /* phi bins */
    fNBins_phi(72),
    fMin_phi(0.0),
    fMax_phi(TMath::Pi() * 2.0),
    fWidth_phi(TMath::Pi() * 2.0 / 72.0),
    fExcludeTOFHole(false),
    /* eta bins */
    fNBins_eta(20),
    fMin_eta(-1.0),
    fMax_eta(1.0),
    fWidth_eta(0.1),
    fNBins_etaPhi(1440),
    fNBins_etaPhiPt(25920),
    fNBins_zEtaPhiPt(1036800),
    /* accumulated values */
    fN1{},
    fNnw1{},
    fSum1Pt{},
    fSum1Ptnw{},
    /* histograms */
    fhN1_vsPt{},
    fhN1_vsEtaPhi{},
    fhSum1Pt_vsEtaPhi{},
    fhN1_vsZEtaPhiPt{},
    fhSum1Pt_vsZEtaPhiPt{},
    /* vs centrality profiles */
    fhN1_vsC{},
    fhSum1Pt_vsC{},
    fhN1nw_vsC{},
    fhSum1Ptnw_vsC{},
    fSpeciesNames{},
    fPreRejectResonances(false),
    fPostRejectResonances(false),
    fhResonanceRoughMasses(NULL),
    fhResonanceMasses(NULL),
    fhDiscardedResonanceMasses(NULL)
{
  for (Int_t ires = 0; ires < 16; ires++)
    fThresholdMult[ires] = 0x0;
}

/// Normal constructor
/// \param name the name for the object instance
AliTwoParticleCorrelationsBase::AliTwoParticleCorrelationsBase(const char* name)
  : TNamed(name, name),
    fOutput(NULL),
    fVertexZ(0.0),
    fIxVertexZ(0),
    fCentrality(0.0),
    fSinglesOnly(kTRUE),
    fUseWeights(kFALSE),
    fUseSimulation(kFALSE),
    /* the arrays with tracks 1 and 2 information */
    fArraySize(10 * 1024),
    fNoOfTracks(0),
    fFlags(nullptr),
    fPt(nullptr),
    fEta(nullptr),
    fPhi(nullptr),
    fPID(nullptr),
    /* correction weights */
    fCorrectionWeights{},
    /* efficiency correction */
    fEfficiencyCorrection{},
    fPairsEfficiency{},
    /* simulation pdfs */
    fTrackPdfs{},
    fTrackCurrentPdf{},
    fSimEventsPerEvent(1),
    /* vertex bins */
    fNBins_vertexZ(40),
    fMin_vertexZ(-10.0),
    fMax_vertexZ(10.0),
    fWidth_vertexZ(0.5),
    /* phi origin shift */
    fNBinsPhiShift(0.0),
    /* pT bins */
    fNBins_pt(18),
    fMin_pt(0.2),
    fMax_pt(2.0),
    fWidth_pt(0.1),
    /* phi bins */
    fNBins_phi(72),
    fMin_phi(0.0),
    fMax_phi(TMath::Pi() * 2.0),
    fWidth_phi(TMath::Pi() * 2.0 / 72.0),
    fExcludeTOFHole(false),
    /* eta bins */
    fNBins_eta(20),
    fMin_eta(-1.0),
    fMax_eta(1.0),
    fWidth_eta(0.1),
    fNBins_etaPhi(1440),
    fNBins_etaPhiPt(25920),
    fNBins_zEtaPhiPt(1036800),
    /* accumulated values */
    fN1{},
    fNnw1{},
    fSum1Pt{},
    fSum1Ptnw{},
    /* histograms */
    fhN1_vsPt{},
    fhN1_vsEtaPhi{},
    fhSum1Pt_vsEtaPhi{},
    fhN1_vsZEtaPhiPt{},
    fhSum1Pt_vsZEtaPhiPt{},
    /* vs centrality profiles */
    fhN1_vsC{},
    fhSum1Pt_vsC{},
    fhN1nw_vsC{},
    fhSum1Ptnw_vsC{},
    fSpeciesNames{},
    fPreRejectResonances(false),
    fPostRejectResonances(false),
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
  delete[] fFlags;
  delete[] fPt;
  delete[] fEta;
  delete[] fPhi;
  delete[] fPID;

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
    if (a->GetEntries() != DPTDPTCORRBINCONFIGPAR+1) {
      delete a;
      return false;
    } else {
      if (TString(a->At(DPTDPTCORRBINCONFIGPAR)->GetName()).EqualTo("TOF")) {
        fExcludeTOFHole = true;
      } else {
        delete a;
        return false;
      }
    }
  }
  sscanf(stritem, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d",
      &fMin_vertexZ, &fMax_vertexZ, &fWidth_vertexZ,
      &min_pt, &max_pt, &width_pt,
      &min_eta, &max_eta, &width_eta, &nBins_phi);
  delete a;

  fMin_pt = min_pt;
  fMax_pt = max_pt;
  fWidth_pt = width_pt;
  fMin_eta = min_eta;
  fMax_eta = max_eta;
  fWidth_eta = width_eta;
  fNBins_phi = nBins_phi;
  return true;
}

/// \brief Establishes the resonances rejection configuration
/// \param confstring string containing the resonances rejection configuration parameters
/// Basically one digit per supported resonance and the digit is the factor in one fourth
/// of the modulus around the resonance mass
void AliTwoParticleCorrelationsBase::ConfigureResonances(const char *confstring) {
  AliInfo(TString::Format("Resonances code string: %s", confstring));

  /* few sanity checks */
  TString str = confstring;
  if (!str.Contains("resonances:"))
    return;
  if (str.Length() > 45) {
    AliFatal(TString::Format("Resonances string suspisciouslly large %d: %s", str.Length(), str.Data()).Data());
  }
  if (!str.Contains(" ")) {
    AliFatal(TString::Format("Resonances string with old format: %s", str.Data()).Data());
  }

  Int_t rescode;
  char prepost[32] = "";
  sscanf(str.Data(), "resonances:%s %d", prepost, &rescode);

  if (TString(prepost).Contains("pre")) {
    fPreRejectResonances = true;
  }
  if (TString(prepost).Contains("post")) {
    fPostRejectResonances = true;
  }
  if (TString(prepost).Contains("none")) {
    fPreRejectResonances = false;
    fPostRejectResonances = false;
  }

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
TString AliTwoParticleCorrelationsBase::GetBinningConfigurationString() const
{
  return TString::Format("%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%d%s",
                         fMin_vertexZ, fMax_vertexZ, fWidth_vertexZ,
                         fMin_pt, fMax_pt, fWidth_pt,
                         fMin_eta, fMax_eta, fWidth_eta, fNBins_phi,
                         fExcludeTOFHole ? ",TOF" : "");
}

/// \brief Build the resonances rejection configuration string
/// \return the configuration string corresponding to the current resonance rejection configuration
TString AliTwoParticleCorrelationsBase::GetResonancesConfigurationString() const {
  TString str = "resonances:";
  if (fPreRejectResonances) {
    str += "pre";
  }
  if (fPostRejectResonances) {
    str += "post";
  }
  if (!fPreRejectResonances && !fPostRejectResonances) {
    str += "none";
  }

  str += " ";

  for (Int_t ires = 0; ires < fgkNoOfResonances; ires++) {
    str += TString::Format("%01d",fThresholdMult[fgkNoOfResonances - 1 - ires]);
  }
  return str;
}

void AliTwoParticleCorrelationsBase::SetSpeciesNames(std::vector<std::string> names)
{
  fSpeciesNames.clear();
  fSpeciesNames.assign(names.begin(), names.end());

  for (auto s : fSpeciesNames) {
    fN1.push_back(0.0);
    fSum1Pt.push_back(0.0);
    fNnw1.push_back(0.0);
    fSum1Ptnw.push_back(0.0);
  }
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
  fNBins_pt = Int_t(0.5 + (fMax_pt - fMin_pt) / fWidth_pt);
  fNBins_eta = Int_t(0.5 + (fMax_eta - fMin_eta) / fWidth_eta);
  fWidth_phi = (fMax_phi - fMin_phi) / fNBins_phi;
  fNBins_etaPhi = fNBins_phi * fNBins_eta;
  fNBins_etaPhiPt = fNBins_etaPhi * fNBins_pt;
  fNBins_zEtaPhiPt = fNBins_vertexZ * fNBins_etaPhiPt;

  /* incorporate configuration parameters */
  fOutput->Add(new TParameter<Int_t>("NoBinsVertexZ",fNBins_vertexZ,'f'));
  fOutput->Add(new TParameter<Int_t>("NoBinsPt", fNBins_pt, 'f'));
  fOutput->Add(new TParameter<Int_t>("NoBinsEta", fNBins_eta, 'f'));
  fOutput->Add(new TParameter<Int_t>("NoBinsPhi", fNBins_phi, 'f'));
  fOutput->Add(new TParameter<Double_t>("MinVertexZ",fMin_vertexZ,'f'));
  fOutput->Add(new TParameter<Double_t>("MaxVertexZ",fMax_vertexZ,'f'));
  fOutput->Add(new TParameter<Double_t>("MinPt", fMin_pt, 'f'));
  fOutput->Add(new TParameter<Double_t>("MaxPt", fMax_pt, 'f'));
  fOutput->Add(new TParameter<Double_t>("MinEta", fMin_eta, 'f'));
  fOutput->Add(new TParameter<Double_t>("MaxEta", fMax_eta, 'f'));
  fOutput->Add(new TParameter<Double_t>("MinPhi", fMin_phi, 'f'));
  fOutput->Add(new TParameter<Double_t>("MaxPhi", fMax_phi, 'f'));
  fOutput->Add(new TParameter<Double_t>("NoBinsPhiShift", fNBinsPhiShift, 'f'));
  fOutput->Add(new TParameter<bool>("ExcludeTOFHole", fExcludeTOFHole, 'f'));

  /* incorporate the resonance rejection configuration parameter */
  Int_t rescode = 0;
  Int_t mult = 1;
  for (Int_t ires = 0; ires < fgkNoOfResonances; ires++) {
    rescode += fThresholdMult[ires] * mult;
    mult *= 10;
  }
  fOutput->Add(new TParameter<Int_t>("ResonancesCode",rescode,'f'));

  /* after the parameters dump the proper phi limits are set according to the phi shift */
  fMax_phi = fMax_phi - fWidth_phi * fNBinsPhiShift;
  fMin_phi = fMin_phi - fWidth_phi * fNBinsPhiShift;

  /* track 1 storage */
  fFlags = new unsigned int[fArraySize];
  fPt = new Float_t[fArraySize];
  fEta = new Float_t[fArraySize];
  fPhi = new Float_t[fArraySize];
  fPID = new int[fArraySize];

  if (fSinglesOnly)  {
    for (unsigned int i = 0; i < fSpeciesNames.size(); ++i) {
      fhN1_vsPt.push_back(new TH1F(Form("n1_%s_vsPt", fSpeciesNames[i].c_str()), "#LT n_{1} #GT;p_{T} (GeV/c);#LT n_{1} #GT", fNBins_pt, fMin_pt, fMax_pt));
      fhN1_vsZEtaPhiPt.push_back(new TH3F(Form("n1_%s_vsZ_vsEtaPhi_vsPt", fSpeciesNames[i].c_str()), "#LT n_{1} #GT;vtx_{z};#eta#times#varphi;p_{T} (GeV/c)",
                                          fNBins_vertexZ, fMin_vertexZ, fMax_vertexZ, fNBins_etaPhi, 0.0, double(fNBins_etaPhi), fNBins_pt, fMin_pt, fMax_pt));
      fhSum1Pt_vsZEtaPhiPt.push_back(new TH3F(Form("sumPt1_%s_vsZ_vsEtaPhi_vsPt", fSpeciesNames[i].c_str()), "#LT #Sigma p_{T} #GT;vtx_{z};#eta#times#varphi;p_{T} (GeV/c)",
                                              fNBins_vertexZ, fMin_vertexZ, fMax_vertexZ, fNBins_etaPhi, 0.0, Double_t(fNBins_etaPhi), fNBins_pt, fMin_pt, fMax_pt));
      fhN1_vsEtaPhi.push_back(new TH2F(Form("n1_%s_vsEtaPhi", fSpeciesNames[i].c_str()), "#LT n_{1} #GT;#eta;#varphi (radian);#LT n_{1} #GT",
                                       fNBins_eta, fMin_eta, fMax_eta, fNBins_phi, fMin_phi, fMax_phi));
      fhSum1Pt_vsEtaPhi.push_back(new TH2F(Form("sumPt_%s_vsEtaPhi", fSpeciesNames[i].c_str()), "#LT #Sigma p_{T} #GT;#eta;#varphi (radian);#LT #Sigma p_{T} #GT (GeV/c)",
                                           fNBins_eta, fMin_eta, fMax_eta, fNBins_phi, fMin_phi, fMax_phi));
    }
    for (unsigned int i = 0; i < fSpeciesNames.size(); ++i) {
      fOutput->Add(fhN1_vsPt[i]);
      fOutput->Add(fhN1_vsZEtaPhiPt[i]);
      fOutput->Add(fhSum1Pt_vsZEtaPhiPt[i]);
      fOutput->Add(fhN1_vsEtaPhi[i]);
      fOutput->Add(fhSum1Pt_vsEtaPhi[i]);
    }
  } else {
    fhResonanceRoughMasses = new TH2F("ResonanceRoughMasses", "Approx invariant mass;res id;Inv mass", fgkNoOfResonances, -0.5, fgkNoOfResonances - 0.5, 1000, 0.0, 5.0);
    fhResonanceMasses = new TH2F("ResonanceMasses", "Invariant mass;res id;Inv mass", fgkNoOfResonances, -0.5, fgkNoOfResonances - 0.5, 1000, 0.0, 5.0);
    fhDiscardedResonanceMasses = new TH2F("DiscardedResonanceMasses", "Discarded invariant mass;res id;Inv mass", fgkNoOfResonances, -0.5, fgkNoOfResonances - 0.5, 1000, 0.0, 5.0);

    fOutput->Add(fhResonanceRoughMasses);
    fOutput->Add(fhResonanceMasses);
    fOutput->Add(fhDiscardedResonanceMasses);

    for (unsigned int i = 0; i < fSpeciesNames.size(); ++i) {
      fhN1_vsEtaPhi.push_back(new TH2F(Form("n1_%s_vsEtaPhi", fSpeciesNames[i].c_str()), "#LT n_{1} #GT;#eta;#varphi (radian);#LT n_{1} #GT",
                                       fNBins_eta, fMin_eta, fMax_eta, fNBins_phi, fMin_phi, fMax_phi));
      fhSum1Pt_vsEtaPhi.push_back(new TH2F(Form("sumPt_%s_vsEtaPhi", fSpeciesNames[i].c_str()), "#LT #Sigma p_{T} #GT;#eta;#varphi (radian);#LT #Sigma p_{T} #GT (GeV/c)",
                                           fNBins_eta, fMin_eta, fMax_eta, fNBins_phi, fMin_phi, fMax_phi));
      fhN1_vsC.push_back(new TProfile(Form("n1_%s_vsM", fSpeciesNames[i].c_str()), "#LT n_{1} #GT (weighted);Centrality (%);#LT n_{1} #GT", 100, 0.0, 100.0));
      fhSum1Pt_vsC.push_back(new TProfile(Form("sumPt_%s_vsM", fSpeciesNames[i].c_str()), "#LT #Sigma p_{T} #GT (weighted);Centrality (%);#LT #Sigma p_{T} #GT (GeV/c)", 100, 0.0, 100.0));
      fhN1nw_vsC.push_back(new TProfile(Form("n1Nw_%s_vsM", fSpeciesNames[i].c_str()), "#LT n_{1} #GT;Centrality (%);#LT n_{1} #GT", 100, 0.0, 100.0));
      fhSum1Ptnw_vsC.push_back(new TProfile(Form("sumPtNw_%s_vsM", fSpeciesNames[i].c_str()), "#LT #Sigma p_{T} #GT;Centrality (%);#LT #Sigma p_{T} #GT (GeV/c)", 100, 0.0, 100.0));
    }

    for (unsigned int i = 0; i < fSpeciesNames.size(); ++i) {
      fOutput->Add(fhN1_vsEtaPhi[i]);
      fOutput->Add(fhSum1Pt_vsEtaPhi[i]);
      fOutput->Add(fhN1_vsC[i]);
      fOutput->Add(fhSum1Pt_vsC[i]);
      fOutput->Add(fhN1nw_vsC[i]);
      fOutput->Add(fhSum1Ptnw_vsC[i]);
    }
  }
  TH1::AddDirectory(oldstatus);

  /* after histogram creation the proper phi upper limit is set according to inclusion or not of the TOF hole */
  if (fExcludeTOFHole) {
    const int nPhiBinsTOFHole = 25;
    fMax_phi = fMax_phi - fWidth_phi * nPhiBinsTOFHole;
  }
}

/// \brief Stores the correction weights
/// Have in mind that there are several situations
/// In the usual case of not identified analysis only weights
/// for HaP and HaM weill be received (different than nullptr)
/// and they miight or might not inccorporate the efficiency correction
/// For identified analyses the received weights are ususally
/// for HaP and HaM as well but they will not incorporate
/// the species efficiency correction
/// Optionally there could be weights for the identified species
/// as that case is also supported
/// \param h3 vector with the weights histograms
/// \return kTRUE if correctly done kFALSE otherwise
Bool_t AliTwoParticleCorrelationsBase::SetWeigths(std::vector<const TH3*> h3)
{
  Bool_t done = kTRUE;

  fCorrectionWeights.reserve(fSpeciesNames.size());
  if (fUseWeights){
    if (h3.size() > 0) {
      bool storedweights = false;
      for (unsigned int i = 0; i < fSpeciesNames.size(); ++i) {
        if (h3[i] != nullptr) {
          /* allocate memory for each track species weights */
          float* buffer = new float[fNBins_vertexZ * fNBins_etaPhi * fNBins_pt];
          fCorrectionWeights[i] = buffer;

          for (int ixZ = 0; ixZ < fNBins_vertexZ; ixZ++) {
            double zval = fMin_vertexZ + fWidth_vertexZ * (ixZ + 0.5);
            for (int ixEtaPhi = 0; ixEtaPhi < fNBins_etaPhi; ixEtaPhi++) {
              for (int ixPt = 0; ixPt < fNBins_pt; ixPt++) {
                double pTval = fMin_pt + fWidth_pt * (ixPt + 0.5);
                buffer[ixZ * fNBins_etaPhi * fNBins_pt + ixEtaPhi * fNBins_pt + ixPt] = h3[i]->GetBinContent(h3[i]->GetXaxis()->FindBin(zval), ixEtaPhi + 1, h3[i]->GetZaxis()->FindBin(pTval));
              }
            }
          }
          storedweights = true;
        }
      }
      if (!storedweights) {
        AliFatal(TString::Format("The weights histogram are not there. ABORTING!!!"));
        done = false;
      }
    } else {
      AliFatal("The weights histograms for different species are not there. ABORTING!!!");
      done = kFALSE;
    }
  } else {
    for (unsigned int isp = 0; isp < fSpeciesNames.size(); ++isp) {
      fCorrectionWeights[isp] = nullptr;
    }
    AliError("Setting weights for a not use weights instance. Ignoring it");
    done = kFALSE;
  }
  return done;
}

/// \brief Stores the efficiency/contamination correction for the different track species
/// If the efficiency correction histograms are not present a efficiency correction value 1.0 is stored.
/// The correction values are always the stored in the corresponding histograms
/// For identified analysis the weight are separated from the efficiency and only the identified efficiencies
/// are stored. Have in mind that the species names always include the hadrons as first two species
/// The order for the particle species is HaP, HaM, PiP, PiM, KaP, KaM, PrP, PrM
/// For unidentified analysis (only hadrons) the efficiency might be incorporated in the weights
/// \param h1 vector with the efficiency histograms
/// \return kTRUE if everything went ok otherwise kFALSE
Bool_t AliTwoParticleCorrelationsBase::SetEfficiencyCorrection(std::vector<const TH1*> h1)
{
  bool allnull = true;
  for (unsigned int pid = 0; pid < fSpeciesNames.size(); ++pid) {
    if (h1[pid] != nullptr) {
      allnull = false;
    }
  }

  /* pre-allocate the number of entries */
  fEfficiencyCorrection.reserve(fSpeciesNames.size());

  if (!allnull) {
    /* setting the efficiency correction for the different species */
    for (unsigned int isp = 0; isp < fSpeciesNames.size(); ++isp) {
      if (h1[isp] != nullptr) {
        /* allocate memory for each species efficiency correction */
        float* buffer = new float[fNBins_pt];
        fEfficiencyCorrection[isp] = buffer;
        for (int ixPt = 0; ixPt < fNBins_pt; ixPt++) {
          int bin = h1[isp]->GetXaxis()->FindBin(fMin_pt + fWidth_pt / 2.0 + ixPt * fWidth_pt);
          buffer[ixPt] = h1[isp]->GetBinContent(bin);
        }
      } else {
        float* buffer = new float[fNBins_pt];
        fEfficiencyCorrection[isp] = buffer;
        for (int ixPt = 0; ixPt < fNBins_pt; ixPt++) {
          buffer[ixPt] = 1.0;
        }
      }
    }
  } else {
    AliInfo("Setting efficiency correction equal to one for all species");
    for (unsigned int isp = 0; isp < fSpeciesNames.size(); ++isp) {
      float* buffer = new float[fNBins_pt];
      fEfficiencyCorrection[isp] = buffer;
      for (int ixPt = 0; ixPt < fNBins_pt; ixPt++) {
        buffer[ixPt] = 1.0;
      }
    }
  }
  AliInfo(Form("Configured efficiency correcton on %s", GetName()));
  for (unsigned int i = 0; i < fSpeciesNames.size(); ++i) {
    printf("Species %s:\n", fSpeciesNames[i].c_str());
    for (int bin = 0; bin < fNBins_pt; bin++) {
      printf("%f, ", fEfficiencyCorrection[i][bin]);
      if ((bin + 1) % 8 == 0)
        printf("\n");
    }
    printf("\n");
  }
  return true;
}

/// \brief Stores the pairs efficiency correction for track species combinations
/// The correction value is stored as the inverse of the efficiency ones
/// \param hn vector with the pair efficiency histograms
/// \return kTRUE if everything went ok otherwise kFALSE
Bool_t AliTwoParticleCorrelationsBase::SetPairEfficiencyCorrection(std::vector<std::vector<const THn*>> hn)
{
  AliInfo("");

  Bool_t done = kTRUE;

  if (hn.size() > 0) {
    /* few sanity checks */
    bool abort = false;
    bool allnull = true;
    if (hn.size() != fSpeciesNames.size()) {
      abort = true;
    } else {
      for (unsigned int i = 0 && allnull; i < fSpeciesNames.size(); ++i) {
        if (hn[i].size() != fSpeciesNames.size()) {
          abort = true;
          break;
        } else {
          for (unsigned int j = 0; j < fSpeciesNames.size(); ++j) {
            if (hn[i][j] != nullptr) {
              allnull = false;
              break;
            }
          }
        }
      }
    }
    if (!abort) {
      for (unsigned int i = 0; i < fSpeciesNames.size(); ++i) {
        std::vector<const THn*> row;
        for (unsigned int j = 0; j < fSpeciesNames.size(); ++j) {
          if (!allnull) {
            if (hn[i][j] != nullptr) {
              if (hn[i][j]->GetNdimensions() != 4) {
                AliFatal(TString::Format("Pair efficiency correction histogram %s-%s expected to have four dimensions. ABORTING!!!", fSpeciesNames[i].c_str(), fSpeciesNames[j].c_str()));
                done = kFALSE;
              } else {
                /* store the pair 11 efficiency histogram */
                row.push_back(hn[i][j]);
              }
            } else {
              AliFatal(TString::Format("The pair efficiency correction histogram %s-%s is a null pointer. ABORTING!!!", fSpeciesNames[i].c_str(), fSpeciesNames[j].c_str()));
              done = kFALSE;
            }
          } else {
            row.push_back(nullptr);
          }
        }
        fPairsEfficiency.push_back(row);
      }
    } else {
      AliFatal("The pair efficiency correction histograms vector is malformed. ABORTING!!!");
      done = kFALSE;
    }
  }
  return done;
}

/// \brief Stores the pdf for the different species depending on z vertex bin
///
/// Pdf are used to generate simulated tracks
/// \param pluspdf array with z vertex dependent pdf for positive tracks
/// \param minuspdf array with z vertex dependent pdf for negative tracks
/// \return kTRUE if everything went ok otherwise kFALSE
Bool_t AliTwoParticleCorrelationsBase::SetSimultationPdfs(std::vector<const TObjArray*> pdf)
{

  bool done = true;
  bool allnull = true;

  for (unsigned int pid = 0; pid < pdf.size(); ++pid) {
    if (pdf[pid] != nullptr) {
      allnull = false;
      break;
    }
  }

  if (fUseSimulation) {
    if (!allnull) {
      for (unsigned int i = 0; i < fSpeciesNames.size(); ++i) {
        if (i < pdf.size() && pdf[i] != nullptr) {
          if (pdf[i]->GetEntriesFast() == fNBins_vertexZ) {
            fTrackPdfs.push_back(pdf[i]);
          } else {
            AliFatal(Form("The track density histogram array number of pdfs %d differ from the number of zvertex bins %d for species %s. ABORTING!!!",
                          pdf[i]->GetEntriesFast(), fNBins_vertexZ, fSpeciesNames[i].c_str()));
            done = kFALSE;
          }
        } else {
          AliFatal(TString::Format("The track density histogram for species %s is a null pointer. ABORTING!!!", fSpeciesNames[i].c_str()));
          done = kFALSE;
        }
      }
    } else {
      AliFatal("The track density histograms for different species are not there. ABORTING!!!");
      done = kFALSE;
    }
  } else {
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
    for (unsigned int i = 0; i < fSpeciesNames.size(); ++i) {
      fTrackCurrentPdf.push_back((TH3F*)fTrackPdfs[i]->At(fIxVertexZ));
    }
  }

  fNoOfTracks = 0;

  /* reset single counters */
  for (unsigned int i = 0; i < fSpeciesNames.size(); ++i) {
    fN1[i] = fSum1Pt[i] = fNnw1[i] = fSum1Ptnw[i] = 0;
  }

  return kTRUE;
}

/// \brief check if the track is accepted according to the correlations cuts
///
/// \param trk the passed track
/// \return kTRUE if the track is accepted kFALSE otherwise
bool AliTwoParticleCorrelationsBase::IsTrackAccepted(AliVTrack* trk)
{
  return IsTrackAccepted(float(trk->Pt()), float(trk->Eta()), float(trk->Phi()));
}

/// \brief check if the true particle is accepted according to the correlations cuts
///
/// \param part the passed particle
/// \return kTRUE if the particle is properly handled kFALSE otherwise
bool AliTwoParticleCorrelationsBase::IsTrackAccepted(AliVParticle* part)
{
  if (part->Charge() != 0) {
    return IsTrackAccepted(float(part->Pt()), float(part->Eta()), float(part->Phi()));
  }
  else {
    return kFALSE;
  }
}

/// \brief checks if the track is accepted according to the configured cuts 
/// \param pT the track \f$ p_T \f$
/// \param eta the track \f$ \eta \f$
/// \param phi the track \f$ \phi \f$
/// \return kTRUE if the track is accepted
/// This is a pre-check so for the time being the pT cut is not checked
/// to allow PID information depending on momentum
bool AliTwoParticleCorrelationsBase::IsTrackAccepted(float, float eta, float ophi)
{
  /* consider a potential phi origin shift */
  float phi = ophi;
  if (!(phi < fMax_phi))
    phi = phi - 2 * TMath::Pi();
  if (phi < fMin_phi) {
    return kFALSE;
  }
  float ixPhi = int((phi - fMin_phi) / fWidth_phi);
  if (ixPhi < 0 || !(ixPhi < fNBins_phi)) {
    return kFALSE;
  }

  if (eta < fMin_eta || fMax_eta < eta) {
    return kFALSE;
  }

  int ixEta = int((eta - fMin_eta) / fWidth_eta);
  if (ixEta < 0 || !(ixEta < fNBins_eta)) {
    return kFALSE;
  }
  return kTRUE;
}


/// \brief process a track and store its parameters if feasible
///
/// If simulation is ordered the actual track is discarded and a new one with the
/// same charge is produced out of the corresponding track pdf
/// \param pid the external track Id
/// \param trk the passed track
/// \return kTRUE if the track is properly handled kFALSE otherwise
bool AliTwoParticleCorrelationsBase::ProcessTrack(int pid, AliVTrack* trk)
{

  if (fUseSimulation) {
    if (pid > 1) {
      AliFatal("Simulation still not prepared for PID");
      return false;
    } else {
      double pT;
      double eta;
      double phi;

      fTrackCurrentPdf[pid]->GetRandom3(eta, phi, pT);

      return ProcessTrack((trk->Charge() > 0) ? 0 : 1, pT, eta, phi);
    }
  } else {
    if (trk->Charge() > 0) {
      return ProcessTrack(2*pid, float(trk->Pt()), float(trk->Eta()), float(trk->Phi()));
    } else if (trk->Charge() < 0) {
      return ProcessTrack(2*pid+1, float(trk->Pt()), float(trk->Eta()), float(trk->Phi()));
    } else {
      return kFALSE;
    }
  }
}

/// \brief process a true particle and store its parameters if feasible
///
/// If simulation is orderd the track is discarded and kFALSE is returned
/// \param pid the external particle Id
/// \param part the passed particle
/// \return kTRUE if the particle is properly handled kFALSE otherwise
bool AliTwoParticleCorrelationsBase::ProcessTrack(int pid, AliVParticle* part)
{

  if (fUseSimulation) {
    return kFALSE;
  } else {
    if (part->Charge() > 0) {
      return ProcessTrack(2*pid, float(part->Pt()), float(part->Eta()), float(part->Phi()));
    } else if (part->Charge() < 0) {
      return ProcessTrack(2*pid+1, float(part->Pt()), float(part->Eta()), float(part->Phi()));
    } else {
      return kFALSE;
    }
  }
}

/// \brief Flag the potential products of conversions and / or resonances
void AliTwoParticleCorrelationsBase::FlagConversionsAndResonances() {
  AliInfo("");
  /* pair with different charges. Both track list are different */
  for (Int_t ix1 = 0; ix1 < fNoOfTracks; ix1++) {
      for (Int_t ix2 = ix1 + 1; ix2 < fNoOfTracks; ix2++) {

      if (fPID[ix1] != fPID[ix2]) {
        for (Int_t ires = 0; ires < fgkNoOfResonances; ires++) {
          if (fThresholdMult[ires] != 0) {
            Float_t mass = checkIfResonance(ires, kTRUE, fPt[ix1], fEta[ix1], fPhi[ix1], fPt[ix2], fEta[ix2], fPhi[ix2]);

            if (0 < mass) {
              fFlags[ix1] |= (0x1 << ires);
              fFlags[ix2] |= (0x1 << ires);
              fhResonanceMasses->Fill(ires, TMath::Sqrt(mass));
            }
          }
        }
      }
      } // ix2
  }     // ix1
}

/// \cond CLASSIMP
ClassImp(AliTwoParticleCorrelationsBase);
/// \endcond
