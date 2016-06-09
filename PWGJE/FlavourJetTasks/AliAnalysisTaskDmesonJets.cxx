/*************************************************************************
* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

// Root
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TVector3.h>
#include <THnSparse.h>
#include <TParticle.h>
#include <TMath.h>
#include <THashList.h>
#include <TFile.h>

// Aliroot general
#include "AliLog.h"
#include "AliEMCALGeometry.h"

// Aliroot HF
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliHFAODMCParticleContainer.h"
#include "AliHFTrackContainer.h"

// Aliroot EMCal jet framework
#include "AliEmcalJetTask.h"
#include "AliEmcalJet.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliEmcalParticle.h"
#include "AliFJWrapper.h"

#include "AliAnalysisTaskDmesonJets.h"

// Definitions of class AliAnalysisTaskDmesonJets::AliDmesonJetInfo

/// Reset all fields to their default values
void AliAnalysisTaskDmesonJets::AliDmesonJetInfo::Reset()
{
  fD.SetPtEtaPhiE(0,0,0,0);
  fSoftPionPt = 0;
  fInvMass2Prong = 0;
  for (auto &jet : fJets) {
    jet.second.fMomentum.SetPtEtaPhiE(0,0,0,0);
    jet.second.fNConstituents = 0;
    jet.second.fNEF = 0;
    jet.second.fMaxChargedPt = 0;
    jet.second.fMaxNeutralPt = 0;
  }
}

/// Prints the content of this object in the standard output.
void AliAnalysisTaskDmesonJets::AliDmesonJetInfo::Print() const
{
  Printf("Printing D Meson Jet object.");
  Printf("D Meson: pT = %.3f, eta = %.3f, phi = %.3f, inv. mass = %.3f", fD.Pt(), fD.Eta(), fD.Phi_0_2pi(), fD.M());
  Printf("Soft pion pT: %.3f. 2-Prong Invariant mass = %.3f", fSoftPionPt, fInvMass2Prong);
  for (auto &jet : fJets) {
    Printf("Jet %s: pT = %.3f, eta = %.3f, phi = %.3f", jet.first.c_str(), jet.second.Pt(), jet.second.Eta(), jet.second.Phi_0_2pi());
    Printf("Jet N Consituents = %d", jet.second.fNConstituents);
  }
}

/// Calculates the parallel fraction
///
/// \return the fraction of the momentum of the particle parallel to the jet over the total jet momentum
Double_t AliAnalysisTaskDmesonJets::AliDmesonJetInfo::GetZ(std::string n) const
{
  std::map<std::string, AliJetInfo>::const_iterator it = fJets.find(n);
  if (it == fJets.end()) return 0;

  Double_t z = 0;

  if ((*it).second.Pt() > 0) {
    TVector3 dvect = fD.Vect();
    TVector3 jvect = (*it).second.fMomentum.Vect();

    Double_t jetMom = jvect * jvect;

    if (jetMom < 1e-6) {
      ::Error("AliAnalysisTaskDmesonJets::AliDmesonJetInfo::GetZ", "Zero jet momentum!");
      z = 0.999;
    }
    else {
      z = (dvect * jvect) / jetMom;
    }

    if (z == 1 || (z > 1 && z - 1 < 1e-3)) z = 0.999; // so that it will contribute to the bin 0.9-1 rather than 1-1.1
  }

  return z;
}

/// Calculates the distance between the D meson and the jet axis
///
/// \param dR reference where the distance will be returned
/// \param deta reference where the eta distance will be returned
/// \param dphi reference where the phi distance will be returned
Double_t AliAnalysisTaskDmesonJets::AliDmesonJetInfo::GetDistance(std::string n, Double_t& deta, Double_t& dphi) const
{
  std::map<std::string, AliJetInfo>::const_iterator it = fJets.find(n);
  if (it == fJets.end()) return 0;

  dphi = TVector2::Phi_mpi_pi(fD.Phi() - (*it).second.Phi());;
  deta = fD.Eta() - (*it).second.Eta();
  return TMath::Sqrt(dphi*dphi + deta*deta);
}

/// Calculates the distance between the D meson and the jet axis
///
/// \param dR reference where the distance will be returned
Double_t AliAnalysisTaskDmesonJets::AliDmesonJetInfo::GetDistance(std::string n) const
{
  Double_t deta = 0;
  Double_t dphi = 0;
  return GetDistance(n, deta, dphi);
}

const AliAnalysisTaskDmesonJets::AliJetInfo* AliAnalysisTaskDmesonJets::AliDmesonJetInfo::GetJet(std::string n) const
{
  std::map<std::string, AliJetInfo>::const_iterator it = fJets.find(n);
  if (it == fJets.end()) {
    ::Error("AliAnalysisTaskDmesonJets::AliDmesonJetInfo::GetJet", "Could not find jet info for the jet definition '%s'!",
       n.c_str());
    return 0;
  }
  return &((*it).second);
}

AliAnalysisTaskDmesonJets::AliJetInfo* AliAnalysisTaskDmesonJets::AliDmesonJetInfo::GetJet(std::string n)
{
  std::map<std::string, AliJetInfo>::iterator it = fJets.find(n);
  if (it == fJets.end()) {
    ::Error("AliAnalysisTaskDmesonJets::AliDmesonJetInfo::GetJet", "Could not find jet info for the jet definition '%s'!",
        n.c_str());
    return 0;
  }
  return &((*it).second);
}

// Definitions of class AliAnalysisTaskDmesonJets::AliJetInfoSummary

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJets::AliJetInfoSummary);
/// \endcond

/// Constructor that uses an AliDmesonJetInfo as its source
///
/// \param source A const reference to a valid AliDmesonJetInfo object
/// \param i      Index of the jet to be copied
AliAnalysisTaskDmesonJets::AliJetInfoSummary::AliJetInfoSummary(const AliDmesonJetInfo& source, std::string n) :
  fPt(0),
  fEta(0),
  fPhi(0),
  fR(0),
  fZ(0)
{
  Set(source, n);
}

/// Reset the current object
void AliAnalysisTaskDmesonJets::AliJetInfoSummary::Reset()
{
  fPt = 0;
  fEta = 0;
  fPhi = 0;
  fR = 0;
  fZ = 0;
}

/// Set the current object using an instance of AliDmesonJetInfo as its source
///
/// \param source A const reference to a valid AliDmesonJetInfo object
/// \param i      Index of the jet to be copied
void AliAnalysisTaskDmesonJets::AliJetInfoSummary::Set(const AliDmesonJetInfo& source, std::string n)
{
  std::map<std::string, AliJetInfo>::const_iterator it = source.fJets.find(n);
  if (it == source.fJets.end()) return;

  fPt = (*it).second.Pt();
  fEta = (*it).second.Eta();
  fPhi = (*it).second.Phi_0_2pi();
  fR = source.GetDistance(n);
  fZ = source.GetZ(n);
}

// Definitions of class AliAnalysisTaskDmesonJets::AliDmesonInfoSummary

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJets::AliDmesonInfoSummary);
/// \endcond

/// Constructor that uses an AliDmesonJetInfo as its source
///
/// \param source A const reference to a valid AliDmesonJetInfo object
AliAnalysisTaskDmesonJets::AliDmesonInfoSummary::AliDmesonInfoSummary(const AliDmesonJetInfo& source) :
  fPt(0),
  fEta(0),
  fPhi(0)
{
  Set(source);
}

/// Set the current object using an instance of AliDmesonJetInfo as its source
///
/// \param source A const reference to a valid AliDmesonJetInfo object
void AliAnalysisTaskDmesonJets::AliDmesonInfoSummary::Set(const AliDmesonJetInfo& source)
{
  fPt = source.fD.Pt();
  fEta = source.fD.Eta();
  fPhi = source.fD.Phi_0_2pi();
}

// Definitions of class AliAnalysisTaskDmesonJets::AliD0InfoSummary

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJets::AliD0InfoSummary);
/// \endcond

/// Constructor that uses an AliDmesonJetInfo as its source
///
/// \param source A const reference to a valid AliDmesonJetInfo object
AliAnalysisTaskDmesonJets::AliD0InfoSummary::AliD0InfoSummary(const AliDmesonJetInfo& source) :
  AliDmesonInfoSummary(source),
  fInvMass(source.fD.M())
{
}

void AliAnalysisTaskDmesonJets::AliD0InfoSummary::Set(const AliDmesonJetInfo& source)
{
  fInvMass = source.fD.M();
  AliDmesonInfoSummary::Set(source);
}

// Definitions of class AliAnalysisTaskDmesonJets::AliDStarInfoSummary

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJets::AliDStarInfoSummary);
/// \endcond

/// Constructor that uses an AliDmesonJetInfo as its source
///
/// \param source A const reference to a valid AliDmesonJetInfo object
AliAnalysisTaskDmesonJets::AliDStarInfoSummary::AliDStarInfoSummary(const AliDmesonJetInfo& source) :
  AliDmesonInfoSummary(source),
  f2ProngInvMass(source.fInvMass2Prong),
  fDeltaInvMass(source.fD.M() - source.fInvMass2Prong)
{
}

void AliAnalysisTaskDmesonJets::AliDStarInfoSummary::Set(const AliDmesonJetInfo& source)
{
  f2ProngInvMass = source.fInvMass2Prong;
  fDeltaInvMass = source.fD.M() - source.fInvMass2Prong;
  AliDmesonInfoSummary::Set(source);
}

// Definitions of class AliAnalysisTaskDmesonJets::AliJetDefinition

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJets::AliHFJetDefinition);
/// \endcond

/// This is the default constructor, used for ROOT I/O purposes.
AliAnalysisTaskDmesonJets::AliHFJetDefinition::AliHFJetDefinition() :
  TObject(),
  fJetType(AliJetContainer::kChargedJet),
  fRadius(0),
  fJetAlgo(AliJetContainer::antikt_algorithm),
  fRecoScheme(AliJetContainer::pt_scheme),
  fAcceptance(AliJetContainer::kUser),
  fMinJetPt(0.),
  fMinJetPhi(0.),
  fMaxJetPhi(0.),
  fMinJetEta(0.),
  fMaxJetEta(0.),
  fMinChargedPt(0.),
  fMaxChargedPt(0.),
  fMinNeutralPt(0.),
  fMaxNeutralPt(0.)
{
}

/// Default constructor
///
/// \param type Jet type (full, charged, neutral)
/// \param r    Jet resolution parameter
/// \param algo Jet algorithm (anit-kt, kt,...)
/// \param reco Jet recombination scheme (pt_scheme, E_scheme,...)
AliAnalysisTaskDmesonJets::AliHFJetDefinition::AliHFJetDefinition(EJetType_t type, Double_t r, EJetAlgo_t algo, ERecoScheme_t reco) :
  TObject(),
  fJetType(type),
  fRadius(r),
  fJetAlgo(algo),
  fRecoScheme(reco),
  fAcceptance(AliJetContainer::kUser),
  fMinJetPt(0.),
  fMinJetPhi(0.),
  fMaxJetPhi(0.),
  fMinJetEta(0.),
  fMaxJetEta(0.),
  fMinChargedPt(0.),
  fMaxChargedPt(0.),
  fMinNeutralPt(0.),
  fMaxNeutralPt(0.)
{
  // By default set detector fiducial acceptance
  switch (type) {
  case AliJetContainer::kFullJet:
  case AliJetContainer::kNeutralJet:
    fAcceptance = AliJetContainer::kEMCALfid;
    break;
  case AliJetContainer::kChargedJet:
    fAcceptance = AliJetContainer::kTPCfid;
    break;
  }
}

/// Copy constructor
///
/// \param source Reference to an AliJetDefinition object to copy from
AliAnalysisTaskDmesonJets::AliHFJetDefinition::AliHFJetDefinition(const AliHFJetDefinition &source) :
  TObject(),
  fJetType(source.fJetType),
  fRadius(source.fRadius),
  fJetAlgo(source.fJetAlgo),
  fRecoScheme(source.fRecoScheme),
  fAcceptance(source.fAcceptance),
  fMinJetPt(source.fMinJetPt),
  fMinJetPhi(source.fMinJetPhi),
  fMaxJetPhi(source.fMaxJetPhi),
  fMinJetEta(source.fMinJetEta),
  fMaxJetEta(source.fMaxJetEta),
  fMinChargedPt(source.fMinChargedPt),
  fMaxChargedPt(source.fMaxChargedPt),
  fMinNeutralPt(source.fMinNeutralPt),
  fMaxNeutralPt(source.fMaxNeutralPt)
{
}

/// Assignment operator
///
/// \param source Reference to an AliJetDefinition object to copy from
AliAnalysisTaskDmesonJets::AliHFJetDefinition& AliAnalysisTaskDmesonJets::AliHFJetDefinition::operator=(const AliHFJetDefinition& source)
{
  new (this) AliHFJetDefinition(source);
  return *this;
}

/// Generate a name for this jet definition
const char* AliAnalysisTaskDmesonJets::AliHFJetDefinition::GetName() const
{
  static TString name;

  name = AliJetContainer::GenerateJetName(fJetType, fJetAlgo, fRecoScheme, fRadius, 0, 0, "Jet");

  return name.Data();
}

/// Decides whether the jet passes the acceptance cut defined in the object
///
/// \param jet Const reference to a AliJetInfo object
/// \return kTRUE if the jet passes the cuts
Bool_t AliAnalysisTaskDmesonJets::AliHFJetDefinition::IsJetInAcceptance(const AliJetInfo& jet) const
{
  if (fMinJetEta < fMaxJetEta && (jet.Eta() < fMinJetEta || jet.Eta() > fMaxJetEta)) return kFALSE;
  if (fMinJetPhi < fMaxJetPhi && (jet.Phi() < fMinJetPhi || jet.Phi() > fMaxJetPhi)) return kFALSE;
  if (jet.Pt() < fMinJetPt) return kFALSE;
  if (jet.fMaxChargedPt < fMinChargedPt || jet.fMaxChargedPt > fMaxChargedPt) return kFALSE;
  if (jet.fMaxNeutralPt < fMinNeutralPt || jet.fMaxNeutralPt > fMaxNeutralPt) return kFALSE;

  return kTRUE;
}

/// Decides whether the jet passes the acceptance cut defined in the object
///
/// \param jet Const reference to a AliJetInfo object
/// \return kTRUE if the jet passes the cuts
Bool_t AliAnalysisTaskDmesonJets::AliHFJetDefinition::IsJetInAcceptance(const AliDmesonJetInfo& dMesonJet, std::string n) const
{
  const AliJetInfo* jet = dMesonJet.GetJet(n);
  if (!jet) return kFALSE;
  return IsJetInAcceptance((*jet));
}

/// Sets the eta/phi acceptance of the jets to the detector boundaries
///
/// \param geom Const pointer to the EMCal/DCal geometry
/// \param run  Run number (needed because certain runs in 2012/2013 had the 1/3 EMCal modules off
void AliAnalysisTaskDmesonJets::AliHFJetDefinition::SetDetectorJetEtaPhiRange(const AliEMCALGeometry* const geom, Int_t run)
{
  Double_t r = 0;
  switch (fAcceptance) {
  case AliJetContainer::kTPCfid:
    r = fRadius;
    // enforce fiducial acceptance
    /* no break */
  case AliJetContainer::kTPC:
    SetJetEtaRange(-0.9 + r, 0.9 - r);
    SetJetPhiRange(0, 0);  // No cut on phi
    break;

  case AliJetContainer::kEMCALfid:
    r = fRadius;
    // enforce fiducial acceptance
    /* no break */
  case AliJetContainer::kEMCAL:
    if (geom) {
      SetJetEtaRange(geom->GetArm1EtaMin() + r, geom->GetArm1EtaMax() - r);

      if(run>=177295 && run<=197470) {//small SM masked in 2012 and 2013
        SetJetPhiRange(1.405 + r,3.135 - r);
      }
      else {
        SetJetPhiRange(geom->GetArm1PhiMin() * TMath::DegToRad() + r, geom->GetEMCALPhiMax() * TMath::DegToRad() - r);
      }
    }
    else {
      AliWarning("Could not get instance of AliEMCALGeometry. Using manual settings for EMCAL year 2011!!");
      SetJetEtaRange(-0.7 + r, 0.7 - r);
      SetJetPhiRange(1.405 + r, 3.135 - r);
    }
    break;

  case AliJetContainer::kDCALfid:
    r = fRadius;
    // enforce fiducial acceptance
    /* no break */
  case AliJetContainer::kDCAL:
    if (geom) {
      SetJetEtaRange(geom->GetArm1EtaMin() + r, geom->GetArm1EtaMax() - r);
      SetJetPhiRange(geom->GetDCALPhiMin() * TMath::DegToRad() + r, geom->GetDCALPhiMax() * TMath::DegToRad() - r);
    }
    else {
      AliWarning("Could not get instance of AliEMCALGeometry. Using manual settings for DCAL year 2015!!");
      SetJetEtaRange(-0.7 + r, 0.7 - r);
      SetJetPhiRange(4.538 + r, 5.727 - r);
    }
    break;

  case AliJetContainer::kUser:
    // Nothing to be done
    break;
  }
}

/// Compares 2 jet definitions.
/// The ordering is based on: jet type, radius, algorithm and recombination scheme, in this order
///
/// \param lhs Reference to the first AliJetDefinition object
/// \param rhs Reference to the second AliJetDefinition object
bool operator<(const AliAnalysisTaskDmesonJets::AliHFJetDefinition& lhs, const AliAnalysisTaskDmesonJets::AliHFJetDefinition& rhs)
{
  if (lhs.fJetType > rhs.fJetType) return false;
  else if (lhs.fJetType < rhs.fJetType) return true;
  else {
    if (lhs.fRadius > rhs.fRadius) return false;
    else if (lhs.fRadius < rhs.fRadius) return true;
    else {
      if (lhs.fJetAlgo > rhs.fJetAlgo) return false;
      else if (lhs.fJetAlgo < rhs.fJetAlgo) return true;
      else {
        if (lhs.fRecoScheme < rhs.fRecoScheme) return true;
        else return false;
      }
    }
  }
}

/// Compares 2 jet definitions.
/// Two analysis engines are considerate equal if they are exactly the same
///
/// \param lhs Reference to the first AliJetDefinition object
/// \param rhs Reference to the second AliJetDefinition object
bool operator==(const AliAnalysisTaskDmesonJets::AliHFJetDefinition& lhs, const AliAnalysisTaskDmesonJets::AliHFJetDefinition& rhs)
{
  if (lhs.fJetType != rhs.fJetType) return false;
  if (lhs.fRadius != rhs.fRadius) return false;
  if (lhs.fJetAlgo != rhs.fJetAlgo) return false;
  if (lhs.fRecoScheme != rhs.fRecoScheme) return false;
  return true;
}

// Definitions of class AliAnalysisTaskDmesonJets::AnalysisEngine

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJets::AnalysisEngine);
/// \endcond

/// This is the default constructor, used for ROOT I/O purposes.
AliAnalysisTaskDmesonJets::AnalysisEngine::AnalysisEngine() :
  TObject(),
  fCandidateType(kD0toKpi),
  fCandidateName(),
  fCandidatePDG(0),
  fNDaughters(0),
  fPDGdaughters(),
  fBranchName(),
  fMCMode(kNoMC),
  fNMassBins(0),
  fMinMass(0),
  fMaxMass(0),
  fRDHFCuts(0),
  fRejectedOrigin(0),
  fAcceptedDecay(0),
  fInhibit(kFALSE),
  fJetDefinitions(),
  fPtBinWidth(0.5),
  fMaxPt(100),
  fDataSlotNumber(-1),
  fTree(0),
  fCurrentDmesonJetInfo(0),
  fCurrentJetInfo(0),
  fCandidateArray(0),
  fMCContainer(0),
  fTrackContainer(0),
  fClusterContainer(0),
  fAodEvent(0),
  fFastJetWrapper(0),
  fHistManager(0)
{
}

/// This is the standard constructor.
///
/// \param type      One of the enum constants of ECandidateType_t
/// \param bkgMode   One of the enum constants of EMCMode_t
/// \param cuts      D meson cuts (if null, it will use standard cuts)
/// \param nMassBins Number of bins in the mass axis
/// \param range     Range of the mass axis (will be centered around the PDG mass)
AliAnalysisTaskDmesonJets::AnalysisEngine::AnalysisEngine(ECandidateType_t type, EMCMode_t MCmode, AliRDHFCuts* cuts, Int_t nMassBins, Double_t range) :
  TObject(),
  fCandidateType(type),
  fCandidateName(),
  fCandidatePDG(0),
  fNDaughters(0),
  fPDGdaughters(),
  fBranchName(),
  fMCMode(MCmode),
  fNMassBins(nMassBins),
  fMinMass(0),
  fMaxMass(0),
  fRDHFCuts(cuts),
  fRejectedOrigin(kUnknownQuark | kFromBottom),
  fAcceptedDecay(kAnyDecay),
  fInhibit(kFALSE),
  fJetDefinitions(),
  fPtBinWidth(0.5),
  fMaxPt(100),
  fDataSlotNumber(-1),
  fTree(0),
  fCurrentDmesonJetInfo(0),
  fCurrentJetInfo(0),
  fCandidateArray(0),
  fMCContainer(0),
  fTrackContainer(0),
  fClusterContainer(0),
  fAodEvent(0),
  fFastJetWrapper(0),
  fHistManager(0)
{
  SetCandidateProperties(range);
}

/// Copy constructor
///
/// \param source Reference to a valid AnalysisEngine to copy from.
AliAnalysisTaskDmesonJets::AnalysisEngine::AnalysisEngine(const AliAnalysisTaskDmesonJets::AnalysisEngine &source) :
  TObject(source),
  fCandidateType(source.fCandidateType),
  fCandidateName(source.fCandidateName),
  fCandidatePDG(source.fCandidatePDG),
  fNDaughters(source.fNDaughters),
  fPDGdaughters(source.fPDGdaughters),
  fBranchName(source.fBranchName),
  fMCMode(source.fMCMode),
  fNMassBins(source.fNMassBins),
  fMinMass(source.fMinMass),
  fMaxMass(source.fMaxMass),
  fRDHFCuts(),
  fRejectedOrigin(source.fRejectedOrigin),
  fAcceptedDecay(source.fAcceptedDecay),
  fInhibit(source.fInhibit),
  fJetDefinitions(source.fJetDefinitions),
  fPtBinWidth(source.fPtBinWidth),
  fMaxPt(source.fMaxPt),
  fDataSlotNumber(-1),
  fTree(0),
  fCurrentDmesonJetInfo(0),
  fCurrentJetInfo(0),
  fCandidateArray(source.fCandidateArray),
  fMCContainer(source.fMCContainer),
  fTrackContainer(source.fTrackContainer),
  fClusterContainer(source.fClusterContainer),
  fAodEvent(source.fAodEvent),
  fFastJetWrapper(source.fFastJetWrapper),
  fHistManager(source.fHistManager)
{
  SetRDHFCuts(source.fRDHFCuts);
}

// Destructor
AliAnalysisTaskDmesonJets::AnalysisEngine::~AnalysisEngine()
{
  delete fRDHFCuts;
}

/// Assignement operator
///
/// \param source Reference to a valid AnalysisEngine to copy from.
AliAnalysisTaskDmesonJets::AnalysisEngine& AliAnalysisTaskDmesonJets::AnalysisEngine::operator=(const AnalysisEngine& source)
{
  new (this) AnalysisEngine(source);
  return *this;
}

/// Checks whether any of the D meson jets is in the acceptance
///
/// \param Const reference to a valid AliDmesonJetInfo object
Bool_t AliAnalysisTaskDmesonJets::AnalysisEngine::IsAnyJetInAcceptance(const AliDmesonJetInfo& dMesonJet) const
{
  for (UInt_t i = 0; i < fJetDefinitions.size(); i++) {
    if (fJetDefinitions[i].IsJetInAcceptance(dMesonJet, fJetDefinitions[i].GetName())) return kTRUE;
  }

  return kFALSE;
}

/// Initialize the analysis engine
void AliAnalysisTaskDmesonJets::AnalysisEngine::Init(const AliEMCALGeometry* const geom, Int_t runNumber)
{
  for (Int_t i = 0; i < fJetDefinitions.size(); i++) {
    fJetDefinitions[i].SetDetectorJetEtaPhiRange(geom, runNumber);
  }
}

/// Sets the D meson candidate properties.
///
/// \param range     Range of the mass axis (will be centered around the PDG mass)
void AliAnalysisTaskDmesonJets::AnalysisEngine::SetCandidateProperties(Double_t range)
{
  switch (fCandidateType) {
  case kD0toKpi:
    fCandidatePDG = 421;
    fCandidateName = "D0";
    fNDaughters = 2;
    fPDGdaughters.Set(fNDaughters);
    fPDGdaughters.Reset();
    fPDGdaughters[0] = 211;  // pi
    fPDGdaughters[1] = 321;  // K
    fBranchName = "D0toKpi";
    fAcceptedDecay = kD0toKpi;
    if (!fRDHFCuts) {
      fRDHFCuts = new AliRDHFCutsD0toKpi();
      fRDHFCuts->SetStandardCutsPP2010();
      fRDHFCuts->SetUsePhysicsSelection(kFALSE);
      fRDHFCuts->SetTriggerClass("","");
    }
    break;
  case kDstartoKpipi:
    fCandidatePDG = 413;
    fCandidateName = "DStar";
    fNDaughters = 3;
    fPDGdaughters.Set(fNDaughters);
    fPDGdaughters.Reset();
    fPDGdaughters[0] = 211; // pi soft
    fPDGdaughters[1] = 211; // pi fromD0
    fPDGdaughters[2] = 321; // K from D0
    fBranchName = "Dstar";
    fAcceptedDecay = kDstartoKpipi;
    if (!fRDHFCuts) {
      fRDHFCuts = new AliRDHFCutsDStartoKpipi();
      fRDHFCuts->SetStandardCutsPP2010();
      fRDHFCuts->SetUsePhysicsSelection(kFALSE);
      fRDHFCuts->SetTriggerClass("","");
    }
    break;
  default:
    ::Error("AliAnalysisTaskDmesonJets::AnalysisEngine::SetCandidateProperties","Candidate %d unknown!", fCandidateType);
  }

  CalculateMassLimits(range, fCandidatePDG, fNMassBins, fMinMass, fMaxMass);
}

/// Adopt the cuts (this class owns the cuts object, which will be destroyed when needed).
///
/// \param Pointer to a AliRDHFCuts object.
void AliAnalysisTaskDmesonJets::AnalysisEngine::AdoptRDHFCuts(AliRDHFCuts* cuts)
{
  if (fRDHFCuts) delete fRDHFCuts;
  fRDHFCuts = cuts;
}

/// Set the cuts (creates a copy, so the original object is not owned by this class).
///
/// \param Pointer to a AliRDHFCuts object.
void AliAnalysisTaskDmesonJets::AnalysisEngine::SetRDHFCuts(AliRDHFCuts* cuts)
{
  if (!cuts) return;
  if (fRDHFCuts) delete fRDHFCuts;
  fRDHFCuts = static_cast<AliRDHFCuts*>(cuts->Clone());
}

/// Generate a name for this analysis parameter set
///
/// \param i  Index of the jet radius array.
const char* AliAnalysisTaskDmesonJets::AnalysisEngine::GetName(const AliHFJetDefinition& jetDef) const
{
  static TString name;

  name = TString::Format("%s_%s", GetName(), jetDef.GetName());

  return name.Data();
}

/// Generate a name for this analysis parameter set
///
/// \param i  Index of the jet radius array.
const char* AliAnalysisTaskDmesonJets::AnalysisEngine::GetName() const
{
  static TString name;

  name = fCandidateName;
  switch (fMCMode) {
  case kBackgroundOnly:
    name += "_kBackgroundOnly";
    break;
  case kSignalOnly:
    name += "_kSignalOnly";
    break;
  case kMCTruth:
    name += "_MCTruth";
    break;
  default:
    break;
  }

  return name.Data();
}

/// Add a new jet definition
/// If the jet definition is already present, it does nothing.
///
/// \param def Reference to a AliJetDefinition object
///
/// \return Pointer to the new jet definition (or to the one that was already present)
AliAnalysisTaskDmesonJets::AliHFJetDefinition* AliAnalysisTaskDmesonJets::AnalysisEngine::AddJetDefinition(const AliAnalysisTaskDmesonJets::AliHFJetDefinition& def)
{
  std::vector<AliHFJetDefinition>::iterator it = FindJetDefinition(def);

  if (it == fJetDefinitions.end() || *it != def) {  // No jet definition was found, adding a new one
    fJetDefinitions.push_back(def);
    ::Info("AliAnalysisTaskDmesonJets::AnalysisEngine::AddJetDefinition", "Jet definition '%s' has been added to analysis engine '%s'."
        "Total number of jet definitions is now %lu.",
        def.GetName(), GetName(), fJetDefinitions.size());
    // For detector level set maximum track pt to 100 GeV/c
    if (fMCMode != kMCTruth) fJetDefinitions[fJetDefinitions.size()-1].SetChargedPtRange(0., 100.);
  }
  else {
    ::Warning("AliAnalysisTaskDmesonJets::AnalysisEngine::AddJetDefinition", "The same jet definition '%s' was already added in analysis engine '%s'.", def.GetName(), GetName());
  }

  return &(*it);
}

/// Add a new jet definition
/// If the jet definition is already present, it does nothing.
///
/// \param type Jet type
/// \param r    Jet radius
/// \param algo Jet algorithm
/// \param reco Recombination scheme
///
/// \return Pointer to the new jet definition (or to the one that was already present)
AliAnalysisTaskDmesonJets::AliHFJetDefinition*
AliAnalysisTaskDmesonJets::AnalysisEngine::AddJetDefinition(EJetType_t type, Double_t r, EJetAlgo_t algo, ERecoScheme_t reco)
{
  AliHFJetDefinition def(type, r, algo, reco);

  return AddJetDefinition(def);
}

/// Look for a jet definition that is equal
///
/// \param def Reference to a jet definition object
///
/// \return An iterator to the jet definition object, if it is found. An iterator to the end if not found.
std::vector<AliAnalysisTaskDmesonJets::AliHFJetDefinition>::iterator AliAnalysisTaskDmesonJets::AnalysisEngine::FindJetDefinition(const AliAnalysisTaskDmesonJets::AliHFJetDefinition& def)
{
  std::vector<AliHFJetDefinition>::iterator it = fJetDefinitions.begin();
  while (it != fJetDefinitions.end() && (*it) != def) it++;
  return it;
}

/// Compares 2 analysis engines.
/// The ordering is based on the candidate type first and then on the MC mode.
///
/// \param lhs Reference to the first AnalysisEngine object
/// \param rhs Reference to the second AnalysisEngine object
bool operator<(const AliAnalysisTaskDmesonJets::AnalysisEngine& lhs, const AliAnalysisTaskDmesonJets::AnalysisEngine& rhs)
{
  if (lhs.fCandidateType > rhs.fCandidateType) return false;
  else if (lhs.fCandidateType < rhs.fCandidateType) return true;
  else {
    if (lhs.fMCMode < rhs.fMCMode) return true;
    else return false;
  }
}

/// Compares 2 analysis engines.
/// Two analysis engines are considerate equal if they have both the same candidate type and MC mode.
///
/// \param lhs Reference to the first AnalysisEngine object
/// \param rhs Reference to the second AnalysisEngine object
bool operator==(const AliAnalysisTaskDmesonJets::AnalysisEngine& lhs, const AliAnalysisTaskDmesonJets::AnalysisEngine& rhs)
{
  if (lhs.fCandidateType != rhs.fCandidateType) return false;
  if (lhs.fMCMode != rhs.fMCMode) return false;
  return true;
}

/// Extract attributes of the D meson (particle level).
///
/// \param part Pointer to a AliAODMCParticle representing the D meson
/// \param DmesonJet Reference to an AliDmesonJetInfo object where the D meson information will be copied
/// \param i Either 0 or 1, for the two possible mass hypothesis assignment (since it is particle level it will return kFALSE for i > 0)
///
/// \return Always kTRUE on success
Bool_t AliAnalysisTaskDmesonJets::AnalysisEngine::ExtractParticleLevelHFAttributes(const AliAODMCParticle* part, AliDmesonJetInfo& DmesonJet)
{
  DmesonJet.fD.SetPtEtaPhiM(part->Pt(), part->Eta(), part->Phi(), part->M());
  return kTRUE;
}

/// Extract attributes of the D meson candidate.
///
/// \param Dcand Pointer to a AliAODRecoDecayHF2Prong representing the D meson candidate
/// \param DmesonJet Reference to an AliDmesonJetInfo object where the D meson candidate information will be copied
/// \param i Either 0 or 1, for the two possible mass hypothesis assignments
///
/// \return kTRUE on success
Bool_t AliAnalysisTaskDmesonJets::AnalysisEngine::ExtractRecoDecayAttributes(const AliAODRecoDecayHF2Prong* Dcand, AliDmesonJetInfo& DmesonJet, UInt_t i)
{
  if (fCandidateType == kD0toKpi) { // D0 candidate
    return ExtractD0Attributes(Dcand, DmesonJet, i);
  }
  else if (fCandidateType == kDstartoKpipi) { // Dstar candidate
    return ExtractDstarAttributes(static_cast<const AliAODRecoCascadeHF*>(Dcand), DmesonJet, i);
  }
  else {
    return kFALSE;
  }
}

/// Extract attributes of the D0 meson candidate.
///
/// \param Dcand Pointer to a AliAODRecoDecayHF2Prong representing the D0 meson candidate
/// \param DmesonJet Reference to an AliDmesonJetInfo object where the D0 meson candidate information will be copied
/// \param i Either 0 or 1, for the two possible mass hypothesis assignments
///
/// \return kTRUE on success
Bool_t AliAnalysisTaskDmesonJets::AnalysisEngine::ExtractD0Attributes(const AliAODRecoDecayHF2Prong* Dcand, AliDmesonJetInfo& DmesonJet, UInt_t i)
{
  Int_t MCtruthPdgCode = 0;

  Double_t invMassD = 0;

  if (fMCMode == kBackgroundOnly || fMCMode == kSignalOnly) {
    Int_t mcLab = Dcand->MatchToMC(fCandidatePDG, fMCContainer->GetArray(), fNDaughters, fPDGdaughters.GetArray());
    if (mcLab >= 0) {
      AliAODMCParticle* aodMcPart = static_cast<AliAODMCParticle*>(fMCContainer->GetArray()->At(mcLab));

      if (aodMcPart) {
        if (fRejectedOrigin && fMCMode == kSignalOnly) {
          EMesonOrigin_t origin = CheckOrigin(aodMcPart, fMCContainer->GetArray());

          if ((origin & fRejectedOrigin) == origin) return kFALSE;
        }
        MCtruthPdgCode = aodMcPart->PdgCode();
      }
    }
  }

  //AliDebug(2,"Checking if D0 meson is selected");
  Int_t isSelected = fRDHFCuts->IsSelected(const_cast<AliAODRecoDecayHF2Prong*>(Dcand), AliRDHFCuts::kAll, fAodEvent);
  if (isSelected == 1) { // selected as a D0
    if (i > 0) return kFALSE; // only one mass hypothesis thanks to PID

    if (fMCMode == kNoMC ||
        (MCtruthPdgCode == fCandidatePDG && fMCMode == kSignalOnly) ||
        (MCtruthPdgCode != fCandidatePDG && fMCMode == kBackgroundOnly)) {
      // both background and signal are requested OR (it is a true D0 AND signal is requested) OR (it is NOT a D0 and background is requested)
      //AliDebug(2,"Selected as D0");
      invMassD = Dcand->InvMassD0();
    }
    else { // conditions above not passed, so return FALSE
      return kFALSE;
    }
  }
  else if (isSelected == 2) { // selected as a D0bar
    if (i > 0) return kFALSE; // only one mass hypothesis thanks to PID

    if (fMCMode == kNoMC ||
        (MCtruthPdgCode == -fCandidatePDG && fMCMode == kSignalOnly) ||
        (MCtruthPdgCode != -fCandidatePDG && fMCMode == kBackgroundOnly)) {
      // both background and signal are requested OR (it is a true D0bar AND signal is requested) OR (it is NOT a D0bar and background is requested)
      //AliDebug(2,"Selected as D0bar");
      invMassD = Dcand->InvMassD0bar();
    }
    else { // conditions above not passed, so return FALSE
      return kFALSE;
    }
  }
  else if (isSelected == 3) { // selected as either a D0bar or a D0 (PID on K and pi undecisive)
    //AliDebug(2,"Selected as either D0 or D0bar");

    // Accept the correct mass hypothesis for signal-only and the wrong one for background-only
    if ((MCtruthPdgCode == fCandidatePDG && fMCMode == kSignalOnly) ||
        (MCtruthPdgCode == -fCandidatePDG && fMCMode == kBackgroundOnly)) {
      if (i > 0) return kFALSE;
      //AliDebug(2, "MC truth is D0");
      invMassD = Dcand->InvMassD0();
    }
    else if ((MCtruthPdgCode == -fCandidatePDG && fMCMode == kSignalOnly) ||
             (MCtruthPdgCode == fCandidatePDG && fMCMode == kBackgroundOnly)) {
      if (i > 0) return kFALSE;
      //AliDebug(2, "MC truth is D0bar");
      invMassD = Dcand->InvMassD0bar();
    }
    else { // (This candidate is neither a D0 nor a D0bar) OR (background-and-signal was requested)

      // Only accept it if background-only OR background-and-signal was requested
      if (fMCMode == kBackgroundOnly || fMCMode == kNoMC) {
        // Select D0 or D0bar depending on the i-parameter
        if (i == 0) {
          //AliDebug(2, "Returning invariant mass with D0 hypothesis");
          invMassD = Dcand->InvMassD0();
        }
        else if (i == 1) {
          //AliDebug(2, "Returning invariant mass with D0bar hypothesis");
          invMassD = Dcand->InvMassD0bar();
        }
        else {  // i > 1
          return kFALSE;
        }
      }
      else { // signal-only was requested but this is not a true D0
        return kFALSE;
      }
    }
  }
  DmesonJet.fD.SetPtEtaPhiM(Dcand->Pt(), Dcand->Eta(), Dcand->Phi(), invMassD);
  return kTRUE;
}

/// Extract attributes of the D* meson candidate.
///
/// \param DstarCand Pointer to a AliAODRecoCascadeHF representing the D* meson candidate
/// \param DmesonJet Reference to an AliDmesonJetInfo object where the D* meson candidate information will be copied
/// \param i Either 0 or 1, for the two possible mass hypothesis assignments (since there is only one mass hypothesis for D*, returns kFALSE for i > 0)
///
/// \return kTRUE on success
Bool_t AliAnalysisTaskDmesonJets::AnalysisEngine::ExtractDstarAttributes(const AliAODRecoCascadeHF* DstarCand, AliDmesonJetInfo& DmesonJet, UInt_t i)
{
  if (i > 0) return kFALSE; // only one mass hypothesis for the D*

  Int_t MCtruthPdgCode = 0;

  Double_t invMassD = 0;

  if (fMCMode == kBackgroundOnly || fMCMode == kSignalOnly) {
    Int_t pdgDgDStartoD0pi[2] = { 421, 211 };  // D0,pi
    Int_t pdgDgD0toKpi[2] = { 321, 211 };      // K, pi

    Int_t mcLab = DstarCand->MatchToMC(fCandidatePDG, 421, pdgDgDStartoD0pi, pdgDgD0toKpi, fMCContainer->GetArray());
    //AliDebug(2, Form("MC label is %d", mcLab));
    if (mcLab >= 0) {
      AliAODMCParticle* aodMcPart = static_cast<AliAODMCParticle*>(fMCContainer->GetArray()->At(mcLab));

      if (aodMcPart) {
        if (fRejectedOrigin && fMCMode == kSignalOnly) {
          EMesonOrigin_t origin = CheckOrigin(aodMcPart, fMCContainer->GetArray());

          if ((origin & fRejectedOrigin) == origin) return kFALSE;
        }

        MCtruthPdgCode = aodMcPart->PdgCode();
        //AliDebug(2, Form("MC truth pdg code is %d",MCtruthPdgCode));
      }
    }
  }

  Int_t absMCtruthPdgCode = TMath::Abs(MCtruthPdgCode);
  if (fMCMode == kNoMC ||
      (absMCtruthPdgCode == 413 && fMCMode == kSignalOnly) ||
      (absMCtruthPdgCode != 413 && fMCMode == kBackgroundOnly)) {
    // both background and signal are requested OR (it is a true D*/D*bar AND signal is requested) OR (it is NOT a D*/D*bar and background is requested)
    invMassD = DstarCand->InvMassDstarKpipi();
    DmesonJet.fSoftPionPt = DstarCand->GetBachelor()->Pt();
    DmesonJet.fInvMass2Prong = DstarCand->InvMassD0();
    DmesonJet.fD.SetPtEtaPhiM(DstarCand->Pt(), DstarCand->Eta(), DstarCand->Phi(), invMassD);
    return kTRUE;
  }
  else { // conditions above not passed, so return FALSE
    return kFALSE;
  }
}

/// Checks the decay channel of a D meson
///
/// \param part Pointer to an AliAODMCParticle object for which decay channel is requested
/// \param mcArray Pointer to a TClonesArray object where to look for particles
///
/// \return One of the enum constants of AliAnalysisTaskDmesonJets::EMesonDecayChannel_t (D0->Kpi or D*->D0pi->Kpipi)
AliAnalysisTaskDmesonJets::EMesonDecayChannel_t AliAnalysisTaskDmesonJets::AnalysisEngine::CheckDecayChannel(const AliAODMCParticle* part, TClonesArray* mcArray)
{
  if (!part) return kDecayOther;
  if (!mcArray) return kDecayOther;

  EMesonDecayChannel_t decay = kDecayOther;

  Int_t absPdgPart = TMath::Abs(part->GetPdgCode());

  if (part->GetNDaughters() == 2) {

    AliAODMCParticle* d1 = static_cast<AliAODMCParticle*>(mcArray->At(part->GetDaughter(0)));
    AliAODMCParticle* d2 = static_cast<AliAODMCParticle*>(mcArray->At(part->GetDaughter(1)));

    if (!d1 || !d2) {
      return decay;
    }

    Int_t absPdg1 = TMath::Abs(d1->GetPdgCode());
    Int_t absPdg2 = TMath::Abs(d2->GetPdgCode());

    if (absPdgPart == 421) { // D0 -> K pi

      if ((absPdg1 == 211 && absPdg2 == 321) || // pi K
          (absPdg1 == 321 && absPdg2 == 211)) { // K pi
        decay = kDecayD0toKpi;
      }
    }

    if (absPdgPart == 413) { // D* -> D0 pi

      if (absPdg1 == 421 && absPdg2 == 211) {  // D0 pi
        Int_t D0decay = CheckDecayChannel(d1, mcArray);
        if (D0decay == kDecayD0toKpi) {
          decay = kDecayDStartoKpipi;
        }
      }

      if (absPdg1 == 211 && absPdg2 == 421) {  // pi D0
        Int_t D0decay = CheckDecayChannel(d2, mcArray);
        if (D0decay == kDecayD0toKpi) {
          decay = kDecayDStartoKpipi;
        }
      }
    }
  }

  return decay;
}

/// Checks the origin of a D meson
///
/// \param part Pointer to an AliAODMCParticle object for which originating quark is required
/// \param mcArray Pointer to a TClonesArray object where to look for particles
///
/// \return One of the enum constants of AliAnalysisTaskDmesonJets::EMesonOrigin_t (unknown quark, bottom or charm)
AliAnalysisTaskDmesonJets::EMesonOrigin_t AliAnalysisTaskDmesonJets::AnalysisEngine::CheckOrigin(const AliAODMCParticle* part, TClonesArray* mcArray)
{
  // Checks whether the mother of the particle comes from a charm or a bottom quark.

  if (!part) return kUnknownQuark;
  if (!mcArray) return kUnknownQuark;

  Int_t pdgGranma = 0;
  Int_t mother = part->GetMother();
  Int_t istep = 0;
  Int_t abspdgGranma = 0;
  Bool_t isFromB = kFALSE;
  Bool_t isQuarkFound = kFALSE;

  while (mother >= 0) {
    istep++;
    AliAODMCParticle* mcGranma = static_cast<AliAODMCParticle*>(mcArray->At(mother));
    if (mcGranma >= 0) {
      pdgGranma = mcGranma->GetPdgCode();
      abspdgGranma = TMath::Abs(pdgGranma);
      if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)) {
        isFromB = kTRUE;
      }

      if (abspdgGranma == 4 || abspdgGranma == 5) isQuarkFound = kTRUE;
      mother = mcGranma->GetMother();
    }
    else {
      ::Error("AliAnalysisTaskDmesonJets::AnalysisParams::CheckOrigin", "Could not retrieve mother particle %d!", mother);
      break;
    }
  }

  if (isQuarkFound) {
    if (isFromB) {
      return kFromBottom;
    }
    else {
      return kFromCharm;
    }
  }
  else {
    return kUnknownQuark;
  }
}

/// Run the analysis
void AliAnalysisTaskDmesonJets::AnalysisEngine::RunAnalysis()
{
  fDmesonJets.clear();

  if (fMCMode == kMCTruth) {
    RunParticleLevelAnalysis();
  }
  else {
    RunDetectorLevelAnalysis();
  }
}

/// Run a detector level analysis
void AliAnalysisTaskDmesonJets::AnalysisEngine::RunDetectorLevelAnalysis()
{
  const Int_t nD = fCandidateArray->GetEntriesFast();

  AliDmesonJetInfo DmesonJet;

  Int_t nAccCharm = 0;
  for (Int_t icharm = 0; icharm < nD; icharm++) {   //loop over D candidates
    Int_t isSelected = 0;

    AliAODRecoDecayHF2Prong* charmCand = static_cast<AliAODRecoDecayHF2Prong*>(fCandidateArray->At(icharm)); // D candidates
    if (!charmCand) continue;

    Int_t nprongs = charmCand->GetNProngs();

    if (fCandidateType == kDstartoKpipi) {
      if (!charmCand->InheritsFrom("AliAODRecoCascadeHF")) {
        ::Error("AliAnalysisTaskDmesonJets::AnalysisParams::RunDetectorLevelAnalysis","Candidate type is D* but object type is wrong (should be AliAODRecoCascadeHF)");
        continue;
      }
    }

    // region of interest + cuts
    if (!fRDHFCuts->IsInFiducialAcceptance(charmCand->Pt(), charmCand->Y(fCandidatePDG))) continue;

    //candidate selected by cuts and PID
    isSelected = fRDHFCuts->IsSelected(charmCand, AliRDHFCuts::kAll, fAodEvent); //selected

    if (!isSelected) continue;

    for (Int_t im = 0; im < 2; im++)  {  // 2 mass hypothesis (when available)
      DmesonJet.Reset();
      if (ExtractRecoDecayAttributes(charmCand, DmesonJet, im)) {
        for (std::vector<AliHFJetDefinition>::iterator itdef = fJetDefinitions.begin(); itdef != fJetDefinitions.end(); itdef++) {
          if (!FindJet(charmCand, DmesonJet, *itdef)) {
            AliWarning(Form("Could not find jet '%s' for D meson '%s': pT = %.3f, eta = %.3f, phi = %.3f",
                (*itdef).GetName(), GetName(), DmesonJet.fD.Pt(), DmesonJet.fD.Eta(), DmesonJet.fD.Phi_0_2pi()));
          }
        }
        fDmesonJets.push_back(DmesonJet);
      }
    }
    nAccCharm++;
  } // end of D cand loop

  TString hname;

  hname = TString::Format("%s/fHistNAcceptedDmesons", GetName());
  fHistManager->FillTH1(hname, nAccCharm);

  hname = TString::Format("%s/fHistNDmesons", GetName());
  fHistManager->FillTH1(hname, nD);
}

/// Find the jet that contains a D meson candidate.
/// The jet finding algorithm is always anti-kt
/// Tracks and clusters are accessed through fTrackContainer and fClusterContainer
///
/// \param Dcand Valid pointer to a D meson candidate object
/// \param DmesonJet Reference to a AliDmesonJetInfo object where the result will be stored
/// \param r Jet radius
///
/// \return kTRUE on success, kFALSE otherwise
Bool_t AliAnalysisTaskDmesonJets::AnalysisEngine::FindJet(AliAODRecoDecayHF2Prong* Dcand, AliDmesonJetInfo& DmesonJet, AliHFJetDefinition& jetDef)
{
  TString hname;

  fFastJetWrapper->Clear();
  fFastJetWrapper->SetR(jetDef.fRadius);
  fFastJetWrapper->SetAlgorithm(AliEmcalJetTask::ConvertToFJAlgo(jetDef.fJetAlgo));
  fFastJetWrapper->SetRecombScheme(AliEmcalJetTask::ConvertToFJRecoScheme(jetDef.fRecoScheme));

  fFastJetWrapper->AddInputVector(DmesonJet.fD.Px(), DmesonJet.fD.Py(), DmesonJet.fD.Pz(), DmesonJet.fD.E(), 0);

  if (fTrackContainer && jetDef.fJetType != AliJetContainer::kNeutralJet) {
    fTrackContainer->SetDMesonCandidate(Dcand);
    hname = TString::Format("%s/%s/fHistTrackRejectionReason", GetName(), jetDef.GetName());
    AddInputVectors(fTrackContainer, 100, static_cast<TH2*>(fHistManager->FindObject(hname)));

    hname = TString::Format("%s/%s/fHistDMesonDaughterNotInJet", GetName(), jetDef.GetName());
    TH1* histDaughterNotInJet = static_cast<TH1*>(fHistManager->FindObject(hname));
    const TObjArray daughterNotInJet = fTrackContainer->GetDaughterList();
    for (Int_t i = 0; i < daughterNotInJet.GetEntriesFast(); i++) {
      AliVParticle* daughter = static_cast<AliVParticle*>(daughterNotInJet.At(i));
      if (!daughter) continue;
      histDaughterNotInJet->Fill(daughter->Pt());
    }
  }

  if (fClusterContainer && jetDef.fJetType != AliJetContainer::kChargedJet) {
    hname = TString::Format("%s/%s/fHistClusterRejectionReason", GetName(), jetDef.GetName());
    AddInputVectors(fClusterContainer, -100, static_cast<TH2*>(fHistManager->FindObject(hname)));
  }

  // run jet finder
  fFastJetWrapper->Run();

  std::vector<fastjet::PseudoJet> jets_incl = fFastJetWrapper->GetInclusiveJets();

  for (UInt_t ijet = 0; ijet < jets_incl.size(); ++ijet) {
    std::vector<fastjet::PseudoJet> constituents(fFastJetWrapper->GetJetConstituents(ijet));

    Bool_t isDmesonJet = kFALSE;

    Double_t maxChPt = 0;
    Double_t maxNePt = 0;
    Double_t totalNeutralPt = 0;

    for (UInt_t ic = 0; ic < constituents.size(); ++ic) {
      if (constituents[ic].user_index() == 0) {
        isDmesonJet = kTRUE;
      }
      else if (constituents[ic].user_index() >= 100) {
        if (constituents[ic].pt() > maxChPt) maxChPt = constituents[ic].pt();
      }
      else if (constituents[ic].user_index() <= -100) {
        totalNeutralPt += constituents[ic].pt();
        if (constituents[ic].pt() > maxNePt) maxChPt = constituents[ic].pt();
      }
    }

    if (isDmesonJet) {
      DmesonJet.fJets[jetDef.GetName()].fMomentum.SetPxPyPzE(jets_incl[ijet].px(), jets_incl[ijet].py(), jets_incl[ijet].pz(), jets_incl[ijet].E());
      DmesonJet.fJets[jetDef.GetName()].fNConstituents = constituents.size();
      DmesonJet.fJets[jetDef.GetName()].fMaxChargedPt = maxChPt;
      DmesonJet.fJets[jetDef.GetName()].fMaxNeutralPt = maxNePt;
      DmesonJet.fJets[jetDef.GetName()].fNEF = totalNeutralPt / jets_incl[ijet].pt();

      return kTRUE;
    }
  }

  return kFALSE;
}

/// Adds all the particles contained in the container into the fastjet wrapper
///
/// \param cont Pointer to a valid AliEmcalContainer object
void AliAnalysisTaskDmesonJets::AnalysisEngine::AddInputVectors(AliEmcalContainer* cont, Int_t offset, TH2* rejectHist)
{
  AliEmcalIterableMomentumContainer itcont = cont->all_momentum();
  for (AliEmcalIterableMomentumContainer::iterator it = itcont.begin(); it != itcont.end(); it++) {
    UInt_t rejectionReason = 0;
    if (!cont->AcceptObject(it.current_index(), rejectionReason)) {
      rejectHist->Fill(cont->GetRejectionReasonBitPosition(rejectionReason), it->first.Pt());
      continue;
    }
    Int_t uid = offset >= 0 ? it.current_index() + offset: -it.current_index() - offset;
    fFastJetWrapper->AddInputVector(it->first.Px(), it->first.Py(), it->first.Pz(), it->first.E(), uid);
  }
}

/// Run a particle level analysis
void AliAnalysisTaskDmesonJets::AnalysisEngine::RunParticleLevelAnalysis()
{
  TString hname;

  fMCContainer->SetSpecialPDG(fCandidatePDG);
  fMCContainer->SetRejectedOriginMap(fRejectedOrigin);
  fMCContainer->SetAcceptedDecayMap(fAcceptedDecay);

  std::map<int, AliDmesonJetInfo> dMesonJets;

  for (auto &jetDef : fJetDefinitions) {

    fFastJetWrapper->Clear();
    fFastJetWrapper->SetR(jetDef.fRadius);
    fFastJetWrapper->SetAlgorithm(AliEmcalJetTask::ConvertToFJAlgo(jetDef.fJetAlgo));
    fFastJetWrapper->SetRecombScheme(AliEmcalJetTask::ConvertToFJRecoScheme(jetDef.fRecoScheme));

    hname = TString::Format("%s/%s/fHistMCParticleRejectionReason", GetName(), jetDef.GetName());
    AddInputVectors(fMCContainer, 100, static_cast<TH2*>(fHistManager->FindObject(hname)));

    fFastJetWrapper->Run();

    std::vector<fastjet::PseudoJet> jets_incl = fFastJetWrapper->GetInclusiveJets();

    for (UInt_t ijet = 0; ijet < jets_incl.size(); ++ijet) {
      std::vector<fastjet::PseudoJet> constituents(fFastJetWrapper->GetJetConstituents(ijet));

      Bool_t isDmesonJet = kFALSE;

      for (UInt_t ic = 0; ic < constituents.size(); ++ic) {
        Int_t iPart = constituents[ic].user_index() - 100;
        AliVParticle* part = fMCContainer->GetParticle(iPart);
        if (!part) {
          ::Error("AliAnalysisTaskDmesonJets::AnalysisEngine::RunParticleLevelAnalysis", "Could not find jet constituent %d!", iPart);
          continue;
        }
        if (TMath::Abs(part->PdgCode()) == fCandidatePDG) {
          std::map<int, AliDmesonJetInfo>::iterator dMesonJetIt = dMesonJets.find(iPart);
          if (dMesonJetIt == dMesonJets.end()) { // This D meson does not exist yet
            std::pair<int, AliDmesonJetInfo> element;
            element.first = iPart;
            element.second = AliDmesonJetInfo();
            dMesonJetIt = dMesonJets.insert(element).first;
            (*dMesonJetIt).second.fD.SetPxPyPzE(part->Px(), part->Py(), part->Pz(), part->E());
          }

          (*dMesonJetIt).second.fJets[jetDef.GetName()].fMomentum.SetPxPyPzE(jets_incl[ijet].px(), jets_incl[ijet].py(), jets_incl[ijet].pz(), jets_incl[ijet].E());
          (*dMesonJetIt).second.fJets[jetDef.GetName()].fNConstituents = constituents.size();
        }
      }
    }
  }

  for (std::map<int, AliDmesonJetInfo>::iterator dMesonJetIt = dMesonJets.begin();
      dMesonJetIt != dMesonJets.end();
      dMesonJetIt++) {
    fDmesonJets.push_back((*dMesonJetIt).second);
  }
}

/// Builds the tree where the output will be posted
///
/// \return Pointer to the new tree
TTree* AliAnalysisTaskDmesonJets::AnalysisEngine::BuildTree(const char* taskName)
{
  TString classname;
  switch (fCandidateType) {
  case kD0toKpi:
    classname = "AliAnalysisTaskDmesonJets::AliD0InfoSummary";
    fCurrentDmesonJetInfo = new AliD0InfoSummary();
    break;
  case kDstartoKpipi:
    classname = "AliAnalysisTaskDmesonJets::AliDStarInfoSummary";
    fCurrentDmesonJetInfo = new AliDStarInfoSummary();
    break;
  }
  TString treeName = TString::Format("%s_%s", taskName, GetName());
  fTree = new TTree(treeName, treeName);
  fTree->Branch("DmesonJet", classname, &fCurrentDmesonJetInfo, 32000, 0);
  fCurrentJetInfo = new AliJetInfoSummary*[fJetDefinitions.size()];
  for (Int_t i = 0; i < fJetDefinitions.size(); i++) {
    fCurrentJetInfo[i] = new AliJetInfoSummary();
    fTree->Branch(fJetDefinitions[i].GetName(), "AliAnalysisTaskDmesonJets::AliJetInfoSummary", &fCurrentJetInfo[i]);
  }

  return fTree;
}

/// Allocate a THnSparse histogram
///
/// \param param Analysis parameters used to properly set some of the axis
void AliAnalysisTaskDmesonJets::AnalysisEngine::BuildHnSparse(UInt_t enabledAxis)
{
  TString hname;

  Int_t nPtBins = TMath::CeilNint(fMaxPt / fPtBinWidth);

  for (auto &jetDef : fJetDefinitions) {

    AliDebug(2,Form("Now working on '%s'", jetDef.GetName()));

    Double_t radius = jetDef.fRadius;

    TString  title[30] = {""};
    Int_t    nbins[30] = {0 };
    Double_t min  [30] = {0.};
    Double_t max  [30] = {0.};
    Int_t    dim       = 0   ;

    title[dim] = "#it{p}_{T,D} (GeV/#it{c})";
    nbins[dim] = nPtBins;
    min[dim] = 0;
    max[dim] = fMaxPt;
    dim++;

    if ((enabledAxis & kPositionD) != 0) {
      title[dim] = "#eta_{D}";
      nbins[dim] = 50;
      min[dim] = -1;
      max[dim] = 1;
      dim++;

      title[dim] = "#phi_{D} (rad)";
      nbins[dim] = 150;
      min[dim] = 0;
      max[dim] = TMath::TwoPi();
      dim++;
    }

    if ((enabledAxis & kInvMass) != 0 && fCandidateType == kDstartoKpipi) {
      title[dim] = "#it{M}_{K#pi#pi} (GeV/#it{c}^{2})";
      nbins[dim] = fNMassBins;
      min[dim] = fMinMass;
      max[dim] = fMaxMass;
      dim++;
    }

    if (fCandidateType == kD0toKpi) {
      title[dim] = "#it{M}_{K#pi} (GeV/#it{c}^{2})";
      nbins[dim] = fNMassBins;
      min[dim] = fMinMass;
      max[dim] = fMaxMass;
      dim++;
    }

    if ((enabledAxis & k2ProngInvMass) != 0 && fCandidateType == kDstartoKpipi) {
      title[dim] = "#it{M}_{K#pi} (GeV/#it{c}^{2})";
      nbins[dim] = fNMassBins;
      CalculateMassLimits(fMaxMass - fMinMass, 421, fNMassBins, min[dim], max[dim]);
      dim++;
    }

    if (fCandidateType == kDstartoKpipi) {
      title[dim] = "#it{M}_{K#pi#pi} - #it{M}_{K#pi} (GeV/#it{c}^{2})";
      nbins[dim] = fNMassBins*6;
      CalculateMassLimits(0.20, 413, nbins[dim], min[dim], max[dim]);

      // subtract mass of D0
      Double_t D0mass = TDatabasePDG::Instance()->GetParticle(421)->Mass();
      min[dim] -= D0mass;
      max[dim] -= D0mass;

      dim++;
    }

    if ((enabledAxis & kSoftPionPt) != 0 && fCandidateType == kDstartoKpipi) {
      title[dim] = "#it{p}_{T,#pi} (GeV/#it{c})";
      nbins[dim] = 100;
      min[dim] = 0;
      max[dim] = 25;
      dim++;
    }

    title[dim] = "#it{z}_{D}";
    nbins[dim] = 110;
    min[dim] = 0;
    max[dim] = 1.10;
    dim++;

    if ((enabledAxis & kDeltaR) != 0) {
      title[dim] = "#Delta R_{D-jet}";
      nbins[dim] = 100;
      min[dim] = 0;
      max[dim] = radius * 1.5;
      dim++;
    }

    if ((enabledAxis & kDeltaEta) != 0) {
      title[dim] = "#eta_{D} - #eta_{jet}";
      nbins[dim] = 100;
      min[dim] = -radius * 1.2;
      max[dim] = radius * 1.2;
      dim++;
    }

    if ((enabledAxis & kDeltaPhi) != 0) {
      title[dim] = "#phi_{D} - #phi_{jet} (rad)";
      nbins[dim] = 100;
      min[dim] = -radius * 1.2;
      max[dim] = radius * 1.2;
      dim++;
    }

    title[dim] = "#it{p}_{T,jet} (GeV/#it{c})";
    nbins[dim] = nPtBins;
    min[dim] = 0;
    max[dim] = fMaxPt;
    dim++;

    if ((enabledAxis & kPositionJet) != 0) {
      title[dim] = "#eta_{jet}";
      nbins[dim] = 50;
      min[dim] = -1;
      max[dim] = 1;
      dim++;

      title[dim] = "#phi_{jet} (rad)";
      nbins[dim] = 150;
      min[dim] = 0;
      max[dim] = TMath::TwoPi();
      dim++;
    }

    if ((enabledAxis & kJetConstituents) != 0) {
      title[dim] = "No. of constituents";
      nbins[dim] = 50;
      min[dim] = -0.5;
      max[dim] = 49.5;
      dim++;
    }

    hname = TString::Format("%s/%s/fDmesonJets", GetName(), jetDef.GetName());
    THnSparse* h = fHistManager->CreateTHnSparse(hname,hname,dim,nbins,min,max);
    for (Int_t j = 0; j < dim; j++) {
      h->GetAxis(j)->SetTitle(title[j]);
    }
  }
}

/// Post the output with D meson jets found in the current event
///
/// \return kTRUE on success
Bool_t AliAnalysisTaskDmesonJets::AnalysisEngine::FillTree(Bool_t applyKinCuts)
{
  TString hname;

  for (Int_t id = 0; id < fDmesonJets.size(); id++) {
    fCurrentDmesonJetInfo->Set(fDmesonJets[id]);
    Int_t accJets = 0;
    for (UInt_t ij = 0; ij < fJetDefinitions.size(); ij++) {
      fCurrentJetInfo[ij]->Reset();
      AliJetInfo* jet = fDmesonJets[id].GetJet(fJetDefinitions[ij].GetName());
      if (!jet) continue;
      if (applyKinCuts && !fJetDefinitions[ij].IsJetInAcceptance(*jet)) {
        hname = TString::Format("%s/%s/fHistRejectedJetPt", GetName(), fJetDefinitions[ij].GetName());
        fHistManager->FillTH1(hname, jet->Pt());
        hname = TString::Format("%s/%s/fHistRejectedJetPhi", GetName(), fJetDefinitions[ij].GetName());
        fHistManager->FillTH1(hname, jet->Phi_0_2pi());
        hname = TString::Format("%s/%s/fHistRejectedJetEta", GetName(), fJetDefinitions[ij].GetName());
        fHistManager->FillTH1(hname, jet->Eta());
        continue;
      }
      fCurrentJetInfo[ij]->Set(fDmesonJets[id], fJetDefinitions[ij].GetName());
      accJets++;
    }
    if (accJets > 0) {
      fTree->Fill();
    }
    else {
      hname = TString::Format("%s/fHistRejectedDMesonPt", GetName());
      fHistManager->FillTH1(hname, fDmesonJets[id].fD.Pt());
      hname = TString::Format("%s/fHistRejectedDMesonPhi", GetName());
      fHistManager->FillTH1(hname, fDmesonJets[id].fD.Phi_0_2pi());
      hname = TString::Format("%s/fHistRejectedDMesonEta", GetName());
      fHistManager->FillTH1(hname, fDmesonJets[id].fD.Eta());
      if (fCandidateType == kD0toKpi) {
        hname = TString::Format("%s/fHistRejectedDMesonInvMass", GetName());
        fHistManager->FillTH1(hname, fDmesonJets[id].fD.M());
      }
      else if (fCandidateType == kDstartoKpipi) {
        hname = TString::Format("%s/fHistRejectedDMeson2ProngInvMass", GetName());
        fHistManager->FillTH1(hname, fDmesonJets[id].fInvMass2Prong);

        hname = TString::Format("%s/fHistRejectedDMesonDeltaInvMass", GetName());
        fHistManager->FillTH1(hname, fDmesonJets[id].fD.M() - fDmesonJets[id].fInvMass2Prong);
      }
    }
  }
  return kTRUE;
}

/// Post the output with D meson jets found in the current event
///
/// \return kTRUE on success
Bool_t AliAnalysisTaskDmesonJets::AnalysisEngine::FillHnSparse(Bool_t applyKinCuts)
{
  TString hname;

  for (Int_t id = 0; id < fDmesonJets.size(); id++) {
    if (!IsAnyJetInAcceptance(fDmesonJets[id])) {
      hname = TString::Format("%s/fHistRejectedDMesonPt", GetName());
      fHistManager->FillTH1(hname, fDmesonJets[id].fD.Pt());
      hname = TString::Format("%s/fHistRejectedDMesonPhi", GetName());
      fHistManager->FillTH1(hname, fDmesonJets[id].fD.Phi_0_2pi());
      hname = TString::Format("%s/fHistRejectedDMesonEta", GetName());
      fHistManager->FillTH1(hname, fDmesonJets[id].fD.Eta());
    }
  }

  for (auto &jetDef : fJetDefinitions) {

    hname = TString::Format("%s/%s/fDmesonJets", GetName(), jetDef.GetName());
    THnSparse* h = static_cast<THnSparse*>(fHistManager->FindObject(hname));

    for (Int_t id = 0; id < fDmesonJets.size(); id++) {
      const AliJetInfo* jet = fDmesonJets[id].GetJet(jetDef.GetName());
      if (!jet) continue;
      if (!jetDef.IsJetInAcceptance(*jet)) {
        hname = TString::Format("%s/%s/fHistRejectedJetPt", GetName(), jetDef.GetName());
        fHistManager->FillTH1(hname, jet->Pt());
        hname = TString::Format("%s/%s/fHistRejectedJetPhi", GetName(), jetDef.GetName());
        fHistManager->FillTH1(hname, jet->Phi_0_2pi());
        hname = TString::Format("%s/%s/fHistRejectedJetEta", GetName(), jetDef.GetName());
        fHistManager->FillTH1(hname, jet->Eta());
        continue;
      }
      FillHnSparse(h, fDmesonJets[id], jetDef.GetName());
    }
  }

  return kTRUE;
}

/// Fill a THnSparse using information from a AliDmesonJetInfo object
///
/// \param h          Valid pointer to a THnSparse object
/// \param DmesonJet  Const reference to an AliDmesonJetInfo object
/// \param n          Jet name
Bool_t AliAnalysisTaskDmesonJets::AnalysisEngine::FillHnSparse(THnSparse* h, const AliDmesonJetInfo& DmesonJet, std::string n)
{
  // Fill the THnSparse histogram.

  Double_t contents[30] = {0.};

  Double_t z = DmesonJet.GetZ(n);
  Double_t deltaPhi = 0;
  Double_t deltaEta = 0;
  Double_t deltaR = DmesonJet.GetDistance(n, deltaEta, deltaPhi);

  std::map<std::string, AliJetInfo>::const_iterator it = DmesonJet.fJets.find(n);
  if (it == DmesonJet.fJets.end()) return kFALSE;

  for (Int_t i = 0; i < h->GetNdimensions(); i++) {
    TString title(h->GetAxis(i)->GetTitle());
    if      (title=="#it{p}_{T,D} (GeV/#it{c})")                     contents[i] = DmesonJet.fD.Pt();
    else if (title=="#eta_{D}")                                      contents[i] = DmesonJet.fD.Eta();
    else if (title=="#phi_{D} (rad)")                                contents[i] = DmesonJet.fD.Phi_0_2pi();
    else if (title=="#it{M}_{K#pi} (GeV/#it{c}^{2})")                contents[i] = DmesonJet.fInvMass2Prong > 0 ? DmesonJet.fInvMass2Prong : DmesonJet.fD.M();
    else if (title=="#it{M}_{K#pi#pi} (GeV/#it{c}^{2})")             contents[i] = DmesonJet.fD.M();
    else if (title=="#it{M}_{K#pi#pi} - #it{M}_{K#pi} (GeV/#it{c}^{2})") contents[i] = DmesonJet.fD.M() - DmesonJet.fInvMass2Prong;
    else if (title=="#it{p}_{T,#pi} (GeV/#it{c})")                   contents[i] = DmesonJet.fSoftPionPt;
    else if (title=="#it{z}_{D}")                                    contents[i] = z;
    else if (title=="#Delta R_{D-jet}")                              contents[i] = deltaR;
    else if (title=="#eta_{D} - #eta_{jet}")                         contents[i] = deltaEta;
    else if (title=="#phi_{D} - #phi_{jet} (rad)")                   contents[i] = deltaPhi;
    else if (title=="#it{p}_{T,jet} (GeV/#it{c})")                   contents[i] = (*it).second.Pt();
    else if (title=="#eta_{jet}")                                    contents[i] = (*it).second.Eta();
    else if (title=="#phi_{jet} (rad)")                              contents[i] = (*it).second.Phi_0_2pi();
    else if (title=="No. of constituents")                           contents[i] = (*it).second.fNConstituents;
    else AliWarning(Form("Unable to fill dimension '%s'!",title.Data()));
  }

  h->Fill(contents);

  return kTRUE;
}

// Definitions of class AliAnalysisTaskDmesonJets

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJets);
/// \endcond

/// This is the default constructor, used for ROOT I/O purposes.
AliAnalysisTaskDmesonJets::AliAnalysisTaskDmesonJets() :
  AliAnalysisTaskEmcalLight(),
  fAnalysisEngines(),
  fEnabledAxis(0),
  fTreeOutput(kFALSE),
  fHistManager(),
  fApplyKinematicCuts(kTRUE),
  fNOutputTrees(0),
  fAodEvent(0),
  fFastJetWrapper(0)
{
  SetMakeGeneralHistograms(kTRUE);
}

/// This is the standard named constructor.
///
/// \param name Name of the task
AliAnalysisTaskDmesonJets::AliAnalysisTaskDmesonJets(const char* name, Int_t nOutputTrees) :
  AliAnalysisTaskEmcalLight(name, kTRUE),
  fAnalysisEngines(),
  fEnabledAxis(k2ProngInvMass),
  fTreeOutput(kFALSE),
  fHistManager(name),
  fApplyKinematicCuts(kTRUE),
  fNOutputTrees(nOutputTrees),
  fAodEvent(0),
  fFastJetWrapper(0)
{
  SetMakeGeneralHistograms(kTRUE);
  for (Int_t i = 0; i < nOutputTrees; i++){
    DefineOutput(2+i, TTree::Class());
  }
}

/// This is the standard destructor.
AliAnalysisTaskDmesonJets::~AliAnalysisTaskDmesonJets()
{
  if (fFastJetWrapper) delete fFastJetWrapper;
}

/// Load D meson cuts from a file.
///
/// \param cutfname Name of the file containing the cut object
/// \param cutsname Name of the object cuts
///
/// \return Pointer to the AliRDHFCuts object if successful, NULL otherwise.
AliRDHFCuts* AliAnalysisTaskDmesonJets::LoadDMesonCutsFromFile(TString cutfname, TString cutsname)
{
  AliRDHFCuts* analysiscuts = 0;
  TFile* filecuts = TFile::Open(cutfname);
  if (!filecuts || filecuts->IsZombie()) {
    ::Warning("AddTaskDmesonJets", "Input file not found: will use std cuts.");
    filecuts = 0;
  }

  if (filecuts) {
    analysiscuts = dynamic_cast<AliRDHFCuts*>(filecuts->Get(cutsname));
    if (!analysiscuts) {
      ::Warning("AddTaskDmesonJetCorr", "Could not find analysis cuts '%s' in '%s'. Using std cuts.", cutsname.Data(), cutfname.Data());
    }
  }

  return analysiscuts;
}

/// Add a new AnalysisEngine object.
///
/// \param type      One of the enum constants of ECandidateType_t
/// \param bkgMode   One of the enum constants of EMCMode_t
/// \param jetradius Radius of the jet
/// \param cuts      Name of the file that container D meson cut object (if null, it will use standard cuts)
///
/// \return Pointer to the AnalysisEngine added to the list.
AliAnalysisTaskDmesonJets::AnalysisEngine* AliAnalysisTaskDmesonJets::AddAnalysisEngine(ECandidateType_t type, EMCMode_t MCmode, EJetType_t jettype, Double_t jetradius, TString cutfname)
{
  AliHFJetDefinition jetDef(jettype, jetradius, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme);
  return AddAnalysisEngine(type, MCmode, jetDef, cutfname);
}

/// Add a new AnalysisEngine object.
///
/// \param type      One of the enum constants of ECandidateType_t
/// \param bkgMode   One of the enum constants of EMCMode_t
/// \param jetradius Radius of the jet
/// \param cuts      Name of the file that container D meson cut object (if null, it will use standard cuts)
///
/// \return Pointer to the AnalysisEngine added to the list.
AliAnalysisTaskDmesonJets::AnalysisEngine* AliAnalysisTaskDmesonJets::AddAnalysisEngine(ECandidateType_t type, EMCMode_t MCmode, const AliHFJetDefinition& jetDef, TString cutfname)
{
  AliRDHFCuts* cuts = 0;

  if (!cutfname.IsNull()) {
    TString cutsname;

    switch (type) {
    case kD0toKpi :
      cutsname = "D0toKpiCuts";
      break;
    case kDstartoKpipi :
      cutsname = "DStartoKpipiCuts";
      break;
    default:
      return 0;
    }

    cuts = LoadDMesonCutsFromFile(cutfname, cutsname);
  }

  AnalysisEngine eng(type, MCmode, cuts);

  std::list<AnalysisEngine>::iterator it = FindAnalysisEngine(eng);

  if (it == fAnalysisEngines.end() || *it != eng) {  // No analysis engine was found, adding a new one
    eng.AddJetDefinition(jetDef);
    it = fAnalysisEngines.insert(it, eng);
    ::Info("AliAnalysisTaskDmesonJets::AddAnalysisEngine", "A new analysis engine '%s' has been added. The total number of analysis engines is %lu.", eng.GetName(jetDef), fAnalysisEngines.size());
  }
  else {
    AnalysisEngine* found_eng = &(*it);
    ::Info("AliAnalysisTaskDmesonJets::AddAnalysisEngine", "An analysis engine '%s' with %lu jet definitions has been found. The total number of analysis engines is %lu. A new jet definition '%s' is being added.", found_eng->GetName(), found_eng->fJetDefinitions.size(), fAnalysisEngines.size(), jetDef.GetName());
    found_eng->AddJetDefinition(jetDef);

    if (cuts && found_eng->fRDHFCuts != 0) {
      ::Warning("AliAnalysisTaskDmesonJets::AddAnalysisEngine", "D meson cuts were already defined for this D meson type. They will be overwritten.");
      found_eng->SetRDHFCuts(cuts);
    }
  }

  return &(*it);
}

std::list<AliAnalysisTaskDmesonJets::AnalysisEngine>::iterator AliAnalysisTaskDmesonJets::FindAnalysisEngine(const AliAnalysisTaskDmesonJets::AnalysisEngine& eng)
{
  std::list<AnalysisEngine>::iterator it = fAnalysisEngines.begin();
  while (it != fAnalysisEngines.end() && (*it) < eng) it++;
  return it;
}

/// Creates the output containers.
void AliAnalysisTaskDmesonJets::UserCreateOutputObjects()
{
  ::Info("UserCreateOutputObjects", "CreateOutputObjects of task %s", GetName());

  AliAnalysisTaskEmcalLight::UserCreateOutputObjects();

  // Define histograms
  // the TList fOutput is already defined in  AliAnalysisTaskEmcalLight::UserCreateOutputObjects()

  TString hname;
  TString htitle;
  TH1* h = 0;
  Int_t treeSlot = 0;

  ::Info("AliAnalysisTaskDmesonJets::UserCreateOutputObjects", "Allocating histograms for task '%s' (%lu analysis engines)", GetName(), fAnalysisEngines.size());
  for (auto &param : fAnalysisEngines) {
    ::Info("AliAnalysisTaskDmesonJets::UserCreateOutputObjects", "Allocating histograms for analysis engine '%s' (%lu jet definitions)", param.GetName(), param.fJetDefinitions.size());

    fHistManager.CreateHistoGroup(param.GetName());

    param.fHistManager = &fHistManager;

    hname = TString::Format("%s/fHistNAcceptedDmesons", param.GetName());
    htitle = hname + ";Number of D accepted meson candidates;counts";
    h = fHistManager.CreateTH1(hname, htitle, 51, -0.5, 50.5);

    hname = TString::Format("%s/fHistNDmesons", param.GetName());
    htitle = hname + ";Number of D meson candidates;counts";
    h = fHistManager.CreateTH1(hname, htitle, 101, -0.5, 100.5);

    hname = TString::Format("%s/fHistNEvents", param.GetName());
    htitle = hname + ";Event status;counts";
    h = fHistManager.CreateTH1(hname, htitle, 2, 0, 2);
    h->GetXaxis()->SetBinLabel(1, "Accepted");
    h->GetXaxis()->SetBinLabel(2, "Rejected");

    hname = TString::Format("%s/fHistEventRejectionReasons", param.GetName());
    htitle = hname + ";Rejection reason;counts";
    h = fHistManager.CreateTH1(hname, htitle, 32, 0, 32);

    hname = TString::Format("%s/fHistRejectedDMesonPt", param.GetName());
    htitle = hname + ";#it{p}_{T,D} (GeV/#it{c});counts";
    fHistManager.CreateTH1(hname, htitle, 150, 0, 150);

    hname = TString::Format("%s/fHistRejectedDMesonEta", param.GetName());
    htitle = hname + ";#it{#eta}_{D};counts";
    fHistManager.CreateTH1(hname, htitle, 100, -2, 2);

    hname = TString::Format("%s/fHistRejectedDMesonPhi", param.GetName());
    htitle = hname + ";#it{#phi}_{D};counts";
    fHistManager.CreateTH1(hname, htitle, 200, 0, TMath::TwoPi());

    if (param.fCandidateType == kD0toKpi) {
      hname = TString::Format("%s/fHistRejectedDMesonInvMass", param.GetName());
      htitle = hname + ";#it{M}_{K#pi} (GeV/#it{c}^{2});counts";
      fHistManager.CreateTH1(hname, htitle, param.fNMassBins, param.fMinMass, param.fMaxMass);
    }
    else if (param.fCandidateType == kDstartoKpipi) {
      Double_t min = 0;
      Double_t max = 0;

      hname = TString::Format("%s/fHistRejectedDMeson2ProngInvMass", param.GetName());
      htitle = hname + ";#it{M}_{K#pi} (GeV/#it{c}^{2});counts";
      CalculateMassLimits(param.fMaxMass - param.fMinMass, 421, param.fNMassBins, min, max);
      fHistManager.CreateTH1(hname, htitle, param.fNMassBins, min, max);

      Double_t D0mass = TDatabasePDG::Instance()->GetParticle(421)->Mass();
      hname = TString::Format("%s/fHistRejectedDMesonDeltaInvMass", param.GetName());
      htitle = hname + ";#it{M}_{K#pi#pi} - #it{M}_{K#pi} (GeV/#it{c}^{2});counts";
      CalculateMassLimits(0.20, 413, param.fNMassBins*6, min, max);
      fHistManager.CreateTH1(hname, htitle, param.fNMassBins*6, min-D0mass, max-D0mass);
    }

    for (std::vector<AliHFJetDefinition>::iterator itdef = param.fJetDefinitions.begin(); itdef != param.fJetDefinitions.end(); itdef++) {
      AliHFJetDefinition* jetDef = &(*itdef);
      ::Info("AliAnalysisTaskDmesonJets::UserCreateOutputObjects", "Allocating histograms for jet definition '%s'", jetDef->GetName());

      fHistManager.CreateHistoGroup(jetDef->GetName(), param.GetName());

      hname = TString::Format("%s/%s/fHistMCParticleRejectionReason", param.GetName(), jetDef->GetName());
      htitle = hname + ";Track rejection reason;#it{p}_{T,track} (GeV/#it{c});counts";
      h = fHistManager.CreateTH2(hname, htitle, 32, 0, 32, 150, 0, 150);
      SetRejectionReasonLabels(h->GetXaxis());

      hname = TString::Format("%s/%s/fHistTrackRejectionReason", param.GetName(), jetDef->GetName());
      htitle = hname + ";Track rejection reason;#it{p}_{T,track} (GeV/#it{c});counts";
      h = fHistManager.CreateTH2(hname, htitle, 32, 0, 32, 150, 0, 150);
      SetRejectionReasonLabels(h->GetXaxis());

      hname = TString::Format("%s/%s/fHistClusterRejectionReason", param.GetName(), jetDef->GetName());
      htitle = hname + ";Cluster rejection reason;#it{p}_{T,cluster} (GeV/#it{c});counts";
      h = fHistManager.CreateTH2(hname, htitle, 32, 0, 32, 150, 0, 150);
      SetRejectionReasonLabels(h->GetXaxis());

      hname = TString::Format("%s/%s/fHistDMesonDaughterNotInJet", param.GetName(), jetDef->GetName());
      htitle = hname + ";#it{p}_{T,track} (GeV/#it{c});counts";
      fHistManager.CreateTH1(hname, htitle, 200, 0, 100);

      hname = TString::Format("%s/%s/fHistRejectedJetPt", param.GetName(), jetDef->GetName());
      htitle = hname + ";#it{p}_{T,jet} (GeV/#it{c});counts";
      fHistManager.CreateTH1(hname, htitle, 150, 0, 150);

      hname = TString::Format("%s/%s/fHistRejectedJetEta", param.GetName(), jetDef->GetName());
      htitle = hname + ";#it{#eta}_{jet};counts";
      fHistManager.CreateTH1(hname, htitle, 100, -2, 2);

      hname = TString::Format("%s/%s/fHistRejectedJetPhi", param.GetName(), jetDef->GetName());
      htitle = hname + ";#it{#phi}_{jet};counts";
      fHistManager.CreateTH1(hname, htitle, 200, 0, TMath::TwoPi());
    }
    if (fTreeOutput) {
      param.BuildTree(GetName());
      if (treeSlot < fNOutputTrees) {
        param.AssignDataSlot(treeSlot+2);
        treeSlot++;
        PostDataFromAnalysisEngine(param);
      }
      else {
        AliError(Form("Number of data output slots %d not sufficient. Tree of analysis engine %s will not be posted!", fNOutputTrees, param.GetName()));
      }
    }
    else {
      param.BuildHnSparse(fEnabledAxis);
    }
  }

  fOutput->Add(fHistManager.GetListOfHistograms());

  PostData(1, fOutput);
}

/// Does some specific initializations for the analysis engines,
/// then calls the base class ExecOnce() method.
void AliAnalysisTaskDmesonJets::ExecOnce()
{
  AliAnalysisTaskEmcalLight::ExecOnce();

  // Load the event
  fAodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);

  fFastJetWrapper = new AliFJWrapper(fName, fTitle);

  fFastJetWrapper->SetAreaType(fastjet::active_area);
  fFastJetWrapper->SetGhostArea(1);

  if (!fAodEvent) {
     AliError(Form("This task need an AOD event! Task '%s' will be disabled!", GetName()));
     return;
  }

  for (auto &params : fAnalysisEngines) {

    params.fAodEvent = fAodEvent;
    params.fFastJetWrapper = fFastJetWrapper;
    params.Init(fGeom, fAodEvent->GetRunNumber());

    if (params.fMCMode != kMCTruth) {
      params.fCandidateArray = dynamic_cast<TClonesArray*>(fAodEvent->GetList()->FindObject(params.fBranchName.Data()));

      if (params.fCandidateArray) {
        if (!params.fCandidateArray->GetClass()->InheritsFrom("AliAODRecoDecayHF2Prong")) {
          ::Error("AliAnalysisTaskDmesonJets::ExecOnce",
              "%s: Objects of type %s in %s are not inherited from AliAODRecoDecayHF2Prong! Task will be disabled!",
              GetName(), params.fCandidateArray->GetClass()->GetName(), params.fCandidateArray->GetName());
          params.fCandidateArray = 0;
          params.fInhibit = kTRUE;
        }
      }
      else {
        ::Error("AliAnalysisTaskDmesonJets::ExecOnce",
            "Could not find candidate array '%s', skipping the event. Analysis engine '%s' will be disabled!",
            params.fBranchName.Data(), params.GetName());
        params.fInhibit = kTRUE;
      }
    }

    if (params.fMCMode != kNoMC) {
      params.fMCContainer = dynamic_cast<AliHFAODMCParticleContainer*>(GetParticleContainer(0));

      if (!params.fMCContainer) params.fMCContainer = dynamic_cast<AliHFAODMCParticleContainer*>(GetParticleContainer(1));

      if (!params.fMCContainer) {
        ::Error("AliAnalysisTaskDmesonJets::ExecOnce",
            "No MC particle container was provided. Analysis engine '%s' will be disabled!",
            params.GetName());
        params.fInhibit = kTRUE;
      }
    }

    if (params.fMCMode != kMCTruth) {
      params.fTrackContainer = dynamic_cast<AliHFTrackContainer*>(GetParticleContainer(0));
      if (!params.fTrackContainer) params.fTrackContainer = dynamic_cast<AliHFTrackContainer*>(GetParticleContainer(1));

      params.fClusterContainer = GetClusterContainer(0);

      if (!params.fTrackContainer && !params.fClusterContainer) {
        ::Error("AliAnalysisTaskDmesonJets::ExecOnce",
            "No track container and no cluster container were provided. Analysis engine '%s' will be disabled!",
            params.GetName());
        params.fInhibit = kTRUE;
      }
    }
  }
}

/// Run the analysis
///
/// \return kTRUE on success
Bool_t AliAnalysisTaskDmesonJets::Run()
{
  if (!fAodEvent) return kFALSE;

  TString hname;

  // fix for temporary bug in ESDfilter
  // the AODs with null vertex pointer didn't pass the PhysSel
  if (!fAodEvent->GetPrimaryVertex() || TMath::Abs(fAodEvent->GetMagneticField()) < 0.001) return kFALSE;

  for (auto &eng : fAnalysisEngines) {

    if (eng.fInhibit) continue;

    //Event selection
    hname = TString::Format("%s/fHistNEvents", eng.GetName());
    Bool_t iseventselected = eng.fRDHFCuts->IsEventSelected(fAodEvent);
    if (!iseventselected) {
      fHistManager.FillTH1(hname, "Rejected");
      hname = TString::Format("%s/fHistEventRejectionReasons", eng.GetName());
      UInt_t bitmap = eng.fRDHFCuts->GetEventRejectionBitMap();
      TString label;
      do {
        label = GetHFEventRejectionReasonLabel(bitmap);
        if (label.IsNull()) break;
        fHistManager.FillTH1(hname, label);
      } while (true);

      continue;
    }

    fHistManager.FillTH1(hname, "Accepted");

    AliDebug(2, "Event selected");

    eng.RunAnalysis();
  }
  return kTRUE;
}

/// Fill the histograms.
///
/// \return Always kTRUE
Bool_t AliAnalysisTaskDmesonJets::FillHistograms()
{
  TString hname;
  for (auto &param : fAnalysisEngines) {

    if (param.fInhibit) continue;

    if (fTreeOutput) {
      param.FillTree(fApplyKinematicCuts);
    }
    else {
      param.FillHnSparse(fApplyKinematicCuts);
    }
  }
  return kTRUE;
}

/// Set the mass limits for the histograms using information from TDatabasePDG.
///
/// \param range   This parameter is used to calculate the mass range as [mass - range/2 ; mass + range/2]
/// \param pdg     PDG code of the candidate
/// \param nbins   Number of bins in the histogram
/// \param minMass Reference to a Double_t where the minimum mass will be stored
/// \param maxMass Reference to a Double_t where the maximum mass will be stored
void AliAnalysisTaskDmesonJets::CalculateMassLimits(Double_t range, Int_t pdg, Int_t nbins, Double_t& minMass, Double_t& maxMass)
{
  TParticlePDG* part = TDatabasePDG::Instance()->GetParticle(TMath::Abs(pdg));
  
  Double_t mass = part->Mass();

  // To make sure that the PDG mass value is not at the edge of a bin
  if (nbins % 2 == 0) {
    minMass = mass - range / 2 - range / nbins / 2;
    maxMass = mass + range / 2 - range / nbins / 2;
  }
  else {
    minMass = mass - range / 2;
    maxMass = mass + range / 2;
  }
}

/// Takes a bitmap and converts the first rejection reason bit to a string; it unsets the first bit.
///
/// \param bitmap Bitmap with one or more bit sets by AliRDHFCuts (only the first one will be considered)
///
/// \return A string that corresponds to the last bit set in the bitmap (a null string if not bit is set)
const char* AliAnalysisTaskDmesonJets::GetHFEventRejectionReasonLabel(UInt_t& bitmap)
{
  static TString label;
  label = "";

  if (bitmap & BIT(AliRDHFCuts::kNotSelTrigger)) {
    label = "NotSelTrigger";
    bitmap &= ~BIT(AliRDHFCuts::kNotSelTrigger);
    return label.Data();
  }
  if (bitmap & BIT(AliRDHFCuts::kNoVertex)) {
    label = "NoVertex";
    bitmap &= ~BIT(AliRDHFCuts::kNoVertex);
    return label.Data();
  }
  if (bitmap & BIT(AliRDHFCuts::kTooFewVtxContrib)) {
    label = "TooFewVtxContrib";
    bitmap &= ~BIT(AliRDHFCuts::kTooFewVtxContrib);
    return label.Data();
  }
  if (bitmap & BIT(AliRDHFCuts::kZVtxOutFid)) {
    label = "ZVtxOutFid";
    bitmap &= ~BIT(AliRDHFCuts::kZVtxOutFid);
    return label.Data();
  }
  if (bitmap & BIT(AliRDHFCuts::kPileup)) {
    label = "Pileup";
    bitmap &= ~BIT(AliRDHFCuts::kPileup);
    return label.Data();
  }
  if (bitmap & BIT(AliRDHFCuts::kOutsideCentrality)) {
    label = "OutsideCentrality";
    bitmap &= ~BIT(AliRDHFCuts::kOutsideCentrality);
    return label.Data();
  }
  if (bitmap & BIT(AliRDHFCuts::kPhysicsSelection)) {
    label = "PhysicsSelection";
    bitmap &= ~BIT(AliRDHFCuts::kPhysicsSelection);
    return label.Data();
  }
  if (bitmap & BIT(AliRDHFCuts::kBadSPDVertex)) {
    label = "BadSPDVertex";
    bitmap &= ~BIT(AliRDHFCuts::kBadSPDVertex);
    return label.Data();
  }
  if (bitmap & BIT(AliRDHFCuts::kZVtxSPDOutFid)) {
    label = "ZVtxSPDOutFid";
    bitmap &= ~BIT(AliRDHFCuts::kZVtxSPDOutFid);
    return label.Data();
  }
  if (bitmap & BIT(AliRDHFCuts::kCentralityFlattening)) {
    label = "CentralityFlattening";
    bitmap &= ~BIT(AliRDHFCuts::kCentralityFlattening);
    return label.Data();
  }
  if (bitmap & BIT(AliRDHFCuts::kBadTrackV0Correl)) {
    label = "BadTrackV0Correl";
    bitmap &= ~BIT(AliRDHFCuts::kBadTrackV0Correl);
    return label.Data();
  }

  return label.Data();
}

/// Post the tree of an analysis engine in the data slot (if the tree exists and the data slot has been assigned)
///
/// \param eng Constant reference to an analysis engine
///
/// \return -1 if unsuccessful, an integer number corresponding to the data slot if successful
Int_t AliAnalysisTaskDmesonJets::PostDataFromAnalysisEngine(const AnalysisEngine& eng)
{
  if (eng.GetDataSlotNumber() >= 0 && eng.GetTree()) {
    PostData(eng.GetDataSlotNumber(), eng.GetTree());
    return eng.GetDataSlotNumber();
  }
  else {
    return -1;
  }
}
