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
#include <TRandom3.h>

// Aliroot general
#include "AliLog.h"
#include "AliEMCALGeometry.h"
#include "AliAnalysisManager.h"
#include "AliVEventHandler.h"

// Aliroot HF
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliHFAODMCParticleContainer.h"
#include "AliHFTrackContainer.h"
#include "AliAnalysisVertexingHF.h"

// Aliroot EMCal jet framework
#include "AliEmcalJetTask.h"
#include "AliEmcalJet.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliEmcalParticle.h"
#include "AliFJWrapper.h"
#include "AliRhoParameter.h"

#include "AliAnalysisTaskDmesonJets.h"

// Definitions of class AliAnalysisTaskDmesonJets::AliJetInfo

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJets::AliJetInfo);
/// \endcond

/// Calculates the distance between this jet and another jet
///
/// \param jet Const reference to a AliJetInfo object
/// \param deta reference where the eta distance will be returned
/// \param dphi reference where the phi distance will be returned
/// \return The distance between this jet and the jet reference provided
Double_t AliAnalysisTaskDmesonJets::AliJetInfo::GetDistance(const AliJetInfo& jet, Double_t& deta, Double_t& dphi) const
{
  dphi = TVector2::Phi_mpi_pi(fMomentum.Phi() - jet.Phi());;
  deta = fMomentum.Eta() - jet.Eta();
  return TMath::Sqrt(dphi*dphi + deta*deta);
}

/// Calculates the distance between this jet and another jet
///
/// \param jet Const reference to a AliJetInfo object
/// \return The distance between this jet and the jet reference provided
Double_t AliAnalysisTaskDmesonJets::AliJetInfo::GetDistance(const AliJetInfo& jet) const
{
  Double_t deta = 0;
  Double_t dphi = 0;
  return GetDistance(jet, deta, dphi);
}

// Definitions of class AliAnalysisTaskDmesonJets::AliDmesonJetInfo

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJets::AliDmesonJetInfo);
/// \endcond

/// Default constructor
AliAnalysisTaskDmesonJets::AliDmesonJetInfo::AliDmesonJetInfo() :
  fDmesonParticle(0),
  fD(),
  fSoftPionPt(0),
  fInvMass2Prong(0),
  fJets(),
  fMCLabel(-1),
  fReconstructed(kFALSE),
  fFirstParton(0),
  fFirstPartonType(0),
  fLastParton(0),
  fLastPartonType(0),
  fSelectionType(0)
{
}

/// Copy constructor
///
/// \param source AliDmesonJetInfo object to copy from
AliAnalysisTaskDmesonJets::AliDmesonJetInfo::AliDmesonJetInfo(const AliDmesonJetInfo &source) :
  fDmesonParticle(source.fDmesonParticle),
  fD(source.fD),
  fSoftPionPt(source.fSoftPionPt),
  fInvMass2Prong(source.fInvMass2Prong),
  fJets(source.fJets),
  fMCLabel(source.fMCLabel),
  fReconstructed(source.fReconstructed),
  fFirstParton(source.fFirstParton),
  fFirstPartonType(source.fFirstPartonType),
  fLastParton(source.fLastParton),
  fLastPartonType(source.fLastPartonType),
  fSelectionType(source.fSelectionType)
{
}

/// Assignment operator
///
/// \param source AliDmesonJetInfo object to copy from
AliAnalysisTaskDmesonJets::AliDmesonJetInfo& AliAnalysisTaskDmesonJets::AliDmesonJetInfo::operator=(const AliDmesonJetInfo& source)
{
  new (this) AliDmesonJetInfo(source);
  return *this;
}

/// Reset all fields to their default values
void AliAnalysisTaskDmesonJets::AliDmesonJetInfo::Reset()
{
  fD.SetPtEtaPhiE(0,0,0,0);
  fSoftPionPt = 0;
  fInvMass2Prong = 0;
  fDmesonParticle = 0;
  fMCLabel = -1;
  fReconstructed = kFALSE;
  fFirstParton = 0;
  fFirstPartonType = 0;
  fLastParton = 0;
  fLastPartonType = 0;
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
/// \param n Name of the jet definition
/// \param deta reference where the eta distance will be returned
/// \param dphi reference where the phi distance will be returned
/// \return The distance between the D meson and the jet axis
Double_t AliAnalysisTaskDmesonJets::AliDmesonJetInfo::GetDistance(std::string n, Double_t& deta, Double_t& dphi) const
{
  std::map<std::string, AliJetInfo>::const_iterator it = fJets.find(n);
  if (it == fJets.end()) return 0;

  return GetDistance((*it).second, deta, dphi);
}

/// Calculates the distance between the D meson and the jet axis
///
/// \param n Name of the jet definition
/// \return The distance between the D meson and the jet axis
Double_t AliAnalysisTaskDmesonJets::AliDmesonJetInfo::GetDistance(std::string n) const
{
  Double_t deta = 0;
  Double_t dphi = 0;
  return GetDistance(n, deta, dphi);
}

/// Calculates the distance between the D meson and the jet axis
///
/// \param jet Const reference to a AliJetInfo object
/// \param deta reference where the eta distance will be returned
/// \param dphi reference where the phi distance will be returned
/// \return The distance between the D meson and the jet axis
Double_t AliAnalysisTaskDmesonJets::AliDmesonJetInfo::GetDistance(const AliJetInfo& jet, Double_t& deta, Double_t& dphi) const
{
  dphi = TVector2::Phi_mpi_pi(fD.Phi() - jet.Phi());;
  deta = fD.Eta() - jet.Eta();
  return TMath::Sqrt(dphi*dphi + deta*deta);
}

/// Calculates the distance between the D meson and the jet axis
///
/// \param jet Const reference to a AliJetInfo object
/// \return The distance between the D meson and the jet axis
Double_t AliAnalysisTaskDmesonJets::AliDmesonJetInfo::GetDistance(const AliJetInfo& jet) const
{
  Double_t deta = 0;
  Double_t dphi = 0;
  return GetDistance(jet, deta, dphi);
}

/// Find jet info object corresponding a jet definition provided as a string
///
/// \param n String containing the jet definition
/// \return Constant pointer to the jet info object, if the jet definition was found. Null pointer otherwise.
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

/// Find jet info object corresponding a jet definition provided as a string
///
/// \param n String containing the jet definition
/// \return Pointer to the jet info object, if the jet definition was found. Null pointer otherwise.
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
  fZ(0),
  fN(0)
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
  fN = 0;
}

/// Set the current object using an instance of AliDmesonJetInfo as its source
///
/// \param source A const reference to a valid AliDmesonJetInfo object
/// \param i      Index of the jet to be copied
void AliAnalysisTaskDmesonJets::AliJetInfoSummary::Set(const AliDmesonJetInfo& source, std::string n)
{
  std::map<std::string, AliJetInfo>::const_iterator it = source.fJets.find(n);
  if (it == source.fJets.end()) return;

  Set((*it).second);

  fR = source.GetDistance(n);
  fZ = source.GetZ(n);
}

/// Set the current object using an instance of AliJetInfo as its source
///
/// \param source A const reference to a valid AliJetInfo object
void AliAnalysisTaskDmesonJets::AliJetInfoSummary::Set(const AliJetInfo& source)
{
  fPt = source.Pt();
  fEta = source.Eta();
  fPhi = source.Phi_0_2pi();
  fN = source.GetNConstituents();
  fR = 0;
  fZ = 0;
}

// Definitions of class AliAnalysisTaskDmesonJets::AliJetInfoPbPbSummary

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJets::AliJetInfoPbPbSummary);
/// \endcond

/// Constructor that uses an AliDmesonJetInfo as its source
///
/// \param source A const reference to a valid AliDmesonJetInfo object
/// \param i      Index of the jet to be copied
AliAnalysisTaskDmesonJets::AliJetInfoPbPbSummary::AliJetInfoPbPbSummary(const AliDmesonJetInfo& source, std::string n) :
  AliJetInfoSummary(),
  fCorrPt(0),
  fArea(0)
{
  Set(source, n);
}

/// Reset the current object
void AliAnalysisTaskDmesonJets::AliJetInfoPbPbSummary::Reset()
{
  AliJetInfoSummary::Reset();
  fCorrPt = 0;
  fArea = 0;
}

/// Set the current object using an instance of AliJetInfo as its source
///
/// \param source A const reference to a valid AliJetInfo object
void AliAnalysisTaskDmesonJets::AliJetInfoPbPbSummary::Set(const AliJetInfo& source)
{
  AliJetInfoSummary::Set(source);
  fArea = source.fArea;
  fCorrPt = source.fCorrPt;
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

/// Reset the object
void AliAnalysisTaskDmesonJets::AliDmesonInfoSummary::Reset()
{
  fPt = 0;
  fEta = 0;
  fPhi = 0;
}

// Definitions of class AliAnalysisTaskDmesonJets::AliDmesonMCInfoSummary

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJets::AliDmesonMCInfoSummary);
/// \endcond

/// Constructor that uses an AliDmesonJetInfo as its source
///
/// \param source A const reference to a valid AliDmesonJetInfo object
AliAnalysisTaskDmesonJets::AliDmesonMCInfoSummary::AliDmesonMCInfoSummary(const AliDmesonJetInfo& source) :
  AliDmesonInfoSummary(source),
  fFirstPartonType(0),
  fFirstPartonPt(0),
  fLastPartonType(0),
  fLastPartonPt(0)
{
  Set(source);
}

/// Set the current object using an instance of AliDmesonJetInfo as its source
///
/// \param source A const reference to a valid AliDmesonJetInfo object
void AliAnalysisTaskDmesonJets::AliDmesonMCInfoSummary::Set(const AliDmesonJetInfo& source)
{
  AliDmesonInfoSummary::Set(source);

  fFirstPartonType = source.fFirstPartonType;
  fLastPartonType = source.fLastPartonType;

  if (source.fFirstParton) {
    fFirstPartonPt = source.fFirstParton->Pt();
  }
  else {
    fFirstPartonPt = 0.;
  }

  if (source.fLastParton) {
    fLastPartonPt = source.fLastParton->Pt();
  }
  else {
    fLastPartonPt = 0.;
  }
}

/// Reset the object
void AliAnalysisTaskDmesonJets::AliDmesonMCInfoSummary::Reset()
{
  AliDmesonInfoSummary::Reset();
  fFirstPartonType = 0,
  fFirstPartonPt = 0.;
  fLastPartonType = 0,
  fLastPartonPt = 0.;
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
  fInvMass(source.fD.M()),
  fSelectionType(0)
{
}

/// Set the current object using an instance of AliDmesonJetInfo as its source
///
/// \param source A const reference to a valid AliDmesonJetInfo object
void AliAnalysisTaskDmesonJets::AliD0InfoSummary::Set(const AliDmesonJetInfo& source)
{
  fInvMass = source.fD.M();
  fSelectionType = source.fSelectionType;
  AliDmesonInfoSummary::Set(source);
}

/// Reset the object
void AliAnalysisTaskDmesonJets::AliD0InfoSummary::Reset()
{
  AliDmesonInfoSummary::Reset();
  fSelectionType = 0;
  fInvMass = 0;
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

/// Set the current object using an instance of AliDmesonJetInfo as its source
///
/// \param source A const reference to a valid AliDmesonJetInfo object
void AliAnalysisTaskDmesonJets::AliDStarInfoSummary::Set(const AliDmesonJetInfo& source)
{
  f2ProngInvMass = source.fInvMass2Prong;
  fDeltaInvMass = source.fD.M() - source.fInvMass2Prong;
  AliDmesonInfoSummary::Set(source);
}

/// Reset the object
void AliAnalysisTaskDmesonJets::AliDStarInfoSummary::Reset()
{
  AliDmesonInfoSummary::Reset();

  f2ProngInvMass = 0;
  fDeltaInvMass = 0;
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
  fMinJetPt(0.),
  fMaxJetPt(500.),
  fMinJetPhi(0.),
  fMaxJetPhi(0.),
  fMinJetEta(-1.),
  fMaxJetEta(1.),
  fMinChargedPt(0.),
  fMaxChargedPt(0.),
  fMinNeutralPt(0.),
  fMaxNeutralPt(0.),
  fRhoName(),
  fRho(0)
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
  fMinJetPt(0.),
  fMaxJetPt(500.),
  fMinJetPhi(0.),
  fMaxJetPhi(0.),
  fMinJetEta(-1.),
  fMaxJetEta(1.),
  fMinChargedPt(0.),
  fMaxChargedPt(0.),
  fMinNeutralPt(0.),
  fMaxNeutralPt(0.),
  fRhoName(),
  fRho(0)
{
}

/// Default constructor
///
/// \param type Jet type (full, charged, neutral)
/// \param r    Jet resolution parameter
/// \param algo Jet algorithm (anit-kt, kt,...)
/// \param reco Jet recombination scheme (pt_scheme, E_scheme,...)
AliAnalysisTaskDmesonJets::AliHFJetDefinition::AliHFJetDefinition(EJetType_t type, Double_t r, EJetAlgo_t algo, ERecoScheme_t reco, TString rhoName) :
  TObject(),
  fJetType(type),
  fRadius(r),
  fJetAlgo(algo),
  fRecoScheme(reco),
  fMinJetPt(0.),
  fMaxJetPt(500.),
  fMinJetPhi(0.),
  fMaxJetPhi(0.),
  fMinJetEta(-1.),
  fMaxJetEta(1.),
  fMinChargedPt(0.),
  fMaxChargedPt(0.),
  fMinNeutralPt(0.),
  fMaxNeutralPt(0.),
  fRhoName(rhoName),
  fRho(0)
{
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
  fMinJetPt(source.fMinJetPt),
  fMaxJetPt(source.fMaxJetPt),
  fMinJetPhi(source.fMinJetPhi),
  fMaxJetPhi(source.fMaxJetPhi),
  fMinJetEta(source.fMinJetEta),
  fMaxJetEta(source.fMaxJetEta),
  fMinChargedPt(source.fMinChargedPt),
  fMaxChargedPt(source.fMaxChargedPt),
  fMinNeutralPt(source.fMinNeutralPt),
  fMaxNeutralPt(source.fMaxNeutralPt),
  fRhoName(source.fRhoName),
  fRho(0)
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
  if (jet.Pt() > fMaxJetPt || jet.Pt() < fMinJetPt) return kFALSE;
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
  fFirstPartons(),
  fLastPartons(),
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
  fRandomGen(0),
  fTrackEfficiency(0),
  fDataSlotNumber(-1),
  fTree(0),
  fCurrentDmesonJetInfo(0),
  fCurrentJetInfo(0),
  fCandidateArray(0),
  fMCContainer(),
  fTrackContainers(),
  fClusterContainers(),
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
  fFirstPartons(),
  fLastPartons(),
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
  fRejectedOrigin(0),
  fAcceptedDecay(kAnyDecay),
  fInhibit(kFALSE),
  fJetDefinitions(),
  fPtBinWidth(0.5),
  fMaxPt(100),
  fRandomGen(0),
  fTrackEfficiency(0),
  fDataSlotNumber(-1),
  fTree(0),
  fCurrentDmesonJetInfo(0),
  fCurrentJetInfo(0),
  fCandidateArray(0),
  fMCContainer(),
  fTrackContainers(),
  fClusterContainers(),
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
  fFirstPartons(source.fFirstPartons),
  fLastPartons(source.fLastPartons),
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
  fRandomGen(source.fRandomGen),
  fTrackEfficiency(source.fTrackEfficiency),
  fDataSlotNumber(-1),
  fTree(0),
  fCurrentDmesonJetInfo(0),
  fCurrentJetInfo(0),
  fCandidateArray(source.fCandidateArray),
  fMCContainer(source.fMCContainer),
  fTrackContainers(source.fTrackContainers),
  fClusterContainers(source.fClusterContainers),
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
void AliAnalysisTaskDmesonJets::AnalysisEngine::Init(const AliEMCALGeometry* const /*geom*/, Int_t /*runNumber*/)
{
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
    fAcceptedDecay = kDecayD0toKpi;
    break;
  case kD0toKpiLikeSign:
    fCandidatePDG = 421;
    fCandidateName = "2ProngLikeSign";
    fNDaughters = 2;
    fPDGdaughters.Set(fNDaughters);
    fPDGdaughters.Reset();
    fPDGdaughters[0] = 211;  // pi
    fPDGdaughters[1] = 321;  // K
    fBranchName = "LikeSign2Prong";
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
    fAcceptedDecay = kDecayDStartoKpipi;
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
  case kWrongPID:
    name += "_WrongPID";
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

/// Extract attributes of the D meson candidate.
///
/// \param Dcand Pointer to a AliAODRecoDecayHF2Prong representing the D meson candidate
/// \param DmesonJet Reference to an AliDmesonJetInfo object where the D meson candidate information will be copied
/// \param i Either 0 or 1, for the two possible mass hypothesis assignments
///
/// \return kTRUE on success
Bool_t AliAnalysisTaskDmesonJets::AnalysisEngine::ExtractRecoDecayAttributes(const AliAODRecoDecayHF2Prong* Dcand, AliDmesonJetInfo& DmesonJet, UInt_t i)
{
  if (fCandidateType == kD0toKpi || fCandidateType == kD0toKpiLikeSign) { // D0 candidate
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
  AliDebug(10,"Checking if D0 meson is selected");
  Int_t isSelected = fRDHFCuts->IsSelected(const_cast<AliAODRecoDecayHF2Prong*>(Dcand), AliRDHFCuts::kAll, fAodEvent);
  if (isSelected == 0) return kFALSE;

  Int_t MCtruthPdgCode = 0;

  Double_t invMassD = 0;

  // If the analysis require knowledge of the MC truth, look for generated D meson matched to reconstructed candidate
  // Checks also the origin, and if it matches the rejected origin mask, return false
  if (fMCMode == kBackgroundOnly || fMCMode == kSignalOnly || fMCMode == kWrongPID) {
    Int_t mcLab = Dcand->MatchToMC(fCandidatePDG, fMCContainer->GetArray(), fNDaughters, fPDGdaughters.GetArray());
    DmesonJet.fMCLabel = mcLab;

    // Retrieve the generated particle (if exists) and its PDG code
    if (mcLab >= 0) {
      AliAODMCParticle* aodMcPart = static_cast<AliAODMCParticle*>(fMCContainer->GetArray()->At(mcLab));

      if (aodMcPart) {
        // Check origin and return false if it matches the rejected origin mask
        if (fRejectedOrigin) {
          auto origin = CheckOrigin(aodMcPart, fMCContainer->GetArray());
          if ((origin.first & fRejectedOrigin) == origin.first) return kFALSE;
        }
        MCtruthPdgCode = aodMcPart->PdgCode();
      }
    }
  }

  if (isSelected == 1) { // selected as a D0
    if (i != 0) return kFALSE; // only one mass hypothesis thanks to PID

    if (fMCMode == kNoMC ||
        (MCtruthPdgCode == fCandidatePDG && fMCMode == kSignalOnly) ||
        (MCtruthPdgCode != fCandidatePDG && fMCMode == kBackgroundOnly) ||
        (MCtruthPdgCode == -fCandidatePDG && fMCMode == kWrongPID)) {
      // both background and signal are requested OR (it is a true D0 AND signal is requested) OR (it is NOT a D0 and background is requested)
      AliDebug(10,"Selected as D0");
      invMassD = Dcand->InvMassD0();
    }
    else { // conditions above not passed, so return FALSE
      return kFALSE;
    }
  }
  else if (isSelected == 2) { // selected as a D0bar
    if (i != 1) return kFALSE; // only one mass hypothesis thanks to PID

    if (fMCMode == kNoMC ||
        (MCtruthPdgCode == -fCandidatePDG && fMCMode == kSignalOnly) ||
        (MCtruthPdgCode != -fCandidatePDG && fMCMode == kBackgroundOnly) ||
        (MCtruthPdgCode == fCandidatePDG && fMCMode == kWrongPID)) {
      // both background and signal are requested OR (it is a true D0bar AND signal is requested) OR (it is NOT a D0bar and background is requested)
      AliDebug(10,"Selected as D0bar");
      invMassD = Dcand->InvMassD0bar();
    }
    else { // conditions above not passed, so return FALSE
      return kFALSE;
    }
  }
  else if (isSelected == 3) { // selected as either a D0bar or a D0 (PID on K and pi undecisive)
    AliDebug(10,"Selected as either D0 or D0bar");

    // Accept the correct mass hypothesis for signal-only and the wrong one for background-only
    if ((MCtruthPdgCode == fCandidatePDG && fMCMode == kSignalOnly) ||
        (MCtruthPdgCode == -fCandidatePDG && (fMCMode == kBackgroundOnly || fMCMode == kWrongPID))) {
      if (i != 0) return kFALSE;
      AliDebug(10, "MC truth is D0");
      invMassD = Dcand->InvMassD0();
    }
    else if ((MCtruthPdgCode == -fCandidatePDG && fMCMode == kSignalOnly) ||
             (MCtruthPdgCode == fCandidatePDG && (fMCMode == kBackgroundOnly || fMCMode == kWrongPID))) {
      if (i != 1) return kFALSE;
      AliDebug(10, "MC truth is D0bar");
      invMassD = Dcand->InvMassD0bar();
    }
    else { // (This candidate is neither a D0 nor a D0bar) OR (background-and-signal was requested)

      // Only accept it if background-only OR background-and-signal was requested
      if (fMCMode == kBackgroundOnly || fMCMode == kNoMC) {
        // Select D0 or D0bar depending on the i-parameter
        if (i == 0) {
          AliDebug(10, "Returning invariant mass with D0 hypothesis");
          invMassD = Dcand->InvMassD0();
        }
        else if (i == 1) {
          AliDebug(10, "Returning invariant mass with D0bar hypothesis");
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
  AliDebug(10,"Checking if D* meson is selected");
  Int_t isSelected = fRDHFCuts->IsSelected(const_cast<AliAODRecoCascadeHF*>(DstarCand), AliRDHFCuts::kAll, fAodEvent);
  if (isSelected == 0) return kFALSE;

  if ((i == 1 && DstarCand->Charge()>0) || (i == 0 && DstarCand->Charge()<0) || i > 1) return kFALSE; // only one mass hypothesis for the D*

  Int_t MCtruthPdgCode = 0;

  Double_t invMassD = 0;

  if (fMCMode == kBackgroundOnly || fMCMode == kSignalOnly) {
    Int_t pdgDgDStartoD0pi[2] = { 421, 211 };  // D0,pi
    Int_t pdgDgD0toKpi[2] = { 321, 211 };      // K, pi

    Int_t mcLab = DstarCand->MatchToMC(fCandidatePDG, 421, pdgDgDStartoD0pi, pdgDgD0toKpi, fMCContainer->GetArray());
    AliDebug(10, Form("MC label is %d", mcLab));
    DmesonJet.fMCLabel = mcLab;
    if (mcLab >= 0) {
      AliAODMCParticle* aodMcPart = static_cast<AliAODMCParticle*>(fMCContainer->GetArray()->At(mcLab));

      if (aodMcPart) {
        if (fRejectedOrigin) {
          auto origin = CheckOrigin(aodMcPart, fMCContainer->GetArray());
          if ((origin.first & fRejectedOrigin) == origin.first) return kFALSE;
        }

        MCtruthPdgCode = aodMcPart->PdgCode();
        AliDebug(10, Form("MC truth pdg code is %d",MCtruthPdgCode));
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
  if (!part) return kUnknownDecay;
  if (!mcArray) return kUnknownDecay;

  EMesonDecayChannel_t decay = kUnknownDecay;

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
      else if (absPdg1 == 211 && absPdg2 == 421) {  // pi D0
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
std::pair<AliAnalysisTaskDmesonJets::EMesonOrigin_t, AliAODMCParticle*> AliAnalysisTaskDmesonJets::AnalysisEngine::CheckOrigin(const AliAODMCParticle* part, TClonesArray* mcArray, Bool_t firstParton)
{
  // Checks whether the mother of the particle comes from a charm or a bottom quark.

  std::pair<AliAnalysisTaskDmesonJets::EMesonOrigin_t, AliAODMCParticle*> result(kUnknownQuark, 0);

  if (!part) return result;
  if (!mcArray) return result;

  Int_t mother = part->GetMother();
  while (mother >= 0) {
    AliAODMCParticle* mcGranma = static_cast<AliAODMCParticle*>(mcArray->At(mother));
    if (mcGranma) {
      Int_t abspdgGranma = TMath::Abs(mcGranma->GetPdgCode());

      if (abspdgGranma == 1) result = {kFromDown, mcGranma};
      if (abspdgGranma == 2) result = {kFromUp, mcGranma};
      if (abspdgGranma == 3) result = {kFromStrange, mcGranma};
      if (abspdgGranma == 4) result = {kFromCharm, mcGranma};
      if (abspdgGranma == 5) result = {kFromBottom, mcGranma};
      if (abspdgGranma == 6) result = {kFromTop, mcGranma};
      if (abspdgGranma == 9 || abspdgGranma == 21) result = {kFromGluon, mcGranma};

      // If looking for the very first parton in the hard scattering, it will continue the loop until it cannot find a mother particle
      if (result.first != kUnknownQuark && !firstParton) return result;

      mother = mcGranma->GetMother();
    }
    else {
      ::Error("AliAnalysisTaskDmesonJets::AnalysisParams::CheckOrigin", "Could not retrieve mother particle %d!", mother);
      break;
    }
  }

  return result;
}

/// Run the analysis
void AliAnalysisTaskDmesonJets::AnalysisEngine::RunAnalysis()
{
  for (auto& jetDef : fJetDefinitions) {
    jetDef.fJets.clear();
  }

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
  // Fill the vertex info of the candidates
  // Needed for reduced delta AOD, where the vertex info has been deleted
  // to reduce the delta AOD file size
  AliAnalysisVertexingHF vHF;

  const Int_t nD = fCandidateArray->GetEntriesFast();

  AliDmesonJetInfo DmesonJet;

  Int_t nAccCharm[3] = {0};
  for (Int_t icharm = 0; icharm < nD; icharm++) {   //loop over D candidates
    AliAODRecoDecayHF2Prong* charmCand = static_cast<AliAODRecoDecayHF2Prong*>(fCandidateArray->At(icharm)); // D candidates
    if (!charmCand) continue;
    if(!(vHF.FillRecoCand(fAodEvent,charmCand))) continue;

    // region of interest + cuts
    if (!fRDHFCuts->IsInFiducialAcceptance(charmCand->Pt(), charmCand->Y(fCandidatePDG))) continue;
    Int_t nMassHypo = 0; // number of mass hypothesis accepted for this D meson
    for (Int_t im = 0; im < 2; im++)  {  // 2 mass hypothesis (when available)
      DmesonJet.Reset();
      DmesonJet.fDmesonParticle = charmCand;
      DmesonJet.fSelectionType = im + 1;
      if (ExtractRecoDecayAttributes(charmCand, DmesonJet, im)) {
        for (auto& def : fJetDefinitions) {
          if (!FindJet(charmCand, DmesonJet, def)) {
            AliWarning(Form("Could not find jet '%s' for D meson '%s': pT = %.3f, eta = %.3f, phi = %.3f",
                def.GetName(), GetName(), DmesonJet.fD.Pt(), DmesonJet.fD.Eta(), DmesonJet.fD.Phi_0_2pi()));
          }
        }
        fDmesonJets[(icharm+1)*(1-(im*2))] = DmesonJet;
        nMassHypo++;
        nAccCharm[im]++;
      }
    }
    if (nMassHypo == 2) {
      nAccCharm[0]--;
      nAccCharm[1]--;
      nAccCharm[2] += 2;
    }
    if (nMassHypo == 2) { // both mass hypothesis accepted
      fDmesonJets[(icharm+1)].fSelectionType = 3;
      fDmesonJets[-(icharm+1)].fSelectionType = 3;
    }
  } // end of D cand loop

  TString hname;

  hname = TString::Format("%s/fHistNTotAcceptedDmesons", GetName());
  fHistManager->FillTH1(hname, "D", nAccCharm[0]);
  fHistManager->FillTH1(hname, "Anti-D", nAccCharm[1]);
  fHistManager->FillTH1(hname, "Both", nAccCharm[2]);

  hname = TString::Format("%s/fHistNAcceptedDmesonsVsNtracks", GetName());
  Int_t ntracks = 0;
  for (auto track_cont : fTrackContainers) ntracks += track_cont->GetNAcceptedTracks();
  fHistManager->FillTH2(hname, ntracks, nAccCharm[0]+nAccCharm[1]+nAccCharm[2]);

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

  Double_t rho = 0;
  if (jetDef.fRho) rho = jetDef.fRho->GetVal();

  fFastJetWrapper->Clear();
  fFastJetWrapper->SetR(jetDef.fRadius);
  fFastJetWrapper->SetAlgorithm(AliEmcalJetTask::ConvertToFJAlgo(jetDef.fJetAlgo));
  fFastJetWrapper->SetRecombScheme(AliEmcalJetTask::ConvertToFJRecoScheme(jetDef.fRecoScheme));

  fFastJetWrapper->AddInputVector(DmesonJet.fD.Px(), DmesonJet.fD.Py(), DmesonJet.fD.Pz(), DmesonJet.fD.E(), 0);

  if (jetDef.fJetType != AliJetContainer::kNeutralJet) {
    for (auto track_cont : fTrackContainers) {
      AliHFTrackContainer* hftrack_cont = dynamic_cast<AliHFTrackContainer*>(track_cont);
      if (hftrack_cont) hftrack_cont->SetDMesonCandidate(Dcand);
      hname = TString::Format("%s/%s/fHistTrackRejectionReason", GetName(), jetDef.GetName());
      AddInputVectors(track_cont, 100, static_cast<TH2*>(fHistManager->FindObject(hname)), fTrackEfficiency);

      if (hftrack_cont) {
        hname = TString::Format("%s/%s/fHistDMesonDaughterNotInJet", GetName(), jetDef.GetName());
        TH1* histDaughterNotInJet = static_cast<TH1*>(fHistManager->FindObject(hname));
        const TObjArray& daughters = hftrack_cont->GetDaughterList();
        for (Int_t i = 0; i < daughters.GetEntriesFast(); i++) {
          AliVParticle* daughter = static_cast<AliVParticle*>(daughters.At(i));
          if (!hftrack_cont->GetArray()->FindObject(daughter)) histDaughterNotInJet->Fill(daughter->Pt());
        }
      }
    }
  }

  if (jetDef.fJetType != AliJetContainer::kChargedJet) {
    for (auto clus_cont : fClusterContainers) {
      hname = TString::Format("%s/%s/fHistClusterRejectionReason", GetName(), jetDef.GetName());
      AddInputVectors(clus_cont, -100, static_cast<TH2*>(fHistManager->FindObject(hname)));
    }
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
      DmesonJet.fJets[jetDef.GetName()].fArea = jets_incl[ijet].area();
      DmesonJet.fJets[jetDef.GetName()].fCorrPt = DmesonJet.fJets[jetDef.GetName()].fMomentum.Pt() - jets_incl[ijet].area() * rho;


      return kTRUE;
    }
  }

  return kFALSE;
}

/// Adds all the particles contained in the container into the fastjet wrapper
///
/// \param cont Pointer to a valid AliEmcalContainer object
void AliAnalysisTaskDmesonJets::AnalysisEngine::AddInputVectors(AliEmcalContainer* cont, Int_t offset, TH2* rejectHist, Double_t eff)
{
  auto itcont = cont->all_momentum();
  for (AliEmcalIterableMomentumContainer::iterator it = itcont.begin(); it != itcont.end(); it++) {
    UInt_t rejectionReason = 0;
    if (!cont->AcceptObject(it.current_index(), rejectionReason)) {
      if (rejectHist) rejectHist->Fill(AliEmcalContainer::GetRejectionReasonBitPosition(rejectionReason), it->first.Pt());
      continue;
    }
    if (fRandomGen && eff > 0 && eff < 1) {
      Double_t rnd = fRandomGen->Rndm();
      if (eff < rnd) {
        if (rejectHist) rejectHist->Fill(6, it->first.Pt());
        continue;
      }
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

  if (!fMCContainer->IsSpecialPDGFound()) return;

  Int_t nAccCharm[3] = {0};

  for (auto &jetDef : fJetDefinitions) {
    Double_t rho = 0;
    if (jetDef.fRho) rho = jetDef.fRho->GetVal();
    hname = TString::Format("%s/%s/fHistNDmesonsVsNconstituents", GetName(), jetDef.GetName());
    TH1* histNDmesonsVsNconstituents = static_cast<TH1*>(fHistManager->FindObject(hname));

    switch (jetDef.fJetType) {
    case AliJetContainer::kFullJet:
      fMCContainer->SetCharge(AliParticleContainer::EChargeCut_t::kNoChargeCut);
      break;
    case AliJetContainer::kChargedJet:
      fMCContainer->SetCharge(AliParticleContainer::EChargeCut_t::kCharged);
      break;
    case AliJetContainer::kNeutralJet:
      fMCContainer->SetCharge(AliParticleContainer::EChargeCut_t::kNeutral);
      break;
    }

    fFastJetWrapper->Clear();
    fFastJetWrapper->SetR(jetDef.fRadius);
    fFastJetWrapper->SetAlgorithm(AliEmcalJetTask::ConvertToFJAlgo(jetDef.fJetAlgo));
    fFastJetWrapper->SetRecombScheme(AliEmcalJetTask::ConvertToFJRecoScheme(jetDef.fRecoScheme));

    hname = TString::Format("%s/%s/fHistMCParticleRejectionReason", GetName(), jetDef.GetName());
    AddInputVectors(fMCContainer, 100, static_cast<TH2*>(fHistManager->FindObject(hname)));

    fFastJetWrapper->Run();

    std::vector<fastjet::PseudoJet> jets_incl = fFastJetWrapper->GetInclusiveJets();

    for (auto jet : jets_incl) {
      Int_t nDmesonsInJet = 0;

      for (auto constituent : jet.constituents()) {
        Int_t iPart = constituent.user_index() - 100;
        AliAODMCParticle* part = fMCContainer->GetMCParticle(iPart);
        if (!part) {
          ::Error("AliAnalysisTaskDmesonJets::AnalysisEngine::RunParticleLevelAnalysis", "Could not find jet constituent %d!", iPart);
          continue;
        }
        if (TMath::Abs(part->PdgCode()) == fCandidatePDG) {
          nDmesonsInJet++;
          std::map<int, AliDmesonJetInfo>::iterator dMesonJetIt = fDmesonJets.find(iPart);
          if (dMesonJetIt == fDmesonJets.end()) { // This D meson does not exist yet
            std::pair<int, AliDmesonJetInfo> element;
            element.first = iPart;
            dMesonJetIt = fDmesonJets.insert(element).first;
            (*dMesonJetIt).second.fD.SetPxPyPzE(part->Px(), part->Py(), part->Pz(), part->E());
            (*dMesonJetIt).second.fDmesonParticle = part;
            (*dMesonJetIt).second.fSelectionType = part->PdgCode() > 0 ? 1 : 2;

            UShort_t p = 0;
            UInt_t rs = 0;

            auto firstParton = CheckOrigin(part, fMCContainer->GetArray(), kTRUE);
            p = 0;
            rs = firstParton.first;
            while (rs >>= 1) { p++; }
            (*dMesonJetIt).second.fFirstPartonType = p;
            (*dMesonJetIt).second.fFirstParton = firstParton.second;

            auto lastParton = CheckOrigin(part, fMCContainer->GetArray(), kFALSE);
            p = 0;
            rs = lastParton.first;
            while (rs >>= 1) { p++; }
            (*dMesonJetIt).second.fLastPartonType = p;
            (*dMesonJetIt).second.fLastParton = lastParton.second;

            if (part->PdgCode() > 0) {
              nAccCharm[0]++;
            }
            else {
              nAccCharm[1]++;
            }
          }

          (*dMesonJetIt).second.fJets[jetDef.GetName()].fMomentum.SetPxPyPzE(jet.px(), jet.py(), jet.pz(), jet.E());
          (*dMesonJetIt).second.fJets[jetDef.GetName()].fNConstituents = jet.constituents().size();
          (*dMesonJetIt).second.fJets[jetDef.GetName()].fArea = jet.area();
          (*dMesonJetIt).second.fJets[jetDef.GetName()].fCorrPt = (*dMesonJetIt).second.fJets[jetDef.GetName()].fMomentum.Pt() - jet.area() * rho;
        } // if constituent is a D meson
      } // for each constituent
      if (nDmesonsInJet > 0) histNDmesonsVsNconstituents->Fill(jet.constituents().size(), nDmesonsInJet);
    } // for each jet
  } // for each jet definition

  if (fDmesonJets.size() != nAccCharm[0]+nAccCharm[1]) AliError(Form("I found %lu mesons (%d)?", fDmesonJets.size(), nAccCharm[0]+nAccCharm[1]));
  hname = TString::Format("%s/fHistNTotAcceptedDmesons", GetName());
  fHistManager->FillTH1(hname, "D", nAccCharm[0]);
  fHistManager->FillTH1(hname, "Anti-D", nAccCharm[1]);
  fHistManager->FillTH1(hname, "Both", nAccCharm[2]);

  hname = TString::Format("%s/fHistNAcceptedDmesonsVsNtracks", GetName());
  fHistManager->FillTH2(hname, fMCContainer->GetNAcceptedParticles(), nAccCharm[0]+nAccCharm[1]+nAccCharm[2]);

  hname = TString::Format("%s/fHistNDmesons", GetName());
  fHistManager->FillTH1(hname, nAccCharm[0]+nAccCharm[1]+nAccCharm[2]); // same as the number of accepted D mesons, since no selection is performed
}

/// Builds the tree where the output will be posted
///
/// \return Pointer to the new tree
TTree* AliAnalysisTaskDmesonJets::AnalysisEngine::BuildTree(const char* taskName)
{
  TString classname;
  if (fMCMode == kMCTruth) {
    classname = "AliAnalysisTaskDmesonJets::AliDmesonMCInfoSummary";
    fCurrentDmesonJetInfo = new AliDmesonMCInfoSummary();
  }
  else {
    switch (fCandidateType) {
    case kD0toKpi:
    case kD0toKpiLikeSign:
      classname = "AliAnalysisTaskDmesonJets::AliD0InfoSummary";
      fCurrentDmesonJetInfo = new AliD0InfoSummary();
      break;
    case kDstartoKpipi:
      classname = "AliAnalysisTaskDmesonJets::AliDStarInfoSummary";
      fCurrentDmesonJetInfo = new AliDStarInfoSummary();
      break;
    }
  }
  TString treeName = TString::Format("%s_%s", taskName, GetName());
  fTree = new TTree(treeName, treeName);
  fTree->Branch("DmesonJet", classname, &fCurrentDmesonJetInfo);
  fCurrentJetInfo = new AliJetInfoSummary*[fJetDefinitions.size()];
  for (Int_t i = 0; i < fJetDefinitions.size(); i++) {
    if (fJetDefinitions[i].fRhoName.IsNull()) {
      fCurrentJetInfo[i] = new AliJetInfoSummary();
      fTree->Branch(fJetDefinitions[i].GetName(), "AliAnalysisTaskDmesonJets::AliJetInfoSummary", &fCurrentJetInfo[i]);
    }
    else {
      fCurrentJetInfo[i] = new AliJetInfoPbPbSummary();
      fTree->Branch(fJetDefinitions[i].GetName(), "AliAnalysisTaskDmesonJets::AliJetInfoPbPbSummary", &fCurrentJetInfo[i]);
    }
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

    if (fCandidateType == kD0toKpi || fCandidateType == kD0toKpiLikeSign) {
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
  fFirstPartons.clear();
  fLastPartons.clear();
  for (auto& dmeson_pair : fDmesonJets) {
    fCurrentDmesonJetInfo->Set(dmeson_pair.second);
    Int_t accJets = 0;
    for (UInt_t ij = 0; ij < fJetDefinitions.size(); ij++) {
      fCurrentJetInfo[ij]->Reset();
      AliJetInfo* jet = dmeson_pair.second.GetJet(fJetDefinitions[ij].GetName());
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
      fCurrentJetInfo[ij]->Set(dmeson_pair.second, fJetDefinitions[ij].GetName());
      accJets++;
    }
    if (accJets > 0) {
      fFirstPartons[dmeson_pair.second.fFirstParton] = dmeson_pair.second.fFirstPartonType;
      fLastPartons[dmeson_pair.second.fLastParton] = dmeson_pair.second.fLastPartonType;

      fTree->Fill();
    }
    else {
      hname = TString::Format("%s/fHistRejectedDMesonPt", GetName());
      fHistManager->FillTH1(hname, dmeson_pair.second.fD.Pt());
      hname = TString::Format("%s/fHistRejectedDMesonPhi", GetName());
      fHistManager->FillTH1(hname, dmeson_pair.second.fD.Phi_0_2pi());
      hname = TString::Format("%s/fHistRejectedDMesonEta", GetName());
      fHistManager->FillTH1(hname, dmeson_pair.second.fD.Eta());
      if (fMCMode != kMCTruth) {
        if (fCandidateType == kD0toKpi || fCandidateType == kD0toKpiLikeSign) {
          hname = TString::Format("%s/fHistRejectedDMesonInvMass", GetName());
          fHistManager->FillTH1(hname, dmeson_pair.second.fD.M());
        }
        else if (fCandidateType == kDstartoKpipi) {
          hname = TString::Format("%s/fHistRejectedDMeson2ProngInvMass", GetName());
          fHistManager->FillTH1(hname, dmeson_pair.second.fInvMass2Prong);

          hname = TString::Format("%s/fHistRejectedDMesonDeltaInvMass", GetName());
          fHistManager->FillTH1(hname, dmeson_pair.second.fD.M() - dmeson_pair.second.fInvMass2Prong);
        }
      }
    }
  }

  hname = TString::Format("%s/fHistFirstPartonPt", GetName());
  TH1* histFirstPartonPt = static_cast<TH1*>(fHistManager->FindObject(hname));
  hname = TString::Format("%s/fHistFirstPartonEta", GetName());
  TH1* histFirstPartonEta = static_cast<TH1*>(fHistManager->FindObject(hname));
  hname = TString::Format("%s/fHistFirstPartonPhi", GetName());
  TH1* histFirstPartonPhi = static_cast<TH1*>(fHistManager->FindObject(hname));
  hname = TString::Format("%s/fHistFirstPartonType", GetName());
  TH1* histFirstPartonType = static_cast<TH1*>(fHistManager->FindObject(hname));

  for (auto parton : fFirstPartons) {
    if (!parton.first) continue;
    histFirstPartonPt->Fill(parton.first->Pt());
    histFirstPartonEta->Fill(parton.first->Eta());
    histFirstPartonPhi->Fill(TVector2::Phi_0_2pi(parton.first->Phi()));
    histFirstPartonType->Fill(parton.second);
  }

  hname = TString::Format("%s/fHistLastPartonPt", GetName());
  TH1* histLastPartonPt = static_cast<TH1*>(fHistManager->FindObject(hname));
  hname = TString::Format("%s/fHistLastPartonEta", GetName());
  TH1* histLastPartonEta = static_cast<TH1*>(fHistManager->FindObject(hname));
  hname = TString::Format("%s/fHistLastPartonPhi", GetName());
  TH1* histLastPartonPhi = static_cast<TH1*>(fHistManager->FindObject(hname));
  hname = TString::Format("%s/fHistLastPartonType", GetName());
  TH1* histLastPartonType = static_cast<TH1*>(fHistManager->FindObject(hname));

  for (auto parton : fLastPartons) {
    if (!parton.first) continue;
    histLastPartonPt->Fill(parton.first->Pt());
    histLastPartonEta->Fill(parton.first->Eta());
    histLastPartonPhi->Fill(TVector2::Phi_0_2pi(parton.first->Phi()));
    histLastPartonType->Fill(parton.second);
  }

  return kTRUE;
}

/// Fills QA histograms. This method is not used by the AliAnalysisTaskDmesonJets task,
/// but can be used by derived tasks that have a custom implementation to fill the output objects.
///
/// \return Always kTRUE
Bool_t AliAnalysisTaskDmesonJets::AnalysisEngine::FillQA(Bool_t applyKinCuts)
{
  TString hname;
  fFirstPartons.clear();
  fLastPartons.clear();
  for (auto& dmeson_pair : fDmesonJets) {
    Int_t accJets = 0;
    for (UInt_t ij = 0; ij < fJetDefinitions.size(); ij++) {
      AliJetInfo* jet = dmeson_pair.second.GetJet(fJetDefinitions[ij].GetName());
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
      accJets++;
    }
    if (accJets > 0) {
      fFirstPartons[dmeson_pair.second.fFirstParton] = dmeson_pair.second.fFirstPartonType;
      fLastPartons[dmeson_pair.second.fLastParton] = dmeson_pair.second.fLastPartonType;
    }
    else {
      hname = TString::Format("%s/fHistRejectedDMesonPt", GetName());
      fHistManager->FillTH1(hname, dmeson_pair.second.fD.Pt());
      hname = TString::Format("%s/fHistRejectedDMesonPhi", GetName());
      fHistManager->FillTH1(hname, dmeson_pair.second.fD.Phi_0_2pi());
      hname = TString::Format("%s/fHistRejectedDMesonEta", GetName());
      fHistManager->FillTH1(hname, dmeson_pair.second.fD.Eta());
      if (fMCMode != kMCTruth) {
        if (fCandidateType == kD0toKpi || fCandidateType == kD0toKpiLikeSign) {
          hname = TString::Format("%s/fHistRejectedDMesonInvMass", GetName());
          fHistManager->FillTH1(hname, dmeson_pair.second.fD.M());
        }
        else if (fCandidateType == kDstartoKpipi) {
          hname = TString::Format("%s/fHistRejectedDMeson2ProngInvMass", GetName());
          fHistManager->FillTH1(hname, dmeson_pair.second.fInvMass2Prong);

          hname = TString::Format("%s/fHistRejectedDMesonDeltaInvMass", GetName());
          fHistManager->FillTH1(hname, dmeson_pair.second.fD.M() - dmeson_pair.second.fInvMass2Prong);
        }
      }
    }
  }

  hname = TString::Format("%s/fHistFirstPartonPt", GetName());
  TH1* histFirstPartonPt = static_cast<TH1*>(fHistManager->FindObject(hname));
  hname = TString::Format("%s/fHistFirstPartonEta", GetName());
  TH1* histFirstPartonEta = static_cast<TH1*>(fHistManager->FindObject(hname));
  hname = TString::Format("%s/fHistFirstPartonPhi", GetName());
  TH1* histFirstPartonPhi = static_cast<TH1*>(fHistManager->FindObject(hname));
  hname = TString::Format("%s/fHistFirstPartonType", GetName());
  TH1* histFirstPartonType = static_cast<TH1*>(fHistManager->FindObject(hname));

  for (auto parton : fFirstPartons) {
    if (!parton.first) continue;
    histFirstPartonPt->Fill(parton.first->Pt());
    histFirstPartonEta->Fill(parton.first->Eta());
    histFirstPartonPhi->Fill(TVector2::Phi_0_2pi(parton.first->Phi()));
    histFirstPartonType->Fill(parton.second);
  }

  hname = TString::Format("%s/fHistLastPartonPt", GetName());
  TH1* histLastPartonPt = static_cast<TH1*>(fHistManager->FindObject(hname));
  hname = TString::Format("%s/fHistLastPartonEta", GetName());
  TH1* histLastPartonEta = static_cast<TH1*>(fHistManager->FindObject(hname));
  hname = TString::Format("%s/fHistLastPartonPhi", GetName());
  TH1* histLastPartonPhi = static_cast<TH1*>(fHistManager->FindObject(hname));
  hname = TString::Format("%s/fHistLastPartonType", GetName());
  TH1* histLastPartonType = static_cast<TH1*>(fHistManager->FindObject(hname));

  for (auto parton : fLastPartons) {
    if (!parton.first) continue;
    histLastPartonPt->Fill(parton.first->Pt());
    histLastPartonEta->Fill(parton.first->Eta());
    histLastPartonPhi->Fill(TVector2::Phi_0_2pi(parton.first->Phi()));
    histLastPartonType->Fill(parton.second);
  }

  return kTRUE;
}

/// Post the output with D meson jets found in the current event
///
/// \return kTRUE on success
Bool_t AliAnalysisTaskDmesonJets::AnalysisEngine::FillHnSparse(Bool_t applyKinCuts)
{
  TString hname;

  for (auto& dmeson_pair : fDmesonJets) {
    if (!IsAnyJetInAcceptance(dmeson_pair.second)) {
      hname = TString::Format("%s/fHistRejectedDMesonPt", GetName());
      fHistManager->FillTH1(hname, dmeson_pair.second.fD.Pt());
      hname = TString::Format("%s/fHistRejectedDMesonPhi", GetName());
      fHistManager->FillTH1(hname, dmeson_pair.second.fD.Phi_0_2pi());
      hname = TString::Format("%s/fHistRejectedDMesonEta", GetName());
      fHistManager->FillTH1(hname, dmeson_pair.second.fD.Eta());
    }
  }

  for (auto &jetDef : fJetDefinitions) {

    hname = TString::Format("%s/%s/fDmesonJets", GetName(), jetDef.GetName());
    THnSparse* h = static_cast<THnSparse*>(fHistManager->FindObject(hname));

    for (auto& dmeson_pair : fDmesonJets) {
      const AliJetInfo* jet = dmeson_pair.second.GetJet(jetDef.GetName());
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
      FillHnSparse(h, dmeson_pair.second, jetDef.GetName());
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
  fOutputType(kTreeOutput),
  fHistManager(),
  fApplyKinematicCuts(kTRUE),
  fNOutputTrees(0),
  fTrackEfficiency(0),
  fMCContainer(0),
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
  fOutputType(kTreeOutput),
  fHistManager(name),
  fApplyKinematicCuts(kTRUE),
  fNOutputTrees(nOutputTrees),
  fTrackEfficiency(0),
  fMCContainer(0),
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
    ::Error("AliAnalysisTaskDmesonJets::LoadDMesonCutsFromFile", "Input file not found: will use std cuts.");
    filecuts = 0;
  }

  if (filecuts) analysiscuts = dynamic_cast<AliRDHFCuts*>(filecuts->Get(cutsname));

  if (!analysiscuts) {
    ::Error("AliAnalysisTaskDmesonJets::LoadDMesonCutsFromFile", "Could not find analysis cuts '%s' in '%s'.", cutsname.Data(), cutfname.Data());
    if (filecuts) {
      filecuts->ls();
    }
  }
  else {
    ::Info("AliAnalysisTaskDmesonJets::LoadDMesonCutsFromFile", "Cuts '%s' loaded from file '%s'", cutsname.Data(), cutfname.Data());
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
AliAnalysisTaskDmesonJets::AnalysisEngine* AliAnalysisTaskDmesonJets::AddAnalysisEngine(ECandidateType_t type, TString cutfname, EMCMode_t MCmode, EJetType_t jettype, Double_t jetradius, TString rhoName)
{
  AliHFJetDefinition jetDef(jettype, jetradius, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, rhoName);
  return AddAnalysisEngine(type, cutfname, MCmode, jetDef, rhoName);
}

/// Add a new AnalysisEngine object.
///
/// \param type      One of the enum constants of ECandidateType_t
/// \param bkgMode   One of the enum constants of EMCMode_t
/// \param jetradius Radius of the jet
/// \param cuts      Name of the file that container D meson cut object (if null, it will use standard cuts)
///
/// \return Pointer to the AnalysisEngine added to the list.
AliAnalysisTaskDmesonJets::AnalysisEngine* AliAnalysisTaskDmesonJets::AddAnalysisEngine(ECandidateType_t type, TString cutfname, EMCMode_t MCmode, const AliHFJetDefinition& jetDef, TString rhoName)
{
  AliRDHFCuts* cuts = 0;

  if (!cutfname.IsNull()) {
    TString cutsname;

    switch (type) {
    case kD0toKpi:
    case kD0toKpiLikeSign:
      cutsname = "D0toKpiCuts";
      break;
    case kDstartoKpipi:
      cutsname = "DStartoKpipiCuts";
      break;
    default:
      return 0;
    }

    cuts = LoadDMesonCutsFromFile(cutfname, cutsname);
    if (cuts) cuts->PrintAll();
  }

  AnalysisEngine eng(type, MCmode, cuts);

  std::list<AnalysisEngine>::iterator it = FindAnalysisEngine(eng);

  if (it == fAnalysisEngines.end() || *it != eng) {  // No analysis engine was found, adding a new one
    eng.AddJetDefinition(jetDef);
    it = fAnalysisEngines.insert(it, eng);
    ::Info("AliAnalysisTaskDmesonJets::AddAnalysisEngine", "A new analysis engine '%s' has been added. The total number of analysis engines is %lu.", eng.GetName(), fAnalysisEngines.size());
  }
  else {
    AnalysisEngine* found_eng = &(*it);
    ::Info("AliAnalysisTaskDmesonJets::AddAnalysisEngine", "An analysis engine '%s' with %lu jet definitions has been found. The total number of analysis engines is %lu. A new jet definition '%s' is being added.", found_eng->GetName(), found_eng->fJetDefinitions.size(), fAnalysisEngines.size(), jetDef.GetName());
    found_eng->AddJetDefinition(jetDef);

    if (cuts) {
      if (found_eng->fRDHFCuts != 0) ::Warning("AliAnalysisTaskDmesonJets::AddAnalysisEngine", "D meson cuts were already defined for this D meson type. They will be overwritten.");
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

  hname = "fHistCharmPt";
  htitle = hname + ";#it{p}_{T,charm} (GeV/#it{c});counts";
  fHistManager.CreateTH1(hname, htitle, 500, 0, 1000);

  hname = "fHistCharmEta";
  htitle = hname + ";#eta_{charm};counts";
  fHistManager.CreateTH1(hname, htitle, 400, -10, 10);

  hname = "fHistCharmPhi";
  htitle = hname + ";#phi_{charm};counts";
  fHistManager.CreateTH1(hname, htitle, 125, 0, TMath::TwoPi());

  hname = "fHistCharmPt_Eta05";
  htitle = hname + ";#it{p}_{T,charm} (GeV/#it{c});counts";
  fHistManager.CreateTH1(hname, htitle, 500, 0, 1000);

  hname = "fHistBottomPt";
  htitle = hname + ";#it{p}_{T,bottom} (GeV/#it{c});counts";
  fHistManager.CreateTH1(hname, htitle, 500, 0, 1000);

  hname = "fHistBottomEta";
  htitle = hname + ";#eta_{bottom};counts";
  fHistManager.CreateTH1(hname, htitle, 400, -10, 10);

  hname = "fHistBottomPhi";
  htitle = hname + ";#phi_{bottom};counts";
  fHistManager.CreateTH1(hname, htitle, 125, 0, TMath::TwoPi());

  hname = "fHistBottomPt_Eta05";
  htitle = hname + ";#it{p}_{T,bottom} (GeV/#it{c});counts";
  fHistManager.CreateTH1(hname, htitle, 500, 0, 1000);

  hname = "fHistHighestPartonPt";
  htitle = hname + ";#it{p}_{T,bottom} (GeV/#it{c});counts";
  fHistManager.CreateTH1(hname, htitle, 500, 0, 1000);

  hname = "fHistHighestPartonType";
  htitle = hname + ";type;counts";
  fHistManager.CreateTH1(hname, htitle, 10, 0, 10);

  hname = "fHistNHeavyQuarks";
  htitle = hname + ";number of heavy-quarks;counts";
  fHistManager.CreateTH1(hname, htitle, 21, -0.5, 20.5);

  ::Info("AliAnalysisTaskDmesonJets::UserCreateOutputObjects", "Allocating histograms for task '%s' (%lu analysis engines)", GetName(), fAnalysisEngines.size());
  for (auto &param : fAnalysisEngines) {
    ::Info("AliAnalysisTaskDmesonJets::UserCreateOutputObjects", "Allocating histograms for analysis engine '%s' (%lu jet definitions)", param.GetName(), param.fJetDefinitions.size());

    param.fHistManager = &fHistManager;

    hname = TString::Format("%s/fHistNAcceptedDmesonsVsNtracks", param.GetName());
    htitle = hname + ";#it{N}_{tracks};#it{N}_{D};events";
    h = fHistManager.CreateTH2(hname, htitle, 251, -0.5, 250.5, 21, -0.5, 20.5);

    hname = TString::Format("%s/fHistNTotAcceptedDmesons", param.GetName());
    htitle = hname + ";;#it{N}_{D}";
    h = fHistManager.CreateTH1(hname, htitle, 3, 0, 3);

    hname = TString::Format("%s/fHistNDmesons", param.GetName());
    htitle = hname + ";#it{N}_{D};events";
    h = fHistManager.CreateTH1(hname, htitle, 501, -0.5, 500.5);

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

    if (param.fMCMode != kMCTruth) {
      if (param.fCandidateType == kD0toKpi || param.fCandidateType == kD0toKpiLikeSign) {
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
    }

    if (param.fMCMode == kMCTruth) {
      hname = TString::Format("%s/fHistFirstPartonPt", param.GetName());
      htitle = hname + ";#it{p}_{T,parton} (GeV/#it{c});counts";
      fHistManager.CreateTH1(hname, htitle, 500, 0, 1000);

      hname = TString::Format("%s/fHistFirstPartonEta", param.GetName());
      htitle = hname + ";#eta_{parton};counts";
      fHistManager.CreateTH1(hname, htitle, 400, -10, 10);

      hname = TString::Format("%s/fHistFirstPartonPhi", param.GetName());
      htitle = hname + ";#phi_{parton};counts";
      fHistManager.CreateTH1(hname, htitle, 125, 0, TMath::TwoPi());

      hname = TString::Format("%s/fHistFirstPartonType", param.GetName());
      htitle = hname + ";type;counts";
      fHistManager.CreateTH1(hname, htitle, 10, 0, 10);

      hname = TString::Format("%s/fHistLastPartonPt", param.GetName());
      htitle = hname + ";#it{p}_{T,parton} (GeV/#it{c});counts";
      fHistManager.CreateTH1(hname, htitle, 500, 0, 1000);

      hname = TString::Format("%s/fHistLastPartonEta", param.GetName());
      htitle = hname + ";#eta_{parton};counts";
      fHistManager.CreateTH1(hname, htitle, 400, -10, 10);

      hname = TString::Format("%s/fHistLastPartonPhi", param.GetName());
      htitle = hname + ";#phi_{parton};counts";
      fHistManager.CreateTH1(hname, htitle, 125, 0, TMath::TwoPi());

      hname = TString::Format("%s/fHistLastPartonType", param.GetName());
      htitle = hname + ";type;counts";
      fHistManager.CreateTH1(hname, htitle, 10, 0, 10);
    }

    for (auto& jetDef : param.fJetDefinitions) {
      ::Info("AliAnalysisTaskDmesonJets::UserCreateOutputObjects", "Allocating histograms for jet definition '%s'", jetDef.GetName());

      if (param.fMCMode == kMCTruth) {
        hname = TString::Format("%s/%s/fHistNDmesonsVsNconstituents", param.GetName(), jetDef.GetName());
        htitle = hname + ";#it{N}_{constituents};#it{N}_{D};counts";
        h = fHistManager.CreateTH2(hname, htitle, 51, -0.5, 50.5, 10, 0.5, 10.5);
      }

      hname = TString::Format("%s/%s/fHistMCParticleRejectionReason", param.GetName(), jetDef.GetName());
      htitle = hname + ";Track rejection reason;#it{p}_{T,track} (GeV/#it{c});counts";
      h = fHistManager.CreateTH2(hname, htitle, 32, 0, 32, 150, 0, 150);
      SetRejectionReasonLabels(h->GetXaxis());

      hname = TString::Format("%s/%s/fHistTrackRejectionReason", param.GetName(), jetDef.GetName());
      htitle = hname + ";Track rejection reason;#it{p}_{T,track} (GeV/#it{c});counts";
      h = fHistManager.CreateTH2(hname, htitle, 32, 0, 32, 150, 0, 150);
      SetRejectionReasonLabels(h->GetXaxis());

      hname = TString::Format("%s/%s/fHistClusterRejectionReason", param.GetName(), jetDef.GetName());
      htitle = hname + ";Cluster rejection reason;#it{p}_{T,cluster} (GeV/#it{c});counts";
      h = fHistManager.CreateTH2(hname, htitle, 32, 0, 32, 150, 0, 150);
      SetRejectionReasonLabels(h->GetXaxis());

      hname = TString::Format("%s/%s/fHistDMesonDaughterNotInJet", param.GetName(), jetDef.GetName());
      htitle = hname + ";#it{p}_{T,track} (GeV/#it{c});counts";
      fHistManager.CreateTH1(hname, htitle, 200, 0, 100);

      hname = TString::Format("%s/%s/fHistRejectedJetPt", param.GetName(), jetDef.GetName());
      htitle = hname + ";#it{p}_{T,jet} (GeV/#it{c});counts";
      fHistManager.CreateTH1(hname, htitle, 150, 0, 150);

      hname = TString::Format("%s/%s/fHistRejectedJetEta", param.GetName(), jetDef.GetName());
      htitle = hname + ";#it{#eta}_{jet};counts";
      fHistManager.CreateTH1(hname, htitle, 100, -2, 2);

      hname = TString::Format("%s/%s/fHistRejectedJetPhi", param.GetName(), jetDef.GetName());
      htitle = hname + ";#it{#phi}_{jet};counts";
      fHistManager.CreateTH1(hname, htitle, 200, 0, TMath::TwoPi());
    }
    switch (fOutputType) {
    case kTreeOutput:
      param.BuildTree(GetName());
      if (treeSlot < fNOutputTrees) {
        param.AssignDataSlot(treeSlot+2);
        treeSlot++;
        PostDataFromAnalysisEngine(param);
      }
      else {
        AliError(Form("Number of data output slots %d not sufficient. Tree of analysis engine %s will not be posted!", fNOutputTrees, param.GetName()));
      }
      break;
    case kTHnOutput:
      param.BuildHnSparse(fEnabledAxis);
      break;
    case kNoOutput:
      break;
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

  // TODO: make this settable
  fFastJetWrapper->SetAreaType(fastjet::active_area_explicit_ghosts);
  fFastJetWrapper->SetGhostArea(0.005);

  if (!fAodEvent) {
     AliError(Form("This task need an AOD event (Task '%s'). Expect troubles...", GetName()));
     //return;
  }

  TRandom* rnd = 0;
  if (fTrackEfficiency > 0 && fTrackEfficiency < 1) rnd = new TRandom3(0);

  for (auto cont_it : fParticleCollArray) {
    AliHFAODMCParticleContainer* part_cont = dynamic_cast<AliHFAODMCParticleContainer*>(cont_it.second);
    if (part_cont) fMCContainer = part_cont;
  }

  for (auto &params : fAnalysisEngines) {

    params.fAodEvent = fAodEvent;
    params.fFastJetWrapper = fFastJetWrapper;
    params.fTrackEfficiency = fTrackEfficiency;
    params.fRandomGen = rnd;

    for (auto &jetdef: params.fJetDefinitions) {
      if (!jetdef.fRhoName.IsNull()) {
        jetdef.fRho = dynamic_cast<AliRhoParameter*>(fInputEvent->FindListObject(jetdef.fRhoName));
        if (!jetdef.fRho) {
          ::Error("AliAnalysisTaskDmesonJets::ExecOnce",
              "%s: Could not find rho object '%s' for engine '%s'",
              jetdef.fRhoName.Data(), GetName(), params.GetName());
        }
      }
    }

    if (!params.fRDHFCuts) {
      if (params.fMCMode == kMCTruth) {
      ::Warning("AliAnalysisTaskDmesonJets::ExecOnce",
          "%s: RDHF cuts not provided for engine '%s'.",
          GetName(), params.GetName());
      }
      else {
        ::Error("AliAnalysisTaskDmesonJets::ExecOnce",
            "%s: RDHF cuts not provided. Engine '%s' disabled.",
            GetName(), params.GetName());
        params.fInhibit = kTRUE;
      }
    }

    params.fMCContainer = fMCContainer;

    for (auto cont_it : fParticleCollArray) {
      AliTrackContainer* track_cont = dynamic_cast<AliTrackContainer*>(cont_it.second);
      if (track_cont) params.fTrackContainers.push_back(track_cont);
    }

    for (auto cont_it : fClusterCollArray) params.fClusterContainers.push_back(cont_it.second);

    if (fAodEvent) params.Init(fGeom, fAodEvent->GetRunNumber());

    if (params.fMCMode != kMCTruth && fAodEvent) {
      params.fCandidateArray = dynamic_cast<TClonesArray*>(fAodEvent->GetList()->FindObject(params.fBranchName.Data()));

      if (params.fCandidateArray) {
        TString className;
        if (params.fCandidateType == kD0toKpi || params.fCandidateType == kD0toKpiLikeSign) {
          className = "AliAODRecoDecayHF2Prong";
        }
        else if (params.fCandidateType == kDstartoKpipi) {
          className = "AliAODRecoCascadeHF";
        }
        if (!params.fCandidateArray->GetClass()->InheritsFrom(className)) {
          ::Error("AliAnalysisTaskDmesonJets::ExecOnce",
              "%s: Objects of type %s in %s are not inherited from %s! Task will be disabled!",
              GetName(), params.fCandidateArray->GetClass()->GetName(), params.fCandidateArray->GetName(), className.Data());
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
      if (!params.fMCContainer) {
        ::Error("AliAnalysisTaskDmesonJets::ExecOnce",
            "No MC particle container was provided. Analysis engine '%s' will be disabled!",
            params.GetName());
        params.fInhibit = kTRUE;
      }
    }

    if (params.fMCMode != kMCTruth) {
      if (params.fTrackContainers.size() == 0 && params.fClusterContainers.size() == 0) {
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
  TString hname;

  // Fix for temporary bug in ESDfilter
  // The AODs with null vertex pointer didn't pass the PhysSel
  // Now adding an entry in the histogram so as to check that this is actually cutting anything out
  if (fAodEvent && (!fAodEvent->GetPrimaryVertex() || TMath::Abs(fAodEvent->GetMagneticField()) < 0.001)) {
    for (auto &eng : fAnalysisEngines) {
        if (eng.fInhibit) continue;
        hname = TString::Format("%s/fHistEventRejectionReasons", eng.GetName());
        fHistManager.FillTH1(hname, "ESDfilterBug");
    }
    return kFALSE;
  }

  if (fAodEvent) {
    Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
    if (matchingAODdeltaAODlevel <= 0) {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      for (auto &eng : fAnalysisEngines) {
          if (eng.fInhibit) continue;
          hname = TString::Format("%s/fHistEventRejectionReasons", eng.GetName());
          fHistManager.FillTH1(hname, "MismatchDeltaAOD");
      }
      return kFALSE;
    }
  }

  for (auto &eng : fAnalysisEngines) {
    eng.fDmesonJets.clear();
    if (eng.fInhibit) continue;

    //Event selection
    hname = TString::Format("%s/fHistNEvents", eng.GetName());
    if (fAodEvent) {
      Bool_t iseventselected = kTRUE;
      if (eng.fRDHFCuts) iseventselected = eng.fRDHFCuts->IsEventSelected(fAodEvent);
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
  for (auto &param : fAnalysisEngines) {
    if (param.fInhibit) continue;

    if (fOutputType == kTreeOutput) {
      param.FillTree(fApplyKinematicCuts);
      PostDataFromAnalysisEngine(param);
    }
    else if (fOutputType == kTHnOutput) {
      param.FillHnSparse(fApplyKinematicCuts);
    }
  }
  if (fMCContainer) FillPartonLevelHistograms();
  return kTRUE;
}

/// Fill histograms with parton-level information
void AliAnalysisTaskDmesonJets::FillPartonLevelHistograms()
{
  auto itcont = fMCContainer->all_momentum();
  Int_t nHQ = 0;
  Double_t highestPartonPt = 0;
  Int_t absPdgHighParton = 0;
  for (auto part : itcont) {
    Int_t absPdgCode = TMath::Abs(part.second->GetPdgCode());

    // Skip all particles that are not either quarks or gluons
    if (absPdgCode > 9 && absPdgCode != 21) continue;

    // Look for highest momentum parton
    if (highestPartonPt < part.first.Pt()) {
      highestPartonPt = part.first.Pt();
      absPdgHighParton = absPdgCode;
    }
    /*
    // Look for the mother PDG code
    Int_t motherIndex = part.second->GetMother();
    AliAODMCParticle *mother = 0;
    Int_t motherPdg = 0;
    Double_t motherPt = 0;
    if (motherIndex >= 0) {
      mother = fMCContainer->GetMCParticle(motherIndex);
      if (motherIndex) {
        motherPdg =  TMath::Abs(mother->GetPdgCode());
        motherPt = mother->Pt();
      }
    }
    */
    if (absPdgCode != 4 && absPdgCode != 5) continue;
    Bool_t notLastInPartonShower = kFALSE;
    for (Int_t idaugh = 0; idaugh < 2; idaugh++){
      Int_t daughterIndex = part.second->GetDaughter(idaugh);
      if (daughterIndex < 0) {
        AliDebug(10, Form("Could not find daughter of heavy quark (pdg=%d, pt=%.3f)!", absPdgCode, part.first.Pt()));
        continue;
      }
      AliAODMCParticle *daughter = fMCContainer->GetMCParticle(daughterIndex);
      if (!daughter) {
        AliDebug(10, Form("Could not find daughter %d of heavy quark (pdg=%d, pt=%.3f)!", daughterIndex, absPdgCode, part.first.Pt()));
        continue;
      }
      Int_t daughterAbsPdgCode = TMath::Abs(daughter->GetPdgCode());
      if (daughterAbsPdgCode <= 9 || daughterAbsPdgCode == 21) notLastInPartonShower = kTRUE; // this parton is not the last parton in the shower
      AliDebug(10, Form("Found daughter with PDG=%d, pt=%.3f", daughterAbsPdgCode, daughter->Pt()));
    }
    if (notLastInPartonShower) continue;

    if (absPdgCode == 4) {
      fHistManager.FillTH1("fHistCharmPt", part.first.Pt());
      fHistManager.FillTH1("fHistCharmEta", part.first.Eta());
      fHistManager.FillTH1("fHistCharmPhi", part.first.Phi_0_2pi());
      if (TMath::Abs(part.first.Eta()) < 0.5) fHistManager.FillTH1("fHistCharmPt_Eta05", part.first.Pt());
    }
    else if (absPdgCode == 5) {
      fHistManager.FillTH1("fHistBottomPt", part.first.Pt());
      fHistManager.FillTH1("fHistBottomEta", part.first.Eta());
      fHistManager.FillTH1("fHistBottomPhi", part.first.Phi_0_2pi());
      if (TMath::Abs(part.first.Eta()) < 0.5) fHistManager.FillTH1("fHistBottomPt_Eta05", part.first.Pt());
    }
    nHQ++;
  }
  fHistManager.FillTH1("fHistNHeavyQuarks", nHQ);
  fHistManager.FillTH1("fHistHighestPartonPt",highestPartonPt);
  Int_t partonType = 0;
  if (absPdgHighParton == 9 || absPdgHighParton == 21) {
    partonType = 7;
  }
  else {
    partonType = absPdgHighParton;
  }
  fHistManager.FillTH1("fHistHighestPartonType",partonType);
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

/// Create an instance of this class and add it to the analysis manager
///
/// \param ntracks name of the track collection
/// \param nclusters name of the calorimeter cluster collection
/// \param nMCpart name of the MC particle collection
/// \param nMaxTrees number of output trees
/// \param suffix additional suffix that can be added at the end of the task name
/// \return pointer to the new AliAnalysisTaskDmesonJets task
AliAnalysisTaskDmesonJets* AliAnalysisTaskDmesonJets::AddTaskDmesonJets(TString ntracks, TString nclusters, TString nMCpart, Int_t nMaxTrees, TString suffix)
{
  // Get the pointer to the existing analysis manager via the static access method
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskDmesonJets", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler) {
    ::Error("AddTaskEmcalJetSpectraQA", "This task requires an input event handler");
    return NULL;
  }

  EDataType_t dataType = kUnknownDataType;

  if (handler->InheritsFrom("AliESDInputHandler")) {
    dataType = kESD;
  }
  else if (handler->InheritsFrom("AliAODInputHandler")) {
    dataType = kAOD;
  }

  // Init the task and do settings
  if (ntracks == "usedefault") {
    if (dataType == kESD) {
      ntracks = "Tracks";
    }
    else if (dataType == kAOD) {
      ntracks = "tracks";
    }
    else {
      ntracks = "";
    }
  }

  if (nclusters == "usedefault") {
    if (dataType == kESD) {
      nclusters = "CaloClusters";
    }
    else if (dataType == kAOD) {
      nclusters = "caloClusters";
    }
    else {
      nclusters = "";
    }
  }

  if (nMCpart == "usedefault") {
    nMCpart = "mcparticles"; // Always needs AliAODMCParticle objects
  }

  TString name("AliAnalysisTaskDmesonJets");
  if (strcmp(suffix, "") != 0) {
    name += TString::Format("_%s", suffix.Data());
  }

  AliAnalysisTaskDmesonJets* jetTask = new AliAnalysisTaskDmesonJets(name, nMaxTrees);

  if (!ntracks.IsNull()) {
    AliHFTrackContainer* trackCont = new AliHFTrackContainer(ntracks);
    jetTask->AdoptParticleContainer(trackCont);
  }

  if (!nMCpart.IsNull()) {
    AliMCParticleContainer* partCont = new AliHFAODMCParticleContainer(nMCpart);
    partCont->SetEtaLimits(-1.5, 1.5);
    partCont->SetPtLimits(0, 1000);
    jetTask->AdoptParticleContainer(partCont);
  }

  jetTask->AddClusterContainer(nclusters.Data());

  // Final settings, pass to manager and set the containers
  mgr->AddTask(jetTask);

  // Create containers for input/output
  AliAnalysisDataContainer* cinput1 = mgr->GetCommonInputContainer();
  TString contname1(name);
  contname1 += "_histos";
  AliAnalysisDataContainer* coutput1 = mgr->CreateContainer(contname1.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s", AliAnalysisManager::GetCommonFileName()));

  mgr->ConnectInput(jetTask, 0, cinput1);
  mgr->ConnectOutput(jetTask, 1, coutput1);

  for (Int_t i = 0; i < nMaxTrees; i++) {
    TString contname = TString::Format("%s_tree_%d", name.Data(), i);
    AliAnalysisDataContainer *coutput = mgr->CreateContainer(contname.Data(),
        TTree::Class(),AliAnalysisManager::kOutputContainer,
        Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput(jetTask, 2+i, coutput);
  }
  return jetTask;
}

