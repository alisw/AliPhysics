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
  fJet.SetPtEtaPhiE(0,0,0,0);
  fJetLeadingPt = 0;
  fJetNConstituents = 0;
  fDaughterDistances.Reset();
}

/// Prints the content of this object in the standard output.
void AliAnalysisTaskDmesonJets::AliDmesonJetInfo::Print() const
{
  Printf("Printing D Meson Jet object.");
  Printf("D Meson: pT = %.3f, eta = %.3f, phi = %.3f, inv. mass = %.3f", fD.Pt(), fD.Eta(), fD.Phi_0_2pi(), fD.M());
  Printf("Soft pion pT: %.3f. 2-Prong Invariant mass = %.3f", fSoftPionPt, fInvMass2Prong);
  Printf("Jet: pT = %.3f, eta = %.3f, phi = %.3f", fJet.Pt(), fJet.Eta(), fJet.Phi_0_2pi());
  Printf("Leading pT = %.3f. Jet N Consituents = %d", fJetLeadingPt, fJetNConstituents);
}

// Definitions of class AliAnalysisTaskDmesonJets::AliJetDefinition

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJets::AliJetDefinition);
/// \endcond

/// This is the default constructor, used for ROOT I/O purposes.
AliAnalysisTaskDmesonJets::AliJetDefinition::AliJetDefinition() :
  TObject(),
  fJetType(AliJetContainer::kChargedJet),
  fRadius(0),
  fJetAlgo(AliJetContainer::antikt_algorithm),
  fRecoScheme(AliJetContainer::pt_scheme),
  fDmesonJets()
{
}

/// Default constructor
///
/// \param type Jet type (full, charged, neutral)
/// \param r    Jet resolution parameter
/// \param algo Jet algorithm (anit-kt, kt,...)
/// \param reco Jet recombination scheme (pt_scheme, E_scheme,...)
AliAnalysisTaskDmesonJets::AliJetDefinition::AliJetDefinition(EJetType_t type, Double_t r, EJetAlgo_t algo, ERecoScheme_t reco) :
  TObject(),
  fJetType(type),
  fRadius(r),
  fJetAlgo(algo),
  fRecoScheme(reco),
  fDmesonJets()
{
}

/// Copy constructor
///
/// \param source Reference to an AliJetDefinition object to copy from
AliAnalysisTaskDmesonJets::AliJetDefinition::AliJetDefinition(const AliJetDefinition &source) :
  TObject(),
  fJetType(source.fJetType),
  fRadius(source.fRadius),
  fJetAlgo(source.fJetAlgo),
  fRecoScheme(source.fRecoScheme),
  fDmesonJets()
{
}

/// Assignment operator
///
/// \param source Reference to an AliJetDefinition object to copy from
AliAnalysisTaskDmesonJets::AliJetDefinition& AliAnalysisTaskDmesonJets::AliJetDefinition::operator=(const AliJetDefinition& source)
{
  new (this) AliJetDefinition(source);
  return *this;
}

/// Generate a name for this jet definition
const char* AliAnalysisTaskDmesonJets::AliJetDefinition::GetName() const
{
  static TString name;

  name = AliJetContainer::GenerateJetName(fJetType, fJetAlgo, fRecoScheme, fRadius, 0, 0, "Jet");

  return name.Data();
}

/// Compares 2 jet definitions.
/// The ordering is based on: jet type, radius, algorithm and recombination scheme, in this order
///
/// \param lhs Reference to the first AliJetDefinition object
/// \param rhs Reference to the second AliJetDefinition object
bool operator<(const AliAnalysisTaskDmesonJets::AliJetDefinition& lhs, const AliAnalysisTaskDmesonJets::AliJetDefinition& rhs)
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
bool operator==(const AliAnalysisTaskDmesonJets::AliJetDefinition& lhs, const AliAnalysisTaskDmesonJets::AliJetDefinition& rhs)
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
  if (fRDHFCuts) delete fRDHFCuts;
}

/// Assignement operator
///
/// \param source Reference to a valid AnalysisEngine to copy from.
AliAnalysisTaskDmesonJets::AnalysisEngine& AliAnalysisTaskDmesonJets::AnalysisEngine::operator=(const AnalysisEngine& source)
{
  new (this) AnalysisEngine(source);
  return *this;
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
const char* AliAnalysisTaskDmesonJets::AnalysisEngine::GetName(const AliJetDefinition& jetDef) const
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
AliAnalysisTaskDmesonJets::AliJetDefinition* AliAnalysisTaskDmesonJets::AnalysisEngine::AddJetDefinition(const AliAnalysisTaskDmesonJets::AliJetDefinition& def)
{
  std::list<AliJetDefinition>::iterator it = FindJetDefinition(def);

  if (it == fJetDefinitions.end() || *it != def) {  // No jet definition was found, adding a new one
    it = fJetDefinitions.insert(it, def);
    ::Info("AliAnalysisTaskDmesonJets::AnalysisEngine::AddJetDefinition", "Jet definition '%s' has been added to analysis engine '%s'."
        "Total number of jet definitions is now %lu.",
        def.GetName(), GetName(), fJetDefinitions.size());
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
AliAnalysisTaskDmesonJets::AliJetDefinition*
AliAnalysisTaskDmesonJets::AnalysisEngine::AddJetDefinition(EJetType_t type, Double_t r, EJetAlgo_t algo, ERecoScheme_t reco)
{
  AliJetDefinition def(type, r, algo, reco);

  return AddJetDefinition(def);
}

/// Look for a jet definition that is equal
///
/// \param def Reference to a jet definition object
///
/// \return An iterator to the jet definition object, if it is found. An iterator to the end if not found.
std::list<AliAnalysisTaskDmesonJets::AliJetDefinition>::iterator AliAnalysisTaskDmesonJets::AnalysisEngine::FindJetDefinition(const AliAnalysisTaskDmesonJets::AliJetDefinition& def)
{
  std::list<AliJetDefinition>::iterator it = fJetDefinitions.begin();
  while (it != fJetDefinitions.end() && (*it) < def) it++;
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
AliAnalysisTaskDmesonJets::EMesonDecayChannel_t AliAnalysisTaskDmesonJets::AnalysisEngine::CheckDecayChannel(AliAODMCParticle* part, TClonesArray* mcArray)
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
AliAnalysisTaskDmesonJets::EMesonOrigin_t AliAnalysisTaskDmesonJets::AnalysisEngine::CheckOrigin(AliAODMCParticle* part, TClonesArray* mcArray)
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
  for (std::list<AliJetDefinition>::iterator itdef = fJetDefinitions.begin(); itdef != fJetDefinitions.end(); itdef++) {
    AliJetDefinition* jetDef = &(*itdef);
    jetDef->fDmesonJets.clear();
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

    DmesonJet.Reset();

    for (Int_t im = 0; im < 2; im++)  {  // 2 mass hypothesis (when available)
      if (ExtractRecoDecayAttributes(charmCand, DmesonJet, im)) {
        for (std::list<AliJetDefinition>::iterator itdef = fJetDefinitions.begin(); itdef != fJetDefinitions.end(); itdef++) {
          AliJetDefinition* jetDef = &(*itdef);
          FindJet(charmCand, DmesonJet, *jetDef);
        }
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
Bool_t AliAnalysisTaskDmesonJets::AnalysisEngine::FindJet(AliAODRecoDecayHF2Prong* Dcand, AliDmesonJetInfo& DmesonJet, AliJetDefinition& jetDef)
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
    Double_t maxPt = 0;

    for (UInt_t ic = 0; ic < constituents.size(); ++ic) {
      if (constituents[ic].user_index() == 0) {
        isDmesonJet = kTRUE;
      }
      if (constituents[ic].pt() > maxPt) {
        maxPt = constituents[ic].pt();
      }
    }

    if (isDmesonJet) {
      DmesonJet.fJet.SetPxPyPzE(jets_incl[ijet].px(), jets_incl[ijet].py(), jets_incl[ijet].pz(), jets_incl[ijet].E());
      DmesonJet.fJetNConstituents = constituents.size();
      DmesonJet.fJetLeadingPt = maxPt;

      jetDef.fDmesonJets.push_back(DmesonJet);
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
  AliTLorentzVector part;
  for (Int_t i = 0; i < cont->GetNEntries(); i++) {
    cont->GetMomentum(part, i);
    if (!cont->AcceptObject(i)) {
      rejectHist->Fill(cont->GetRejectionReasonBitPosition(), part.Pt());
      continue;
    }
    Int_t uid = offset >= 0 ? i : -i;
    uid += offset;
    fFastJetWrapper->AddInputVector(part.Px(), part.Py(), part.Pz(), part.E(), uid);
  }
}

/// Run a particle level analysis
void AliAnalysisTaskDmesonJets::AnalysisEngine::RunParticleLevelAnalysis()
{
  TString hname;

  fMCContainer->SetSpecialPDG(fCandidatePDG);
  fMCContainer->SetRejectedOriginMap(fRejectedOrigin);
  fMCContainer->SetAcceptedDecayMap(fAcceptedDecay);

  AliDmesonJetInfo DmesonJet;

  for (std::list<AliJetDefinition>::iterator itdef = fJetDefinitions.begin(); itdef != fJetDefinitions.end(); itdef++) {
    AliJetDefinition* jetDef = &(*itdef);

    fFastJetWrapper->Clear();
    fFastJetWrapper->SetR(jetDef->fRadius);
    fFastJetWrapper->SetAlgorithm(AliEmcalJetTask::ConvertToFJAlgo(jetDef->fJetAlgo));
    fFastJetWrapper->SetRecombScheme(AliEmcalJetTask::ConvertToFJRecoScheme(jetDef->fRecoScheme));

    hname = TString::Format("%s/%s/fHistClusterRejectionReason", GetName(), jetDef->GetName());
    AddInputVectors(fMCContainer, 100, static_cast<TH2*>(fHistManager->FindObject(hname)));

    fFastJetWrapper->Run();

    std::vector<fastjet::PseudoJet> jets_incl = fFastJetWrapper->GetInclusiveJets();

    for (UInt_t ijet = 0; ijet < jets_incl.size(); ++ijet) {
      DmesonJet.Reset();
      std::vector<fastjet::PseudoJet> constituents(fFastJetWrapper->GetJetConstituents(ijet));

      Bool_t isDmesonJet = kFALSE;

      for (UInt_t ic = 0; ic < constituents.size(); ++ic) {
        Int_t iPart = constituents[ic].user_index() - 100;
        AliVParticle* part = fMCContainer->GetParticle(iPart);
        if (!part) {
          ::Error("AliAnalysisTaskDmesonJets::AnalysisEngine::RunParticleLevelAnalysis", "Could not find jet constituent %d!", iPart);
          continue;
        }
        if (part->PdgCode() == fCandidatePDG) {
          DmesonJet.fD.SetPxPyPzE(part->Px(), part->Py(), part->Pz(), part->E());
          isDmesonJet = kTRUE;
          break;
        }
      }

      if (isDmesonJet) {
        DmesonJet.fJet.SetPxPyPzE(jets_incl[ijet].px(), jets_incl[ijet].py(), jets_incl[ijet].pz(), jets_incl[ijet].E());
      }
      jetDef->fDmesonJets.push_back(DmesonJet);
    }
  }
}

// Definitions of class AliAnalysisTaskDmesonJets

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJets);
/// \endcond

/// This is the default constructor, used for ROOT I/O purposes.
AliAnalysisTaskDmesonJets::AliAnalysisTaskDmesonJets() :
  AliAnalysisTaskEmcal(),
  fAnalysisEngines(),
  fEnabledAxis(0),
  fHistManager(),
  fAodEvent(0),
  fFastJetWrapper(0)
{
  SetMakeGeneralHistograms(kTRUE);
}

/// This is the standard named constructor.
///
/// \param name Name of the task
AliAnalysisTaskDmesonJets::AliAnalysisTaskDmesonJets(const char* name) :
  AliAnalysisTaskEmcal(name, kTRUE),
  fAnalysisEngines(),
  fEnabledAxis(k2ProngInvMass),
  fHistManager(name),
  fAodEvent(0),
  fFastJetWrapper(0)
{
  SetMakeGeneralHistograms(kTRUE);
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
  AliJetDefinition jetDef(jettype, jetradius, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme);
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
AliAnalysisTaskDmesonJets::AnalysisEngine* AliAnalysisTaskDmesonJets::AddAnalysisEngine(ECandidateType_t type, EMCMode_t MCmode, const AliJetDefinition& jetDef, TString cutfname)
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

  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  // Define histograms
  // the TList fOutput is already defined in  AliAnalysisTaskEmcal::UserCreateOutputObjects()

  TString hname;
  TString htitle;
  TH1* h = 0;

  ::Info("AliAnalysisTaskDmesonJets::UserCreateOutputObjects", "Allocating histograms for task '%s' (%lu analysis engines)", GetName(), fAnalysisEngines.size());
  for (std::list<AnalysisEngine>::iterator it = fAnalysisEngines.begin(); it != fAnalysisEngines.end(); it++) {
    AnalysisEngine* param = &(*it);
    ::Info("AliAnalysisTaskDmesonJets::UserCreateOutputObjects", "Allocating histograms for analysis engine '%s' (%lu jet definitions)", param->GetName(), param->fJetDefinitions.size());

    fHistManager.CreateHistoGroup(param->GetName());

    param->fHistManager = &fHistManager;

    hname = TString::Format("%s/fHistNAcceptedDmesons", param->GetName());
    htitle = hname + ";Number of D accepted meson candidates;counts";
    h = fHistManager.CreateTH1(hname, htitle, 51, -0.5, 50.5);

    hname = TString::Format("%s/fHistNDmesons", param->GetName());
    htitle = hname + ";Number of D meson candidates;counts";
    h = fHistManager.CreateTH1(hname, htitle, 101, -0.5, 100.5);

    hname = TString::Format("%s/fHistNEvents", param->GetName());
    htitle = hname + ";Event status;counts";
    h = fHistManager.CreateTH1(hname, htitle, 2, 0, 2);
    h->GetXaxis()->SetBinLabel(1, "Accepted");
    h->GetXaxis()->SetBinLabel(2, "Rejected");

    hname = TString::Format("%s/fHistEventRejectionReasons", param->GetName());
    htitle = hname + ";Rejection reason;counts";
    h = fHistManager.CreateTH1(hname, htitle, 32, 0, 32);

    for (std::list<AliJetDefinition>::iterator itdef = param->fJetDefinitions.begin(); itdef != param->fJetDefinitions.end(); itdef++) {
      AliJetDefinition* jetDef = &(*itdef);
      ::Info("AliAnalysisTaskDmesonJets::UserCreateOutputObjects", "Allocating histograms for jet definition '%s'", jetDef->GetName());

      fHistManager.CreateHistoGroup(jetDef->GetName(), param->GetName());

      hname = TString::Format("%s/%s/fHistMCParticleRejectionReason", param->GetName(), jetDef->GetName());
      htitle = hname + ";Track rejection reason;#it{p}_{T,track} (GeV/#it{c});counts";
      h = fHistManager.CreateTH2(hname, htitle, 32, 0, 32, 150, 0, 150);
      SetRejectionReasonLabels(h->GetXaxis());

      hname = TString::Format("%s/%s/fHistTrackRejectionReason", param->GetName(), jetDef->GetName());
      htitle = hname + ";Track rejection reason;#it{p}_{T,track} (GeV/#it{c});counts";
      h = fHistManager.CreateTH2(hname, htitle, 32, 0, 32, 150, 0, 150);
      SetRejectionReasonLabels(h->GetXaxis());

      hname = TString::Format("%s/%s/fHistClusterRejectionReason", param->GetName(), jetDef->GetName());
      htitle = hname + ";Cluster rejection reason;#it{p}_{T,cluster} (GeV/#it{c});counts";
      h = fHistManager.CreateTH2(hname, htitle, 32, 0, 32, 150, 0, 150);
      SetRejectionReasonLabels(h->GetXaxis());

      hname = TString::Format("%s/%s/fHistDMesonDaughterNotInJet", param->GetName(), jetDef->GetName());
      htitle = hname + ";#it{p}_{T,track} (GeV/#it{c});counts";
      fHistManager.CreateTH1(hname, htitle, 400, 0, 100);
      SetRejectionReasonLabels(h->GetXaxis());
    }
    AllocateTHnSparse(*param);
  }

  fOutput->Add(fHistManager.GetListOfHistograms());

  PostData(1, fOutput);
}

/// Allocate a THnSparse histogram
///
/// \param param Analysis parameters used to properly set some of the axis
void AliAnalysisTaskDmesonJets::AllocateTHnSparse(const AnalysisEngine& param)
{
  TString hname;

  for (std::list<AliJetDefinition>::const_iterator it = param.fJetDefinitions.begin(); it != param.fJetDefinitions.end(); it++) {
    const AliJetDefinition* jetDef = &(*it);

    Printf("Now working on '%s'", jetDef->GetName());

    Double_t radius = jetDef->fRadius;

    TString  title[30] = {""};
    Int_t    nbins[30] = {0 };
    Double_t min  [30] = {0.};
    Double_t max  [30] = {0.};
    Int_t    dim       = 0   ;

    title[dim] = "#it{p}_{T,D} (GeV/#it{c})";
    nbins[dim] = fNbins;
    min[dim] = 0;
    max[dim] = 100;
    dim++;

    if ((fEnabledAxis & kPositionD) != 0) {
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

    if ((fEnabledAxis & kInvMass) != 0 && param.fCandidateType == kDstartoKpipi) {
      title[dim] = "#it{M}_{K#pi#pi} (GeV/#it{c}^{2})";
      nbins[dim] = param.fNMassBins;
      min[dim] = param.fMinMass;
      max[dim] = param.fMaxMass;
      dim++;
    }

    if (param.fCandidateType == kD0toKpi) {
      title[dim] = "#it{M}_{K#pi} (GeV/#it{c}^{2})";
      nbins[dim] = param.fNMassBins;
      min[dim] = param.fMinMass;
      max[dim] = param.fMaxMass;
      dim++;
    }

    if ((fEnabledAxis & k2ProngInvMass) != 0 && param.fCandidateType == kDstartoKpipi) {
      title[dim] = "#it{M}_{K#pi} (GeV/#it{c}^{2})";
      nbins[dim] = param.fNMassBins;
      CalculateMassLimits(param.fMaxMass - param.fMinMass, 421, param.fNMassBins, min[dim], max[dim]);
      dim++;
    }

    if (param.fCandidateType == kDstartoKpipi) {
      title[dim] = "#it{M}_{K#pi#pi} - #it{M}_{K#pi} (GeV/#it{c}^{2})";
      nbins[dim] = param.fNMassBins*6;
      CalculateMassLimits(0.20, 413, nbins[dim], min[dim], max[dim]);

      // subtract mass of D0
      Double_t D0mass = TDatabasePDG::Instance()->GetParticle(TMath::Abs(421))->Mass();
      min[dim] -= D0mass;
      max[dim] -= D0mass;

      dim++;
    }

    if ((fEnabledAxis & kSoftPionPt) != 0 && param.fCandidateType == kDstartoKpipi) {
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

    if ((fEnabledAxis & kDeltaR) != 0) {
      title[dim] = "#Delta R_{D-jet}";
      nbins[dim] = 100;
      min[dim] = 0;
      max[dim] = radius * 1.5;
      dim++;
    }

    if ((fEnabledAxis & kDeltaEta) != 0) {
      title[dim] = "#eta_{D} - #eta_{jet}";
      nbins[dim] = 100;
      min[dim] = -radius * 1.2;
      max[dim] = radius * 1.2;
      dim++;
    }

    if ((fEnabledAxis & kDeltaPhi) != 0) {
      title[dim] = "#phi_{D} - #phi_{jet} (rad)";
      nbins[dim] = 100;
      min[dim] = -radius * 1.2;
      max[dim] = radius * 1.2;
      dim++;
    }

    title[dim] = "#it{p}_{T,jet} (GeV/#it{c})";
    nbins[dim] = fNbins;
    min[dim] = fMinBinPt;
    max[dim] = fMaxBinPt;
    dim++;

    if ((fEnabledAxis & kPositionJet) != 0) {
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

    if ((fEnabledAxis & kLeadingPt) != 0) {
      title[dim] = "#it{p}_{T,particle}^{leading} (GeV/#it{c})";
      nbins[dim] = 120;
      min[dim] = 0;
      max[dim] = 120;
      dim++;
    }

    if ((fEnabledAxis & kJetConstituents) != 0) {
      title[dim] = "No. of constituents";
      nbins[dim] = 50;
      min[dim] = -0.5;
      max[dim] = 49.5;
      dim++;
    }

    if ((fEnabledAxis & kDaughterDistances) != 0) {
      for (Int_t j = 0; j < param.fNDaughters; j++) {
        title[dim] = Form("#Delta R_{d%d-jet}", j);
        nbins[dim] = 100;
        min[dim] = 0;
        max[dim] = 4;
        dim++;
      }
    }

    hname = TString::Format("%s/%s/fDmesonJets", param.GetName(), jetDef->GetName());
    THnSparse* h = fHistManager.CreateTHnSparse(hname,hname,dim,nbins,min,max);
    for (Int_t j = 0; j < dim; j++) {
      h->GetAxis(j)->SetTitle(title[j]);
    }
  }
}

/// Does some specific initializations for the analysis engines,
/// then calls the base class ExecOnce() method.
void AliAnalysisTaskDmesonJets::ExecOnce()
{
  AliAnalysisTaskEmcal::ExecOnce();

  // Load the event
  fAodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);

  fFastJetWrapper = new AliFJWrapper(fName, fTitle);

  fFastJetWrapper->SetAreaType(fastjet::active_area);
  fFastJetWrapper->SetGhostArea(1);

  if (!fAodEvent) {
     AliError(Form("This task need an AOD event! Task '%s' will be disabled!", GetName()));
     return;
  }

  for (std::list<AnalysisEngine>::iterator it = fAnalysisEngines.begin(); it != fAnalysisEngines.end(); it++) {
    AnalysisEngine* params = &(*it);

    params->fAodEvent = fAodEvent;
    params->fFastJetWrapper = fFastJetWrapper;

    if (params->fMCMode != kMCTruth) {
      params->fCandidateArray = dynamic_cast<TClonesArray*>(fAodEvent->GetList()->FindObject(params->fBranchName.Data()));

      if (params->fCandidateArray) {
        if (!params->fCandidateArray->GetClass()->InheritsFrom("AliAODRecoDecayHF2Prong")) {
          ::Error("AliAnalysisTaskDmesonJets::ExecOnce",
              "%s: Objects of type %s in %s are not inherited from AliAODRecoDecayHF2Prong! Task will be disabled!",
              GetName(), params->fCandidateArray->GetClass()->GetName(), params->fCandidateArray->GetName());
          params->fCandidateArray = 0;
          params->fInhibit = kTRUE;
        }
      }
      else {
        ::Error("AliAnalysisTaskDmesonJets::ExecOnce",
            "Could not find candidate array '%s', skipping the event. Analysis engine '%s' will be disabled!",
            params->fBranchName.Data(), params->GetName());
        params->fInhibit = kTRUE;
      }
    }

    if (params->fMCMode != kNoMC) {
      params->fMCContainer = dynamic_cast<AliHFAODMCParticleContainer*>(GetParticleContainer(0));

      if (!params->fMCContainer) params->fMCContainer = dynamic_cast<AliHFAODMCParticleContainer*>(GetParticleContainer(1));

      if (!params->fMCContainer) {
        ::Error("AliAnalysisTaskDmesonJets::ExecOnce",
            "No MC particle container was provided. Analysis engine '%s' will be disabled!",
            params->GetName());
        params->fInhibit = kTRUE;
      }
    }

    if (params->fMCMode != kMCTruth) {
      params->fTrackContainer = dynamic_cast<AliHFTrackContainer*>(GetParticleContainer(0));
      if (!params->fTrackContainer) params->fTrackContainer = dynamic_cast<AliHFTrackContainer*>(GetParticleContainer(1));

      params->fClusterContainer = GetClusterContainer(0);

      if (!params->fTrackContainer && !params->fClusterContainer) {
        ::Error("AliAnalysisTaskDmesonJets::ExecOnce",
            "No track container and no cluster container were provided. Analysis engine '%s' will be disabled!",
            params->GetName());
        params->fInhibit = kTRUE;
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

  for (std::list<AnalysisEngine>::iterator it = fAnalysisEngines.begin(); it != fAnalysisEngines.end(); it++) {
    AnalysisEngine* eng = &(*it);

    if (eng->fInhibit) continue;

    //Event selection
    hname = TString::Format("%s/fHistNEvents", eng->GetName());
    Bool_t iseventselected = eng->fRDHFCuts->IsEventSelected(fAodEvent);
    if (!iseventselected) {
      fHistManager.FillTH1(hname, "Rejected");
      hname = TString::Format("%s/fHistEventRejectionReasons", eng->GetName());
      UInt_t bitmap = eng->fRDHFCuts->GetEventRejectionBitMap();
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

    eng->RunAnalysis();
  }
  return kTRUE;
}

/// Fill the histograms.
///
/// \return Always kTRUE
Bool_t AliAnalysisTaskDmesonJets::FillHistograms()
{
  TString hname;
  Int_t ip = 0;
  for (std::list<AnalysisEngine>::iterator it = fAnalysisEngines.begin(); it != fAnalysisEngines.end(); it++) {
    AnalysisEngine* param = &(*it);

    if (param->fInhibit) continue;


    for (std::list<AliJetDefinition>::iterator itdef = param->fJetDefinitions.begin(); itdef != param->fJetDefinitions.end(); itdef++) {
      AliJetDefinition* jetDef = &(*itdef);

      Double_t radius = jetDef->fRadius;

      hname = TString::Format("%s/%s/fDmesonJets", param->GetName(), jetDef->GetName());
      THnSparse* h = static_cast<THnSparse*>(fHistManager.FindObject(hname));

      for (Int_t ij = 0; ij < jetDef->fDmesonJets.size(); ij++) {
        FillTHnSparse(h, jetDef->fDmesonJets[ij]);
      }
    }
    ip++;
  }
  return kTRUE;
}

/// Fill a THnSparse using information from a AliDmesonJetInfo object
///
/// \param h          Valid pointer to a THnSparse object
/// \param DmesonJet  Const reference to an AliDmesonJetInfo object
void AliAnalysisTaskDmesonJets::FillTHnSparse(THnSparse* h, const AliDmesonJetInfo& DmesonJet)
{
  // Fill the THnSparse histogram.

  Double_t contents[30] = {0.};

  Double_t z = 1.;
  Double_t deltaR = 1.;
  Double_t deltaPhi = 1.;
  Double_t deltaEta = 1.;

  if (DmesonJet.fJet.Pt() > 0) {
    TVector3 dvect = DmesonJet.fD.Vect();
    TVector3 jvect = DmesonJet.fJet.Vect();
    
    Double_t jetMom = jvect * jvect;

    if (jetMom < 1e-6) {
      ::Error("AliAnalysisTaskDmesonJets::FillTHnSparse", "Zero jet momentum!");
      z = 0.999;
    }
    else {
      z = (dvect * jvect) / jetMom;
    }

    if (z == 1 || (z > 1 && z - 1 < 1e-3)) z = 0.999; // so that it will contribute to the bin 0.9-1 rather than 1-1.1
    
    deltaPhi = TVector2::Phi_mpi_pi(DmesonJet.fD.Phi() - DmesonJet.fJet.Phi());
    deltaEta = DmesonJet.fD.Eta() - DmesonJet.fJet.Eta();
    
    deltaR = TMath::Sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);
  }
  
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
    else if (title=="#it{p}_{T,jet} (GeV/#it{c})")                   contents[i] = DmesonJet.fJet.Pt();
    else if (title=="#eta_{jet}")                                    contents[i] = DmesonJet.fJet.Eta();
    else if (title=="#phi_{jet} (rad)")                              contents[i] = DmesonJet.fJet.Phi_0_2pi();
    else if (title=="#it{p}_{T,particle}^{leading} (GeV/#it{c})")    contents[i] = DmesonJet.fJetLeadingPt;
    else if (title=="No. of constituents")                           contents[i] = DmesonJet.fJetNConstituents;
    else if (title=="#Delta R_{d0-jet}")                             contents[i] = DmesonJet.fDaughterDistances[0];
    else if (title=="#Delta R_{d1-jet}")                             contents[i] = DmesonJet.fDaughterDistances[1];
    else if (title=="#Delta R_{d2-jet}")                             contents[i] = DmesonJet.fDaughterDistances[2];
    else AliWarning(Form("Unable to fill dimension '%s'!",title.Data()));
  }

  h->Fill(contents);
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
