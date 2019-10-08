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
//copy paste the task of Salvatore with some additions for substructure
//C++
#include <sstream>
#include <array>

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

#include "AliAnalysisTaskDmesonJetsSub.h"

AliAnalysisTaskDmesonJetsSub::AliEventNotFound::AliEventNotFound(const std::string& class_name, const std::string& method_name) :
  std::exception(),
  fClassName(class_name),
  fAccessMethodName(method_name)
{
  std::stringstream what_str;
  what_str << "ALICE event not found in class '" <<  fClassName << "' using method '" << method_name << "'.";
  fWhat = what_str.str();
}

#if !(defined(__CINT__) || defined(__MAKECINT__))
const char* AliAnalysisTaskDmesonJetsSub::AliEventNotFound::what() const noexcept
{
  return fWhat.c_str();
}
#endif


// Definitions of class AliAnalysisTaskDmesonJetsSub::AliEventInfoSummary

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJetsSub::AliEventInfoSummary);
/// \endcond

/// Constructor that sets the object with the provided event information
///
/// \param cent Event centrality
/// \param ep Event plane
AliAnalysisTaskDmesonJetsSub::AliEventInfoSummary::AliEventInfoSummary(EventInfo event) :
  fWeight(1),
  fPtHard(0)
{
  Set(event);
}

/// Reset the object
void AliAnalysisTaskDmesonJetsSub::AliEventInfoSummary::Reset()
{
  fWeight = 1;
  fPtHard = 0;
}

/// Set the object with the provided event information
///
/// \param cent Event centrality
/// \param ep Event plane
void AliAnalysisTaskDmesonJetsSub::AliEventInfoSummary::Set(EventInfo event)
{
  fWeight = event.fWeight;
  fPtHard = event.fPtHard;
}


// Definitions of class AliAnalysisTaskDmesonJetsSub::AliJetInfo

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJetsSub::AliJetInfo);
/// \endcond

/// Calculates the distance between this jet and another jet
///
/// \param jet Const reference to a AliJetInfo object
/// \param deta reference where the eta distance will be returned
/// \param dphi reference where the phi distance will be returned
/// \return The distance between this jet and the jet reference provided
Double_t AliAnalysisTaskDmesonJetsSub::AliJetInfo::GetDistance(const AliJetInfo& jet, Double_t& deta, Double_t& dphi) const
{
  dphi = TVector2::Phi_mpi_pi(fMomentum.Phi() - jet.Phi());;
  deta = fMomentum.Eta() - jet.Eta();
  return TMath::Sqrt(dphi*dphi + deta*deta);
}

/// Calculates the distance between this jet and another jet
///
/// \param jet Const reference to a AliJetInfo object
/// \return The distance between this jet and the jet reference provided
Double_t AliAnalysisTaskDmesonJetsSub::AliJetInfo::GetDistance(const AliJetInfo& jet) const
{
  Double_t deta = 0;
  Double_t dphi = 0;
  return GetDistance(jet, deta, dphi);
}

// Definitions of class AliAnalysisTaskDmesonJetsSub::AliDmesonJetInfo

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJetsSub::AliDmesonJetInfo);
/// \endcond

/// Default constructor
AliAnalysisTaskDmesonJetsSub::AliDmesonJetInfo::AliDmesonJetInfo() :
  fDmesonParticle(0),
  fD(),
  fSoftPionPt(0),
  fInvMass2Prong(0),
  fJets(),
  fMCLabel(-1),
  fReconstructed(kFALSE),
  fParton(0),
  fPartonType(0),
  fAncestor(0),
  fD0D0bar(kFALSE),
  fSelectionType(0),
  fEvent(nullptr)
{
}

/// Copy constructor
///
/// \param source AliDmesonJetInfo object to copy from
AliAnalysisTaskDmesonJetsSub::AliDmesonJetInfo::AliDmesonJetInfo(const AliDmesonJetInfo &source) :
  fDmesonParticle(source.fDmesonParticle),
  fD(source.fD),
  fSoftPionPt(source.fSoftPionPt),
  fInvMass2Prong(source.fInvMass2Prong),
  fJets(source.fJets),
  fMCLabel(source.fMCLabel),
  fReconstructed(source.fReconstructed),
  fParton(source.fParton),
  fPartonType(source.fPartonType),
  fAncestor(source.fAncestor),
  fD0D0bar(source.fD0D0bar),
  fSelectionType(source.fSelectionType),
  fEvent(source.fEvent)
{
}

/// Assignment operator
///
/// \param source AliDmesonJetInfo object to copy from
AliAnalysisTaskDmesonJetsSub::AliDmesonJetInfo& AliAnalysisTaskDmesonJetsSub::AliDmesonJetInfo::operator=(const AliDmesonJetInfo& source)
{
  new (this) AliDmesonJetInfo(source);
  return *this;
}

/// Reset all fields to their default values
void AliAnalysisTaskDmesonJetsSub::AliDmesonJetInfo::Reset()
{
  fD.SetPtEtaPhiE(0,0,0,0);
  fSoftPionPt = 0;
  fInvMass2Prong = 0;
  fDmesonParticle = 0;
  fMCLabel = -1;
  fReconstructed = kFALSE;
  fParton = 0;
  fPartonType = 0;
  fAncestor = 0;
  fD0D0bar = kFALSE;
  for (auto &jet : fJets) {
    jet.second.fMomentum.SetPtEtaPhiE(0,0,0,0);
    jet.second.fNConstituents = 0;
    jet.second.fNEF = 0;
    jet.second.fMaxChargedPt = 0;
    jet.second.fMaxNeutralPt = 0;
  }
}

/// Prints the content of this object in the standard output.
void AliAnalysisTaskDmesonJetsSub::AliDmesonJetInfo::Print() const
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
Double_t AliAnalysisTaskDmesonJetsSub::AliDmesonJetInfo::GetZ(std::string n) const
{
  std::map<std::string, AliJetInfo>::const_iterator it = fJets.find(n);
  if (it == fJets.end()) return 0;

  Double_t z = 0;

  if ((*it).second.Pt() > 0) {
    TVector3 dvect = fD.Vect();
    TVector3 jvect = (*it).second.fMomentum.Vect();

    Double_t jetMom = jvect * jvect;

    if (jetMom < 1e-6) {
      ::Error("AliAnalysisTaskDmesonJetsSub::AliDmesonJetInfo::GetZ", "Zero jet momentum!");
      z = 0.999;
    }
    else {
      z = (dvect * jvect) / jetMom;
    }

    if (z == 1 || (z > 1 && z - 1 < 1e-3)) z = 0.999; // so that it will contribute to the bin 0.9-1 rather than 1-1.1
  }

  return z;
}

/// Calculates the parallel fraction
///
/// \return the fraction of the momentum of the particle parallel to the jet over the total jet momentum
Double_t AliAnalysisTaskDmesonJetsSub::AliDmesonJetInfo::GetCorrZ(std::string n) const
{
  std::map<std::string, AliJetInfo>::const_iterator it = fJets.find(n);
  if (it == fJets.end()) return 0;

  Double_t z = 0;

  if ((*it).second.Pt() > 0) {
    TVector3 dvect = fD.Vect();
    TVector3 jvect = (*it).second.fMomentum.Vect();
    // If the corr pt is < 0, assign 0.
    Double_t corrpt = (*it).second.fCorrPt > 0 ? (*it).second.fCorrPt : 0.;
    jvect.SetPerp(corrpt);

    Double_t jetMom = jvect * jvect;

    if (jetMom < 1e-6) {
      z = 1.0;
    }
    else {
      z = (dvect * jvect) / jetMom;
    }
  }

  return z;
}

/// Calculates the distance between the D meson and the jet axis
///
/// \param n Name of the jet definition
/// \param deta reference where the eta distance will be returned
/// \param dphi reference where the phi distance will be returned
/// \return The distance between the D meson and the jet axis
Double_t AliAnalysisTaskDmesonJetsSub::AliDmesonJetInfo::GetDistance(std::string n, Double_t& deta, Double_t& dphi) const
{
  std::map<std::string, AliJetInfo>::const_iterator it = fJets.find(n);
  if (it == fJets.end()) return 0;

  return GetDistance((*it).second, deta, dphi);
}

/// Calculates the distance between the D meson and the jet axis
///
/// \param n Name of the jet definition
/// \return The distance between the D meson and the jet axis
Double_t AliAnalysisTaskDmesonJetsSub::AliDmesonJetInfo::GetDistance(std::string n) const
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
Double_t AliAnalysisTaskDmesonJetsSub::AliDmesonJetInfo::GetDistance(const AliJetInfo& jet, Double_t& deta, Double_t& dphi) const
{
  dphi = TVector2::Phi_mpi_pi(fD.Phi() - jet.Phi());;
  deta = fD.Eta() - jet.Eta();
  return TMath::Sqrt(dphi*dphi + deta*deta);
}

/// Calculates the distance between the D meson and the jet axis
///
/// \param jet Const reference to a AliJetInfo object
/// \return The distance between the D meson and the jet axis
Double_t AliAnalysisTaskDmesonJetsSub::AliDmesonJetInfo::GetDistance(const AliJetInfo& jet) const
{
  Double_t deta = 0;
  Double_t dphi = 0;
  return GetDistance(jet, deta, dphi);
}

/// Find jet info object corresponding a jet definition provided as a string
///
/// \param n String containing the jet definition
/// \return Constant pointer to the jet info object, if the jet definition was found. Null pointer otherwise.
const AliAnalysisTaskDmesonJetsSub::AliJetInfo* AliAnalysisTaskDmesonJetsSub::AliDmesonJetInfo::GetJet(std::string n) const
{
  std::map<std::string, AliJetInfo>::const_iterator it = fJets.find(n);
  if (it == fJets.end()) {
    ::Error("AliAnalysisTaskDmesonJetsSub::AliDmesonJetInfo::GetJet", "Could not find jet info for the jet definition '%s'!",
       n.c_str());
    return 0;
  }
  return &((*it).second);
}

/// Find jet info object corresponding a jet definition provided as a string
///
/// \param n String containing the jet definition
/// \return Pointer to the jet info object, if the jet definition was found. Null pointer otherwise.
AliAnalysisTaskDmesonJetsSub::AliJetInfo* AliAnalysisTaskDmesonJetsSub::AliDmesonJetInfo::GetJet(std::string n)
{
  std::map<std::string, AliJetInfo>::iterator it = fJets.find(n);
  if (it == fJets.end()) {
    ::Error("AliAnalysisTaskDmesonJetsSub::AliDmesonJetInfo::GetJet", "Could not find jet info for the jet definition '%s'!",
        n.c_str());
    return 0;
  }
  return &((*it).second);
}

// Definitions of class AliAnalysisTaskDmesonJetsSub::AliJetInfoSummary

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJetsSub::AliJetInfoSummary);
/// \endcond

/// Constructor that uses an AliDmesonJetInfo as its source
///
/// \param source A const reference to a valid AliDmesonJetInfo object
/// \param i      Index of the jet to be copied
AliAnalysisTaskDmesonJetsSub::AliJetInfoSummary::AliJetInfoSummary(const AliDmesonJetInfo& source, std::string n) :
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
void AliAnalysisTaskDmesonJetsSub::AliJetInfoSummary::Reset()
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
void AliAnalysisTaskDmesonJetsSub::AliJetInfoSummary::Set(const AliDmesonJetInfo& source, std::string n)
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
void AliAnalysisTaskDmesonJetsSub::AliJetInfoSummary::Set(const AliJetInfo& source)
{
  fPt = source.Pt();
  fEta = source.Eta();
  fPhi = source.Phi_0_2pi();
  fN = source.GetNConstituents();
  fR = 0;
  fZ = 0;
}

// Definitions of class AliAnalysisTaskDmesonJetsSub::AliJetInfoPbPbSummary

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJetsSub::AliJetInfoPbPbSummary);
/// \endcond

/// Constructor that uses an AliDmesonJetInfo as its source
///
/// \param source A const reference to a valid AliDmesonJetInfo object
/// \param i      Index of the jet to be copied
AliAnalysisTaskDmesonJetsSub::AliJetInfoPbPbSummary::AliJetInfoPbPbSummary(const AliDmesonJetInfo& source, std::string n) :
  AliJetInfoSummary(),
  fCorrPt(0),
  fCorrZ(0),
  fArea(0)
{
  Set(source, n);
}

/// Reset the current object
void AliAnalysisTaskDmesonJetsSub::AliJetInfoPbPbSummary::Reset()
{
  AliJetInfoSummary::Reset();
  fCorrPt = 0;
  fCorrZ = 0;
  fArea = 0;
}

/// Set the current object using an instance of AliJetInfo as its source
///
/// \param source A const reference to a valid AliJetInfo object
void AliAnalysisTaskDmesonJetsSub::AliJetInfoPbPbSummary::Set(const AliJetInfo& source)
{
  AliJetInfoSummary::Set(source);
  fArea = source.fArea;
  fCorrPt = source.fCorrPt;
}

/// Set the current object using an instance of AliDmesonJetInfo as its source
///
/// \param source A const reference to a valid AliDmesonJetInfo object
/// \param i      Index of the jet to be copied
void AliAnalysisTaskDmesonJetsSub::AliJetInfoPbPbSummary::Set(const AliDmesonJetInfo& source, std::string n)
{
  AliJetInfoSummary::Set(source, n);
  fCorrZ = source.GetCorrZ(n);
}

// Definitions of class AliAnalysisTaskDmesonJetsSub::AliDmesonInfoSummary

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJetsSub::AliDmesonInfoSummary);
/// \endcond

/// Constructor that uses an AliDmesonJetInfo as its source
///
/// \param source A const reference to a valid AliDmesonJetInfo object
AliAnalysisTaskDmesonJetsSub::AliDmesonInfoSummary::AliDmesonInfoSummary(const AliDmesonJetInfo& source) :
  fPt(0),
  fEta(0),
  fPhi(0)
{
  Set(source);
}

/// Set the current object using an instance of AliDmesonJetInfo as its source
///
/// \param source A const reference to a valid AliDmesonJetInfo object
void AliAnalysisTaskDmesonJetsSub::AliDmesonInfoSummary::Set(const AliDmesonJetInfo& source)
{
  fPt = source.fD.Pt();
  fEta = source.fD.Eta();
  fPhi = source.fD.Phi_0_2pi();
}

/// Reset the object
void AliAnalysisTaskDmesonJetsSub::AliDmesonInfoSummary::Reset()
{
  fPt = 0;
  fEta = 0;
  fPhi = 0;
}

// Definitions of class AliAnalysisTaskDmesonJetsSub::AliDmesonMCInfoSummary

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJetsSub::AliDmesonMCInfoSummary);
/// \endcond

/// Constructor that uses an AliDmesonJetInfo as its source
///
/// \param source A const reference to a valid AliDmesonJetInfo object
AliAnalysisTaskDmesonJetsSub::AliDmesonMCInfoSummary::AliDmesonMCInfoSummary(const AliDmesonJetInfo& source) :
  AliDmesonInfoSummary(source),
  fPartonType(0),
  fPartonPt(0),
  fAncestorPDG(0)
{
  Set(source);
}

/// Set the current object using an instance of AliDmesonJetInfo as its source
///
/// \param source A const reference to a valid AliDmesonJetInfo object
void AliAnalysisTaskDmesonJetsSub::AliDmesonMCInfoSummary::Set(const AliDmesonJetInfo& source)
{
  AliDmesonInfoSummary::Set(source);

  fPartonType = source.fPartonType;

  if (source.fParton) {
    fPartonPt = source.fParton->Pt();
  }
  else {
    fPartonPt = 0.;
  }

  if (source.fAncestor) {
    fAncestorPDG = (UShort_t)((UInt_t)(TMath::Abs(source.fAncestor->GetPdgCode())));
  }
}

/// Reset the object
void AliAnalysisTaskDmesonJetsSub::AliDmesonMCInfoSummary::Reset()
{
  AliDmesonInfoSummary::Reset();
  fPartonType = 0,
  fPartonPt = 0.;
  fAncestorPDG = 0;
}

// Definitions of class AliAnalysisTaskDmesonJetsSub::AliD0InfoSummary

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJetsSub::AliD0InfoSummary);
/// \endcond

/// Constructor that uses an AliDmesonJetInfo as its source
///
/// \param source A const reference to a valid AliDmesonJetInfo object
AliAnalysisTaskDmesonJetsSub::AliD0InfoSummary::AliD0InfoSummary(const AliDmesonJetInfo& source) :
  AliDmesonInfoSummary(source),
  fInvMass(source.fD.M()),
  fSelectionType(0)
{
}

/// Set the current object using an instance of AliDmesonJetInfo as its source
///
/// \param source A const reference to a valid AliDmesonJetInfo object
void AliAnalysisTaskDmesonJetsSub::AliD0InfoSummary::Set(const AliDmesonJetInfo& source)
{
  fInvMass = source.fD.M();
  fSelectionType = source.GetSelectionTypeSummary();
  AliDmesonInfoSummary::Set(source);
}

/// Reset the object
void AliAnalysisTaskDmesonJetsSub::AliD0InfoSummary::Reset()
{
  AliDmesonInfoSummary::Reset();
  fSelectionType = 0;
  fInvMass = 0;
}

// Definitions of class AliAnalysisTaskDmesonJetsSub::AliD0ExtendedInfoSummary

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJetsSub::AliD0ExtendedInfoSummary);
/// \endcond

/// Constructor that uses an AliDmesonJetInfo as its source
///
/// \param source A const reference to a valid AliDmesonJetInfo object
AliAnalysisTaskDmesonJetsSub::AliD0ExtendedInfoSummary::AliD0ExtendedInfoSummary(const AliDmesonJetInfo& source) :
  AliD0InfoSummary(source),
  fDCA(0),
  fCosThetaStar(0),
  fd0K(0),
  fd0Pi(0),
  fd0d0(0),
  fCosPointing(0),
  fMaxNormd0(0)
{
  Set(source);
}

/// Set the current object using an instance of AliDmesonJetInfo as its source
///
/// \param source A const reference to a valid AliDmesonJetInfo object
void AliAnalysisTaskDmesonJetsSub::AliD0ExtendedInfoSummary::Set(const AliDmesonJetInfo& source)
{
  AliD0InfoSummary::Set(source);

  AliAODRecoDecayHF2Prong* recoDecay = dynamic_cast<AliAODRecoDecayHF2Prong*>(source.fDmesonParticle);
  if (recoDecay) {
    fDCA = recoDecay->GetDCA();
    if (source.fSelectionType == 1) { // D0
      fCosThetaStar = recoDecay->CosThetaStarD0();
      fPtK = recoDecay->PtProng(0);
      fPtPi = recoDecay->PtProng(1);
      fd0K = recoDecay->Getd0Prong(0);
      fd0Pi = recoDecay->Getd0Prong(1);
    }
    else { //D0bar
      fCosThetaStar = recoDecay->CosThetaStarD0bar();
      fPtK = recoDecay->PtProng(1);
      fPtPi = recoDecay->PtProng(0);
      fd0K = recoDecay->Getd0Prong(1);
      fd0Pi = recoDecay->Getd0Prong(0);
    }

    fMaxNormd0 = 0.;
    // Based on Int_t AliRDHFCutsD0toKpi::IsSelected(TObject* obj,Int_t selectionLevel,AliAODEvent* aod)
    // Line 480 and following
    if (source.fEvent) {
      for (Int_t ipr=0; ipr < 2; ipr++) {
        Double_t diffIP = 0., errdiffIP = 0.;
        recoDecay->Getd0MeasMinusExpProng(ipr, source.fEvent->GetMagneticField(), diffIP, errdiffIP);
        Double_t normdd0 = 0.;
        if (errdiffIP > 0.) {
          normdd0 = diffIP / errdiffIP;
        }
        else {
          if (diffIP == 0) {
            normdd0 = 0;
          }
          else {
            normdd0 = diffIP > 0 ? 9999. : -9999.;
          }
        }
        if (TMath::Abs(normdd0) > TMath::Abs(fMaxNormd0)) {
          fMaxNormd0 = normdd0;
        }
      }
    }
    else {
      throw AliAnalysisTaskDmesonJetsSub::AliEventNotFound("AliAnalysisTaskDmesonJetsSub::AliDmesonJetInfo", "fEvent");
    }

    fd0d0 = recoDecay->Prodd0d0();
    fCosPointing = recoDecay->CosPointingAngle();
  }
}

/// Reset the object
void AliAnalysisTaskDmesonJetsSub::AliD0ExtendedInfoSummary::Reset()
{
  AliD0InfoSummary::Reset();
  fDCA = 0;
  fCosThetaStar = 0;
  fd0K = 0;
  fd0Pi = 0;
  fd0d0 = 0;
  fCosPointing = 0;
  fMaxNormd0 = 0;
}

// Definitions of class AliAnalysisTaskDmesonJetsSub::AliDStarInfoSummary

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJetsSub::AliDStarInfoSummary);
/// \endcond

/// Constructor that uses an AliDmesonJetInfo as its source
///
/// \param source A const reference to a valid AliDmesonJetInfo object
AliAnalysisTaskDmesonJetsSub::AliDStarInfoSummary::AliDStarInfoSummary(const AliDmesonJetInfo& source) :
  AliDmesonInfoSummary(source),
  f2ProngInvMass(source.fInvMass2Prong),
  fDeltaInvMass(source.fD.M() - source.fInvMass2Prong)
{
}

/// Set the current object using an instance of AliDmesonJetInfo as its source
///
/// \param source A const reference to a valid AliDmesonJetInfo object
void AliAnalysisTaskDmesonJetsSub::AliDStarInfoSummary::Set(const AliDmesonJetInfo& source)
{
  f2ProngInvMass = source.fInvMass2Prong;
  fDeltaInvMass = source.fD.M() - source.fInvMass2Prong;
  AliDmesonInfoSummary::Set(source);
}

/// Reset the object
void AliAnalysisTaskDmesonJetsSub::AliDStarInfoSummary::Reset()
{
  AliDmesonInfoSummary::Reset();

  f2ProngInvMass = 0;
  fDeltaInvMass = 0;
}

// Definitions of class AliAnalysisTaskDmesonJetsSub::OutputHandler

/// Constructor
AliAnalysisTaskDmesonJetsSub::OutputHandler::OutputHandler() :
  fCandidateType(kD0toKpi),
  fMCMode(kNoMC),
  fNMassBins(0),
  fMinMass(0),
  fMaxMass(0),
  fJetDefinitions(nullptr),
  fPtBinWidth(0.5),
  fMaxPt(100),
  fD0Extended(kFALSE),
  fEventInfo(nullptr),
  fDmesonJets(nullptr),
  fHistManager(nullptr),
  fName()
{
}

/// Constructor
AliAnalysisTaskDmesonJetsSub::OutputHandler::OutputHandler(AnalysisEngine* eng) :
  fCandidateType(eng->fCandidateType),
  fMCMode(eng->fMCMode),
  fNMassBins(eng->fNMassBins),
  fMinMass(eng->fMinMass),
  fMaxMass(eng->fMaxMass),
  fJetDefinitions(&eng->fJetDefinitions),
  fPtBinWidth(eng->fPtBinWidth),
  fMaxPt(eng->fMaxPt),
  fD0Extended(eng->fD0Extended),
  fEventInfo(&eng->fEventInfo),
  fDmesonJets(&eng->fDmesonJets),
  fHistManager(eng->fHistManager),
  fName(eng->GetName())
{
}

/// Fills QA histograms. This method is not used by the AliAnalysisTaskDmesonJetsSub task,
/// but can be used by derived tasks that have a custom implementation to fill the output objects.
///
/// \return Always kTRUE
Bool_t AliAnalysisTaskDmesonJetsSub::OutputHandler::FillOutput(Bool_t applyKinCuts)
{
  TString hname;

  TH1* histAncestor = nullptr;
  TH1* histPrompt = nullptr;

  if (fMCMode == kSignalOnly || fMCMode == kMCTruth) {
    hname = TString::Format("%s/fHistPrompt", fName.Data());
    histPrompt = static_cast<TH1*>(fHistManager->FindObject(hname));

    hname = TString::Format("%s/fHistAncestor", fName.Data());
    histAncestor = static_cast<TH1*>(fHistManager->FindObject(hname));
  }

  std::map<AliAODMCParticle*, Short_t> partons ; // set of the partons in the shower that produced each D meson
  for (auto& dmeson_pair : *fDmesonJets) {
    Int_t accJets = 0;
    for (UInt_t ij = 0; ij < fJetDefinitions->size(); ij++) {
      AliJetInfo* jet = dmeson_pair.second.GetJet(fJetDefinitions->at(ij).GetName());
      if (!jet) continue;
      if (applyKinCuts && !fJetDefinitions->at(ij).IsJetInAcceptance(*jet)) {
        hname = TString::Format("%s/%s/fHistRejectedJetPt", fName.Data(), fJetDefinitions->at(ij).GetName());
        fHistManager->FillTH1(hname, jet->Pt());
        hname = TString::Format("%s/%s/fHistRejectedJetPhi", fName.Data(), fJetDefinitions->at(ij).GetName());
        fHistManager->FillTH1(hname, jet->Phi_0_2pi());
        hname = TString::Format("%s/%s/fHistRejectedJetEta", fName.Data(), fJetDefinitions->at(ij).GetName());
        fHistManager->FillTH1(hname, jet->Eta());
        continue;
      }
      accJets++;
    }
    if (accJets > 0) {
      if (histPrompt) {
        if (dmeson_pair.second.fParton) {
          partons[dmeson_pair.second.fParton] = dmeson_pair.second.fPartonType;
          UInt_t absPdgParton = TMath::Abs(dmeson_pair.second.fParton->GetPdgCode());
          if (absPdgParton == 4) {
            histPrompt->Fill("Prompt", 1);
          }
          else if (absPdgParton == 5) {
            histPrompt->Fill("Non-Prompt", 1);
          }
          else {
            histPrompt->Fill("Unknown", 1);
          }
        }
        else {
          histPrompt->Fill("Unknown", 1);
        }
      }

      if (histAncestor) {
        if (dmeson_pair.second.fAncestor) {
          UInt_t absPdgAncestor = TMath::Abs(dmeson_pair.second.fAncestor->GetPdgCode());
          if (absPdgAncestor == 4) {
            histAncestor->Fill("Charm", 1);
          }
          else if (absPdgAncestor == 5) {
            histAncestor->Fill("Bottom", 1);
          }
          else if (absPdgAncestor == 2212) {
            histAncestor->Fill("Proton", 1);
          }
          else {
            histAncestor->Fill("Unknown", 1);
          }
        }
        else {
          histAncestor->Fill("Unknown", 1);
        }
      }
    }
    else {
      hname = TString::Format("%s/fHistRejectedDMesonPt", fName.Data());
      fHistManager->FillTH1(hname, dmeson_pair.second.fD.Pt());
      hname = TString::Format("%s/fHistRejectedDMesonPhi", fName.Data());
      fHistManager->FillTH1(hname, dmeson_pair.second.fD.Phi_0_2pi());
      hname = TString::Format("%s/fHistRejectedDMesonEta", fName.Data());
      fHistManager->FillTH1(hname, dmeson_pair.second.fD.Eta());
      if (fMCMode != kMCTruth) {
        if (fCandidateType == kD0toKpi || fCandidateType == kD0toKpiLikeSign) {
          hname = TString::Format("%s/fHistRejectedDMesonInvMass", fName.Data());
          fHistManager->FillTH1(hname, dmeson_pair.second.fD.M());
        }
        else if (fCandidateType == kDstartoKpipi) {
          hname = TString::Format("%s/fHistRejectedDMeson2ProngInvMass", fName.Data());
          fHistManager->FillTH1(hname, dmeson_pair.second.fInvMass2Prong);

          hname = TString::Format("%s/fHistRejectedDMesonDeltaInvMass", fName.Data());
          fHistManager->FillTH1(hname, dmeson_pair.second.fD.M() - dmeson_pair.second.fInvMass2Prong);
        }
      }
    }
  }

  if (fMCMode == kSignalOnly || fMCMode == kMCTruth) {
    hname = TString::Format("%s/fHistPartonPt", fName.Data());
    TH1* histPartonPt = static_cast<TH1*>(fHistManager->FindObject(hname));
    hname = TString::Format("%s/fHistPartonEta", fName.Data());
    TH1* histPartonEta = static_cast<TH1*>(fHistManager->FindObject(hname));
    hname = TString::Format("%s/fHistPartonPhi", fName.Data());
    TH1* histPartonPhi = static_cast<TH1*>(fHistManager->FindObject(hname));
    hname = TString::Format("%s/fHistPartonType", fName.Data());
    TH1* histPartonType = static_cast<TH1*>(fHistManager->FindObject(hname));

    for (auto parton : partons) {
      if (!parton.first) continue;
      histPartonPt->Fill(parton.first->Pt());
      histPartonEta->Fill(parton.first->Eta());
      histPartonPhi->Fill(TVector2::Phi_0_2pi(parton.first->Phi()));
      histPartonType->Fill(parton.second);
    }
  }

  return kTRUE;
}

// Definitions of class AliAnalysisTaskDmesonJetsSub::OutputHandlerTHnSparse

/// Constructor
AliAnalysisTaskDmesonJetsSub::OutputHandlerTHnSparse::OutputHandlerTHnSparse() :
  OutputHandler(),
  fEnabledAxis(0)
{
}

/// Constructor
AliAnalysisTaskDmesonJetsSub::OutputHandlerTHnSparse::OutputHandlerTHnSparse(AnalysisEngine* eng) :
  OutputHandler(eng),
  fEnabledAxis(0)
{
}


/// Allocate a THnSparse histogram
///
/// \param param Analysis parameters used to properly set some of the axis
void AliAnalysisTaskDmesonJetsSub::OutputHandlerTHnSparse::BuildOutputObject(const char* /*taskName*/)
{
  TString hname;

  Int_t nPtBins = TMath::CeilNint(fMaxPt / fPtBinWidth);

  for (auto &jetDef : *fJetDefinitions) {

    AliDebugGeneralStream("AliAnalysisTaskDmesonJetsSub::OutputHandlerTHnSparse::BuildOutputObject", 2) << "Now working on '" << jetDef.GetName() << "'" << std::endl;

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

    if ((fEnabledAxis & kInvMass) != 0 && fCandidateType == kDstartoKpipi) {
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

    if ((fEnabledAxis & k2ProngInvMass) != 0 && fCandidateType == kDstartoKpipi) {
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

    if ((fEnabledAxis & kSoftPionPt) != 0 && fCandidateType == kDstartoKpipi) {
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
    nbins[dim] = nPtBins;
    min[dim] = 0;
    max[dim] = fMaxPt;
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

    if ((fEnabledAxis & kJetConstituents) != 0) {
      title[dim] = "No. of constituents";
      nbins[dim] = 50;
      min[dim] = -0.5;
      max[dim] = 49.5;
      dim++;
    }

    hname = TString::Format("%s/%s/fDmesonJets", fName.Data(), jetDef.GetName());
    THnSparse* h = fHistManager->CreateTHnSparse(hname,hname,dim,nbins,min,max);
    for (Int_t j = 0; j < dim; j++) {
      h->GetAxis(j)->SetTitle(title[j]);
    }
  }
}

/// Checks whether any of the D meson jets is in the acceptance
///
/// \param Const reference to a valid AliDmesonJetInfo object
Bool_t AliAnalysisTaskDmesonJetsSub::OutputHandlerTHnSparse::IsAnyJetInAcceptance(const AliDmesonJetInfo& dMesonJet) const
{
  for (UInt_t i = 0; i < fJetDefinitions->size(); i++) {
    if (fJetDefinitions->at(i).IsJetInAcceptance(dMesonJet, fJetDefinitions->at(i).GetName())) return kTRUE;
  }

  return kFALSE;
}

/// Post the output with D meson jets found in the current event
///
/// \return kTRUE on success
Bool_t AliAnalysisTaskDmesonJetsSub::OutputHandlerTHnSparse::FillOutput(Bool_t applyKinCuts)
{
  TString hname;

  for (auto& dmeson_pair : *fDmesonJets) {
    if (!IsAnyJetInAcceptance(dmeson_pair.second)) {
      hname = TString::Format("%s/fHistRejectedDMesonPt", fName.Data());
      fHistManager->FillTH1(hname, dmeson_pair.second.fD.Pt());
      hname = TString::Format("%s/fHistRejectedDMesonPhi", fName.Data());
      fHistManager->FillTH1(hname, dmeson_pair.second.fD.Phi_0_2pi());
      hname = TString::Format("%s/fHistRejectedDMesonEta", fName.Data());
      fHistManager->FillTH1(hname, dmeson_pair.second.fD.Eta());
    }
  }

  for (auto &jetDef : *fJetDefinitions) {

    hname = TString::Format("%s/%s/fDmesonJets", fName.Data(), jetDef.GetName());
    THnSparse* h = static_cast<THnSparse*>(fHistManager->FindObject(hname));

    for (auto& dmeson_pair : *fDmesonJets) {
      const AliJetInfo* jet = dmeson_pair.second.GetJet(jetDef.GetName());
      if (!jet) continue;
      if (!jetDef.IsJetInAcceptance(*jet)) {
        hname = TString::Format("%s/%s/fHistRejectedJetPt", fName.Data(), jetDef.GetName());
        fHistManager->FillTH1(hname, jet->Pt());
        hname = TString::Format("%s/%s/fHistRejectedJetPhi", fName.Data(), jetDef.GetName());
        fHistManager->FillTH1(hname, jet->Phi_0_2pi());
        hname = TString::Format("%s/%s/fHistRejectedJetEta", fName.Data(), jetDef.GetName());
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
Bool_t AliAnalysisTaskDmesonJetsSub::OutputHandlerTHnSparse::FillHnSparse(THnSparse* h, const AliDmesonJetInfo& DmesonJet, std::string n)
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
    else {
      AliWarningGeneralStream("AliAnalysisTaskDmesonJetsSub::OutputHandlerTHnSparse::FillHnSparse") << "Unable to fill dimension '" << title.Data() << "'!" << std::endl;
    }
  }

  h->Fill(contents);

  return kTRUE;
}

// Definitions of class AliAnalysisTaskDmesonJetsSub::OutputHandlerTTree

/// Constructor
AliAnalysisTaskDmesonJetsSub::OutputHandlerTTree::OutputHandlerTTree() :
  OutputHandler(),
  fDataSlotNumber(-1),
  fTree(0),
  fCurrentDmesonJetInfo(0),
  fCurrentJetInfo(0)
{
}

/// Constructor
AliAnalysisTaskDmesonJetsSub::OutputHandlerTTree::OutputHandlerTTree(AnalysisEngine* eng) :
  OutputHandler(eng),
  fDataSlotNumber(-1),
  fTree(0),
  fCurrentDmesonJetInfo(0),
  fCurrentJetInfo(0)
{
}

/// Builds the tree where the output will be posted
///
/// \return Pointer to the new tree
void AliAnalysisTaskDmesonJetsSub::OutputHandlerTTree::BuildOutputObject(const char* taskName)
{
  TString classname;
  if (fMCMode == kMCTruth) {
    classname = "AliAnalysisTaskDmesonJetsSub::AliDmesonMCInfoSummary";
    fCurrentDmesonJetInfo = new AliDmesonMCInfoSummary();
  }
  else {
    switch (fCandidateType) {
    case kD0toKpi:
    case kD0toKpiLikeSign:
      if (fD0Extended) {
        classname = "AliAnalysisTaskDmesonJetsSub::AliD0ExtendedInfoSummary";
        fCurrentDmesonJetInfo = new AliD0ExtendedInfoSummary();
      }
      else {
        classname = "AliAnalysisTaskDmesonJetsSub::AliD0InfoSummary";
        fCurrentDmesonJetInfo = new AliD0InfoSummary();
      }
      break;
    case kDstartoKpipi:
      classname = "AliAnalysisTaskDmesonJetsSub::AliDStarInfoSummary";
      fCurrentDmesonJetInfo = new AliDStarInfoSummary();
      break;
    }
  }
  TString treeName = TString::Format("%s_%s", taskName, fName.Data());
  fTree = new TTree(treeName, treeName);
  fTree->Branch("DmesonJet", classname, &fCurrentDmesonJetInfo);
  fCurrentJetInfo = new AliJetInfoSummary*[fJetDefinitions->size()];
  for (Int_t i = 0; i < fJetDefinitions->size(); i++) {
    if (fJetDefinitions->at(i).fRhoName.IsNull()) {
      fCurrentJetInfo[i] = new AliJetInfoSummary();
      fTree->Branch(fJetDefinitions->at(i).GetName(), "AliAnalysisTaskDmesonJetsSub::AliJetInfoSummary", &fCurrentJetInfo[i]);
    }
    else {
      fCurrentJetInfo[i] = new AliJetInfoPbPbSummary();
      fTree->Branch(fJetDefinitions->at(i).GetName(), "AliAnalysisTaskDmesonJetsSub::AliJetInfoPbPbSummary", &fCurrentJetInfo[i]);
    }
  }
}

/// Post the output with D meson jets found in the current event
///
/// \return kTRUE on success
Bool_t AliAnalysisTaskDmesonJetsSub::OutputHandlerTTree::FillOutput(Bool_t applyKinCuts)
{
  TString hname;

  std::map<AliAODMCParticle*, Short_t> partons ; //!<! set of the partons in the shower that produced each D meson

  TH1* histAncestor = nullptr;
  TH1* histPrompt = nullptr;

  if (fMCMode == kSignalOnly || fMCMode == kMCTruth) {
    hname = TString::Format("%s/fHistPrompt", fName.Data());
    histPrompt = static_cast<TH1*>(fHistManager->FindObject(hname));

    hname = TString::Format("%s/fHistAncestor", fName.Data());
    histAncestor = static_cast<TH1*>(fHistManager->FindObject(hname));
  }

  for (auto& dmeson_pair : *fDmesonJets) {
    fCurrentDmesonJetInfo->Set(dmeson_pair.second);
    Int_t accJets = 0;
    for (UInt_t ij = 0; ij < fJetDefinitions->size(); ij++) {
      fCurrentJetInfo[ij]->Reset();
      AliJetInfo* jet = dmeson_pair.second.GetJet(fJetDefinitions->at(ij).GetName());
      if (!jet) continue;
      if (applyKinCuts && !fJetDefinitions->at(ij).IsJetInAcceptance(*jet)) {
        hname = TString::Format("%s/%s/fHistRejectedJetPt", fName.Data(), fJetDefinitions->at(ij).GetName());
        fHistManager->FillTH1(hname, jet->Pt());
        hname = TString::Format("%s/%s/fHistRejectedJetPhi", fName.Data(), fJetDefinitions->at(ij).GetName());
        fHistManager->FillTH1(hname, jet->Phi_0_2pi());
        hname = TString::Format("%s/%s/fHistRejectedJetEta", fName.Data(), fJetDefinitions->at(ij).GetName());
        fHistManager->FillTH1(hname, jet->Eta());
        continue;
      }
      fCurrentJetInfo[ij]->Set(dmeson_pair.second, fJetDefinitions->at(ij).GetName());
      accJets++;
    }
    if (accJets > 0) {
      if (histPrompt) {
        if (dmeson_pair.second.fParton) {
          partons[dmeson_pair.second.fParton] = dmeson_pair.second.fPartonType;
          UInt_t absPdgParton = TMath::Abs(dmeson_pair.second.fParton->GetPdgCode());
          if (absPdgParton == 4) {
            histPrompt->Fill("Prompt", 1);
          }
          else if (absPdgParton == 5) {
            histPrompt->Fill("Non-Prompt", 1);
          }
          else {
            histPrompt->Fill("Unknown", 1);
          }
        }
        else {
          histPrompt->Fill("Unknown", 1);
        }
      }

      if (histAncestor) {
        if (dmeson_pair.second.fAncestor) {
          UInt_t absPdgAncestor = TMath::Abs(dmeson_pair.second.fAncestor->GetPdgCode());
          if (absPdgAncestor == 4) {
            histAncestor->Fill("Charm", 1);
          }
          else if (absPdgAncestor == 5) {
            histAncestor->Fill("Bottom", 1);
          }
          else if (absPdgAncestor == 2212) {
            histAncestor->Fill("Proton", 1);
          }
          else {
            histAncestor->Fill("Unknown", 1);
          }
        }
        else {
          histAncestor->Fill("Unknown", 1);
        }
      }

      fTree->Fill();
    }
    else {
      hname = TString::Format("%s/fHistRejectedDMesonPt", fName.Data());
      fHistManager->FillTH1(hname, dmeson_pair.second.fD.Pt());
      hname = TString::Format("%s/fHistRejectedDMesonPhi", fName.Data());
      fHistManager->FillTH1(hname, dmeson_pair.second.fD.Phi_0_2pi());
      hname = TString::Format("%s/fHistRejectedDMesonEta", fName.Data());
      fHistManager->FillTH1(hname, dmeson_pair.second.fD.Eta());
      if (fMCMode != kMCTruth) {
        if (fCandidateType == kD0toKpi || fCandidateType == kD0toKpiLikeSign) {
          hname = TString::Format("%s/fHistRejectedDMesonInvMass", fName.Data());
          fHistManager->FillTH1(hname, dmeson_pair.second.fD.M());
        }
        else if (fCandidateType == kDstartoKpipi) {
          hname = TString::Format("%s/fHistRejectedDMeson2ProngInvMass", fName.Data());
          fHistManager->FillTH1(hname, dmeson_pair.second.fInvMass2Prong);

          hname = TString::Format("%s/fHistRejectedDMesonDeltaInvMass", fName.Data());
          fHistManager->FillTH1(hname, dmeson_pair.second.fD.M() - dmeson_pair.second.fInvMass2Prong);
        }
      }
    }
  }

  if (fMCMode == kSignalOnly || fMCMode == kMCTruth) {
    hname = TString::Format("%s/fHistPartonPt", fName.Data());
    TH1* histPartonPt = static_cast<TH1*>(fHistManager->FindObject(hname));
    hname = TString::Format("%s/fHistPartonEta", fName.Data());
    TH1* histPartonEta = static_cast<TH1*>(fHistManager->FindObject(hname));
    hname = TString::Format("%s/fHistPartonPhi", fName.Data());
    TH1* histPartonPhi = static_cast<TH1*>(fHistManager->FindObject(hname));
    hname = TString::Format("%s/fHistPartonType", fName.Data());
    TH1* histPartonType = static_cast<TH1*>(fHistManager->FindObject(hname));

    for (auto parton : partons) {
      if (!parton.first) continue;
      histPartonPt->Fill(parton.first->Pt());
      histPartonEta->Fill(parton.first->Eta());
      histPartonPhi->Fill(TVector2::Phi_0_2pi(parton.first->Phi()));
      histPartonType->Fill(parton.second);
    }
  }

  return kTRUE;
}

// Definitions of class AliAnalysisTaskDmesonJetsSub::OutputHandlerTTreeExtended

/// Constructor
AliAnalysisTaskDmesonJetsSub::OutputHandlerTTreeExtendedBase::OutputHandlerTTreeExtendedBase() :
  OutputHandler(),
  fDataSlotNumber(-1),
  fTree(0),
  fEventClassName(),
  fDMesonClassName(),
  fJetClassName()
{
}

/// Constructor
AliAnalysisTaskDmesonJetsSub::OutputHandlerTTreeExtendedBase::OutputHandlerTTreeExtendedBase(AnalysisEngine* eng) :
  OutputHandler(eng),
  fDataSlotNumber(-1),
  fTree(0),
  fEventClassName(),
  fDMesonClassName(),
  fJetClassName()
{
}

/// Builds the tree where the output will be posted
///
/// \return Pointer to the new tree
AliAnalysisTaskDmesonJetsSub::OutputHandlerTTreeExtendedBase* AliAnalysisTaskDmesonJetsSub::OutputHandlerTTreeExtendedBase::GenerateOutputHandler(AnalysisEngine* eng)
{
  TString event_class_name = "AliAnalysisTaskDmesonJetsSub::AliEventInfoSummary";
  TString d_meson_class_name;

  TString jet_class_name = "std::vector<AliAnalysisTaskDmesonJetsSub::AliJetInfoSummary>";
  Bool_t RhoJet = kFALSE;
  for (auto jetDef : eng->GetJetDefinitions()) {
    if (!jetDef.fRhoName.IsNull()) {
      RhoJet = kTRUE;
      jet_class_name = "std::vector<AliAnalysisTaskDmesonJetsSub::AliJetInfoPbPbSummary>";
    }
  }

  OutputHandlerTTreeExtendedBase* result = nullptr;
  if (eng->GetMCMode() == kMCTruth) {
    d_meson_class_name = "std::vector<AliAnalysisTaskDmesonJetsSub::AliDmesonMCInfoSummary>";
    if (RhoJet) {
      result = new AliAnalysisTaskDmesonJetsSub::OutputHandlerTTreeExtended<AliEventInfoSummary, AliDmesonMCInfoSummary, AliJetInfoPbPbSummary>(eng);
    }
    else {
      result = new AliAnalysisTaskDmesonJetsSub::OutputHandlerTTreeExtended<AliEventInfoSummary, AliDmesonMCInfoSummary, AliJetInfoSummary>(eng);
    }
  }
  else {
    switch (eng->GetCandidateType()) {
    case kD0toKpi:
    case kD0toKpiLikeSign:
      if (eng->IsD0Extended()) {
        d_meson_class_name = "std::vector<AliAnalysisTaskDmesonJetsSub::AliD0ExtendedInfoSummary>";
        if (RhoJet) {
          result = new AliAnalysisTaskDmesonJetsSub::OutputHandlerTTreeExtended<AliEventInfoSummary, AliD0ExtendedInfoSummary, AliJetInfoPbPbSummary>(eng);
        }
        else {
          result = new AliAnalysisTaskDmesonJetsSub::OutputHandlerTTreeExtended<AliEventInfoSummary, AliD0ExtendedInfoSummary, AliJetInfoSummary>(eng);
        }
      }
      else {
        d_meson_class_name = "AliAnalysisTaskDmesonJetsSub::AliD0InfoSummary";
        if (RhoJet) {
          result = new AliAnalysisTaskDmesonJetsSub::OutputHandlerTTreeExtended<AliEventInfoSummary, AliD0InfoSummary, AliJetInfoPbPbSummary>(eng);
        }
        else {
          result = new AliAnalysisTaskDmesonJetsSub::OutputHandlerTTreeExtended<AliEventInfoSummary, AliD0InfoSummary, AliJetInfoSummary>(eng);
        }
      }
      break;
    case kDstartoKpipi:
      d_meson_class_name = "AliAnalysisTaskDmesonJetsSub::AliDStarInfoSummary";
      if (RhoJet) {
        result = new AliAnalysisTaskDmesonJetsSub::OutputHandlerTTreeExtended<AliEventInfoSummary, AliDStarInfoSummary, AliJetInfoPbPbSummary>(eng);
      }
      else {
        result = new AliAnalysisTaskDmesonJetsSub::OutputHandlerTTreeExtended<AliEventInfoSummary, AliDStarInfoSummary, AliJetInfoSummary>(eng);
      }
      break;
    }
  }

  result->fEventClassName = event_class_name;
  result->fDMesonClassName = d_meson_class_name;
  result->fJetClassName = jet_class_name;

  return result;
}

// Definitions of class AliAnalysisTaskDmesonJetsSub::OutputHandlerTTree

/// Constructor
template<class EVENTTYPE, class DMESONTYPE, class JETTYPE>
AliAnalysisTaskDmesonJetsSub::OutputHandlerTTreeExtended<EVENTTYPE, DMESONTYPE, JETTYPE>::OutputHandlerTTreeExtended() :
  OutputHandlerTTreeExtendedBase(),
  fCurrentEventInfo(),
  fCurrentDmesonInfo(),
  fCurrentJetInfo()
{
}

/// Constructor
template<class EVENTTYPE, class DMESONTYPE, class JETTYPE>
AliAnalysisTaskDmesonJetsSub::OutputHandlerTTreeExtended<EVENTTYPE, DMESONTYPE, JETTYPE>::OutputHandlerTTreeExtended(AnalysisEngine* eng) :
  OutputHandlerTTreeExtendedBase(eng),
  fCurrentEventInfo(),
  fCurrentDmesonInfo(),
  fCurrentJetInfo()
{
}

/// Builds the tree where the output will be posted
///
/// \return Pointer to the new tree
template<class EVENTTYPE, class DMESONTYPE, class JETTYPE>
void AliAnalysisTaskDmesonJetsSub::OutputHandlerTTreeExtended<EVENTTYPE, DMESONTYPE, JETTYPE>::BuildOutputObject(const char* taskName)
{
  TString treeName = TString::Format("%s_%s", taskName, fName.Data());
  fTree = new TTree(treeName, treeName);
  fTree->Branch("Event", fEventClassName, &fCurrentEventInfo);
  fTree->Branch("Dmesons", fDMesonClassName, &fCurrentDmesonInfo);

  for (auto jetDef : *fJetDefinitions) {
    fCurrentJetInfo[jetDef.GetName()] = std::vector<JETTYPE>();
    fTree->Branch(jetDef.GetName(), fJetClassName, &fCurrentJetInfo[jetDef.GetName()]);
  }
}

/// Post the output with D meson jets found in the current event
///
/// \return kTRUE on success
template<class EVENTTYPE, class DMESONTYPE, class JETTYPE>
Bool_t AliAnalysisTaskDmesonJetsSub::OutputHandlerTTreeExtended<EVENTTYPE, DMESONTYPE, JETTYPE>::FillOutput(Bool_t applyKinCuts)
{
  if (fDmesonJets->empty()) return kFALSE;

  TString hname;

  std::map<AliAODMCParticle*, Short_t> partons ; //!<! set of the partons in the shower that produced each D meson

  TH1* histAncestor = nullptr;
  TH1* histPrompt = nullptr;

  if (fMCMode == kSignalOnly || fMCMode == kMCTruth) {
    hname = TString::Format("%s/fHistPrompt", fName.Data());
    histPrompt = static_cast<TH1*>(fHistManager->FindObject(hname));

    hname = TString::Format("%s/fHistAncestor", fName.Data());
    histAncestor = static_cast<TH1*>(fHistManager->FindObject(hname));
  }

  fCurrentEventInfo.Set(*fEventInfo);

  fCurrentDmesonInfo.clear();
  for (auto& jetInfo : fCurrentJetInfo) {
    jetInfo.second.clear();
  }

  for (auto& dmeson_pair : *fDmesonJets) {
    DMESONTYPE dmeson_tree;
    dmeson_tree.Set(dmeson_pair.second);
    Int_t accJets = 0;
    for (auto jetDef : *fJetDefinitions) {
      JETTYPE jet_tree;
      AliJetInfo* jet = dmeson_pair.second.GetJet(jetDef.GetName());
      if (jet) {
        if (applyKinCuts && !jetDef.IsJetInAcceptance(*jet)) {
          hname = TString::Format("%s/%s/fHistRejectedJetPt", fName.Data(), jetDef.GetName());
          fHistManager->FillTH1(hname, jet->Pt());
          hname = TString::Format("%s/%s/fHistRejectedJetPhi", fName.Data(), jetDef.GetName());
          fHistManager->FillTH1(hname, jet->Phi_0_2pi());
          hname = TString::Format("%s/%s/fHistRejectedJetEta", fName.Data(), jetDef.GetName());
          fHistManager->FillTH1(hname, jet->Eta());
        }
        else {
          jet_tree.Set(dmeson_pair.second, jetDef.GetName());
          accJets++;
        }
      }
      fCurrentJetInfo[jetDef.GetName()].push_back(jet_tree);
    }
    if (accJets > 0) {
      if (histPrompt) {
        if (dmeson_pair.second.fParton) {
          partons[dmeson_pair.second.fParton] = dmeson_pair.second.fPartonType;
          UInt_t absPdgParton = TMath::Abs(dmeson_pair.second.fParton->GetPdgCode());
          if (absPdgParton == 4) {
            histPrompt->Fill("Prompt", 1);
          }
          else if (absPdgParton == 5) {
            histPrompt->Fill("Non-Prompt", 1);
          }
          else {
            histPrompt->Fill("Unknown", 1);
          }
        }
        else {
          histPrompt->Fill("Unknown", 1);
        }
      }

      if (histAncestor) {
        if (dmeson_pair.second.fAncestor) {
          UInt_t absPdgAncestor = TMath::Abs(dmeson_pair.second.fAncestor->GetPdgCode());
          if (absPdgAncestor == 4) {
            histAncestor->Fill("Charm", 1);
          }
          else if (absPdgAncestor == 5) {
            histAncestor->Fill("Bottom", 1);
          }
          else if (absPdgAncestor == 2212) {
            histAncestor->Fill("Proton", 1);
          }
          else {
            histAncestor->Fill("Unknown", 1);
          }
        }
        else {
          histAncestor->Fill("Unknown", 1);
        }
      }

      fCurrentDmesonInfo.push_back(dmeson_tree);
    }
    else {
      hname = TString::Format("%s/fHistRejectedDMesonPt", fName.Data());
      fHistManager->FillTH1(hname, dmeson_pair.second.fD.Pt());
      hname = TString::Format("%s/fHistRejectedDMesonPhi", fName.Data());
      fHistManager->FillTH1(hname, dmeson_pair.second.fD.Phi_0_2pi());
      hname = TString::Format("%s/fHistRejectedDMesonEta", fName.Data());
      fHistManager->FillTH1(hname, dmeson_pair.second.fD.Eta());
      if (fMCMode != kMCTruth) {
        if (fCandidateType == kD0toKpi || fCandidateType == kD0toKpiLikeSign) {
          hname = TString::Format("%s/fHistRejectedDMesonInvMass", fName.Data());
          fHistManager->FillTH1(hname, dmeson_pair.second.fD.M());
        }
        else if (fCandidateType == kDstartoKpipi) {
          hname = TString::Format("%s/fHistRejectedDMeson2ProngInvMass", fName.Data());
          fHistManager->FillTH1(hname, dmeson_pair.second.fInvMass2Prong);

          hname = TString::Format("%s/fHistRejectedDMesonDeltaInvMass", fName.Data());
          fHistManager->FillTH1(hname, dmeson_pair.second.fD.M() - dmeson_pair.second.fInvMass2Prong);
        }
      }
      for (auto& jetInfo : fCurrentJetInfo) {
        jetInfo.second.pop_back();
      }
    }
  }

  if (!fCurrentDmesonInfo.empty()) fTree->Fill();

  if (fMCMode == kSignalOnly || fMCMode == kMCTruth) {
    hname = TString::Format("%s/fHistPartonPt", fName.Data());
    TH1* histPartonPt = static_cast<TH1*>(fHistManager->FindObject(hname));
    hname = TString::Format("%s/fHistPartonEta", fName.Data());
    TH1* histPartonEta = static_cast<TH1*>(fHistManager->FindObject(hname));
    hname = TString::Format("%s/fHistPartonPhi", fName.Data());
    TH1* histPartonPhi = static_cast<TH1*>(fHistManager->FindObject(hname));
    hname = TString::Format("%s/fHistPartonType", fName.Data());
    TH1* histPartonType = static_cast<TH1*>(fHistManager->FindObject(hname));

    for (auto parton : partons) {
      if (!parton.first) continue;
      histPartonPt->Fill(parton.first->Pt());
      histPartonEta->Fill(parton.first->Eta());
      histPartonPhi->Fill(TVector2::Phi_0_2pi(parton.first->Phi()));
      histPartonType->Fill(parton.second);
    }
  }

  return kTRUE;
}

// Definitions of class AliAnalysisTaskDmesonJetsSub::AliJetDefinition

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJetsSub::AliHFJetDefinition);
/// \endcond

/// This is the default constructor, used for ROOT I/O purposes.
AliAnalysisTaskDmesonJetsSub::AliHFJetDefinition::AliHFJetDefinition() :
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
AliAnalysisTaskDmesonJetsSub::AliHFJetDefinition::AliHFJetDefinition(EJetType_t type, Double_t r, EJetAlgo_t algo, ERecoScheme_t reco) :
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
AliAnalysisTaskDmesonJetsSub::AliHFJetDefinition::AliHFJetDefinition(EJetType_t type, Double_t r, EJetAlgo_t algo, ERecoScheme_t reco, TString rhoName) :
  TObject(),
  fJetType(type),
  fRadius(r),
  fJetAlgo(algo),
  fRecoScheme(reco),
  fMinJetPt(0.),
  fMaxJetPt(0.),
  fMinJetPhi(0.),
  fMaxJetPhi(0.),
  fMinJetEta(0.),
  fMaxJetEta(0.),
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
AliAnalysisTaskDmesonJetsSub::AliHFJetDefinition::AliHFJetDefinition(const AliHFJetDefinition &source) :
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
AliAnalysisTaskDmesonJetsSub::AliHFJetDefinition& AliAnalysisTaskDmesonJetsSub::AliHFJetDefinition::operator=(const AliHFJetDefinition& source)
{
  new (this) AliHFJetDefinition(source);
  return *this;
}

/// Generate a name for this jet definition
const char* AliAnalysisTaskDmesonJetsSub::AliHFJetDefinition::GetName() const
{
  static TString name;

  name = AliJetContainer::GenerateJetName(fJetType, fJetAlgo, fRecoScheme, fRadius, 0, 0, "Jet");

  return name.Data();
}

/// Decides whether the jet passes the acceptance cut defined in the object
///
/// \param jet Const reference to a AliJetInfo object
/// \return kTRUE if the jet passes the cuts
Bool_t AliAnalysisTaskDmesonJetsSub::AliHFJetDefinition::IsJetInAcceptance(const AliJetInfo& jet) const
{
  if (fMinJetEta < fMaxJetEta && (jet.Eta() < fMinJetEta || jet.Eta() > fMaxJetEta)) return kFALSE;
  if (fMinJetPhi < fMaxJetPhi && (jet.Phi() < fMinJetPhi || jet.Phi() > fMaxJetPhi)) return kFALSE;
  if (fMinJetPt < fMaxJetPt && (jet.Pt() > fMaxJetPt || jet.Pt() < fMinJetPt)) return kFALSE;
  if (fMinChargedPt < fMaxChargedPt && (jet.fMaxChargedPt < fMinChargedPt || jet.fMaxChargedPt > fMaxChargedPt)) return kFALSE;
  if (fMinNeutralPt < fMaxNeutralPt && (jet.fMaxNeutralPt < fMinNeutralPt || jet.fMaxNeutralPt > fMaxNeutralPt)) return kFALSE;

  return kTRUE;
}

/// Decides whether the jet passes the acceptance cut defined in the object
///
/// \param jet Const reference to a AliJetInfo object
/// \return kTRUE if the jet passes the cuts
Bool_t AliAnalysisTaskDmesonJetsSub::AliHFJetDefinition::IsJetInAcceptance(const AliDmesonJetInfo& dMesonJet, std::string n) const
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
bool operator<(const AliAnalysisTaskDmesonJetsSub::AliHFJetDefinition& lhs, const AliAnalysisTaskDmesonJetsSub::AliHFJetDefinition& rhs)
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
bool operator==(const AliAnalysisTaskDmesonJetsSub::AliHFJetDefinition& lhs, const AliAnalysisTaskDmesonJetsSub::AliHFJetDefinition& rhs)
{
  if (lhs.fJetType != rhs.fJetType) return false;
  if (lhs.fRadius != rhs.fRadius) return false;
  if (lhs.fJetAlgo != rhs.fJetAlgo) return false;
  if (lhs.fRecoScheme != rhs.fRecoScheme) return false;
  return true;
}

// Definitions of class AliAnalysisTaskDmesonJetsSub::AnalysisEngine

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJetsSub::AnalysisEngine);
/// \endcond

/// This is the default constructor, used for ROOT I/O purposes.
AliAnalysisTaskDmesonJetsSub::AnalysisEngine::AnalysisEngine() :
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
  fD0Extended(kFALSE),
  fOutputHandler(nullptr),
  fRandomGen(0),
  fTrackEfficiency(0),
  fRejectISR(kFALSE),
  fDmesonJets(),
  fCandidateArray(0),
  fMCContainer(),
  fTrackContainers(),
  fClusterContainers(),
  fAodEvent(0),
  fFastJetWrapper(0),
  fHistManager(0),
  fEventInfo(),
  fName()
{
}

/// This is the standard constructor.
///
/// \param type      One of the enum constants of ECandidateType_t
/// \param bkgMode   One of the enum constants of EMCMode_t
/// \param cuts      D meson cuts (if null, it will use standard cuts)
/// \param nMassBins Number of bins in the mass axis
/// \param range     Range of the mass axis (will be centered around the PDG mass)
AliAnalysisTaskDmesonJetsSub::AnalysisEngine::AnalysisEngine(ECandidateType_t type, EMCMode_t MCmode, AliRDHFCuts* cuts, Int_t nMassBins, Double_t range) :
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
  fRejectedOrigin(0),
  fAcceptedDecay(kAnyDecay),
  fInhibit(kFALSE),
  fJetDefinitions(),
  fPtBinWidth(0.5),
  fMaxPt(100),
  fD0Extended(kFALSE),
  fOutputHandler(nullptr),
  fRandomGen(0),
  fTrackEfficiency(0),
  fDmesonJets(),
  fCandidateArray(0),
  fMCContainer(),
  fTrackContainers(),
  fClusterContainers(),
  fAodEvent(0),
  fFastJetWrapper(0),
  fHistManager(0),
  fEventInfo(),
  fName()
{
  SetCandidateProperties(range);
}

/// Copy constructor
///
/// \param source Reference to a valid AnalysisEngine to copy from.
AliAnalysisTaskDmesonJetsSub::AnalysisEngine::AnalysisEngine(const AliAnalysisTaskDmesonJetsSub::AnalysisEngine &source) :
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
  fD0Extended(source.fD0Extended),
  fRandomGen(source.fRandomGen),
  fTrackEfficiency(source.fTrackEfficiency),
  fDmesonJets(),
  fCandidateArray(source.fCandidateArray),
  fMCContainer(source.fMCContainer),
  fTrackContainers(source.fTrackContainers),
  fClusterContainers(source.fClusterContainers),
  fAodEvent(source.fAodEvent),
  fFastJetWrapper(source.fFastJetWrapper),
  fHistManager(source.fHistManager),
  fEventInfo(),
  fName()
{
  SetRDHFCuts(source.fRDHFCuts);
}

// Destructor
AliAnalysisTaskDmesonJetsSub::AnalysisEngine::~AnalysisEngine()
{
  delete fRDHFCuts;
}

/// Assignement operator
///
/// \param source Reference to a valid AnalysisEngine to copy from.
AliAnalysisTaskDmesonJetsSub::AnalysisEngine& AliAnalysisTaskDmesonJetsSub::AnalysisEngine::operator=(const AnalysisEngine& source)
{
  new (this) AnalysisEngine(source);
  return *this;
}

/// Initialize the analysis engine
void AliAnalysisTaskDmesonJetsSub::AnalysisEngine::Init(const AliEMCALGeometry* const /*geom*/, Int_t /*runNumber*/)
{
}

/// Sets the D meson candidate properties.
///
/// \param range     Range of the mass axis (will be centered around the PDG mass)
void AliAnalysisTaskDmesonJetsSub::AnalysisEngine::SetCandidateProperties(Double_t range)
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
    ::Error("AliAnalysisTaskDmesonJetsSub::AnalysisEngine::SetCandidateProperties","Candidate %d unknown!", fCandidateType);
  }

  CalculateMassLimits(range, fCandidatePDG, fNMassBins, fMinMass, fMaxMass);
}

/// Adopt the cuts (this class owns the cuts object, which will be destroyed when needed).
///
/// \param Pointer to a AliRDHFCuts object.
void AliAnalysisTaskDmesonJetsSub::AnalysisEngine::AdoptRDHFCuts(AliRDHFCuts* cuts)
{
  if (fRDHFCuts) delete fRDHFCuts;
  fRDHFCuts = cuts;
}

/// Set the cuts (creates a copy, so the original object is not owned by this class).
///
/// \param Pointer to a AliRDHFCuts object.
void AliAnalysisTaskDmesonJetsSub::AnalysisEngine::SetRDHFCuts(AliRDHFCuts* cuts)
{
  if (!cuts) return;
  if (fRDHFCuts) delete fRDHFCuts;
  fRDHFCuts = static_cast<AliRDHFCuts*>(cuts->Clone());
}

/// Generate a name for this analysis parameter set
///
/// \param i  Index of the jet radius array.
const char* AliAnalysisTaskDmesonJetsSub::AnalysisEngine::GetName(const AliHFJetDefinition& jetDef) const
{
  static TString name;

  name = TString::Format("%s_%s", GetName(), jetDef.GetName());

  return name.Data();
}

/// Generate a name for this analysis parameter set
///
/// \param i  Index of the jet radius array.
const char* AliAnalysisTaskDmesonJetsSub::AnalysisEngine::GetName() const
{
  fName = fCandidateName;
  switch (fMCMode) {
  case kBackgroundOnly:
    fName += "_BackgroundOnly";
    break;
  case kSignalOnly:
    fName += "_SignalOnly";
    break;
  case kMCTruth:
    fName += "_MCTruth";
    break;
  case kD0Reflection:
    fName += "_D0Reflection";
    break;
  case kOnlyWrongPIDAccepted:
    fName += "_OnlyWrongPIDAccepted";
    break;
  default:
    break;
  }

  if (fRDHFCuts) fName += TString::Format("_%s", fRDHFCuts->GetName());

  return fName.Data();
}

/// Add a new jet definition
/// If the jet definition is already present, it does nothing.
///
/// \param def Reference to a AliJetDefinition object
///
/// \return Pointer to the new jet definition (or to the one that was already present)
AliAnalysisTaskDmesonJetsSub::AliHFJetDefinition* AliAnalysisTaskDmesonJetsSub::AnalysisEngine::AddJetDefinition(const AliAnalysisTaskDmesonJetsSub::AliHFJetDefinition& def)
{
  std::vector<AliHFJetDefinition>::iterator it = FindJetDefinition(def);

  if (it == fJetDefinitions.end() || *it != def) {  // No jet definition was found, adding a new one
    fJetDefinitions.push_back(def);
    ::Info("AliAnalysisTaskDmesonJetsSub::AnalysisEngine::AddJetDefinition", "Jet definition '%s' has been added to analysis engine '%s'."
        "Total number of jet definitions is now %lu.",
        def.GetName(), GetName(), fJetDefinitions.size());
    // For detector level set maximum track pt to 100 GeV/c
    if (fMCMode != kMCTruth) fJetDefinitions[fJetDefinitions.size()-1].SetChargedPtRange(0., 100.);
  }
  else {
    ::Warning("AliAnalysisTaskDmesonJetsSub::AnalysisEngine::AddJetDefinition", "The same jet definition '%s' was already added in analysis engine '%s'.", def.GetName(), GetName());
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
AliAnalysisTaskDmesonJetsSub::AliHFJetDefinition*
AliAnalysisTaskDmesonJetsSub::AnalysisEngine::AddJetDefinition(EJetType_t type, Double_t r, EJetAlgo_t algo, ERecoScheme_t reco)
{
  AliHFJetDefinition def(type, r, algo, reco);

  return AddJetDefinition(def);
}

/// Look for a jet definition that is equal
///
/// \param def Reference to a jet definition object
///
/// \return An iterator to the jet definition object, if it is found. An iterator to the end if not found.
std::vector<AliAnalysisTaskDmesonJetsSub::AliHFJetDefinition>::iterator AliAnalysisTaskDmesonJetsSub::AnalysisEngine::FindJetDefinition(const AliAnalysisTaskDmesonJetsSub::AliHFJetDefinition& def)
{
  std::vector<AliHFJetDefinition>::iterator it = fJetDefinitions.begin();
  while (it != fJetDefinitions.end() && (*it) != def) it++;
  return it;
}

/// Set the jet phi range of all jet definitions
/// \param min Lower bound
/// \param max Upper bound
void AliAnalysisTaskDmesonJetsSub::AnalysisEngine::SetJetPhiRange(Double_t min, Double_t max)
{
  for (auto &jetdef : fJetDefinitions) jetdef.SetJetPhiRange(min, max);
}

/// Set the jet eta range of all jet definitions
/// \param min Lower bound
/// \param max Upper bound
void AliAnalysisTaskDmesonJetsSub::AnalysisEngine::SetJetEtaRange(Double_t min, Double_t max)
{
  for (auto &jetdef : fJetDefinitions) jetdef.SetJetEtaRange(min, max);
}

/// Set the jet pt range of all jet definitions
/// \param min Lower bound
/// \param max Upper bound
void AliAnalysisTaskDmesonJetsSub::AnalysisEngine::SetJetPtRange(Double_t min, Double_t max)
{
  for (auto &jetdef : fJetDefinitions) jetdef.SetJetPtRange(min, max);
}

/// Set the jet leading charged constituent pt range of all jet definitions
/// \param min Lower bound
/// \param max Upper bound
void AliAnalysisTaskDmesonJetsSub::AnalysisEngine::SetChargedPtRange(Double_t min, Double_t max)
{
  for (auto &jetdef : fJetDefinitions) jetdef.SetChargedPtRange(min, max);
}

/// Set the jet leading neutral constituent pt range range of all jet definitions
/// \param min Lower bound
/// \param max Upper bound
void AliAnalysisTaskDmesonJetsSub::AnalysisEngine::SetNeutralPtRange(Double_t min, Double_t max)
{
  for (auto &jetdef : fJetDefinitions) jetdef.SetNeutralPtRange(min, max);
}

/// Compares 2 analysis engines.
/// The ordering is based on the candidate type first and then on the MC mode.
///
/// \param lhs Reference to the first AnalysisEngine object
/// \param rhs Reference to the second AnalysisEngine object
bool operator<(const AliAnalysisTaskDmesonJetsSub::AnalysisEngine& lhs, const AliAnalysisTaskDmesonJetsSub::AnalysisEngine& rhs)
{
  if (lhs.fCandidateType < rhs.fCandidateType) {
    return true;
  }
  else if (lhs.fCandidateType > rhs.fCandidateType) {
    return false;
  }
  else if (lhs.fMCMode < rhs.fMCMode) {
    return true;
  }
  else if (lhs.fMCMode > rhs.fMCMode) {
    return false;
  }
  else if (lhs.fRDHFCuts && !rhs.fRDHFCuts) {
    return true;
  }
  else if (lhs.fRDHFCuts && rhs.fRDHFCuts && strcmp(lhs.fRDHFCuts->GetName(), rhs.fRDHFCuts->GetName()) < 0) {
    return true;
  }
  else {
    return false;
  }
}

/// Compares 2 analysis engines.
/// Two analysis engines are considerate equal if they have both the same candidate type and MC mode.
///
/// \param lhs Reference to the first AnalysisEngine object
/// \param rhs Reference to the second AnalysisEngine object
bool operator==(const AliAnalysisTaskDmesonJetsSub::AnalysisEngine& lhs, const AliAnalysisTaskDmesonJetsSub::AnalysisEngine& rhs)
{
  if (lhs.fCandidateType != rhs.fCandidateType) return false;
  if (lhs.fMCMode != rhs.fMCMode) return false;
  if (lhs.fRDHFCuts == nullptr && rhs.fRDHFCuts != nullptr) return false;
  if (lhs.fRDHFCuts != nullptr && rhs.fRDHFCuts == nullptr) return false;
  if (lhs.fRDHFCuts && rhs.fRDHFCuts && strcmp(lhs.fRDHFCuts->GetName(), rhs.fRDHFCuts->GetName()) != 0) return false;
  return true;
}

/// Extract attributes of the D meson candidate.
///
/// \param Dcand Pointer to a AliAODRecoDecayHF2Prong representing the D meson candidate
/// \param DmesonJet Reference to an AliDmesonJetInfo object where the D meson candidate information will be copied
/// \param i Either 0 or 1, for the two possible mass hypothesis assignments
///
/// \return kTRUE on success
Bool_t AliAnalysisTaskDmesonJetsSub::AnalysisEngine::ExtractRecoDecayAttributes(const AliAODRecoDecayHF2Prong* Dcand, AliDmesonJetInfo& DmesonJet, UInt_t i)
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

Bool_t AliAnalysisTaskDmesonJetsSub::AnalysisEngine::ExtractEfficiencies(const AliAODRecoDecayHF2Prong* Dcand, AliDmesonJetInfo& DmesonJet, AliHFJetDefinition& jetDef,UInt_t i)
{
  if (fCandidateType == kD0toKpi || fCandidateType == kD0toKpiLikeSign) { // D0 candidate
   
    return ExtractD0Efficiencies(Dcand, DmesonJet, jetDef, i);
  }
    else {
    return kFALSE;
  }
}

Bool_t AliAnalysisTaskDmesonJetsSub::AnalysisEngine::ExtractD0Efficiencies(const AliAODRecoDecayHF2Prong* Dcand, AliDmesonJetInfo& DmesonJet,AliHFJetDefinition& jetDef, UInt_t i)
{
  AliDebug(10,"Checking if D0 meson is selected");
  Int_t isSelected = fRDHFCuts->IsSelected(const_cast<AliAODRecoDecayHF2Prong*>(Dcand), AliRDHFCuts::kAll, fAodEvent);
  if (isSelected == 0) return kFALSE;
  TString hname1;
  TString hname2;
 
  Int_t MCtruthPdgCode = 0;
 
  Int_t myflag=0;
  Double_t jeteta=0;
  Double_t jetpt=0;
  hname1 = TString::Format("%s/EfficiencyMatchesPrompt", fName.Data());
  TH2* EfficiencyMatchesPrompt = static_cast<TH2*>(fHistManager->FindObject(hname1));
  AliAODMCParticle* aodMcPart; 
  
  
  hname2 = TString::Format("%s/EfficiencyMatchesNonPrompt", fName.Data());
  TH2* EfficiencyMatchesNonPrompt = static_cast<TH2*>(fHistManager->FindObject(hname2));
  
  // If the analysis require knowledge of the MC truth, look for generated D meson matched to reconstructed candidate
  // Checks also the origin, and if it matches the rejected origin mask, return false
 Double_t jetPtdet = DmesonJet.fJets[jetDef.GetName()].fMomentum.Pt();
 Double_t jetEtadet=  DmesonJet.fJets[jetDef.GetName()].fMomentum.Eta();
 if(TMath::Abs(jetEtadet)>0.5) return kFALSE;
 if(jetPtdet>100) return kFALSE;

 
   if (fMCMode != kNoMC) {
    Int_t mcLab = Dcand->MatchToMC(fCandidatePDG, fMCContainer->GetArray(), fNDaughters, fPDGdaughters.GetArray());
    DmesonJet.fMCLabel = mcLab;

    // Retrieve the generated particle (if exists) and its PDG code
    if (mcLab >= 0) {
      aodMcPart = static_cast<AliAODMCParticle*>(fMCContainer->GetArray()->At(mcLab));

      if (aodMcPart) {
        // Check origin and return false if it matches the rejected origin mask
        if (fRejectedOrigin) {
          auto origin = IsPromptCharm(aodMcPart, fMCContainer->GetArray());
          if ((origin.first & fRejectedOrigin) == origin.first) return kFALSE;
        }
        MCtruthPdgCode = aodMcPart->PdgCode();
      }
    }
  }
  

 if(fMCMode==kSignalOnly){
     fMCContainer->SetCharge(AliParticleContainer::EChargeCut_t::kCharged);
    fFastJetWrapper->Clear();
    fFastJetWrapper->SetR(jetDef.fRadius);
    fFastJetWrapper->SetAlgorithm(AliEmcalJetTask::ConvertToFJAlgo(jetDef.fJetAlgo));
    fFastJetWrapper->SetRecombScheme(AliEmcalJetTask::ConvertToFJRecoScheme(jetDef.fRecoScheme));
    AddInputVectors(fMCContainer, 100);
    fFastJetWrapper->Run();
    std::vector<fastjet::PseudoJet> jets_incl = sorted_by_pt(fFastJetWrapper->GetInclusiveJets());

    for (auto jet : jets_incl) {
      std::vector<fastjet::PseudoJet> constituents = jet.constituents();
      // if(constituents.size()<2) continue;
          for (auto constituent : jet.constituents()) {
             Int_t iPart = constituent.user_index() - 100;
             if (constituent.perp() < 1e-6) continue; // reject ghost particles
        AliAODMCParticle* part = fMCContainer->GetMCParticle(iPart);
        if (!part) {
          ::Error("AliAnalysisTaskDmesonJetsSub::AnalysisEngine::RunParticleLevelAnalysis", "Could not find jet constituent %d!", iPart);
          continue;
        }
	if(part==aodMcPart){myflag=1;
	       jetpt=jet.perp();
	       jeteta=jet.eta();
	  break;}

	  }}
    if(myflag==0) return kFALSE;
    if(TMath::Abs(jeteta)>0.5) return kFALSE;
  
          Double_t maxFiducialY=-0.2/15*aodMcPart->Pt()*aodMcPart->Pt()+1.9/15*aodMcPart->Pt()+0.5;
	  Double_t minFiducialY = 0.2/15*aodMcPart->Pt()*aodMcPart->Pt()-1.9/15*aodMcPart->Pt()-0.5;	
   
    if (isSelected == 1) { // selected as a D0
    if (i != 0) return kFALSE; // only one mass hypothesis thanks to PID
   
      if(MCtruthPdgCode == fCandidatePDG){
         auto origin = IsPromptCharm(aodMcPart, fMCContainer->GetArray());

	 if(origin.first == kFromCharm) {
	  	if(aodMcPart->Pt()<=5) if(aodMcPart->Y()>=minFiducialY && aodMcPart->Y()<=maxFiducialY) EfficiencyMatchesPrompt->Fill(aodMcPart->Pt(),jetpt);
	if(aodMcPart->Pt()>5) if(TMath::Abs(aodMcPart->Y())<0.8)EfficiencyMatchesPrompt->Fill(aodMcPart->Pt(),jetpt);}
          if(origin.first == kFromBottom) {
	    	if(aodMcPart->Pt()<=5) if(aodMcPart->Y()>=minFiducialY && aodMcPart->Y()<=maxFiducialY)EfficiencyMatchesNonPrompt->Fill(aodMcPart->Pt(),jetpt);
	if(aodMcPart->Pt()>5) if(TMath::Abs(aodMcPart->Y())<0.8)EfficiencyMatchesNonPrompt->Fill(aodMcPart->Pt(),jetpt);}
      }

  }
  else if (isSelected == 2) { // selected as a D0bar
    if (i != 1) return kFALSE; // only one mass hypothesis thanks to PID

      
      if(MCtruthPdgCode == -fCandidatePDG){
      auto origin = IsPromptCharm(aodMcPart, fMCContainer->GetArray());
      if(origin.first == kFromCharm){
 
      	  if(aodMcPart->Pt()<=5) if(aodMcPart->Y()>=minFiducialY && aodMcPart->Y()<=maxFiducialY)EfficiencyMatchesPrompt->Fill(aodMcPart->Pt(),jetpt);
	  if(aodMcPart->Pt()>5) if(TMath::Abs(aodMcPart->Y())<0.8)EfficiencyMatchesPrompt->Fill(aodMcPart->Pt(),jetpt);}
        if(origin.first == kFromBottom){
	  if(aodMcPart->Pt()<=5) if(aodMcPart->Y()>=minFiducialY && aodMcPart->Y()<=maxFiducialY)EfficiencyMatchesNonPrompt->Fill(aodMcPart->Pt(),jetpt);
	  if(aodMcPart->Pt()>5) if(TMath::Abs(aodMcPart->Y())<0.8)EfficiencyMatchesNonPrompt->Fill(aodMcPart->Pt(),jetpt);
	  }
  }

  }

  else if (isSelected == 3) { // selected as either a D0bar or a D0 (PID on K and pi undecisive)
  

    // Accept the correct mass hypothesis for signal-only and the wrong one for background-only
    if (MCtruthPdgCode == fCandidatePDG){
      
      if (i == 0){
      auto origin = IsPromptCharm(aodMcPart, fMCContainer->GetArray());
      if(origin.first == kFromCharm){
	if(aodMcPart->Pt()<=5) if(aodMcPart->Y()>=minFiducialY && aodMcPart->Y()<=maxFiducialY)EfficiencyMatchesPrompt->Fill(aodMcPart->Pt(),jetpt);
	if(aodMcPart->Pt()>5) if(TMath::Abs(aodMcPart->Y())<0.8)EfficiencyMatchesPrompt->Fill(aodMcPart->Pt(),jetpt);
      }
        if(origin.first == kFromBottom){
	if(aodMcPart->Pt()<=5) if(aodMcPart->Y()>=minFiducialY && aodMcPart->Y()<=maxFiducialY)EfficiencyMatchesNonPrompt->Fill(aodMcPart->Pt(),jetpt);
	if(aodMcPart->Pt()>5) if(TMath::Abs(aodMcPart->Y())<0.8)EfficiencyMatchesNonPrompt->Fill(aodMcPart->Pt(),jetpt);


	}}
    }
    else if (MCtruthPdgCode == -fCandidatePDG){
      if (i == 1){

       auto origin = IsPromptCharm(aodMcPart, fMCContainer->GetArray());
      if(origin.first == kFromCharm){
		if(aodMcPart->Pt()<=5) if(aodMcPart->Y()>=minFiducialY && aodMcPart->Y()<=maxFiducialY)EfficiencyMatchesPrompt->Fill(aodMcPart->Pt(),jetpt);
	if(aodMcPart->Pt()>5) if(TMath::Abs(aodMcPart->Y())<0.8)EfficiencyMatchesPrompt->Fill(aodMcPart->Pt(),jetpt);}
        if(origin.first == kFromBottom){
	  	if(aodMcPart->Pt()<=5) if(aodMcPart->Y()>=minFiducialY && aodMcPart->Y()<=maxFiducialY)EfficiencyMatchesNonPrompt->Fill(aodMcPart->Pt(),jetpt);
	if(aodMcPart->Pt()>5) if(TMath::Abs(aodMcPart->Y())<0.8)EfficiencyMatchesNonPrompt->Fill(aodMcPart->Pt(),jetpt);}}
      
    }
  }










 }
  return kTRUE;
}


Bool_t AliAnalysisTaskDmesonJetsSub::AnalysisEngine::GetEfficiencyDenominator(AliHFJetDefinition& jetDef)
{
  
  TString hname;
  TString hname1;
  TString hname2;
  
  Double_t jeteta=0;
  Double_t jetpt=0;
  Int_t TheTrueCode = 0;

  fMCContainer->SetSpecialPDG(fCandidatePDG);
  fMCContainer->SetRejectedOriginMap(fRejectedOrigin);
  fMCContainer->SetAcceptedDecayMap(fAcceptedDecay);
  fMCContainer->SetRejectISR(fRejectISR);
  fMCContainer->SetSpecialPDG(fCandidatePDG);
  fMCContainer->SetCharge(AliParticleContainer::EChargeCut_t::kCharged);
  
     if (!fMCContainer->IsSpecialPDGFound()) return kFALSE;
  
  hname1 = TString::Format("%s/EfficiencyGeneratorPrompt", fName.Data());
  TH2* EfficiencyGeneratorPrompt = static_cast<TH2*>(fHistManager->FindObject(hname1));
 
  
  
  hname2 = TString::Format("%s/EfficiencyGeneratorNonPrompt", fName.Data());
  TH2* EfficiencyGeneratorNonPrompt = static_cast<TH2*>(fHistManager->FindObject(hname2));
  if(fMCMode==kSignalOnly){
   
    fFastJetWrapper->Clear();
    fFastJetWrapper->SetR(jetDef.fRadius);
    fFastJetWrapper->SetAlgorithm(AliEmcalJetTask::ConvertToFJAlgo(jetDef.fJetAlgo));
    fFastJetWrapper->SetRecombScheme(AliEmcalJetTask::ConvertToFJRecoScheme(jetDef.fRecoScheme));
     hname = TString::Format("%s/%s/fHistMCParticleRejectionReason", GetName(), jetDef.GetName());
     AddInputVectors(fMCContainer, 100, static_cast<TH2*>(fHistManager->FindObject(hname))); 
    
    fFastJetWrapper->Run();
    std::vector<fastjet::PseudoJet> jets_incl =  sorted_by_pt(fFastJetWrapper->GetInclusiveJets());
   
 
     
      for (auto jet : jets_incl) {
      jetpt=0;
      jeteta=0;
        std::vector<fastjet::PseudoJet> constituents = jet.constituents();
    
            for (auto constituent : jet.constituents()) {
             Int_t iPart = constituent.user_index() - 100;
             if (constituent.perp() < 1e-6) continue; // reject ghost particles
             AliAODMCParticle* part = fMCContainer->GetMCParticle(iPart);
             if (!part) {
	       ::Error("AliAnalysisTaskDmesonJetsSub::AnalysisEngine::RunParticleLevelAnalysis", "Could not find jet constituent %d!", iPart);
              continue;
                         }

	     TheTrueCode=part->PdgCode();
	     //  Int_t mother=part->GetMother();
	     // AliAODMCParticle* mypart = fMCContainer->GetMCParticle(mother);
	     if(TMath::Abs(TheTrueCode)==fCandidatePDG){
	  
      	  jetpt=jet.perp();
	  jeteta=jet.eta();
          if(TMath::Abs(jeteta)>0.5) continue;
	  
          auto origin = IsPromptCharm(part, fMCContainer->GetArray());
	  Double_t maxFiducialY=-0.2/15*part->Pt()*part->Pt()+1.9/15*part->Pt()+0.5;
	  Double_t minFiducialY = 0.2/15*part->Pt()*part->Pt()-1.9/15*part->Pt()-0.5;	
      if(origin.first == kFromCharm){
	if(part->Pt()>5) if(TMath::Abs(part->Y())<0.8) EfficiencyGeneratorPrompt->Fill(part->Pt(),jetpt);
        if(part->Pt()<=5) if(part->Y()>=minFiducialY && part->Y()<=maxFiducialY) EfficiencyGeneratorPrompt->Fill(part->Pt(),jetpt);
      }
        if(origin.first == kFromBottom){ 
	  if(part->Pt()>5)if(TMath::Abs(part->Y())<0.8) EfficiencyGeneratorNonPrompt->Fill(part->Pt(),jetpt);
	  if(part->Pt()<=5) if(part->Y()>=minFiducialY && part->Y()<=maxFiducialY) EfficiencyGeneratorNonPrompt->Fill(part->Pt(),jetpt);
	}
       
	     }

  
   
	    }}
     
  }
 
  return kTRUE;
}

Bool_t AliAnalysisTaskDmesonJetsSub::AnalysisEngine::GetEfficiencyDenominatorOneByOne(AliHFJetDefinition& jetDef)
{
  
  TString hname;
  TString hname1;
  TString hname2;
 
  Double_t jeteta=0;
  Double_t jetpt=0;
  Int_t TheTrueCode = 0;
 
  vector<int> dlabel;
  fMCContainer->SetSpecialPDG(fCandidatePDG);
  fMCContainer->SetRejectedOriginMap(fRejectedOrigin);
  fMCContainer->SetAcceptedDecayMap(fAcceptedDecay);
  fMCContainer->SetRejectISR(fRejectISR);
  fMCContainer->SetSpecialPDG(fCandidatePDG);
  fMCContainer->SetSpecialIndex(-10);
  fMCContainer->SetCharge(AliParticleContainer::EChargeCut_t::kCharged);
     if (!fMCContainer->IsSpecialPDGFound()) return kFALSE;
   
    
  hname1 = TString::Format("%s/EfficiencyGeneratorPrompt", fName.Data());
  TH2* EfficiencyGeneratorPrompt = static_cast<TH2*>(fHistManager->FindObject(hname1));
 
  
  
  hname2 = TString::Format("%s/EfficiencyGeneratorNonPrompt", fName.Data());
  TH2* EfficiencyGeneratorNonPrompt = static_cast<TH2*>(fHistManager->FindObject(hname2));
  if(fMCMode==kSignalOnly){
   
    //here I loop over the container and count the D mesons and store their indexes
    auto cont = fMCContainer->all();
    for (auto it = cont.begin(); it != cont.end(); ++it) {
    UInt_t rejectionReason = 0;
    if((*it)->PdgCode()==fCandidatePDG){
    
      dlabel.push_back(it.current_index());}
      
 
     if (!fMCContainer->AcceptObject(it.current_index(), rejectionReason)) {
      
      continue;
     }

     
     }
   
     // then, for each D meson I replace only its decays (not other D decays) and I  only keep for the jet finding the given D meson
    for(Int_t j=0;j<dlabel.size();j++){
     
          fMCContainer->SetSpecialPDG(-10);
          fMCContainer->SetSpecialIndex(dlabel[j]);
	 
    fFastJetWrapper->Clear();
    fFastJetWrapper->SetR(jetDef.fRadius);
    fFastJetWrapper->SetAlgorithm(AliEmcalJetTask::ConvertToFJAlgo(jetDef.fJetAlgo));
    fFastJetWrapper->SetRecombScheme(AliEmcalJetTask::ConvertToFJRecoScheme(jetDef.fRecoScheme));
     hname = TString::Format("%s/%s/fHistMCParticleRejectionReason", GetName(), jetDef.GetName());
     AddInputVectors(fMCContainer, 100, static_cast<TH2*>(fHistManager->FindObject(hname))); 
     fFastJetWrapper->Run();
    std::vector<fastjet::PseudoJet> jets_incl =  sorted_by_pt(fFastJetWrapper->GetInclusiveJets());
   
 
     
      for (auto jet : jets_incl) {
      jetpt=0;
      jeteta=0;
        std::vector<fastjet::PseudoJet> constituents = jet.constituents();
    
            for (auto constituent : jet.constituents()) {
             Int_t iPart = constituent.user_index() - 100;
             if (constituent.perp() < 1e-6) continue; // reject ghost particles
             AliAODMCParticle* part = fMCContainer->GetMCParticle(iPart);
             if (!part) {
	       ::Error("AliAnalysisTaskDmesonJetsSub::AnalysisEngine::RunParticleLevelAnalysis", "Could not find jet constituent %d!", iPart);
              continue;
                         }

	     TheTrueCode=part->PdgCode();
	    
	     if(TMath::Abs(TheTrueCode)==fCandidatePDG){
	  
      	  jetpt=jet.perp();
	  jeteta=jet.eta();
          if(TMath::Abs(jeteta)>0.5) continue;
	  
          auto origin = IsPromptCharm(part, fMCContainer->GetArray());
          Double_t maxFiducialY=-0.2/15*part->Pt()*part->Pt()+1.9/15*part->Pt()+0.5;
	  Double_t minFiducialY = 0.2/15*part->Pt()*part->Pt()-1.9/15*part->Pt()-0.5;	
      if(origin.first == kFromCharm){
         	if(part->Pt()>5) if(TMath::Abs(part->Y())<0.8) EfficiencyGeneratorPrompt->Fill(part->Pt(),jetpt);
		if(part->Pt()<=5) if(part->Y()>=minFiducialY && part->Y()<=maxFiducialY) EfficiencyGeneratorPrompt->Fill(part->Pt(),jetpt);}
        if(origin.first == kFromBottom){ 
	 if(part->Pt()>5) if(TMath::Abs(part->Y())<0.8) EfficiencyGeneratorNonPrompt->Fill(part->Pt(),jetpt);
		if(part->Pt()<=5) if(part->Y()>=minFiducialY && part->Y()<=maxFiducialY) EfficiencyGeneratorNonPrompt->Fill(part->Pt(),jetpt);

	    }
      
	     }

  
   
	    }}
     
    }
    dlabel.clear();
  }
 
  return kTRUE;
}




/// Extract attributes of the D0 meson candidate.
///
/// \param Dcand Pointer to a AliAODRecoDecayHF2Prong representing the D0 meson candidate
/// \param DmesonJet Reference to an AliDmesonJetInfo object where the D0 meson candidate information will be copied
/// \param i Either 0 or 1, for the two possible mass hypothesis assignments
///
/// \return kTRUE on success
Bool_t AliAnalysisTaskDmesonJetsSub::AnalysisEngine::ExtractD0Attributes(const AliAODRecoDecayHF2Prong* Dcand, AliDmesonJetInfo& DmesonJet, UInt_t i)
{
  AliDebug(10,"Checking if D0 meson is selected");
  Int_t isSelected = fRDHFCuts->IsSelected(const_cast<AliAODRecoDecayHF2Prong*>(Dcand), AliRDHFCuts::kAll, fAodEvent);
  if (isSelected == 0) return kFALSE;
  TString hname;
  TString hname2;
  Int_t MCtruthPdgCode = 0;
  
  Double_t invMassD = 0;

  AliAODMCParticle* aodMcPart; 

  
  // If the analysis require knowledge of the MC truth, look for generated D meson matched to reconstructed candidate
  // Checks also the origin, and if it matches the rejected origin mask, return false
  if (fMCMode != kNoMC) {
    Int_t mcLab = Dcand->MatchToMC(fCandidatePDG, fMCContainer->GetArray(), fNDaughters, fPDGdaughters.GetArray());
    DmesonJet.fMCLabel = mcLab;

    // Retrieve the generated particle (if exists) and its PDG code
    if (mcLab >= 0) {
      aodMcPart = static_cast<AliAODMCParticle*>(fMCContainer->GetArray()->At(mcLab));

      if (aodMcPart) {
        // Check origin and return false if it matches the rejected origin mask
        if (fRejectedOrigin) {
          auto origin = IsPromptCharm(aodMcPart, fMCContainer->GetArray());
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
        (MCtruthPdgCode == -fCandidatePDG && (fMCMode == kD0Reflection || fMCMode == kOnlyWrongPIDAccepted))) {
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
        (MCtruthPdgCode == fCandidatePDG && (fMCMode == kD0Reflection || fMCMode == kOnlyWrongPIDAccepted))) {
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
        (MCtruthPdgCode == -fCandidatePDG && (fMCMode == kBackgroundOnly || fMCMode == kD0Reflection))) {
      if (i != 0) return kFALSE;
      AliDebug(10, "MC truth is D0");
      invMassD = Dcand->InvMassD0();
    }
    else if ((MCtruthPdgCode == -fCandidatePDG && fMCMode == kSignalOnly) ||
             (MCtruthPdgCode == fCandidatePDG && (fMCMode == kBackgroundOnly || fMCMode == kD0Reflection))) {
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
Bool_t AliAnalysisTaskDmesonJetsSub::AnalysisEngine::ExtractDstarAttributes(const AliAODRecoCascadeHF* DstarCand, AliDmesonJetInfo& DmesonJet, UInt_t i)
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
          auto origin = IsPromptCharm(aodMcPart, fMCContainer->GetArray());
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
/// \return One of the enum constants of AliAnalysisTaskDmesonJetsSub::EMesonDecayChannel_t (D0->Kpi or D*->D0pi->Kpipi)
AliAnalysisTaskDmesonJetsSub::EMesonDecayChannel_t AliAnalysisTaskDmesonJetsSub::AnalysisEngine::CheckDecayChannel(const AliAODMCParticle* part, TClonesArray* mcArray)
{
  if (!part) return kUnknownDecay;
  if (!mcArray) return kUnknownDecay;

  EMesonDecayChannel_t decay = kUnknownDecay;

  Int_t absPdgPart = TMath::Abs(part->GetPdgCode());

  if (part->GetNDaughters() == 2) {

    AliAODMCParticle* d1 = static_cast<AliAODMCParticle*>(mcArray->At(part->GetDaughterLabel(0)));
    AliAODMCParticle* d2 = static_cast<AliAODMCParticle*>(mcArray->At(part->GetDaughterLabel(1)));

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

/// Checks whether a particle is the result of the hadronization of a charm quark or a bottom quark
///
/// \param part Pointer to an AliAODMCParticle object for which originating quark is required
/// \param mcArray Pointer to a TClonesArray object where to look for particles
///
/// \return A pair: first is either kFromCharm or kFromBottom; second is the pointer to the quark
std::pair<AliAnalysisTaskDmesonJetsSub::EMesonOrigin_t, AliAODMCParticle*> AliAnalysisTaskDmesonJetsSub::AnalysisEngine::IsPromptCharm(const AliAODMCParticle* part, TClonesArray* mcArray)
{
  std::pair<AliAnalysisTaskDmesonJetsSub::EMesonOrigin_t, AliAODMCParticle*> result(kUnknownQuark, 0);

  if (!part) return result;
  if (!mcArray) return result;

  static std::set<UInt_t> partons = { 4, 5 };

  AliAODMCParticle* parton = FindParticleOrigin(part, mcArray, kFindLast, partons);
  if (parton) {
    result.second = parton;
    UInt_t absPdgParton = TMath::Abs(parton->GetPdgCode());
    if (absPdgParton == 4) result.first = kFromCharm;
    else if (absPdgParton == 5) result.first = kFromBottom;
  }

  return result;
}

/// Finds a particle in the fragmentation tree of a final state particle
///
/// \param part Pointer to an AliAODMCParticle object for which originating quark is required
/// \param mcArray Pointer to a TClonesArray object where to look for particles
/// \param mode See documentation of the enum type EFindParticleOriginMode_t
///
/// \return A pointer to the MC particle found in the fragmentation tree

AliAODMCParticle* AliAnalysisTaskDmesonJetsSub::AnalysisEngine::FindParticleOrigin(const AliAODMCParticle* part, TClonesArray* mcArray, EFindParticleOriginMode_t mode)
{
  static std::set<UInt_t> pdgSet;

  return FindParticleOrigin(part, mcArray, mode, pdgSet);
}

/// Finds a particle in the fragmentation tree of a final state particle
///
/// \param part Pointer to an AliAODMCParticle object for which originating quark is required
/// \param mcArray Pointer to a TClonesArray object where to look for particles
/// \param mode See documentation of the enum type EFindParticleOriginMode_t
/// \param pdgSet A set of PDG codes that are being searched
///
/// \return A pointer to the MC particle found in the fragmentation tree

AliAODMCParticle* AliAnalysisTaskDmesonJetsSub::AnalysisEngine::FindParticleOrigin(const AliAODMCParticle* part, TClonesArray* mcArray, EFindParticleOriginMode_t mode, const std::set<UInt_t>& pdgSet)
{
  AliAODMCParticle* result = nullptr;

  Int_t mother = part->GetMother();
  while (mother >= 0) {
    AliAODMCParticle* mcGranma = static_cast<AliAODMCParticle*>(mcArray->At(mother));
    if (mcGranma) {
      UInt_t abspdgGranma = TMath::Abs(mcGranma->GetPdgCode());

      // If the current particle is one of the particle types that is being searched assign it to the result pointer
      if (pdgSet.empty() || pdgSet.count(abspdgGranma) > 0) {
        result = mcGranma;
        // If the last particle in the fragmentation tree (first when going reverse) was requested then stop the loop
        if (mode == kFindLast) break;
      }
      if (mother == mcGranma->GetMother()) { // avoid infinite loop!
        AliWarningClassStream() << "Particle " << mother << " (PDG=" << mcGranma->PdgCode() << ") is the mother of itself!?" << std::endl;
        break;
      }
      mother = mcGranma->GetMother();
    }
    else {
      AliErrorClassStream() << "Could not retrieve mother particle " << mother << "!" << std::endl;
      break;
    }
  }

  return result;
}

/// Run the analysis
void AliAnalysisTaskDmesonJetsSub::AnalysisEngine::RunAnalysis()
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
void AliAnalysisTaskDmesonJetsSub::AnalysisEngine::RunDetectorLevelAnalysis()
{
  // Fill the vertex info of the candidates
  // Needed for reduced delta AOD, where the vertex info has been deleted
  // to reduce the delta AOD file size
  AliAnalysisVertexingHF vHF;

  const Int_t nD = fCandidateArray->GetEntriesFast();

  AliDmesonJetInfo DmesonJet;
  DmesonJet.fEvent = this->fAodEvent;

  std::map<AliHFJetDefinition*,Double_t> maxJetPt;
  for (auto& def : fJetDefinitions) maxJetPt[&def] = 0;
  Double_t maxDPt = 0;

  std::array<int, 3> nAccCharm = {0};
  std::array<std::array<int, 3>, 5> nAccCharmPt = {{{0}}};


   //fill the mc efficiency//
  for (auto& def : fJetDefinitions)GetEfficiencyDenominator(def);

  
  for (Int_t icharm = 0; icharm < nD; icharm++) {   //loop over D candidates
    AliAODRecoDecayHF2Prong* charmCand = static_cast<AliAODRecoDecayHF2Prong*>(fCandidateArray->At(icharm)); // D candidates
    if (!charmCand) continue;
    if(!(vHF.FillRecoCand(fAodEvent,charmCand))) continue;


   
   
    //region of interest + cuts
    if (!fRDHFCuts->IsInFiducialAcceptance(charmCand->Pt(), charmCand->Y(fCandidatePDG))) continue;
    Int_t nMassHypo = 0; // number of mass hypothesis accepted for this D meson
    if (charmCand->Pt() > maxDPt) maxDPt = charmCand->Pt();
    for (Int_t im = 0; im < 2; im++)  {  // 2 mass hypothesis (when available)
      DmesonJet.Reset();
      DmesonJet.fDmesonParticle = charmCand;
      DmesonJet.fSelectionType = im + 1;
      if (ExtractRecoDecayAttributes(charmCand, DmesonJet, im)) {
        for (auto& def : fJetDefinitions) {
          if (FindJet(charmCand, DmesonJet, def,im)) {
            Double_t jetPt = DmesonJet.fJets[def.GetName()].fMomentum.Pt();
	    ExtractEfficiencies(charmCand,DmesonJet,def,im);
            if (jetPt > maxJetPt[&def]) maxJetPt[&def] = jetPt;
          }
          else {
            AliWarning(Form("Could not find jet '%s' for D meson '%s': pT = %.3f, eta = %.3f, phi = %.3f",
                def.GetName(), GetName(), DmesonJet.fD.Pt(), DmesonJet.fD.Eta(), DmesonJet.fD.Phi_0_2pi()));
          }
        }
        fDmesonJets[(icharm+1)*(1-(im*2))] = DmesonJet;
        nMassHypo++;
        nAccCharm[im]++;

        for (int i = 0; i < nAccCharmPt.size(); i++) {
          if (charmCand->Pt() < i) break;
          nAccCharmPt[i][im]++;
        }
      }
    }
    if (nMassHypo == 2) { // both mass hypothesis accepted
      nAccCharm[0]--;
      nAccCharm[1]--;
      nAccCharm[2]++;

      for (int i = 0; i < nAccCharmPt.size(); i++) {
        if (charmCand->Pt() < i) break;
        nAccCharmPt[i][0]--;
        nAccCharmPt[i][1]--;
        nAccCharmPt[i][2]++;
      }

      fDmesonJets[(icharm+1)].fD0D0bar = kTRUE;
      fDmesonJets[-(icharm+1)].fD0D0bar = kTRUE;
    }
  } // end of D cand loop

  TString hname;

  Int_t ntracks = 0;

  for (auto track_cont : fTrackContainers) {
    AliHFTrackContainer* hftrack_cont = dynamic_cast<AliHFTrackContainer*>(track_cont);
    if (hftrack_cont) hftrack_cont->SetDMesonCandidate(nullptr);
    ntracks += track_cont->GetNAcceptEntries();
  }

  for (auto& def : fJetDefinitions) {
    if (!def.fRho) continue;
    hname = TString::Format("%s/%s/fHistRhoVsLeadJetPt", GetName(), def.GetName());
    fHistManager->FillTH2(hname, maxJetPt[&def], def.fRho->GetVal());

    hname = TString::Format("%s/%s/fHistRhoVsLeadDPt", GetName(), def.GetName());
    fHistManager->FillTH2(hname, maxDPt, def.fRho->GetVal());

    hname = TString::Format("%s/%s/fHistRhoVsCent", GetName(), def.GetName());
    fHistManager->FillTH2(hname, fEventInfo.fCent, def.fRho->GetVal());

    hname = TString::Format("%s/%s/fHistLeadJetPtVsCent", GetName(), def.GetName());
    fHistManager->FillTH2(hname, fEventInfo.fCent, maxJetPt[&def]);

    hname = TString::Format("%s/%s/fHistLeadDPtVsCent", GetName(), def.GetName());
    fHistManager->FillTH2(hname, fEventInfo.fCent, maxDPt);

    hname = TString::Format("%s/%s/fHistRhoVsNTracks", GetName(), def.GetName());
    fHistManager->FillTH2(hname, ntracks, def.fRho->GetVal());

    hname = TString::Format("%s/%s/fHistLeadJetPtVsNTracks", GetName(), def.GetName());
    fHistManager->FillTH2(hname, ntracks, maxJetPt[&def]);

    hname = TString::Format("%s/%s/fHistLeadDPtVsNTracks", GetName(), def.GetName());
    fHistManager->FillTH2(hname, ntracks, maxDPt);
  }

  hname = TString::Format("%s/fHistNTotAcceptedDmesons", GetName());
  fHistManager->FillTH1(hname, "D", nAccCharm[0]);
  fHistManager->FillTH1(hname, "Anti-D", nAccCharm[1]);
  fHistManager->FillTH1(hname, "Both", nAccCharm[2]);

  hname = TString::Format("%s/fHistNAcceptedDmesonsVsNtracks", GetName());
  fHistManager->FillTH2(hname, ntracks, nAccCharm[0]+nAccCharm[1]+nAccCharm[2]);

  for (int i = 0; i < nAccCharmPt.size(); i++) {
    hname = TString::Format("%s/fHistNTotAcceptedDmesonsPt%d", GetName(), i);
    fHistManager->FillTH1(hname, "D", nAccCharmPt[i][0]);
    fHistManager->FillTH1(hname, "Anti-D", nAccCharmPt[i][1]);
    fHistManager->FillTH1(hname, "Both", nAccCharmPt[i][2]);

    hname = TString::Format("%s/fHistNAcceptedDmesonsPt%d", GetName(), i);
    fHistManager->FillTH1(hname, nAccCharmPt[i][0]+nAccCharmPt[i][1]+nAccCharmPt[i][2]);
  }

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
Bool_t AliAnalysisTaskDmesonJetsSub::AnalysisEngine::FindJet(AliAODRecoDecayHF2Prong* Dcand, AliDmesonJetInfo& DmesonJet, AliHFJetDefinition& jetDef, Int_t numcand)
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
    Int_t nConst = 1;

    for (UInt_t ic = 0; ic < constituents.size(); ++ic) {
      if (constituents[ic].user_index() == 0) {
        isDmesonJet = kTRUE;
      }
      else if (constituents[ic].user_index() >= 100) {
        if (constituents[ic].pt() > maxChPt) maxChPt = constituents[ic].pt();
        nConst++;
      }
      else if (constituents[ic].user_index() <= -100) {
        totalNeutralPt += constituents[ic].pt();
        if (constituents[ic].pt() > maxNePt) maxChPt = constituents[ic].pt();
        nConst++;
      }
    }

    if (isDmesonJet) {
      DmesonJet.fJets[jetDef.GetName()].fMomentum.SetPxPyPzE(jets_incl[ijet].px(), jets_incl[ijet].py(), jets_incl[ijet].pz(), jets_incl[ijet].E());
      DmesonJet.fJets[jetDef.GetName()].fNConstituents = nConst;
      DmesonJet.fJets[jetDef.GetName()].fMaxChargedPt = maxChPt;
      DmesonJet.fJets[jetDef.GetName()].fMaxNeutralPt = maxNePt;
      DmesonJet.fJets[jetDef.GetName()].fNEF = totalNeutralPt / jets_incl[ijet].pt();
      DmesonJet.fJets[jetDef.GetName()].fArea = jets_incl[ijet].area();
      DmesonJet.fJets[jetDef.GetName()].fCorrPt = DmesonJet.fJets[jetDef.GetName()].fMomentum.Pt() - jets_incl[ijet].area() * rho;
      IterativeDeclustering(ijet,1,jetDef, DmesonJet.fD.M());
      return kTRUE;
    }
    if(!isDmesonJet && numcand!=1) IterativeDeclustering(ijet,0,jetDef,0.); 
  }

  return kFALSE;
}

/// Adds all the particles contained in the container into the fastjet wrapper
///
/// \param cont Pointer to a valid AliEmcalContainer object
void AliAnalysisTaskDmesonJetsSub::AnalysisEngine::AddInputVectors(AliEmcalContainer* cont, Int_t offset, TH2* rejectHist, Double_t eff)
{
  auto itcont = cont->all_momentum();
  for (AliEmcalIterableMomentumContainer::iterator it = itcont.begin(); it != itcont.end(); it++) {
    UInt_t rejectionReason = 0;
    //only physical primaries are accepted for the jet finding
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


void AliAnalysisTaskDmesonJetsSub::AnalysisEngine::IterativeDeclustering(Int_t ijet,Double_t type,AliHFJetDefinition& jetDef, Double_t invmass)
{
  
   double nall = 0;        
   double zg = 0.;         
   double flagSubjet=0;
   double xconstperp=0;
   TString hname;
   TString hname2;
   fastjet::JetAlgorithm jet_algo(fastjet::cambridge_algorithm);
   double jet_radius_ca = 1.0;
   fastjet::JetDefinition jet_def(jet_algo, jet_radius_ca,static_cast<fastjet::RecombinationScheme>(0), fastjet::Best);
   hname = TString::Format("%s/%s/LundIterative", GetName(), jetDef.GetName());
   THnSparse* h = static_cast<THnSparse*>(fHistManager->FindObject(hname)); 
   hname2 = TString::Format("%s/%s/AngleDifference", GetName(), jetDef.GetName());
   TH2* hdiffangle = static_cast<TH2*>(fHistManager->FindObject(hname2)); 	    
      try{
      std::vector<fastjet::PseudoJet> particles(fFastJetWrapper->GetJetConstituents(ijet));
      fastjet::ClusterSequence cs_ca(particles, jet_def);
      std::vector<fastjet::PseudoJet> output_jets = cs_ca.inclusive_jets(0);
      output_jets = sorted_by_pt(output_jets);
         

      fastjet::PseudoJet jj = output_jets[0];
      fastjet::PseudoJet j1; 
      fastjet::PseudoJet j2;  
      
      while(jj.has_parents(j1,j2)){
         nall = nall + 1;
         if(j1.perp() < j2.perp()) std::swap(j1,j2);
         flagSubjet=0;
          vector < fastjet::PseudoJet > constitj1 = sorted_by_pt(j1.constituents());
          if(type==1){ 
	  for(Int_t j=0;j<constitj1.size();j++){
                if(constitj1[j].user_index()==0){
		  xconstperp=constitj1[j].perp();
		  flagSubjet=1; }}}

	 
         double delta_R = j1.delta_R(j2);
	 double delta_Raxis=j2.delta_R(output_jets[0]);
         zg = j2.perp()/(j1.perp()+j2.perp());   
         double yh=j1.e()+j2.e();            
         double y = log(1.0/delta_R);
         double lnpt_rel = log(j2.perp()*delta_R);

	 double lundEntries[10] = {y, lnpt_rel, output_jets[0].perp(), nall, type, flagSubjet, xconstperp, invmass,yh,TMath::Abs(output_jets[0].eta())};
         h->Fill(lundEntries);
	 hdiffangle->Fill(delta_R, delta_Raxis);
                jj=j1;
      }

      if(nall==0){ double lundEntrieszero[10]={0,0,output_jets[0].perp(),0,type,0,0,invmass,0,TMath::Abs(output_jets[0].eta())};
	h->Fill(lundEntrieszero);}
      
      } catch (fastjet::Error) { /*return -1;*/ }
                       
}
         								 
/// Run a particle level analysis
void AliAnalysisTaskDmesonJetsSub::AnalysisEngine::RunParticleLevelAnalysis()
{
  TString hname;

  fMCContainer->SetSpecialPDG(fCandidatePDG);
  fMCContainer->SetRejectedOriginMap(fRejectedOrigin);
  fMCContainer->SetAcceptedDecayMap(fAcceptedDecay);
  fMCContainer->SetRejectISR(fRejectISR);

  if (!fMCContainer->IsSpecialPDGFound()) return;

  std::array<int,2> nAccCharm = {0};
  std::array<std::array<int, 2>, 5> nAccCharmPt = {{{0}}};

  std::map<AliHFJetDefinition*, Double_t> maxJetPt;
  Double_t maxDPt = 0;

  for (auto &jetDef : fJetDefinitions) {
    maxJetPt[&jetDef] = 0;
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
        if (constituent.perp() < 1e-6) continue; // reject ghost particles
        AliAODMCParticle* part = fMCContainer->GetMCParticle(iPart);
        if (!part) {
          ::Error("AliAnalysisTaskDmesonJetsSub::AnalysisEngine::RunParticleLevelAnalysis", "Could not find jet constituent %d!", iPart);
          continue;
        }
        if (TMath::Abs(part->PdgCode()) == fCandidatePDG) {
          nDmesonsInJet++;
          std::map<int, AliDmesonJetInfo>::iterator dMesonJetIt = fDmesonJets.find(iPart);
          if (dMesonJetIt == fDmesonJets.end()) { // This D meson does not exist yet
            if (part->Pt() > maxDPt) maxDPt = part->Pt();
            std::pair<int, AliDmesonJetInfo> element;
            element.first = iPart;
            dMesonJetIt = fDmesonJets.insert(element).first;
            (*dMesonJetIt).second.fD.SetPxPyPzE(part->Px(), part->Py(), part->Pz(), part->E());
            (*dMesonJetIt).second.fDmesonParticle = part;
            (*dMesonJetIt).second.fSelectionType = part->PdgCode() > 0 ? 1 : 2;

            UShort_t p = 0;
            UInt_t rs = 0;

            auto origin = IsPromptCharm(part, fMCContainer->GetArray());
            p = 0;
            rs = origin.first;
            while (rs >>= 1) { p++; }
            (*dMesonJetIt).second.fPartonType = p;
            (*dMesonJetIt).second.fParton = origin.second;

            (*dMesonJetIt).second.fAncestor = FindParticleOrigin(part, fMCContainer->GetArray(), kFindFirst);

            Int_t im = -1;
            if (part->PdgCode() > 0) {  // D0
              im = 0;
            }
            else { // D0bar
              im = 1;
            }

            nAccCharm[im]++;
            for (int i = 0; i < nAccCharmPt.size(); i++) {
              if (part->Pt() < i) break;
              nAccCharmPt[i][im]++;
            }
          }

          (*dMesonJetIt).second.fJets[jetDef.GetName()].fMomentum.SetPxPyPzE(jet.px(), jet.py(), jet.pz(), jet.E());
          (*dMesonJetIt).second.fJets[jetDef.GetName()].fNConstituents = jet.constituents().size();
          (*dMesonJetIt).second.fJets[jetDef.GetName()].fArea = jet.area();
          (*dMesonJetIt).second.fJets[jetDef.GetName()].fCorrPt = (*dMesonJetIt).second.fJets[jetDef.GetName()].fMomentum.Pt() - jet.area() * rho;
          if (jet.perp() > maxJetPt[&jetDef]) maxJetPt[&jetDef] = jet.perp();
        } // if constituent is a D meson
      } // for each constituent
      if (nDmesonsInJet > 0) histNDmesonsVsNconstituents->Fill(jet.constituents().size(), nDmesonsInJet);
    } // for each jet
  } // for each jet definition

  Int_t npart = fMCContainer->GetNAcceptedParticles();

  for (auto& def : fJetDefinitions) {
    if (!def.fRho) continue;
    hname = TString::Format("%s/%s/fHistRhoVsLeadJetPt", GetName(), def.GetName());
    fHistManager->FillTH2(hname, maxJetPt[&def], def.fRho->GetVal());

    hname = TString::Format("%s/%s/fHistRhoVsLeadDPt", GetName(), def.GetName());
    fHistManager->FillTH2(hname, maxDPt, def.fRho->GetVal());

    hname = TString::Format("%s/%s/fHistRhoVsCent", GetName(), def.GetName());
    fHistManager->FillTH2(hname, fEventInfo.fCent, def.fRho->GetVal());

    hname = TString::Format("%s/%s/fHistLeadJetPtVsCent", GetName(), def.GetName());
    fHistManager->FillTH2(hname, fEventInfo.fCent, maxJetPt[&def]);

    hname = TString::Format("%s/%s/fHistLeadDPtVsCent", GetName(), def.GetName());
    fHistManager->FillTH2(hname, fEventInfo.fCent, maxDPt);

    hname = TString::Format("%s/%s/fHistRhoVsNTracks", GetName(), def.GetName());
    fHistManager->FillTH2(hname, npart, def.fRho->GetVal());

    hname = TString::Format("%s/%s/fHistLeadJetPtVsNTracks", GetName(), def.GetName());
    fHistManager->FillTH2(hname, npart, maxJetPt[&def]);

    hname = TString::Format("%s/%s/fHistLeadDPtVsNTracks", GetName(), def.GetName());
    fHistManager->FillTH2(hname, npart, maxDPt);
  }

  if (fDmesonJets.size() != nAccCharm[0]+nAccCharm[1]) AliError(Form("I found %lu mesons (%d)?", fDmesonJets.size(), nAccCharm[0]+nAccCharm[1]));
  hname = TString::Format("%s/fHistNTotAcceptedDmesons", GetName());
  fHistManager->FillTH1(hname, "D", nAccCharm[0]);
  fHistManager->FillTH1(hname, "Anti-D", nAccCharm[1]);

  hname = TString::Format("%s/fHistNAcceptedDmesonsVsNtracks", GetName());
  fHistManager->FillTH2(hname, npart, nAccCharm[0]+nAccCharm[1]);

  for (int i = 0; i < nAccCharmPt.size(); i++) {
    hname = TString::Format("%s/fHistNTotAcceptedDmesonsPt%d", GetName(), i);
    fHistManager->FillTH1(hname, "D", nAccCharmPt[i][0]);
    fHistManager->FillTH1(hname, "Anti-D", nAccCharmPt[i][1]);

    hname = TString::Format("%s/fHistNAcceptedDmesonsPt%d", GetName(), i);
    fHistManager->FillTH1(hname, nAccCharmPt[i][0]+nAccCharmPt[i][1]);
  }

  hname = TString::Format("%s/fHistNDmesons", GetName());
  fHistManager->FillTH1(hname, nAccCharm[0]+nAccCharm[1]); // same as the number of accepted D mesons, since no selection is performed
}



// Definitions of class AliAnalysisTaskDmesonJetsSub

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJetsSub);
/// \endcond

/// This is the default constructor, used for ROOT I/O purposes.
AliAnalysisTaskDmesonJetsSub::AliAnalysisTaskDmesonJetsSub() :
  AliAnalysisTaskEmcalLight(),
  fAnalysisEngines(),
  fEnabledAxis(0),
  fOutputType(kTreeOutput),
  fHistManager(),
  fApplyKinematicCuts(kTRUE),
  fNOutputTrees(0),
  fTrackEfficiency(0),
  fRejectISR(kFALSE),
  fJetAreaType(fastjet::active_area),
  fJetGhostArea(0.005),
  fMCContainer(0),
  fAodEvent(0),
  fFastJetWrapper(0)
{
  SetMakeGeneralHistograms(kTRUE);
}

/// This is the standard named constructor.
///
/// \param name Name of the task
AliAnalysisTaskDmesonJetsSub::AliAnalysisTaskDmesonJetsSub(const char* name, Int_t nOutputTrees) :
  AliAnalysisTaskEmcalLight(name, kTRUE),
  fAnalysisEngines(),
  fEnabledAxis(k2ProngInvMass),
  fOutputType(kTreeOutput),
  fHistManager(name),
  fApplyKinematicCuts(kTRUE),
  fNOutputTrees(nOutputTrees),
  fTrackEfficiency(0),
  fRejectISR(kFALSE),
  fJetAreaType(fastjet::active_area),
  fJetGhostArea(0.005),
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
AliAnalysisTaskDmesonJetsSub::~AliAnalysisTaskDmesonJetsSub()
{
  if (fFastJetWrapper) delete fFastJetWrapper;
}

/// Load D meson cuts from a file.
///
/// \param cutfname Name of the file containing the cut object
/// \param cutsname Name of the object cuts
///
/// \return Pointer to the AliRDHFCuts object if successful, NULL otherwise.
AliRDHFCuts* AliAnalysisTaskDmesonJetsSub::LoadDMesonCutsFromFile(TString cutfname, TString cutsname)
{
  AliRDHFCuts* analysiscuts = 0;
  TFile* filecuts = TFile::Open(cutfname);
  if (!filecuts || filecuts->IsZombie()) {
    ::Error("AliAnalysisTaskDmesonJetsSub::LoadDMesonCutsFromFile", "Input file not found: will use std cuts.");
    filecuts = 0;
  }

  if (filecuts) analysiscuts = dynamic_cast<AliRDHFCuts*>(filecuts->Get(cutsname));

  if (!analysiscuts) {
    ::Error("AliAnalysisTaskDmesonJetsSub::LoadDMesonCutsFromFile", "Could not find analysis cuts '%s' in '%s'.", cutsname.Data(), cutfname.Data());
    if (filecuts) {
      filecuts->ls();
    }
  }
  else {
    ::Info("AliAnalysisTaskDmesonJetsSub::LoadDMesonCutsFromFile", "Cuts '%s' loaded from file '%s'", cutsname.Data(), cutfname.Data());
  }

  return analysiscuts;
}

/// Add a new AnalysisEngine object.
///
/// \param type      One of the enum constants of ECandidateType_t
/// \param cutfname  Name of the file that contains the D meson cut object
/// \param cuttype   Type of RDHF cuts
/// \param MCmode    One of the enum constants of EMCMode_t
/// \param jettype   Jet type
/// \param jetradius Radius of the jet
/// \param rhoName   Name of the rho object for the subtraction of the jet average background
///
/// \return Pointer to the AnalysisEngine added to the list.
AliAnalysisTaskDmesonJetsSub::AnalysisEngine* AliAnalysisTaskDmesonJetsSub::AddAnalysisEngine(ECandidateType_t type, TString cutfname, TString cuttype, EMCMode_t MCmode, EJetType_t jettype, Double_t jetradius, TString rhoName)
{
  AliHFJetDefinition jetDef(jettype, jetradius, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, rhoName);
  return AddAnalysisEngine(type, cutfname, cuttype, MCmode, jetDef, rhoName);
}

/// Add a new AnalysisEngine object.
///
/// \param type      One of the enum constants of ECandidateType_t
/// \param cutfname  Name of the file that contains the D meson cut object
/// \param cuttype   Type of RDHF cuts
/// \param MCmode    One of the enum constants of EMCMode_t
/// \param jetDef    Jet definition
/// \param rhoName   Name of the rho object for the subtraction of the jet average background
///
/// \return Pointer to the AnalysisEngine added to the list.
AliAnalysisTaskDmesonJetsSub::AnalysisEngine* AliAnalysisTaskDmesonJetsSub::AddAnalysisEngine(ECandidateType_t type, TString cutfname, TString cuttype, EMCMode_t MCmode, const AliHFJetDefinition& jetDef, TString rhoName)
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

    if (!cuttype.IsNull()) {
      cutsname += TString::Format("_%s", cuttype.Data());
    }

    cuts = LoadDMesonCutsFromFile(cutfname, cutsname);
    if (cuts) cuts->PrintAll();
  }

  AnalysisEngine eng(type, MCmode, cuts);

  std::list<AnalysisEngine>::iterator it = FindAnalysisEngine(eng);

  if (it == fAnalysisEngines.end() || *it != eng) {  // No analysis engine was found, adding a new one
    eng.AddJetDefinition(jetDef);
    it = fAnalysisEngines.insert(it, eng);
    ::Info("AliAnalysisTaskDmesonJetsSub::AddAnalysisEngine", "A new analysis engine '%s' has been added. The total number of analysis engines is %lu.", eng.GetName(), fAnalysisEngines.size());
  }
  else {
    AnalysisEngine* found_eng = &(*it);
    ::Info("AliAnalysisTaskDmesonJetsSub::AddAnalysisEngine", "An analysis engine '%s' with %lu jet definitions has been found. The total number of analysis engines is %lu. A new jet definition '%s' is being added.", found_eng->GetName(), found_eng->fJetDefinitions.size(), fAnalysisEngines.size(), jetDef.GetName());
    found_eng->AddJetDefinition(jetDef);

    if (cuts) {
      if (found_eng->fRDHFCuts != 0) ::Warning("AliAnalysisTaskDmesonJetsSub::AddAnalysisEngine", "D meson cuts were already defined for this D meson type. They will be overwritten.");
      found_eng->SetRDHFCuts(cuts);
    }
  }

  return &(*it);
}

std::list<AliAnalysisTaskDmesonJetsSub::AnalysisEngine>::iterator AliAnalysisTaskDmesonJetsSub::FindAnalysisEngine(const AliAnalysisTaskDmesonJetsSub::AnalysisEngine& eng)
{
  std::list<AnalysisEngine>::iterator it = fAnalysisEngines.begin();
  while (it != fAnalysisEngines.end() && (*it) != eng) it++;
  return it;
}

/// Creates the output containers.
void AliAnalysisTaskDmesonJetsSub::UserCreateOutputObjects()
{
  ::Info("UserCreateOutputObjects", "CreateOutputObjects of task %s", GetName());

  AliAnalysisTaskEmcalLight::UserCreateOutputObjects();

  // Define histograms
  // the TList fOutput is already defined in  AliAnalysisTaskEmcalLight::UserCreateOutputObjects()

  TString hname;
  TString htitle;
  TH1* h = 0;
  Int_t treeSlot = 0;

  Int_t maxTracks = 6000;
  Double_t maxRho = 500;
  if (fForceBeamType == kpp) {
    maxRho = 50;
    maxTracks = 200;
  }
  else if (fForceBeamType == kpA) {
    maxRho = 200;
    maxTracks = 500;
  }

      Int_t dimx   = 10;
      Int_t nbinsx[10]   = {50,100,10,20,2,2,200,150,100,9};
      Double_t minx[10] =  {0,-10,0,0,0,0,0,1.6,0,0};
      Double_t maxx[10]  = {5,10,100,20,2,2,100,2.3,100,0.9};
      TString titlex[10]={"log(1/deltaR)","log(zteta)","jet pt","n","type","flagSubjet","ptD","invmass","frac","abs(eta)"};



  
  hname = "fHistCharmPt";
  htitle = hname + ";#it{p}_{T,charm} (GeV/#it{c});counts";
  fHistManager.CreateTH1(hname, htitle, 500, 0, 1000);

  hname = "fHistCharmEta";
  htitle = hname + ";#eta_{charm};counts";
  fHistManager.CreateTH1(hname, htitle, 400, -10, 10);

  hname = "fHistCharmPhi";
  htitle = hname + ";#phi_{charm};counts";
  fHistManager.CreateTH1(hname, htitle, 125, 0, TMath::TwoPi());

  hname = "fHistBottomPt";
  htitle = hname + ";#it{p}_{T,bottom} (GeV/#it{c});counts";
  fHistManager.CreateTH1(hname, htitle, 500, 0, 1000);

  hname = "fHistBottomEta";
  htitle = hname + ";#eta_{bottom};counts";
  fHistManager.CreateTH1(hname, htitle, 400, -10, 10);

  hname = "fHistBottomPhi";
  htitle = hname + ";#phi_{bottom};counts";
  fHistManager.CreateTH1(hname, htitle, 125, 0, TMath::TwoPi());

  hname = "fHistHighestPartonPt";
  htitle = hname + ";#it{p}_{T,bottom} (GeV/#it{c});counts";
  fHistManager.CreateTH1(hname, htitle, 500, 0, 1000);

  hname = "fHistHighestPartonType";
  htitle = hname + ";type;counts";
  fHistManager.CreateTH1(hname, htitle, 10, 0, 10);

  hname = "fHistNHeavyQuarks";
  htitle = hname + ";number of heavy-quarks;counts";
  fHistManager.CreateTH1(hname, htitle, 21, -0.5, 20.5);

  ::Info("AliAnalysisTaskDmesonJetsSub::UserCreateOutputObjects", "Allocating histograms for task '%s' (%lu analysis engines)", GetName(), fAnalysisEngines.size());
  for (auto &param : fAnalysisEngines) {
    ::Info("AliAnalysisTaskDmesonJetsSub::UserCreateOutputObjects", "Allocating histograms for analysis engine '%s' (%lu jet definitions)", param.GetName(), param.fJetDefinitions.size());

    param.fHistManager = &fHistManager;

    hname = TString::Format("%s/fHistNAcceptedDmesonsVsNtracks", param.GetName());
    htitle = hname + ";#it{N}_{tracks};#it{N}_{D};events";
    h = fHistManager.CreateTH2(hname, htitle, 251, -0.5, 250.5, 21, -0.5, 20.5);

    hname = TString::Format("%s/fHistNTotAcceptedDmesons", param.GetName());
    htitle = hname + ";;#it{N}_{D}";
    h = fHistManager.CreateTH1(hname, htitle, 3, 0, 3);

    for (int i = 0 ; i < 5; i++) {
      hname = TString::Format("%s/fHistNAcceptedDmesonsPt%d", param.GetName(), i);
      htitle = hname + ";#it{N}_{D};events";
      h = fHistManager.CreateTH1(hname, htitle, 21, -0.5, 20.5);

      hname = TString::Format("%s/fHistNTotAcceptedDmesonsPt%d", param.GetName(), i);
      htitle = hname + ";;#it{N}_{D}";
      h = fHistManager.CreateTH1(hname, htitle, 3, 0, 3);
    }

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

     hname = TString::Format("%s/EfficiencyMatchesPrompt",param.GetName());
      htitle = hname + ";D meson matches prompt";
      fHistManager.CreateTH2(hname,htitle,100,0,50,20,0,100);

      hname = TString::Format("%s/EfficiencyGeneratorPrompt",param.GetName());
      htitle = hname + ";D meson part level prompt";
      fHistManager.CreateTH2(hname,htitle,100,0,50,20,0,100);

      hname = TString::Format("%s/EfficiencyMatchesNonPrompt",param.GetName());
      htitle = hname + ";D meson matches non prompt";
      fHistManager.CreateTH2(hname,htitle,100,0,50,20,0,100);

      hname = TString::Format("%s/EfficiencyGeneratorNonPrompt",param.GetName());
      htitle = hname + ";D meson part level non prompt";
      fHistManager.CreateTH2(hname,htitle,100,0,50,20,0,100);

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
      hname = TString::Format("%s/fHistPartonPt", param.GetName());
      htitle = hname + ";#it{p}_{T,parton} (GeV/#it{c});counts";
      fHistManager.CreateTH1(hname, htitle, 500, 0, 1000);

      hname = TString::Format("%s/fHistPartonEta", param.GetName());
      htitle = hname + ";#eta_{parton};counts";
      fHistManager.CreateTH1(hname, htitle, 400, -10, 10);

      hname = TString::Format("%s/fHistPartonPhi", param.GetName());
      htitle = hname + ";#phi_{parton};counts";
      fHistManager.CreateTH1(hname, htitle, 125, 0, TMath::TwoPi());

      hname = TString::Format("%s/fHistPartonType", param.GetName());
      htitle = hname + ";type;counts";
      fHistManager.CreateTH1(hname, htitle, 10, 0, 10);

      hname = TString::Format("%s/fHistPrompt", param.GetName());
      htitle = hname + ";Type;counts";
      h = fHistManager.CreateTH1(hname, htitle, 3, 0, 3);
      h->GetXaxis()->SetBinLabel(1, "Unknown");
      h->GetXaxis()->SetBinLabel(2, "Prompt");
      h->GetXaxis()->SetBinLabel(3, "Non-Prompt");

      hname = TString::Format("%s/fHistAncestor", param.GetName());
      htitle = hname + ";Ancestor;counts";
      h = fHistManager.CreateTH1(hname, htitle, 4, 0, 4);
      h->GetXaxis()->SetBinLabel(1, "Unknown");
      h->GetXaxis()->SetBinLabel(2, "Charm");
      h->GetXaxis()->SetBinLabel(3, "Bottom");
      h->GetXaxis()->SetBinLabel(4, "Proton");
    }

    for (auto& jetDef : param.fJetDefinitions) {
      ::Info("AliAnalysisTaskDmesonJetsSub::UserCreateOutputObjects", "Allocating histograms for jet definition '%s'", jetDef.GetName());

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

      if (!jetDef.fRhoName.IsNull()) {
        hname = TString::Format("%s/%s/fHistRhoVsLeadJetPt", param.GetName(), jetDef.GetName());
        htitle = hname + ";#it{p}_{T,jet} (GeV/#it{c});#rho (GeV/#it{c} #times rad^{-1});counts";
        fHistManager.CreateTH2(hname, htitle, 300, 0, 150, 1000, 0, maxRho);

        hname = TString::Format("%s/%s/fHistRhoVsLeadDPt", param.GetName(), jetDef.GetName());
        htitle = hname + ";#it{p}_{T,D} (GeV/#it{c});#rho (GeV/#it{c} #times rad^{-1});counts";
        fHistManager.CreateTH2(hname, htitle, 300, 0, 150, 1000, 0, maxRho);

        hname = TString::Format("%s/%s/fHistRhoVsCent", param.GetName(), jetDef.GetName());
        htitle = hname + ";Centrality (%);#rho (GeV/#it{c} #times rad^{-1});counts";
        fHistManager.CreateTH2(hname, htitle, 100, 0, 100, 1000, 0, maxRho);

        hname = TString::Format("%s/%s/fHistLeadJetPtVsCent", param.GetName(), jetDef.GetName());
        htitle = hname + ";Centrality (%);#it{p}_{T,jet} (GeV/#it{c});counts";
        fHistManager.CreateTH2(hname, htitle, 100, 0, 100, 300, 0, 150);

        hname = TString::Format("%s/%s/fHistLeadDPtVsCent", param.GetName(), jetDef.GetName());
        htitle = hname + ";Centrality (%);#it{p}_{T,D} (GeV/#it{c});counts";
        fHistManager.CreateTH2(hname, htitle, 100, 0, 100, 300, 0, 150);

        hname = TString::Format("%s/%s/fHistRhoVsNTracks", param.GetName(), jetDef.GetName());
        htitle = hname + ";no. of tracks;#rho (GeV/#it{c} #times rad^{-1});counts";
        fHistManager.CreateTH2(hname, htitle, 200, 0, maxTracks, 1000, 0, maxRho);

        hname = TString::Format("%s/%s/fHistLeadJetPtVsNTracks", param.GetName(), jetDef.GetName());
        htitle = hname + ";no. of tracks;#it{p}_{T,jet} (GeV/#it{c});counts";
        fHistManager.CreateTH2(hname, htitle, 200, 0, maxTracks, 300, 0, 150);

        hname = TString::Format("%s/%s/fHistLeadDPtVsNTracks", param.GetName(), jetDef.GetName());
        htitle = hname + ";no. of tracks;#it{p}_{T,D} (GeV/#it{c});counts";
        fHistManager.CreateTH2(hname, htitle, 200, 0, maxTracks, 300, 0, 150);
      }
      hname = TString::Format("%s/%s/LundIterative",param.GetName(),jetDef.GetName());
      cout<<"at the begining"<<hname<<endl;
      THnSparse* h = fHistManager.CreateTHnSparse(hname,hname,dimx,nbinsx,minx,maxx);
      for (Int_t j = 0; j < dimx; j++) {
      h->GetAxis(j)->SetTitle(titlex[j]);}

     
      
      hname = TString::Format("%s/%s/AngleDifference",param.GetName(),jetDef.GetName());
      htitle = hname + ";angle iterative declustering;angle to axis";
      fHistManager.CreateTH2(hname,htitle,100,0.,2*3.1416,100,0,2*3.1416);
      
    }
    switch (fOutputType) {
    case kTreeOutput:
    {
      OutputHandlerTTree* tree_handler = new OutputHandlerTTree(&param);
      param.fOutputHandler = tree_handler;
      tree_handler->BuildOutputObject(GetName());
      if (treeSlot < fNOutputTrees) {
        tree_handler->AssignDataSlot(treeSlot+2);
        treeSlot++;
        PostDataFromAnalysisEngine(tree_handler);
      }
      else {
        AliError(Form("Number of data output slots %d not sufficient. Tree of analysis engine %s will not be posted!", fNOutputTrees, param.GetName()));
      }
    }
    break;
    case kTHnOutput:
    {
      OutputHandlerTHnSparse* thnsparse_handler = new OutputHandlerTHnSparse(&param);
      param.fOutputHandler = thnsparse_handler;
      thnsparse_handler->SetEnabledAxis(fEnabledAxis);
      thnsparse_handler->BuildOutputObject(GetName());
    }
    break;
    case kOnlyQAOutput:
    {
      OutputHandler* qa_handler = new OutputHandler(&param);
      param.fOutputHandler = qa_handler;
    }
    break;
    case kNoOutput:
      break;
    case kTreeExtendedOutput:
    {
      OutputHandlerTTreeExtendedBase* tree_handler = OutputHandlerTTreeExtendedBase::GenerateOutputHandler(&param);
      param.fOutputHandler = tree_handler;
      tree_handler->BuildOutputObject(GetName());
      if (treeSlot < fNOutputTrees) {
        tree_handler->AssignDataSlot(treeSlot+2);
        treeSlot++;
        PostDataFromAnalysisEngine(tree_handler);
      }
      else {
        AliError(Form("Number of data output slots %d not sufficient. Tree of analysis engine %s will not be posted!", fNOutputTrees, param.GetName()));
      }
    }
    break;
    }
  }

  fOutput->Add(fHistManager.GetListOfHistograms());

  PostData(1, fOutput);
}

/// Does some specific initializations for the analysis engines,
/// then calls the base class ExecOnce() method.
void AliAnalysisTaskDmesonJetsSub::ExecOnce()
{
  AliAnalysisTaskEmcalLight::ExecOnce();

  // Load the event
  fAodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);

  fFastJetWrapper = new AliFJWrapper(fName, fTitle);

  fFastJetWrapper->SetAreaType((fastjet::AreaType)fJetAreaType);
  fFastJetWrapper->SetGhostArea(fJetGhostArea);

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
    params.fRejectISR = fRejectISR;
    params.fRandomGen = rnd;

    for (auto &jetdef: params.fJetDefinitions) {
      if (!jetdef.fRhoName.IsNull()) {
        jetdef.fRho = dynamic_cast<AliRhoParameter*>(fInputEvent->FindListObject(jetdef.fRhoName));
        if (!jetdef.fRho) {
          ::Error("AliAnalysisTaskDmesonJetsSub::ExecOnce",
              "%s: Could not find rho object '%s' for engine '%s'",
              GetName(), jetdef.fRhoName.Data(), params.GetName());
        }
      }
    }

    if (!params.fRDHFCuts) {
      if (params.fMCMode == kMCTruth) {
      ::Warning("AliAnalysisTaskDmesonJetsSub::ExecOnce",
          "%s: RDHF cuts not provided for engine '%s'.",
          GetName(), params.GetName());
      }
      else {
        ::Error("AliAnalysisTaskDmesonJetsSub::ExecOnce",
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
          ::Error("AliAnalysisTaskDmesonJetsSub::ExecOnce",
              "%s: Objects of type %s in %s are not inherited from %s! Task will be disabled!",
              GetName(), params.fCandidateArray->GetClass()->GetName(), params.fCandidateArray->GetName(), className.Data()); // @suppress("Ambiguous problem")
          params.fCandidateArray = 0;
          params.fInhibit = kTRUE;
        }
      }
      else {
        ::Error("AliAnalysisTaskDmesonJetsSub::ExecOnce",
            "Could not find candidate array '%s', skipping the event. Analysis engine '%s' will be disabled!",
            params.fBranchName.Data(), params.GetName());
        params.fInhibit = kTRUE;
      }
    }

    if (params.fMCMode != kNoMC) {
      if (!params.fMCContainer) {
        ::Error("AliAnalysisTaskDmesonJetsSub::ExecOnce",
            "No MC particle container was provided. Analysis engine '%s' will be disabled!",
            params.GetName());
        params.fInhibit = kTRUE;
      }
    }

    if (params.fMCMode != kMCTruth) {
      if (params.fTrackContainers.size() == 0 && params.fClusterContainers.size() == 0) {
        ::Error("AliAnalysisTaskDmesonJetsSub::ExecOnce",
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
Bool_t AliAnalysisTaskDmesonJetsSub::Run()
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

    eng.fEventInfo = EventInfo(fCent, fEPV0, fEventWeight, fPtHard);

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
Bool_t AliAnalysisTaskDmesonJetsSub::FillHistograms()
{
  for (auto &param : fAnalysisEngines) {
    if (param.fInhibit) continue;
    param.fOutputHandler->FillOutput(fApplyKinematicCuts);
    PostDataFromAnalysisEngine(param.fOutputHandler);
  }
  if (fMCContainer) FillPartonLevelHistograms();
  return kTRUE;
}

/// Fill histograms with parton-level information
void AliAnalysisTaskDmesonJetsSub::FillPartonLevelHistograms()
{
  auto itcont = fMCContainer->all_momentum();
  Int_t nHQ = 0;
  Double_t highestPartonPt = 0;
  Int_t PdgHighParton = 0;
  for (auto it = itcont.begin(); it != itcont.end(); it++) {
    auto part = *it;
    if (part.first.Pt() == 0) continue;

    Int_t PdgCode = part.second->GetPdgCode();

    // Skip all particles that are not either quarks or gluons
    if ((PdgCode < -9 || PdgCode > 9) && PdgCode != 21  && PdgCode != -21) continue;

    AliDebugStream(5) << "Parton " << it.current_index() <<
        " with pdg=" << PdgCode <<
        ", px=" << part.first.Px() <<
        ", py=" << part.first.Py() <<
        ", pz=" << part.first.Pz() <<
        ", n daughters = " << part.second->GetNDaughters() <<
        std::endl;

    // Skip partons that do not have any children
    // Unclear how this can happen, it would seem that this parton were not fragmented by the generator
    if (part.second->GetNDaughters() == 0) continue;

    // Look for highest momentum parton
    if (highestPartonPt < part.first.Pt()) {
      highestPartonPt = part.first.Pt();
      PdgHighParton = PdgCode;
    }

    // Skip partons that are not HF
    if (PdgCode != 4 && PdgCode != 5 && PdgCode != -4 && PdgCode != -5) continue;

    Bool_t lastInPartonShower = kTRUE;
    Bool_t hadronDaughter = kFALSE;
    for (Int_t daughterIndex = part.second->GetDaughterFirst(); daughterIndex <= part.second->GetDaughterLast(); daughterIndex++){
      if (daughterIndex < 0) {
        AliDebugStream(5) << "Could not find daughter index!" << std::endl;
        continue;
      }
      AliAODMCParticle *daughter = fMCContainer->GetMCParticle(daughterIndex);
      if (!daughter) {
        AliDebugStream(5) << "Could not find particle with index " << daughterIndex << "!" << std::endl;
        continue;
      }
      Int_t daughterPdgCode = daughter->GetPdgCode();
      if (daughter->GetMother() != it.current_index()) {
        AliDebugStream(5) << "Particle " << daughterIndex << " with pdg=" << daughterPdgCode <<
            ", px=" << daughter->Px() <<
            ", py=" << daughter->Py() <<
            ", pz=" << daughter->Pz() <<
            ", is not a daughter of " << it.current_index() <<
            "!" << std::endl;
        continue;
      }

      AliDebugStream(5) << "Found daughter " << daughterIndex <<
        " with pdg=" << daughterPdgCode <<
        ", px=" << daughter->Px() <<
        ", py=" << daughter->Py() <<
        ", pz=" << daughter->Pz() <<
        std::endl;
      // Codes between 81 and 100 are for internal MC code use, they may be intermediate states used e.g. in hadronizaion models
      if (daughterPdgCode == PdgCode) lastInPartonShower = kFALSE; // this parton is not the last parton in the shower
      if (TMath::Abs(daughterPdgCode) >= 111 || (daughterPdgCode >= 81 && daughterPdgCode <= 100)) hadronDaughter = kTRUE;
    }
    if (hadronDaughter) {
      AliDebugStream(5) << "This particle has at least a hadron among its daughters!" << std::endl;
      if (!lastInPartonShower) AliDebugStream(2) << "Odly, quark " << it.current_index() << " with PDG " << PdgCode << " (pt = " << part.first.Pt() << ", eta = " <<  part.first.Eta() << ") is not the last in the parton shower but at least a hadron found among its daughters?!" << std::endl;
    }
    else {
      AliDebugStream(5) << "This particle does not have hadrons among its daughters!" << std::endl;
      if (lastInPartonShower) AliDebugStream(2) << "Odly, quark " << it.current_index() << " with PDG " << PdgCode << " (pt = " << part.first.Pt() << ", eta = " <<  part.first.Eta() << ") is the last in the parton shower but no hadron found among its daughters?!" << std::endl;
      continue;
    }

    if (PdgCode == 4 || PdgCode == -4) {
      fHistManager.FillTH1("fHistCharmPt", part.first.Pt());
      fHistManager.FillTH1("fHistCharmEta", part.first.Eta());
      fHistManager.FillTH1("fHistCharmPhi", part.first.Phi_0_2pi());
    }
    else if (PdgCode == 5 || PdgCode == -5) {
      fHistManager.FillTH1("fHistBottomPt", part.first.Pt());
      fHistManager.FillTH1("fHistBottomEta", part.first.Eta());
      fHistManager.FillTH1("fHistBottomPhi", part.first.Phi_0_2pi());
    }
    nHQ++;
  }
  fHistManager.FillTH1("fHistNHeavyQuarks", nHQ);
  fHistManager.FillTH1("fHistHighestPartonPt",highestPartonPt);
  Int_t partonType = 0;
  Int_t absPdgHighParton = TMath::Abs(PdgHighParton);
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
void AliAnalysisTaskDmesonJetsSub::CalculateMassLimits(Double_t range, Int_t pdg, Int_t nbins, Double_t& minMass, Double_t& maxMass)
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
const char* AliAnalysisTaskDmesonJetsSub::GetHFEventRejectionReasonLabel(UInt_t& bitmap)
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
Int_t AliAnalysisTaskDmesonJetsSub::PostDataFromAnalysisEngine(OutputHandler const* handler)
{
  if (handler->GetDataSlotNumber() >= 0 && handler->GetOutputObject()) {
    PostData(handler->GetDataSlotNumber(), handler->GetOutputObject());
    return handler->GetDataSlotNumber();
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
/// \return pointer to the new AliAnalysisTaskDmesonJetsSub task
AliAnalysisTaskDmesonJetsSub* AliAnalysisTaskDmesonJetsSub::AddTaskDmesonJetsSub(TString ntracks, TString nclusters, TString nMCpart, Int_t nMaxTrees, TString suffix)
{
  // Get the pointer to the existing analysis manager via the static access method
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskDmesonJetsSub", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler) {
    ::Error("AddTaskDmesonJetsSub", "This task requires an input event handler");
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

  TString name("AliAnalysisTaskDmesonJetsSub");
  if (strcmp(suffix, "") != 0) {
    name += TString::Format("_%s", suffix.Data());
  }

  AliAnalysisTaskDmesonJetsSub* jetTask = new AliAnalysisTaskDmesonJetsSub(name, nMaxTrees);

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
