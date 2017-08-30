/// \class AliHFAODMCParticleContainer
/// \brief Select MC particles based on specific prescriptions of HF analysis
///
/// This class derives from AliParticleContainer. It allows
/// to select MC particles based on specific prescriptions of HF analysis.
/// In particular it will select physical primary particles,
/// exceot for specific D mesons, where the D meson itself
/// will be selected whereas its daughter rejected.
///
/// \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
/// \date Feb 8, 2016

#ifndef ALIHFAODMCPARTICLECONTAINER_H
#define ALIHFAODMCPARTICLECONTAINER_H

/**************************************************************************
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

#include "AliMCParticleContainer.h"
#include "AliAnalysisTaskDmesonJets.h"

class AliHFAODMCParticleContainer : public AliMCParticleContainer {
 public:
  AliHFAODMCParticleContainer();
  AliHFAODMCParticleContainer(const char *name);

  void   SetSpecialPDG(Int_t pdg)          { fSpecialPDG = pdg            ; }

  void   SetRejectQuarkNotFound(Bool_t c)  { fRejectedOrigin = c ?  fRejectedOrigin | AliAnalysisTaskDmesonJets::kUnknownQuark : fRejectedOrigin & ~AliAnalysisTaskDmesonJets::kUnknownQuark; }
  Bool_t GetRejectQuarkNotFound() const    { return (fRejectedOrigin & AliAnalysisTaskDmesonJets::kUnknownQuark) != 0 ; }

  void   SetRejectDfromB(Bool_t c)         { fRejectedOrigin = c ?  fRejectedOrigin | AliAnalysisTaskDmesonJets::kFromBottom : fRejectedOrigin & ~AliAnalysisTaskDmesonJets::kFromBottom; }
  Bool_t GetRejectDfromB() const           { return (fRejectedOrigin & AliAnalysisTaskDmesonJets::kFromBottom) != 0 ; }

  void   SetKeepOnlyDfromB(Bool_t c)       { fRejectedOrigin = c ?  fRejectedOrigin | AliAnalysisTaskDmesonJets::kFromCharm : fRejectedOrigin & ~AliAnalysisTaskDmesonJets::kFromCharm; }
  Bool_t GetKeepOnlyDfromB() const         { return (fRejectedOrigin & AliAnalysisTaskDmesonJets::kFromCharm) != 0 ; }

  void   SetKeepOnlyD0toKpi()              { fAcceptedDecay = AliAnalysisTaskDmesonJets::kDecayD0toKpi     ; }
  void   SetKeepOnlyDStartoKpipi()         { fAcceptedDecay = AliAnalysisTaskDmesonJets::kDecayDStartoKpipi; }

  void   SetRejectedOriginMap(UInt_t m)    { fRejectedOrigin = m; }
  void   SetAcceptedDecayMap(UInt_t m)     { fAcceptedDecay  = m; }

  void   SetRejectISR(Bool_t b)            { fRejectISR      = b; }

  void   SelectCharmtoD0toKpi();
  void   SelectCharmtoDStartoKpipi();
  
  virtual Bool_t AcceptMCParticle(const AliAODMCParticle *vp, UInt_t &rejectionReason) const;
  virtual Bool_t AcceptMCParticle(Int_t i, UInt_t &rejectionReason) const;

  Bool_t IsSpecialPDGFound() const;

  void SetHistOrigin(TH1* h) { fHistOrigin = h; }

 protected:
  Bool_t          IsSpecialPDGDaughter(const AliAODMCParticle* part) const;
  Bool_t          IsSpecialPDG(const AliAODMCParticle* part, TH1* histOrigin = 0) const;

  Int_t           fSpecialPDG;             ///<  include particles with this PDG code even if they are not primary particles (and exclude their daughters)
  UInt_t          fRejectedOrigin;         ///<  Bit mask with D meson origins that are rejected
  UInt_t          fAcceptedDecay;          ///<  Bit mask with D meson decays that are accepted
  Bool_t          fRejectISR;              ///<  Reject initial state radiation
  TH1*            fHistOrigin;             //!<! Book-keeping histogram with origin of special PDG particles

 private:
  AliHFAODMCParticleContainer(const AliHFAODMCParticleContainer&);            // not implemented
  AliHFAODMCParticleContainer &operator=(const AliHFAODMCParticleContainer&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliHFAODMCParticleContainer, 3);
  /// \endcond
};
#endif
