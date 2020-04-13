/**************************************************************************************
 * Copyright (C) 2016, Copyright Holders of the ALICE Collaboration                   *
 * All rights reserved.                                                               *
 *                                                                                    *
 * Redistribution and use in source and binary forms, with or without                 *
 * modification, are permitted provided that the following conditions are met:        *
 *     * Redistributions of source code must retain the above copyright               *
 *       notice, this list of conditions and the following disclaimer.                *
 *     * Redistributions in binary form must reproduce the above copyright            *
 *       notice, this list of conditions and the following disclaimer in the          *
 *       documentation and/or other materials provided with the distribution.         *
 *     * Neither the name of the <organization> nor the                               *
 *       names of its contributors may be used to endorse or promote products         *
 *       derived from this software without specific prior written permission.        *
 *                                                                                    *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND    *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED      *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE             *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY                *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES         *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;       *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND        *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT         *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 **************************************************************************************/
#ifndef ALIEMCALJETTASK_H
#define ALIEMCALJETTASK_H

class TClonesArray;
class TObjArray;
class AliVEvent;
class AliEmcalJetUtility;

#include "TF1.h"
#include "TRandom3.h"

#include <AliLog.h>

#include "AliAnalysisTaskEmcal.h"
#include "AliFJWrapper.h"
#include "FJ_includes.h"
#include "AliEmcalJet.h"
#include "AliJetContainer.h"
#if !(defined(__CINT__) || defined(__MAKECINT__))
#include "AliEmcalContainerIndexMap.h"
#endif

namespace fastjet {
  class PseudoJet;
}

/**
 * @class AliEmcalJetTask
 * @brief General jet finder task implementing a wrapper for FastJet
 * @ingroup PWGJEBASE
 * @author Constantin Lozides <cloizides@lbl.gov>, Lawrence Berkeley National Laboratory
 * @author Marta Verweij
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 *
 * This class implements a wrapper for the FastJet jet finder. It allows
 * to set a jet definition (jet algorithm, recombination scheme) and the
 * list of jet constituents. Jet constituents are provided via multiple instances
 * of AliParticleContainer and AliClusterContainer. These classes are delegated for
 * applying cuts and filtering constituents that are then fed to the jet finder.
 * This task will further filter constituents based on whether the jet was
 * defined as being charged, neutral or full. The jet finding is delegated to
 * the class AliFJWrapper which implements an interface to FastJet.
 *
 * The FastJet contrib utilities are available via the AliEmcalJetUtility base class
 * and its derived classes. Utilities can be added via the AddUtility(AliEmcalJetUtility*) method.
 * All the utilities added in the list will be executed. Users can implement new utilities
 * deriving a new class from AliEmcalJetUtility to interface functionalities of the FastJet contribs.
 */
class AliEmcalJetTask : public AliAnalysisTaskEmcal {
 public:

  typedef AliJetContainer::EJetType_t EJetType_t;
  typedef AliJetContainer::EJetAlgo_t EJetAlgo_t;
  typedef AliJetContainer::ERecoScheme_t ERecoScheme_t;

#if !defined(__CINT__) && !defined(__MAKECINT__)
  typedef fastjet::JetAlgorithm FJJetAlgo;
  typedef fastjet::RecombinationScheme FJRecoScheme;
#endif

  AliEmcalJetTask();
  AliEmcalJetTask(const char *name);
  virtual ~AliEmcalJetTask();

  Bool_t Run();

  void                   SetGhostArea(Double_t gharea)              { if (IsLocked()) return; fGhostArea        = gharea; }
  void                   SetJetsName(const char *n)                 { if (IsLocked()) return; fJetsTag          = n     ; }
  void                   SetJetsTag(const char *n)                  { if (IsLocked()) return; fJetsTag          = n     ; }
  void                   SetJetEtaRange(Double_t emi, Double_t ema) { if (IsLocked()) return; fJetEtaMin        = emi   ; fJetEtaMax = ema; }
  void                   SetJetPhiRange(Double_t pmi, Double_t pma) { if (IsLocked()) return; fJetPhiMin        = pmi   ; fJetPhiMax = pma; }
  void                   SetJetAlgo(EJetAlgo_t a)                   { if (IsLocked()) return; fJetAlgo          = a     ; }
  void                   SetJetType(EJetType_t t)                   { if (IsLocked()) return; fJetType          = t     ; }
  void                   SetLocked()                                { fLocked = kTRUE;}
  void                   SetMinJetArea(Double_t a)                  { if (IsLocked()) return; fMinJetArea       = a     ; }
  void                   SetMinJetPt(Double_t j)                    { if (IsLocked()) return; fMinJetPt         = j     ; }
  void                   SetRecombScheme(ERecoScheme_t scheme)      { if (IsLocked()) return; fRecombScheme     = scheme; }
  void                   SetTrackEfficiency(Double_t t)             { if (IsLocked()) return; fTrackEfficiency  = t     ; }
  void                   SetQoverPtShift(Double_t shift)            { if (IsLocked()) return; fQoverPtShift     = shift; fApplyQoverPtShift = true; }
  void                   SetTrackEfficiencyOnlyForEmbedding(Bool_t b) { if (IsLocked()) return; fTrackEfficiencyOnlyForEmbedding = b     ; }
  void                   SetEnableAliBasicParticleCompatibility(Bool_t b) { if (IsLocked()) return; fEnableAliBasicParticleCompatibility = b; }
  void                   SetLegacyMode(Bool_t mode)                 { if (IsLocked()) return; fLegacyMode       = mode  ; }
  void                   SetFillGhost(Bool_t b=kTRUE)               { if (IsLocked()) return; fFillGhost        = b     ; }
  void                   SetRadius(Double_t r)                      { if (IsLocked()) return; fRadius           = r     ; }

  void                   SetEtaRange(Double_t emi, Double_t ema);
  void                   SetMinJetClusPt(Double_t min);
  void                   SetMinJetClusE(Double_t min);
  void                   SetMinJetTrackPt(Double_t min);
  void                   SetPhiRange(Double_t pmi, Double_t pma);

  AliEmcalJetUtility*    AddUtility(AliEmcalJetUtility* utility);

  Double_t               GetGhostArea()                   { return fGhostArea         ; }
  const char*            GetJetsName()                    { return fJetsName.Data()   ; }
  const char*            GetJetsTag()                     { return fJetsTag.Data()    ; }
  Double_t               GetJetEtaMin()                   { return fJetEtaMin         ; }
  Double_t               GetJetEtaMax()                   { return fJetEtaMax         ; }
  Double_t               GetJetPhiMin()                   { return fJetPhiMin         ; }
  Double_t               GetJetPhiMax()                   { return fJetPhiMax         ; }
  UInt_t                 GetJetType()                     { return fJetType           ; }
  UInt_t                 GetJetAlgo()                     { return fJetAlgo           ; }
  Bool_t                 GetLegacyMode()                  { return fLegacyMode        ; }
  Double_t               GetMinJetArea()                  { return fMinJetArea        ; }
  Double_t               GetMinJetPt()                    { return fMinJetPt          ; }
  Int_t                  GetMinMCLabel()                  { return fMinMCLabel        ; }
  Double_t               GetRadius()                      { return fRadius            ; }
  Int_t                  GetRecombScheme()                { return fRecombScheme      ; }
  Double_t               GetTrackEfficiency()             { return fTrackEfficiency   ; }
  Bool_t                 GetTrackEfficiencyOnlyForEmbedding() { return fTrackEfficiencyOnlyForEmbedding; }

  TClonesArray*          GetJets()                        { return fJets              ; }
  TObjArray*             GetUtilities()                   { return fUtilities         ; }

  void                   FillJetConstituents(AliEmcalJet *jet, std::vector<fastjet::PseudoJet>& constituents,
                                             std::vector<fastjet::PseudoJet>& constituents_sub, Int_t flag = 0, TString particlesSubName = "");

  UInt_t                 FindJetAcceptanceType(Double_t eta, Double_t phi, Double_t r);
  
  void                   LoadTrackEfficiencyFunction(const std::string & path, const std::string & name);

  Bool_t                 IsLocked() const;
  void                   SelectCollisionCandidates(UInt_t offlineTriggerMask = AliVEvent::kMB);
  void                   SetType(Int_t t);

  /**
   * @brief Switch for whether to fill the AliEmcalJetConstituent objects (clusters and particles/tracks)
   *
   * By default jet constituent object will be filled to the jet.
   *
   * @param doFill Switch for filling jet consituent object
   */
  void                   SetFillJetConsituents(Bool_t doFill) { fFillConstituents = doFill; }

  static AliEmcalJetTask* AddTaskEmcalJet(
      const TString nTracks                      = "usedefault",
      const TString nClusters                    = "usedefault",
      const AliJetContainer::EJetAlgo_t jetAlgo  = AliJetContainer::antikt_algorithm,
      const Double_t radius                      = 0.4,
      const AliJetContainer::EJetType_t jetType  = AliJetContainer::kFullJet,
      const Double_t minTrPt                     = 0.15,
      const Double_t minClPt                     = 0.30,
      const Double_t ghostArea                   = 0.005,
      const AliJetContainer::ERecoScheme_t reco  = AliJetContainer::pt_scheme,
      const TString tag                          = "Jet",
      const Double_t minJetPt                    = 0.,
      const Bool_t lockTask                      = kTRUE,
      const Bool_t bFillGhosts                   = kFALSE,
      const char *suffix                         = ""
    );

#if !defined(__CINT__) && !defined(__MAKECINT__)
  static FJJetAlgo       ConvertToFJAlgo(EJetAlgo_t algo);
  static FJRecoScheme    ConvertToFJRecoScheme(ERecoScheme_t reco);
#endif

 protected:

  Int_t                  FindJets();
  void                   FillJetBranch();
  void                   ExecOnce();
  void                   InitEvent();
  void                   InitUtilities();
  void                   PrepareUtilities();
  void                   ExecuteUtilities(AliEmcalJet* jet, Int_t ij);
  void                   TerminateUtilities();
  Bool_t                 GetSortedArray(Int_t indexes[], std::vector<fastjet::PseudoJet> array) const;
  Bool_t                 IsJetInEmcal(Double_t eta, Double_t phi, Double_t r);
  Bool_t                 IsJetInDcal(Double_t eta, Double_t phi, Double_t r);
  Bool_t                 IsJetInDcalOnly(Double_t eta, Double_t phi, Double_t r);
  Bool_t                 IsJetInPhos(Double_t eta, Double_t phi, Double_t r);

  TString                fJetsTag;                ///< tag of jet collection (usually = "Jets")

  EJetType_t             fJetType;                ///< jet type (full, charged, neutral)
  EJetAlgo_t             fJetAlgo;                ///< jet algorithm (kt, akt, etc)
  ERecoScheme_t          fRecombScheme;           ///< recombination scheme used by fastjet
  Double_t               fRadius;                 ///< jet radius
  Double_t               fMinJetArea;             ///< min area to keep jet in output
  Double_t               fMinJetPt;               ///< min jet pt to keep jet in output
  Double_t               fJetPhiMin;              ///< minimum phi to keep jet in output
  Double_t               fJetPhiMax;              ///< maximum phi to keep jet in output
  Double_t               fJetEtaMin;              ///< minimum eta to keep jet in output
  Double_t               fJetEtaMax;              ///< maximum eta to keep jet in output
  Double_t               fGhostArea;              ///< ghost area
  Double_t               fTrackEfficiency;        ///< artificial tracking inefficiency (0...1)
  Double_t               fQoverPtShift;           ///< artificial q/pt shift
  TObjArray             *fUtilities;              ///< jet utilities (gen subtractor, constituent subtractor etc.)
  Bool_t                 fTrackEfficiencyOnlyForEmbedding; ///< Apply aritificial tracking inefficiency only for embedded tracks
  TF1                   *fTrackEfficiencyFunction;///< Function that describes the artificial tracking efficiency to be applied on top of the nominal tracking efficiency, as a function of track pT
  Bool_t                 fApplyArtificialTrackingEfficiency; ///< Flag to apply artificial tracking efficiency
  Bool_t                 fApplyQoverPtShift;      ///< Apply Q/pt shift
  TRandom3               fRandom;                 //!<! Random number generator for artificial tracking efficiency
  Bool_t                 fLocked;                 ///< true if lock is set
  Bool_t	               fFillConstituents;		 ///< If true jet consituents will be filled to the AliEmcalJet

  TString                fJetsName;               //!<!name of jet collection
  Bool_t                 fIsInit;                 //!<!=true if already initialized
  Bool_t                 fIsPSelSet;              //!<!=true if physics selection was set
  Bool_t                 fIsEmcPart;              //!<!=true if emcal particles are given as input (for clusters)
  Bool_t                 fEnableAliBasicParticleCompatibility; ///< Flag to allow compatibility with AliBasicParticle constituents
  Bool_t                 fLegacyMode;             //!<!=true to enable FJ 2.x behavior
  Bool_t                 fFillGhost;              ///< =true ghost particles will be filled in AliEmcalJet obj

  TClonesArray          *fJets;                   //!<!jet collection
  AliFJWrapper           fFastJetWrapper;         //!<!fastjet wrapper

  static const Int_t     fgkConstIndexShift;      //!<!contituent index shift

#if !(defined(__CINT__) || defined(__MAKECINT__))
  // Handle mapping between index and containers
  AliEmcalContainerIndexMap <AliClusterContainer, AliVCluster> fClusterContainerIndexMap;    //!<! Mapping between index and cluster containers
  AliEmcalContainerIndexMap <AliParticleContainer, AliVParticle> fParticleContainerIndexMap; //!<! Mapping between index and particle containers
#endif

 private:
  AliEmcalJetTask(const AliEmcalJetTask&);            // not implemented
  AliEmcalJetTask &operator=(const AliEmcalJetTask&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliEmcalJetTask, 30);
  /// \endcond
};
#endif
