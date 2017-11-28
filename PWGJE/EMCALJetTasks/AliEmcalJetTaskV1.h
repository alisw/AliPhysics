/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#ifndef ALIEMCALJETTASKV1_H
#define ALIEMCALJETTASKV1_H

#include "AliAnalysisTaskEmcal.h"
#include "AliJetContainer.h"

class TClonesArray;

namespace fastjet {
  class PseudoJet;
}

namespace PWGJE{

namespace EMCALJetTasks {

class AliEmcalJetFinderKernel;

/**
 * @class AliEmcalJetTaskV1
 * @brief General jet finder task implementing a wrapper for FastJet
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
class AliEmcalJetTaskV1 : public AliAnalysisTaskEmcal {
 public:
  AliEmcalJetTaskV1();
  AliEmcalJetTaskV1(const char *name);
  virtual ~AliEmcalJetTaskV1();


  TClonesArray            *GetJets()                      { return fJets; }
  AliEmcalJetFinderKernel *GetJetFinder() const           { return fJetFinder; }
  const TString           &GetJetsName() const            { return fJetsName; }
  const TString           &GetJetsTag() const             { return fJetsTag; }

  void SetJetsTag(const char *tag) { fJetsTag = tag; }

  static AliEmcalJetTaskV1* AddTaskEmcalJet(
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
      const Bool_t bFillGhosts                   = kFALSE
    );

 protected:

  Bool_t                 Run();
  void                   ExecOnce();
  void                   RunChanged(Int_t runnumber);

 private:
  AliEmcalJetFinderKernel  *fJetFinder;              ///< jet finder
  TClonesArray             *fJets;                   //!<!jet collection
  TString                   fJetsTag;                ///< tag of jet collection (usually = "Jets")
  TString                   fJetsName;               //!<!name of jet collection


 private:
  AliEmcalJetTaskV1(const AliEmcalJetTaskV1&);            // not implemented
  AliEmcalJetTaskV1 &operator=(const AliEmcalJetTaskV1&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliEmcalJetTaskV1, 1);
  /// \endcond
};

}
}

#endif
