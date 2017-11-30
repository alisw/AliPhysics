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
#ifndef ALIANALYSISTASKEMCALJETSUBSTRUCTURETREE_H
#define ALIANALYSISTASKEMCALJETSUBSTRUCTURETREE_H

#include "AliAnalysisTaskEmcalJet.h"
#include <exception>
#include <TString.h>
#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>

class THistManager;
class TTree;
class AliClusterContainer;
class AliEmcalJet;
class AliParticleContainer;
class AliTrackContainer;

namespace EmcalTriggerJets {

/**
 * @struct AliNSubjettinessResults
 * @brief Results of the n-subjettiness algorithm
 * @ingroup PWGJETASKS
 */
struct AliNSubjettinessParameters {
  Double_t fOneSubjettiness;      ///< 1-subjettiness
  Double_t fTwoSubjettiness;      ///< 2-subjettiness
};

/**
 * @struct AliNSubjettinessDefiniion
 * @brief Definition of settings for the n-subjettiness algorithm
 * @ingroup PWGJETASKS
 */
struct AliNSubjettinessDefinition{
  Double_t fBeta;                 ///< beta
  Double_t fRadius;               ///< radius
};

/**
 * @struct AliSoftDropParameters
 * @brief Structure for results from the soft drop algorithm
 * @ingroup PWGJETASKS
 */
struct AliSoftDropParameters {
  Double_t fZg;             ///< Groomed jet z
  Double_t fMg;             ///< Groomed jet mass
  Double_t fRg;             ///< Groomed jet radius
  Double_t fPtg;            ///< Groomed jet pt
  Double_t fMug;            ///< Mass Drop parameter
  Int_t fNDropped;          ///< Number of dropped subjets
};

/**
 * @struct AliSoftdropDefinition
 * @brief Definition for the algorithm obtaining the softdrop parameters
 * @ingroup PWGJETASKS
 */
struct AliSoftdropDefinition {
  Double_t fZ;                              ///< Cut on z
  Double_t fBeta;                           ///< Cut on Beta
  fastjet::JetAlgorithm fRecluserAlgo;      ///< Reclusterization algorithm
};

struct AliJetSubstructureSettings {
  AliSoftdropDefinition fSoftdropSettings;
  AliNSubjettinessDefinition fSubjettinessSettings;
};

struct AliJetSubstructureData {
  AliSoftDropParameters fSoftDrop;
  AliNSubjettinessParameters fNsubjettiness;
};

/**
 * @class AliAnalysisTaskEmcalJetSubstructureTree
 * @brief Tree with jet substructure information
 * @ingroup PWGJETASKS
 *
 */
class AliAnalysisTaskEmcalJetSubstructureTree : public AliAnalysisTaskEmcalJet {
public:
  class ReclusterizerException : public std::exception {
  public:
    ReclusterizerException() : std::exception() {}
    virtual ~ReclusterizerException() throw() {}

    virtual const char *what() const throw() { return "Error in reclusterizing in fastjet"; }
  };
  class SubstructureException : public std::exception {
  public:
    SubstructureException() : std::exception() {}
    virtual ~SubstructureException() throw() {}

    virtual const char *what() const throw() { return "Error in builing substructure observable"; }
  };
  class SoftDropException : public std::exception {
  public:
    SoftDropException() : std::exception() {}
    virtual ~SoftDropException() throw() {}

    virtual const char *what() const throw() { return "No associated softdrop structure found for jet - softdrop algorithm failing"; }
  };
  enum Reclusterizer_t {
    kCAAlgo = 0,
    kKTAlgo = 1,
    kAKTAlgo = 2
  };
  enum JetTreeEntry {
    kTRadius = 0,
    kTWeight = 1,
    kTPtJetRec = 2,
    kTPtJetSim = 3,
    kTEJetRec = 4,
    kTEJetSim = 5,
    kTRhoPtRec = 6,
    kTRhoPtSim = 7,
    kTRhoMassRec = 8,
    kTRhoMassSim = 9,
    kTAreaRec = 10,
    kTAreaSim = 11,
    kTNEFRec = 12,
    kTNEFSim = 13,
    kTMassRec = 14,
    kTMassSim = 15,
    kTZgMeasured = 16,
    kTZgTrue = 17,
    kTRgMeasured = 18,
    kTRgTrue = 19,
    kTMgMeasured = 20,
    kTMgTrue = 21,
    kTPtgMeasured = 22,
    kTPtgTrue = 23,
    kTMugMeasured = 24,
    kTMugTrue = 25,
    kTOneNSubjettinessMeasured = 26,
    kTOneNSubjettinessTrue = 27,
    kTTwoNSubjettinessMeasured = 28,
    kTTwoNSubjettinessTrue = 29,
    kTAngularityMeasured = 30,
    kTAngularityTrue = 31,
    kTPtDMeasured = 32,
    kTPtDTrue = 33,
    kTNCharged = 34,
    kTNNeutral = 35,
    kTNConstTrue = 36,
    kTNDroppedMeasured = 37,
    kTNDroppedTrue = 38,
    kTNVar = 39
  };

	AliAnalysisTaskEmcalJetSubstructureTree();
	AliAnalysisTaskEmcalJetSubstructureTree(const char *name);
	virtual ~AliAnalysisTaskEmcalJetSubstructureTree();

	void SetTriggerBits(UInt_t triggersel) { fTriggerSelectionBits = triggersel; }
	void SetTriggerString(TString triggerstring) { fTriggerSelectionString = triggerstring; }
	void SetUseDownscaleWeight(Bool_t usedownscale) { fUseDownscaleWeight = usedownscale; }

	void SetSoftdropDefiniion(Double_t zcut, Double_t betacut, Reclusterizer_t reclusterizer) {
	  fSDZCut = zcut;
	  fSDBetaCut = betacut;
	  fReclusterizer = reclusterizer;
	}

  void SetFillPartLevelBranches(Bool_t doFill) { fFillPart = doFill; }
  void SetFillRhoBranches(Bool_t doFill) { fFillRho = doFill; }
  void SetFillMassBranches(Bool_t doFill) { fFillMass = doFill; }
  void SetFillSoftdropBranches(Bool_t doFill) { fFillSoftDrop = doFill; }
  void SetFillNSubjettinessBranches(Bool_t doFill) { fFillNSub = doFill; }
  void SetFillSubstructureBranches(Bool_t doFill) { fFillStructGlob = doFill; }

	static AliAnalysisTaskEmcalJetSubstructureTree *AddEmcalJetSubstructureTreeMaker(Bool_t isMC, Bool_t isData, Double_t jetradius, AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recombinationScheme, const char *name);

protected:
	virtual void UserCreateOutputObjects();
	virtual bool Run();
	virtual void RunChanged(Int_t newrun);

	AliJetSubstructureData MakeJetSubstructure(const AliEmcalJet &jet, double jetradius, const AliParticleContainer *tracks, const AliClusterContainer *clusters, const AliJetSubstructureSettings &settings) const;

	AliSoftDropParameters MakeSoftDropParameters(const fastjet::PseudoJet &jet, const AliSoftdropDefinition &cut) const;

	AliNSubjettinessParameters MakeNsubjettinessParameters(const fastjet::PseudoJet &jet, const AliNSubjettinessDefinition &cut) const;

	Double_t MakeAngularity(const AliEmcalJet &jet, const AliParticleContainer *tracks, const AliClusterContainer *clusters) const;

	Double_t MakePtD(const AliEmcalJet &jet, const AliParticleContainer *const particles, const AliClusterContainer *const clusters) const;

	void FillTree(double r, double weight, const AliEmcalJet *datajet, const AliEmcalJet *mcjet, AliSoftDropParameters *dataSoftdrop, AliSoftDropParameters *mcsoftdrop, AliNSubjettinessParameters *dataSubjettiness, AliNSubjettinessParameters *mcSubjettiness, Double_t *angularity, Double_t *ptd, Double_t *rhoparameters);

	void DoConstituentQA(const AliEmcalJet *jet, const AliParticleContainer *tracks, const AliClusterContainer *clusters);

  void LinkOutputBranch(const TString &branchname, Double_t *datalocation);

  bool IsPartBranch(const TString &branchname) const;
  bool IsRhoBranch(const TString &branchname) const;
  bool IsMassBranch(const TString &branchname) const;
  bool IsSoftdropBranch(const TString &branchname) const;
  bool IsNSubjettinessBranch(const TString &branchname) const;
  bool IsStructbranch(const TString &branchname) const;

private:
	TTree                       *fJetSubstructureTree;        //!<! Tree with jet substructure information
	Double_t                     fJetTreeData[kTNVar];        ///< Variable storage for the jet tree
	THistManager                *fQAHistos;                   //!<! QA histos

	Double_t                     fSDZCut;                     ///< Soft drop z-cut
	Double_t                     fSDBetaCut;                  ///< Soft drop beta cut
	Reclusterizer_t              fReclusterizer;              ///< Reclusterizer method

	UInt_t                       fTriggerSelectionBits;       ///< Trigger selection bits
  TString                      fTriggerSelectionString;     ///< Trigger selection string
  Bool_t                       fUseDownscaleWeight;         ///< Use 1/downscale as weight

  // Fill levels for tree (save disk space when information is not needed)
  Bool_t                       fFillPart;                   ///< Fill particle level information
  Bool_t                       fFillRho;                    ///< Fill rho parameters
  Bool_t                       fFillMass;                   ///< Fill jet mass
  Bool_t                       fFillSoftDrop;               ///< Fill soft drop parameters
  Bool_t                       fFillNSub;                   ///< Fill N-subjettiness
  Bool_t                       fFillStructGlob;             ///< Fill other substructure variables

	/// \cond CLASSIMP
	ClassDef(AliAnalysisTaskEmcalJetSubstructureTree, 1);
	/// \endcond
};

} /* namespace EmcalTriggerJets */

#endif /* ALIANALYSISTASKEMCALJETSUBSTRUCTURETREE_H */
