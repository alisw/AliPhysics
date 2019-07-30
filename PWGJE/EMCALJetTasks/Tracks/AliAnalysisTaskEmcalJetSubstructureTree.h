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
#include "AliAnalysisEmcalTriggerSelectionHelper.h"
#include <exception>
#include <vector>
#include <string>
#include <TString.h>
#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>

class TH1;
class THistManager;
class TTree;
class AliClusterContainer;
class AliEmcalJet;
class AliParticleContainer;
class AliTrackContainer;

#define EXPERIMENTAL_JETCONSTITUENTS

namespace PWGJE {

namespace EMCALJetTasks {

/**
 * @struct AliNSubjettinessResults
 * @brief Results of the n-subjettiness algorithm
 * @ingroup PWGJETASKS
 */
struct AliNSubjettinessParameters {
  Double_t fOneSubjettiness;      ///< 1-subjettiness
  Double_t fTwoSubjettiness;      ///< 2-subjettiness

  void LinkJetTreeBranches(TTree *jettree, const char *tag);
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
  Double_t fDeltaR;         ///< Delta_r of the branches at the last splitting
  Double_t fMug;            ///< Mass Drop parameter
  Int_t fNDropped;          ///< Number of dropped subjets

  void LinkJetTreeBranches(TTree *jettree, const char *tag);
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
 * @struct AliJetKineParameters
 * @brief Jet kinematic parameters
 * @ingroup PWGJETASKS
 */
struct AliJetKineParameters {
  Double_t fPt;                              ///< Jet Pt
  Double_t fE;                               ///< Jet Energy
  Double_t fMass;                            ///< Jet Mass
  Double_t fEta;                             ///< Jet Eta 
  Double_t fPhi;                             ///< Jet Phi
  Double_t fArea;                            ///< Jet Area
  Double_t fNEF;                             ///< Jet Neutral Energy Fraction
  Int_t    fNCharged;                        ///< Number of charged constituents
  Int_t    fNNeutral;                        ///< Number of neutral constituents
  Double_t fZLeading;                        ///< z of the leading constituent
  Double_t fZLeadingCharged;                 ///< z of the leading charged constituent
  Double_t fZLeadingNeutral;                 ///< z of the leading neutral constituent

  void LinkJetTreeBranches(TTree *jettree, const char *tag);
};

/**
 * @struct AliJetStructureParameters
 * @brief Global jet substructure paramters
 * @ingroup PWGJETASKS
 */
struct AliJetStructureParameters {
  Double_t fAngularity;                       ///< Angularity
  Double_t fPtD;                              ///< Pt dispersion

  void LinkJetTreeBranches(TTree *jettree, const char *tag);
};

struct AliJetTreeGlobalParameters {
  Double_t fJetRadius;                        ///< jet radius
  Double_t fEventWeight;                      ///< event weight (downscale factor)
  Int_t    fTriggerClusterIndex;              ///< Index of the trigger cluster (0 - CENT, 1 - CENTNOTRD)
  Double_t fRhoParamters[4];                  ///< Rho parameters

  void LinkJetTreeBranches(TTree *jettree, bool fillRho);
};

/**
 * @class AliAnalysisTaskEmcalJetSubstructureTree
 * @brief Tree with jet substructure information
 * @ingroup PWGJETASKS
 *
 */
class AliAnalysisTaskEmcalJetSubstructureTree : public AliAnalysisTaskEmcalJet, public AliAnalysisEmcalTriggerSelectionHelperImpl {
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

  enum JetRecType_t {
    kDetLevel = 0,
    kPartLevel = 1
  };

	AliAnalysisTaskEmcalJetSubstructureTree();
	AliAnalysisTaskEmcalJetSubstructureTree(const char *name);
	virtual ~AliAnalysisTaskEmcalJetSubstructureTree();

	void SetTriggerBits(UInt_t triggersel) { fTriggerSelectionBits = triggersel; }
	void SetTriggerString(TString triggerstring) { fTriggerSelectionString = triggerstring; }
	void SetUseDownscaleWeight(Bool_t usedownscale) { fUseDownscaleWeight = usedownscale; }
  void SetGlobalTriggerDecisionContainerName(const char *name) { fNameTriggerDecisionContainer = name; }
  void SetUseTriggerSelectionOnData(bool doUse) { fUseTriggerSelectionForData = doUse; }

	void SetSoftdropDefiniion(Double_t zcut, Double_t betacut, Reclusterizer_t reclusterizer) {
	  fSDZCut = zcut;
	  fSDBetaCut = betacut;
	  fReclusterizer = reclusterizer;
	}

  void SetFillPartLevelBranches(Bool_t doFill) { fFillPart = doFill; }
  void SetFillRhoBranches(Bool_t doFill) { fFillRho = doFill; }
  void SetFillSoftdropBranches(Bool_t doFill) { fFillSoftDrop = doFill; }
  void SetFillNSubjettinessBranches(Bool_t doFill) { fFillNSub = doFill; }
  void SetFillSubstructureBranches(Bool_t doFill) { fFillStructGlob = doFill; }

  void SetUseChargedConstituents(Bool_t doUse) { fUseChargedConstituents = doUse; }
  void SetUseNeutralConstituents(Bool_t doUse) { fUseNeutralConstituents = doUse; }

  void SetHasRecEvent(Bool_t hasrec) { fHasRecEvent = hasrec; }
  void SetHasTrueEvent(Bool_t hastrue) { fHasTrueEvent = hastrue; }

	static AliAnalysisTaskEmcalJetSubstructureTree *AddEmcalJetSubstructureTreeMaker(Bool_t isMC, Bool_t isData, Double_t jetradius, AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recombinationScheme, Bool_t useDCAL, const char *name);

protected:
	virtual void UserCreateOutputObjects();
	virtual bool Run();
	virtual void RunChanged(Int_t newrun);
  virtual void UserExecOnce();
  virtual Bool_t IsTriggerSelected();

	AliJetSubstructureData MakeJetSubstructure(const AliEmcalJet &jet, double jetradius, const AliParticleContainer *tracks, const AliClusterContainer *clusters, const AliJetSubstructureSettings &settings) const;

	AliSoftDropParameters MakeSoftDropParameters(const fastjet::PseudoJet &jet, const AliSoftdropDefinition &cut) const;

	AliNSubjettinessParameters MakeNsubjettinessParameters(const fastjet::PseudoJet &jet, const AliNSubjettinessDefinition &cut) const;
  
  AliJetKineParameters MakeJetKineParameters(const AliEmcalJet &jet, JetRecType_t rectype, const AliParticleContainer *const particles, const AliClusterContainer *const clusters) const;

	Double_t MakeAngularity(const AliEmcalJet &jet, const AliParticleContainer *tracks, const AliClusterContainer *clusters) const;

	Double_t MakePtD(const AliEmcalJet &jet, const AliParticleContainer *const particles, const AliClusterContainer *const clusters) const;

  void FillLuminosity();
  
	void DoConstituentQA(const AliEmcalJet *jet, const AliParticleContainer *tracks, const AliClusterContainer *clusters);

  bool SelectJet(const AliEmcalJet &jet, const AliParticleContainer *particles) const;

private:
	TTree                       *fJetSubstructureTree;        //!<! Tree with jet substructure information
  AliJetTreeGlobalParameters  *fGlobalTreeParams;           //!<! Global jet tree parameters (same for all jets in event)
  AliSoftDropParameters       *fSoftDropMeasured;           //!<! Data field for measured soft drop parameters in jet tree
  AliSoftDropParameters       *fSoftDropTrue;               //!<! Data field for true soft drop parameters in jet tree
  AliNSubjettinessParameters  *fNSubMeasured;               //!<! Data field for measured n-subjettiness parameters in jet tree
  AliNSubjettinessParameters  *fNSubTrue;                   //!<! Data field for true n-subjettiness parameters in jet tree
  AliJetKineParameters        *fKineRec;                    //!<! Detector level jet kinematics
  AliJetKineParameters        *fKineSim;                    //!<! Particle level jet kinematics
  AliJetStructureParameters   *fJetStructureMeasured;       //!<! Measured jet substructure parameters
  AliJetStructureParameters   *fJetStructureTrue;           //!<! True jet substructure paramteres
	THistManager                *fQAHistos;                   //!<! QA histos
  TH1                         *fLumiMonitor;                //!<! Luminosity monitor

	Double_t                     fSDZCut;                     ///< Soft drop z-cut
	Double_t                     fSDBetaCut;                  ///< Soft drop beta cut
	Reclusterizer_t              fReclusterizer;              ///< Reclusterizer method

  Bool_t                       fHasRecEvent;                ///< Has reconstructed event (for trigger selection)
  Bool_t                       fHasTrueEvent;               ///< Has Monte-Carlo truth (for trigger selection)
	UInt_t                       fTriggerSelectionBits;       ///< Trigger selection bits
  TString                      fTriggerSelectionString;     ///< Trigger selection string
  TString                      fNameTriggerDecisionContainer; ///< Global trigger decision container
  Bool_t                       fUseTriggerSelectionForData; ///< Use trigger selection on data (require trigger patch in addition to trigger selection string)
  Bool_t                       fUseDownscaleWeight;         ///< Use 1/downscale as weight
  Bool_t                       fUseChargedConstituents;     ///< Use charged constituents 
  Bool_t                       fUseNeutralConstituents;     ///< Use neutral constituents

  // Fill levels for tree (save disk space when information is not needed)
  Bool_t                       fFillPart;                   ///< Fill particle level information
  Bool_t                       fFillRho;                    ///< Fill rho parameters
  Bool_t                       fFillSoftDrop;               ///< Fill soft drop parameters
  Bool_t                       fFillNSub;                   ///< Fill N-subjettiness
  Bool_t                       fFillStructGlob;             ///< Fill other substructure variables

	ClassDef(AliAnalysisTaskEmcalJetSubstructureTree, 1);
};

/**
 * @brief Helper function linking struct members to branches in the jet substructure tree
 * 
 * @param jettree Jet tree to be linked
 * @param data Data field to be linked
 * @param branchname Name of the branch in the jet tree
 * @param type Variable data type
 */
void LinkBranch(TTree *jettree, void *data, const char *branchname, const char *type);

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIANALYSISTASKEMCALJETSUBSTRUCTURETREE_H */
