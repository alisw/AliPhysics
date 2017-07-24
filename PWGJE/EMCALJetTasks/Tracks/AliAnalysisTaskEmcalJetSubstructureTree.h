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

    const char *what() const throw() { return "Error in reclusterizing in fastjet"; }
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
    kTAreaRec = 4,
    kTAreaSim = 5,
    kTNEFRec = 6,
    kTNEFSim = 7,
    kTMassRec = 8,
    kTMassSim = 9,
    kTZgMeasured = 10,
    kTZgTrue = 11,
    kTRgMeasured = 12,
    kTRgTrue = 13,
    kTMgMeasured = 14,
    kTMgTrue = 15,
    kTPtgMeasured = 16,
    kTPtgTrue = 17,
    kTOneNSubjettinessMeasured = 18,
    kTOneNSubjettinessTrue = 19,
    kTTwoNSubjettinessMeasured = 20,
    kTTwoNSubjettinessTrue = 21,
    kTNCharged = 22,
    kTNNeutral = 23,
    kTNConstTrue = 24,
    kTNDroppedMeasured = 25,
    kTNDroppedTrue = 26,
    kTNVar = 27
  };
	AliAnalysisTaskEmcalJetSubstructureTree();
	AliAnalysisTaskEmcalJetSubstructureTree(const char *name);
	virtual ~AliAnalysisTaskEmcalJetSubstructureTree();

	void SetTriggerBits(UInt_t triggersel) { fTriggerSelectionBits = triggersel; }
	void SetTriggerString(TString triggerstring) { fTriggerSelectionString = triggerstring; }

	void SetSoftdropDefiniion(Double_t zcut, Double_t betacut, Reclusterizer_t reclusterizer) {
	  fSDZCut = zcut;
	  fSDBetaCut = betacut;
	  fReclusterizer = reclusterizer;
	}

	static AliAnalysisTaskEmcalJetSubstructureTree *AddEmcalJetSubstructureTreeMaker(Bool_t isMC, Bool_t isData, Double_t jetradius, const char *name);

protected:
	virtual void UserCreateOutputObjects();
	virtual bool Run();

	AliJetSubstructureData MakeJetSubstructure(const AliEmcalJet &jet, double jetradius, const AliParticleContainer *tracks, const AliClusterContainer *clusters, const AliJetSubstructureSettings &settings) const;

	AliSoftDropParameters MakeSoftDropParameters(const fastjet::PseudoJet &jet, const AliSoftdropDefinition &cut) const;

	AliNSubjettinessParameters MakeNsubjettinessParameters(const fastjet::PseudoJet &jet, const AliNSubjettinessDefinition &cut) const;

	void FillTree(double r, double weight, const AliEmcalJet *datajet, const AliEmcalJet *mcjet, AliSoftDropParameters *dataSoftdrop, AliSoftDropParameters *mcsoftdrop, AliNSubjettinessParameters *dataSubjettiness, AliNSubjettinessParameters *mcSubjettiness);

private:
	TTree                       *fJetSubstructureTree;        //!<! Tree with jet substructure information
	Double_t                     fJetTreeData[kTNVar];        ///< Variable storage for the jet tree

	Double_t                     fSDZCut;                     ///< Soft drop z-cut
	Double_t                     fSDBetaCut;                  ///< Soft drop beta cut
	Reclusterizer_t              fReclusterizer;              ///< Reclusterizer method

	UInt_t                       fTriggerSelectionBits;       ///< Trigger selection bits
	TString                      fTriggerSelectionString;     ///< Trigger selection string


	/// \cond CLASSIMP
	ClassDef(AliAnalysisTaskEmcalJetSubstructureTree, 1);
	/// \endcond
};

} /* namespace EmcalTriggerJets */

#endif /* ALIANALYSISTASKEMCALJETSUBSTRUCTURETREE_H */
