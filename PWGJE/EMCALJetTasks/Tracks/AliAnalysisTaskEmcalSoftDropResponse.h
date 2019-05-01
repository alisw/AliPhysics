/************************************************************************************
 * Copyright (C) 2019, Copyright Holders of the ALICE Collaboration                 *
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
#ifndef ALIANALYSISTASKEMCALSOFTDROPRESPONSE_H
#define ALIANALYSISTASKEMCALSOFTDROPRESPONSE_H

#include <vector>
#include "AliAnalysisTaskEmcalJet.h"

class RooUnfoldResponse;
class TBinning;
class TH2;
class TRandom;

namespace PWGJE{

namespace EMCALJetTasks {

class AliAnalysisTaskEmcalSoftDropResponse : public AliAnalysisTaskEmcalJet {
public:
  enum EBinningMode_t {
    kSDModeINT7,
    kSDModeEJ1,
    kSDModeEJ2,
  };
  enum EReclusterizer_t {
    kCAAlgo = 0,
    kKTAlgo = 1,
    kAKTAlgo = 2
  };

  AliAnalysisTaskEmcalSoftDropResponse();
  AliAnalysisTaskEmcalSoftDropResponse(const char *name);
  virtual ~AliAnalysisTaskEmcalSoftDropResponse();

  void SetBinningMode(EBinningMode_t binmode) { fBinningMode = binmode; }
  void SetCustomPartLevelPtBinning(TBinning *binning) { fPartLevelPtBinning = binning; }
  void SetCustomDetLevelPtBinning(TBinning *binning) { fDetLevelPtBinning = binning; }
  void SetFractionResponseClosure(Double_t fracClosure) { fFractionResponseClosure = fracClosure; }
  void SetBeta(double beta) { fBeta = beta; }
  void SetZcut(double zcut) { fZcut = zcut; }
  void SetReclusterizingAlgorithm(EReclusterizer_t reclust) { fReclusterizer = reclust; }
  void SetSampleFraction(Double_t samplefraction) { fSampleFraction = samplefraction; }
  void SetUseChargedConstituents(bool doUse) { fUseChargedConstituents = doUse; }
  void SetUseNeutralConstituents(bool doUse) { fUseNeutralConstituents = doUse; }
  void SetNameMCParticleContainer(const char *name) { fNameMCParticles = name; }

  static AliAnalysisTaskEmcalSoftDropResponse *AddTaskEmcalSoftDropResponse(Double_t jetradius, AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recombinationScheme, const char *namepartcont, const char *trigger);

protected:
  virtual void UserCreateOutputObjects();
  virtual Bool_t CheckMCOutliers();
  virtual bool Run();

  TBinning *GetDefaultPartLevelPtBinning() const;
  TBinning *GetDefaultDetLevelPtBinning() const;
  TBinning *GetZgBinning() const;

  std::vector<double> MakeSoftdrop(const AliEmcalJet &jet, double jetradius, const AliParticleContainer *tracks, const AliClusterContainer *clusters) const;

private:

  EBinningMode_t                fBinningMode;               ///< Binning adapted to trigger
  Double_t                      fFractionResponseClosure;   ///< Fraction of events used for the response matrix in the closure test
  Double_t                      fZcut;                      ///< Zcut (softdrop definition)
  Double_t                      fBeta;                      ///< Beta (softdrop definition)
  EReclusterizer_t              fReclusterizer;             ///< Reclusterizing algorithm
  Double32_t                    fSampleFraction;            ///< Fraction of statistics used for the analysis
  Bool_t                        fUseChargedConstituents;    ///< Use charged constituents for softdrop
  Bool_t                        fUseNeutralConstituents;    ///< Use neutral constituents for softdrop
  TString                       fNameMCParticles;           ///< Name of the MC particle container
  TRandom                       *fSampleSplitter;           ///< Sample splitter
  TRandom                       *fSampleTrimmer;            ///< Sample trimmer
  TBinning                      *fPartLevelPtBinning;       ///< Particle level pt binning
  TBinning                      *fDetLevelPtBinning;        ///< Detector level pt binning
  RooUnfoldResponse             *fZgResponse;               //!<! RooUnfold response object
  RooUnfoldResponse             *fZgResponseClosure;        //!<! RooUnfold response for the closure test
  TH2                           *fZgPartLevel;              //!<! Zg vs. pt at particle level
  TH2                           *fZgDetLevel;               //!<! Zg vs. pt at detector level
  TH2                           *fZgPartLevelTruncated;     //!<! Zg vs. pt at particle level after truncation at
  TH2                           *fZgPartLevelClosureNoResp; //!<! Zg vs. pt at particle level for closure test (jets not used for response)
  TH2                           *fZgDetLevelClosureNoResp;  //!<! Zg vs. pt at detector level for closure test (jets not used for response)
  TH2                           *fZgPartLevelClosureResp;   //!<! Zg vs. pt at particle level for closure test (jets used for response)
  TH2                           *fZgDetLevelClosureResp;    //!<! Zg vs. pt at detector level for closure test (jets used for response)

  AliAnalysisTaskEmcalSoftDropResponse(const AliAnalysisTaskEmcalSoftDropResponse &);
  AliAnalysisTaskEmcalSoftDropResponse &operator=(const AliAnalysisTaskEmcalSoftDropResponse &);
  
  ClassDef(AliAnalysisTaskEmcalSoftDropResponse, 1);

};

}

}


#endif