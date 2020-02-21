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
#ifndef ALIANALYSISTASKEMCALSOFTDROPDATA_H
#define ALIANALYSISTASKEMCALSOFTDROPDATA_H

#include <AliAnalysisTaskEmcalJet.h>
#include <string>
#include <vector>

class TBinning;
class THistManager;

namespace PWGJE{

namespace EMCALJetTasks {

class AliAnalysisTaskEmcalSoftDropData : public AliAnalysisTaskEmcalJet {
public:
  enum EReclusterizer_t {
    kCAAlgo = 0,
    kKTAlgo = 1,
    kAKTAlgo = 2
  };

  AliAnalysisTaskEmcalSoftDropData();
  AliAnalysisTaskEmcalSoftDropData(EMCAL_STRINGVIEW name);
  virtual ~AliAnalysisTaskEmcalSoftDropData();

  void SetCustomPtBinning(TBinning *binning) { fPtBinning = binning; }
  void SetBeta(double beta) { fBeta = beta; }
  void SetZcut(double zcut) { fZcut = zcut; }
  void SetReclusterizingAlgorithm(EReclusterizer_t reclust) { fReclusterizer = reclust; }
  void SetUseChargedConstituents(bool doUse) { fUseChargedConstituents = doUse; }
  void SetUseNeutralConstituents(bool doUse) { fUseNeutralConstituents = doUse; }
  void SetSelectTrigger(UInt_t triggerbits, const char *triggerstring) { fTriggerBits = triggerbits; fTriggerString = triggerstring; }
  void SetUseDownscaleWeight(Bool_t doUse) { fUseDownscaleWeight = doUse; }

  static AliAnalysisTaskEmcalSoftDropData *AddTaskEmcalSoftDropData(Double_t jetradius, AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recombinationScheme, EMCAL_STRINGVIEW trigger);

protected:
  virtual void UserCreateOutputObjects();
  virtual void RunChanged(Int_t newrun);
  virtual Bool_t IsTriggerSelected();
  virtual Bool_t Run();

  TBinning *GetDefaultPtBinning() const;
  TBinning *GetZgBinning() const;
  TBinning *GetRgBinning(double R) const;

  Double_t GetDownscaleWeight() const;
  std::vector<double> MakeSoftdrop(const AliEmcalJet &jet, double jetradius, const AliParticleContainer *tracks, const AliClusterContainer *clusters);

private:
  UInt_t                        fTriggerBits;               ///< Trigger selection bits
  std::string                   fTriggerString;             ///< Trigger selection string
  Bool_t                        fUseDownscaleWeight;        ///< Usage of downscale weights
  Double_t                      fBeta;                      ///< Beta
  Double_t                      fZcut;                      ///< Zcut
  EReclusterizer_t              fReclusterizer;             ///< Reclusterizer
  Bool_t                        fUseChargedConstituents;    ///< Use also charged constituents
  Bool_t                        fUseNeutralConstituents;    ///< Use also neutral constituents
  THistManager                  *fHistos;                   //!<! Histogram handler
  TBinning                      *fPtBinning;                ///< Detector level pt binning

  ClassDef(AliAnalysisTaskEmcalSoftDropData, 1)
};

}

}

#endif