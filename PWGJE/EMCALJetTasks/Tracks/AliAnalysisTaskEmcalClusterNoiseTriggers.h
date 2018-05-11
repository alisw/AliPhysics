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
#ifndef __ALIANALYSISTASKMEMCALCLUSTERNOISETRIGGERS_H__
#define __ALIANALYSISTASKMEMCALCLUSTERNOISETRIGGERS_H__

#include <vector>
#include "AliAnalysisTaskEmcalTriggerBase.h"
#include <TArrayI.h>

class AliEMCALTriggerPatchInfo;

namespace EMCalTriggerPtAnalysis {

class AliAnalysisTaskEmcalClusterNoiseTriggers : public AliAnalysisTaskEmcalTriggerBase {
public:
  enum ETriggerThreshold_t {
    kThresholdEG1 = 0,
    kThresholdEG2,
    kThresholdDG1,
    kThresholdDG2,
    kThresholdEJ1,
    kThresholdEJ2,
    kThresholdDJ1,
    kThresholdDJ2
  };
  enum EPatchType_t {
    kEGApatches,
    kDGApatches,
    kEJEpatches,
    kDJEpatches
  };
  AliAnalysisTaskEmcalClusterNoiseTriggers();
  AliAnalysisTaskEmcalClusterNoiseTriggers(const char *name);
  virtual ~AliAnalysisTaskEmcalClusterNoiseTriggers() {};

  void SetEnableSumw2(Bool_t doEnable) { fEnableSumw2 = doEnable; }
  void SetOnlineThreshold(ETriggerThreshold_t trigger, Int_t value) { fOnlineThresholds[static_cast<int>(trigger)] = value; }

  static AliAnalysisTaskEmcalClusterNoiseTriggers *AddTaskEmcalClusterNoiseTriggers(const char *suffix);

protected:
  virtual void CreateUserObjects() {}
  virtual void CreateUserHistos();
  virtual bool Run();
  virtual void UserFillHistosAfterEventSelection();
  std::vector<const AliEMCALTriggerPatchInfo *> SelectAllPatchesByType(const TClonesArray &patches, EPatchType_t patchtype) const;
  std::vector<const AliEMCALTriggerPatchInfo *> SelectFiredPatchesByTrigger(const TClonesArray &patches, ETriggerThreshold_t trigger) const;
  std::vector<std::string> GetAcceptedTriggerClusters(const char *triggerstring) const;
  int GetNumberNonOverlappingPatchAreas(const std::vector<const AliEMCALTriggerPatchInfo *> &diredpatches) const;
  bool HasOverlap(const AliEMCALTriggerPatchInfo &ref, const AliEMCALTriggerPatchInfo &test) const;
  bool InRange(int test, int includemin, int includemax) const { return test >= includemin && test <= includemax; }

private:
  Bool_t                fEnableSumw2;         ///< Enable sum of weights
  TArrayI               fOnlineThresholds;    ///< Online thresholds

  AliAnalysisTaskEmcalClusterNoiseTriggers(const AliAnalysisTaskEmcalClusterNoiseTriggers &);
  AliAnalysisTaskEmcalClusterNoiseTriggers &operator=(const AliAnalysisTaskEmcalClusterNoiseTriggers &);

  ClassDef(AliAnalysisTaskEmcalClusterNoiseTriggers, 1);
};

}
#endif