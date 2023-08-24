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
#ifndef __ALIANALYSISTASKMEMCALRECALCPATCHESREF_H__
#define __ALIANALYSISTASKMEMCALRECALCPATCHESREF_H__

#include <vector>
#include "AliAnalysisTaskEmcalTriggerBase.h"
#include <TArrayI.h>

class AliEMCALTriggerPatchInfo;

namespace PWGJE {

namespace EMCALJetTasks {

class AliAnalysisTaskEmcalRecalcPatchesRef : public AliAnalysisTaskEmcalTriggerBase {
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
  AliAnalysisTaskEmcalRecalcPatchesRef();
  AliAnalysisTaskEmcalRecalcPatchesRef(const char *name);
  virtual ~AliAnalysisTaskEmcalRecalcPatchesRef() {};

  void SetEnableSumw2(Bool_t doEnable) { fEnableSumw2 = doEnable; }
  void SetOnlineThreshold(ETriggerThreshold_t trigger, Int_t value) { fOnlineThresholds[static_cast<int>(trigger)] = value; }
  void SetFillTHnSparse(bool doFill) { fFillTHnSparse = doFill; }
  void SetSwapPatches(Bool_t doSwap) { fSwapPatches = doSwap; }

  void AddRequiredTriggerOverlap(const char *trigger);
  void AddExcludedTriggerOverlap(const char *trigger);

  /**
   * @brief Set centrality selection.
   *
   * Note: Needs multiplicity task to run in front
   *
   * @param[in] min Min. value of the centrality interval
   * @param[in] max Max. value of the centrality interval
   */
  void SetCentralityRange(double min, double max) { fCentralityRange.SetLimits(min,max); fRequestCentrality = true; }

  static AliAnalysisTaskEmcalRecalcPatchesRef *AddTaskEmcalRecalcPatches(const char *suffix);

protected:
  virtual void CreateUserObjects() {}
  virtual void CreateUserHistos();
  virtual bool IsUserEventSelected();
  virtual bool Run();
  virtual void RunChanged(Int_t newrun);
  virtual void UserFillHistosAfterEventSelection();
  void LoadTriggerThresholdsFromCDB();
  std::vector<const AliEMCALTriggerPatchInfo *> SelectAllPatchesByType(const TClonesArray &patches, EPatchType_t patchtype) const;
  std::vector<const AliEMCALTriggerPatchInfo *> SelectFiredPatchesByTrigger(const TClonesArray &patches, ETriggerThreshold_t trigger) const;
  std::vector<std::string> GetAcceptedTriggerClusters(const char *triggerstring) const;
  int GetNumberNonOverlappingPatchAreas(const std::vector<const AliEMCALTriggerPatchInfo *> &diredpatches) const;
  bool HasOverlap(const AliEMCALTriggerPatchInfo &ref, const AliEMCALTriggerPatchInfo &test) const;
  bool InRange(int test, int includemin, int includemax) const { return test >= includemin && test <= includemax; }

private:
  Bool_t                fEnableSumw2;         ///< Enable sum of weights
  TArrayI               fOnlineThresholds;    ///< Online thresholds
  Bool_t                fSwapPatches;         ///< Look explicitly for the wrong patches
  TObjArray             fRequiredOverlaps;    ///< Add option to require overlap with certain triggers
  TObjArray             fExcludedOverlaps;    ///< Add option to exclude overlap with certain triggers
  AliCutValueRange<double> fCentralityRange;  ///< Range of accepted event centralities
  Bool_t                fUseRecalcPatches;    ///< Switch between offline (FEE) and recalc (L1) patches
  Bool_t                fRequestCentrality;   ///< Switch for request of centrality selection
  Bool_t                fFillTHnSparse;       ///< Switch for filling THnSparse
  Double_t              fEventCentrality;     //!<! Event centrality

  AliAnalysisTaskEmcalRecalcPatchesRef(const AliAnalysisTaskEmcalRecalcPatchesRef &);
  AliAnalysisTaskEmcalRecalcPatchesRef &operator=(const AliAnalysisTaskEmcalRecalcPatchesRef &);

  ClassDef(AliAnalysisTaskEmcalRecalcPatchesRef, 1);
};

}

}
#endif