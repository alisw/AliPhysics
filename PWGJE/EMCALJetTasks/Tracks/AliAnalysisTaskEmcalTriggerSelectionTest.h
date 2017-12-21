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
#ifndef PWGJE_EMCALJETTASKS_TRACKS_ALIANALYSISTASKEMCALTRIGGERSELECTIONTEST_H
#define PWGJE_EMCALJETTASKS_TRACKS_ALIANALYSISTASKEMCALTRIGGERSELECTIONTEST_H

#include <AliAnalysisTaskEmcal.h>
#include <vector>

namespace PWG{
namespace EMCAL {
class AliEmcalTriggerDecisionContainer;
}
}

namespace PWGJE {
namespace EMCALJetTasks {
namespace Test {

class AliAnalysisTaskEmcalTriggerSelectionTest : public AliAnalysisTaskEmcal{
public:
  AliAnalysisTaskEmcalTriggerSelectionTest();
  AliAnalysisTaskEmcalTriggerSelectionTest(const char *name);
  virtual ~AliAnalysisTaskEmcalTriggerSelectionTest();

  static AliAnalysisTaskEmcalTriggerSelectionTest *AddTaskEmcalTriggerSelectionTest(TString suffix);

protected:

  virtual void UserCreateOutputObjects();
  virtual bool Run();
  virtual bool FillHistograms();
  virtual void UserExecOnce();

  void FillHistosForTrigger(const char *triggername);

private:
  PWG::EMCAL::AliEmcalTriggerDecisionContainer    *fTriggerDecision;        ///< Trigger decision container (from event)
  THistManager                                    *fHistos;                 ///< Histogram handler
  std::vector<TString>                            fSelectedTriggers;        //!<! selected triggers in event

  AliAnalysisTaskEmcalTriggerSelectionTest(const AliAnalysisTaskEmcalTriggerSelectionTest &);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalTriggerSelectionTest, 1);
  /// \endcond
};

} /* namespace Test */
} /* namespace EMCALJetTasks */
} /* namespace PWGJE */

#endif /* PWGJE_EMCALJETTASKS_TRACKS_ALIANALYSISTASKEMCALTRIGGERSELECTIONTEST_H */
