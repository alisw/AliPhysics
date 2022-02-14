/************************************************************************************
 * Copyright (C) 2018, Copyright Holders of the ALICE Collaboration                 *
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
#ifndef ALIANALYSISTASKEMCALEG1CORRELATION_H
#define ALIANALYSISTASKEMCALEG1CORRELATION_H

#include "AliAnalysisTaskSE.h"

#include <string>
#include <vector>

class TH2;
class TList;

namespace PWG{
namespace EMCAL {
class Triggerinfo;
}
}

namespace PWGJE {
  
namespace EMCALJetTasks {

class AliAnalysisTaskEmcalEG1Correlation : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskEmcalEG1Correlation();
  AliAnalysisTaskEmcalEG1Correlation(const char *taskname);
  virtual ~AliAnalysisTaskEmcalEG1Correlation();

  std::vector<std::string> GetSupportedTriggers() const;
  bool IsTriggerSelectedTS(const std::string &trigger, const std::vector<PWG::EMCAL::Triggerinfo> &triggers) const;
  bool IsTriggerSelectedPS(const std::string &trigger, const std::vector<PWG::EMCAL::Triggerinfo> &triggers) const;

  static AliAnalysisTaskEmcalEG1Correlation *AddTaskEmcalEG1Correlation(const char *taskname);

protected:
  void UserCreateOutputObjects();
  void UserExec(Option_t *);

private:
  TList                       *fOutput;                                 //!<!  Container for output histograms
  TH2                         *fCorrelationHistNoPhysSel;               //!<!  Correlation histogram without physics selection applied
  TH2                         *fCorrelationHistWithPhysSel;             //!<!  Correlation histogram with physics selection applied

  AliAnalysisTaskEmcalEG1Correlation(const AliAnalysisTaskEmcalEG1Correlation &);
  AliAnalysisTaskEmcalEG1Correlation &operator=(const AliAnalysisTaskEmcalEG1Correlation &);

  ClassDef(AliAnalysisTaskEmcalEG1Correlation, 1);
};

}

}
#endif
