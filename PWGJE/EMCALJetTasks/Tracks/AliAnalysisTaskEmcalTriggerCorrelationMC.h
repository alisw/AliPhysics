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
#ifndef ALIANALYSISTASKEMCALTRIGGERCORRELATIONMC_H
#define ALIANALYSISTASKEMCALTRIGGERCORRELATIONMC_H

#include "AliAnalysisTaskEmcal.h"
#include <TString.h>

class TH2;

namespace PWGJE {
  
namespace EMCALJetTasks {

class AliAnalysisTaskEmcalTriggerCorrelationMC : public AliAnalysisTaskEmcal {
public:
  AliAnalysisTaskEmcalTriggerCorrelationMC();
  AliAnalysisTaskEmcalTriggerCorrelationMC(const char *name);
  virtual ~AliAnalysisTaskEmcalTriggerCorrelationMC();

  void SetNameTriggerDecisionContainer(const char *name) { fNameTriggerDecisionContainer = name; }

  static AliAnalysisTaskEmcalTriggerCorrelationMC *AddTaskEmcalTriggerCorrelationMC(const char *suffix);

protected:

  virtual void UserCreateOutputObjects();
  virtual bool Run();

private:

  TH2                             *fTriggerCorrelationHist;             //!<! Trigger correlation histogram
  TString                         fNameTriggerDecisionContainer;        ///< Name of the trigger decision container
  
  AliAnalysisTaskEmcalTriggerCorrelationMC(const AliAnalysisTaskEmcalTriggerCorrelationMC &);
  AliAnalysisTaskEmcalTriggerCorrelationMC &operator=(const AliAnalysisTaskEmcalTriggerCorrelationMC &);

  ClassDef(AliAnalysisTaskEmcalTriggerCorrelationMC, 1);
};

}

}

#endif