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
#ifndef ALIANALYSISTASKEMCALTRIGGERCORRELATION_H
#define ALIANALYSISTASKEMCALTRIGGERCORRELATION_H

#include "AliAnalysisTaskEmcalTriggerBase.h"
#include "AliCutValueRange.h"
#include "TString.h"

namespace PWGJE {

namespace EMCALJetTasks {

class AliAnalysisTaskEmcalTriggerCorrelation : public AliAnalysisTaskEmcalTriggerBase {
public:
  AliAnalysisTaskEmcalTriggerCorrelation();
  AliAnalysisTaskEmcalTriggerCorrelation(const char *name);
  virtual ~AliAnalysisTaskEmcalTriggerCorrelation();

  static AliAnalysisTaskEmcalTriggerCorrelation *AddTaskTriggerCorrelation(const char *name);

  void SetCentralityRange(double mincent, double maxcent) { fCentralityRange.SetLimits(mincent, maxcent);  fRequestCentrality = true; }
  void SetCentralityEstimator(const char *estimator) { fCentralityEstimator = estimator; }

protected:
  virtual void CreateUserHistos() {}
  virtual void CreateUserObjects() {}
  virtual bool IsUserEventSelected();

private:
  Bool_t                    fRequestCentrality;           ///< Request centrality
  Double_t                  fEventCentrality;             //!<! Event centrality: Tranisent worker variable
  TString                   fCentralityEstimator;         ///< Centrality estimator (default: V0M)
  AliCutValueRange<double>  fCentralityRange;             ///< Selected centrality range

  AliAnalysisTaskEmcalTriggerCorrelation(const AliAnalysisTaskEmcalTriggerCorrelation &);
  AliAnalysisTaskEmcalTriggerCorrelation &operator=(const AliAnalysisTaskEmcalTriggerCorrelation &);

  ClassDef(AliAnalysisTaskEmcalTriggerCorrelation, 1);
};

}

}

#endif