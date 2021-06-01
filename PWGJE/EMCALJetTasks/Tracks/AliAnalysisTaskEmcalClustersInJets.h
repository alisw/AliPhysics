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
#ifndef ALIANALYSISTASKEMCALCLUSTERSINJETS_H
#define ALIANALYSISTASKEMCALCLUSTERSINJETS_H

#include "AliAnalysisTaskEmcalJet.h"
#include <TObjArray.h>
#include <TString.h>

class THistManager;

namespace PWGJE {

namespace EMCALJetTasks {

class AliAnalysisTaskEmcalClustersInJets : public AliAnalysisTaskEmcalJet {
public:
  AliAnalysisTaskEmcalClustersInJets();
  AliAnalysisTaskEmcalClustersInJets(const char *name);
  virtual ~AliAnalysisTaskEmcalClustersInJets();

  void AddNameJetContainer(const char *name);
  void SetNameClusterContainer(const char *name) { fNameClusterContainer = name; }
  void SetNameTriggerClass(const char *name) { fNameTriggerClass = name; }
  
  static AliAnalysisTaskEmcalClustersInJets *AddTaskEmcalClustersInJets(AliJetContainer::EJetType_t jettype, const char *trigger);

protected:
  virtual void UserCreateOutputObjects();
  virtual bool Run();

private:
  THistManager                          *fHistos;               //!<!   Histogram manager
  TObjArray                             fNamesJetContainers;    ///<    Names of the jet containers
  TString                               fNameClusterContainer;  ///<    Name of the cluster container
  TString                               fNameTriggerClass;      ///<    Name of the trigger class to be selected

  ClassDef(AliAnalysisTaskEmcalClustersInJets, 1)
};

}

}
#endif
