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
#ifndef __ALIANALYSISTASKEMCALTRIGGERNORMALIZATION_H__
#define __ALIANALYSISTASKEMCALTRIGGERNORMALIZATION_H__

#include "AliAnalysisTaskEmcal.h"
#include <exception>
#include <string>
#include <vector>

class THistManager;

namespace PWGJE {

namespace EMCALJetTasks {

class AliAnalysisTaskEmcalTriggerNormalization : public AliAnalysisTaskEmcal {
public:
  AliAnalysisTaskEmcalTriggerNormalization();
  AliAnalysisTaskEmcalTriggerNormalization(const char *name);
  virtual ~AliAnalysisTaskEmcalTriggerNormalization() {}

  void SetTriggerCluster(const char *triggercluster) { fTriggerCluster = triggercluster; } 
  void AddMBTriggerClass(const char *triggerclass) { fMBTriggerClasses.emplace_back(triggerclass); }

  static AliAnalysisTaskEmcalTriggerNormalization *AddTaskTriggerNormalization(const char *name);

  class TriggerClusterNotSetException : public std::exception {
  public:
    TriggerClusterNotSetException() {}
    virtual ~TriggerClusterNotSetException() throw() {}

    virtual const char *what() const throw() { return "Trigger cluster not set."; }
  };

  class MBTriggerNotSetException : public std::exception {
  public:
    MBTriggerNotSetException() {}
    virtual ~MBTriggerNotSetException() throw() {}

    virtual const char *what() const throw() { return "No min. bias trigger defined."; }
  };

  class CentralityNotSetException : public std::exception {
  public:
    CentralityNotSetException() {}
    virtual ~CentralityNotSetException() throw() {}

    virtual const char *what() const throw() { return "Centrality estimator not available"; }
  };

protected:
  virtual void UserCreateOutputObjects();
  virtual bool Run();
  virtual void RunChanged(int newrun);
  std::string MatchTrigger(EMCAL_STRINGVIEW triggerstring, const std::vector<std::string> &triggers, EMCAL_STRINGVIEW triggercluster) const;

private:
  THistManager             *fHistos;                ///< List of histograms
  std::string               fTriggerCluster;        ///< Trigger cluster       
  std::vector<std::string>  fMBTriggerClasses;      ///< List of valid min. bias trigger classes


  AliAnalysisTaskEmcalTriggerNormalization(const AliAnalysisTaskEmcalTriggerNormalization &);
  AliAnalysisTaskEmcalTriggerNormalization &operator=(const AliAnalysisTaskEmcalTriggerNormalization &);

  ClassDef(AliAnalysisTaskEmcalTriggerNormalization, 1);
};

}

}


#endif