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
#ifndef ALIANALYSISTASKEMCALJETENERGYSCALE_H
#define ALIANALYSISTASKEMCALJETENERGYSCALE_H

#if !(defined __CINT__ || defined __MAKECINT__)
#if __cplusplus >= 201103L
// In c++11 mode we will rely on c++11 keywords also in header files 
// (i.e. final, override, default, delete)
#define USECXX11HEADERS
#endif
#endif

#include <TString.h>
#include "AliAnalysisTaskEmcalJet.h"
#include "AliJetContainer.h"

class THistManager;

namespace EmcalTriggerJets {

class AliAnalysisTaskEmcalJetEnergyScale : public AliAnalysisTaskEmcalJet {
public:
  AliAnalysisTaskEmcalJetEnergyScale();
  AliAnalysisTaskEmcalJetEnergyScale(const char *name);
  virtual ~AliAnalysisTaskEmcalJetEnergyScale();

  void SetNameDetJetContainer(const char *name)  { fNameDetectorJets = name; }
  void SetNamePartJetContainer(const char *name) { fNameParticleJets = name; }
  void SetTriggerName(const char *name)          { fTriggerSelectionString = name; }

  static AliAnalysisTaskEmcalJetEnergyScale *AddTaskJetEnergyScale(
    AliJetContainer::EJetType_t       jetType,
    Double_t                          radius,
    Bool_t                            useDCAL,
    const char *                      trigger
  );

protected:
  virtual void UserCreateOutputObjects();
  virtual Bool_t Run(); 
  bool IsSelectEmcalTriggers(const TString &triggerstring) const;

private:
  THistManager                *fHistos;                       //!<! Histogram collection
  TString                     fNameDetectorJets;              ///< Name of the data jet container
  TString                     fNameParticleJets;              ///< Name of the MC jet container
  TString                     fTriggerSelectionString;        ///< Trigger selection string
  TString                     fNameTriggerDecisionContainer;  ///< Global trigger decision container

  AliAnalysisTaskEmcalJetEnergyScale(const AliAnalysisTaskEmcalJetEnergyScale &);
  AliAnalysisTaskEmcalJetEnergyScale &operator=(const AliAnalysisTaskEmcalJetEnergyScale &);

  ClassDef(AliAnalysisTaskEmcalJetEnergyScale, 1);
};

}
#endif // ALIANALYSISTASKEMCALJETENERGYSCALE_H