/************************************************************************************
 * Copyright (C) 2020, Copyright Holders of the ALICE Collaboration                 *
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
#ifndef __ALIANALYSISTASKEMCALQOVERPTSHIFT_H__
#define __ALIANALYSISTASKEMCALQOVERPTSHIFT_H__

#include "AliAnalysisTaskEmcalLight.h"
#include <TString.h>

class THistManager;

namespace PWGJE
{
  
namespace EMCALJetTasks
{

class AliAnalysisTaskEmcalQoverPtShift : public AliAnalysisTaskEmcalLight
{
public:
  AliAnalysisTaskEmcalQoverPtShift();
  AliAnalysisTaskEmcalQoverPtShift(const char *name);
  virtual ~AliAnalysisTaskEmcalQoverPtShift();

  void SetQOverPtShift(Double_t shift) { fQOverPtShift = shift; }
  void SetTriggerSelection(UInt_t triggerbits, const char * triggerstring) { fTriggerBits = triggerbits; fTriggerString = triggerstring; }

  static AliAnalysisTaskEmcalQoverPtShift *AddTaskQOverPtShift(const char *trigger, double shift);

public:
  virtual void UserCreateOutputObjects();
  virtual Bool_t IsTriggerSelected();
  virtual Bool_t Run();

private:
  THistManager *fHistos;          //!<! Histogram manager
  Double_t fQOverPtShift;         ///< Q/pt shift applied in the task
  UInt_t fTriggerBits;            ///< Trigger selection bits
  TString fTriggerString;         ///< Trigger selection string

  ClassDef(AliAnalysisTaskEmcalQoverPtShift, 1);
};

} // namespace EMCALJetTasks
  
} // namespace PWGJE


#endif // __ALIANALYSISTASKEMCALQOVERPTSHIFT_H__