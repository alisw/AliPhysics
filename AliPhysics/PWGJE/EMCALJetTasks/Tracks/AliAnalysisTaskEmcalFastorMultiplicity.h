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
#ifndef __ALIANALYSISTASKEMCALFASTORMULTIPLICITY_H__
#define __ALIANALYSISTASKEMCALFASTORMULTIPLICITY_H__

#include "AliAnalysisTaskEmcal.h"
#include <string>

class THistManager;

namespace PWGJE
{
namespace EMCALJetTasks
{

class AliAnalysisTaskEmcalFastorMultiplicity : public AliAnalysisTaskEmcal
{
public:
  AliAnalysisTaskEmcalFastorMultiplicity();
  AliAnalysisTaskEmcalFastorMultiplicity(const char *name);
  virtual ~AliAnalysisTaskEmcalFastorMultiplicity();

  void SetTriggerClass(const char *name) { fTriggerClass = name; }

  static AliAnalysisTaskEmcalFastorMultiplicity *AddTaskEmcalFastorMultiplicity(const char *name);

protected:
  void UserCreateOutputObjects();
  virtual bool IsTriggerSelected();
  bool Run();

private:
  THistManager *fHistos;        //!<! Histogram manager
  std::string fTriggerClass;    ///< Trigger class name

  AliAnalysisTaskEmcalFastorMultiplicity(const AliAnalysisTaskEmcalFastorMultiplicity &);
  AliAnalysisTaskEmcalFastorMultiplicity &operator=(const AliAnalysisTaskEmcalFastorMultiplicity &);

  ClassDef(AliAnalysisTaskEmcalFastorMultiplicity, 1);
};

} // namespace EMCALJetTasks

} // namespace PWGJE

#endif // __ALIANALYSISTASKEMCALFASTORMULTIPLICITY_H__