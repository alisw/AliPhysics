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
#ifndef __ALIEMCALMCPARTONINFOCREATOR_H__
#define __ALIEMCALMCPARTONINFOCREATOR_H__

#include "AliAnalysisTaskSE.h"
#include <TString.h>

namespace PWG {

namespace EMCAL {

class AliEmcalMCPartonInfo;

class AliEmcalMCPartonInfoCreator : public AliAnalysisTaskSE {
public:
  AliEmcalMCPartonInfoCreator();
  AliEmcalMCPartonInfoCreator(const char *name);
  virtual ~AliEmcalMCPartonInfoCreator();

  void SetNamePartonInfo(const char *name) { fNamePartonInfo = name; }

  static AliEmcalMCPartonInfoCreator *AddTaskMCPartonInfoCreator(const char *taskname);

protected:
  virtual void UserCreateOutputObjects() {}
  virtual void UserExec(Option_t *option);

private:
  AliEmcalMCPartonInfo *fMCPartonInfo;          //!<! Direct parton info object
  Bool_t fInitialized;                          //!<! Check whether output object is initialized
  TString fNamePartonInfo;                      ///< Name of the MC parton info object

  AliEmcalMCPartonInfoCreator(const AliEmcalMCPartonInfoCreator &);
  AliEmcalMCPartonInfoCreator &operator=(const AliEmcalMCPartonInfoCreator &);

  ClassDef(AliEmcalMCPartonInfoCreator, 1);
};


}

}


#endif