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
#ifndef __ALIEMCALTRIGGERSTRINGDECODER_H__
#define __ALIEMCALTRIGGERSTRINGDECODER_H__
#include <string>
#include <vector>
#include <TObject.h>
#include "AliEmcalStringView.h"

namespace PWG {

namespace EMCAL {

/**
 * @class Triggerinfo
 * @brief Decoded structure of a trigger string
 * 
 * A trigger class string consists of 4 characteristic information
 * - Name of the trigger input
 * - Bunch crossing type
 * - Type of the past-future protection algorithm
 * - Name of the trigger cluster
 * 
 * For easy determinaton of the various information within a trigger
 * string the struct provides access to the different information
 * within the trigger string.
 */
class Triggerinfo : public TObject {
public:
  Triggerinfo(): TObject(), fTriggerClass(), fBunchCrossing(), fPastFutureProtection(), fTriggerCluster() {}
  Triggerinfo(const std::string &triggerclass, const std::string &bc, const std::string &pf, const std::string &clust): 
    TObject(), 
    fTriggerClass(triggerclass), fBunchCrossing(bc), fPastFutureProtection(pf), fTriggerCluster(clust) {}
  virtual ~Triggerinfo() {}

  /**
   * @brief Reconstruct trigger string from information in the Triggerinfo object
   * @return std::string Trigger class representation of the trigger info
   */
  std::string ExpandClassName() const; 

  /**
   * @brief Check if the trigger info corresponds to a certain trigger input class
   * 
   * @param triggerclass Name of the trigger class without C
   * @return true Trigger class matches
   * @return false Trigger class does not match
   */
  bool IsTriggerClass(EMCAL_STRINGVIEW triggerclass) const;

  const std::string &Triggerclass() const { return fTriggerClass; }
  const std::string &BunchCrossing() const { return fBunchCrossing; }
  const std::string &PastFutureProtection() const { return fPastFutureProtection; }
  const std::string &Triggercluster() const { return fTriggerCluster; }

  /**
   * @brief Decoding trigger string
   * 
   * For easy access of the various components of a trigger string
   * the trigger string is decoded into Tiggerinfo structs. Each
   * entry in the struct corresponds to one trigger class present
   * in the trigger string
   * 
   * @param triggerstring Valid trigger class string
   * @return std::vector<Triggerinfo> Trigger info objects for all trigger classes found in the trigger string
   */
  static std::vector<PWG::EMCAL::Triggerinfo> DecodeTriggerString(EMCAL_STRINGVIEW triggerstring);
private:
  std::string fTriggerClass;              ///< Trigger class
  std::string fBunchCrossing;             ///< Bunch crossing type
  std::string fPastFutureProtection;      ///< Type of the past-future protection
  std::string fTriggerCluster;            ///< Trigger cluster
  ClassDef(Triggerinfo, 1);
};


}
}
#endif
