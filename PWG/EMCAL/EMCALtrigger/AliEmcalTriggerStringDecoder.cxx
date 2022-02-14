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
#include <sstream>
#include "AliEmcalTriggerStringDecoder.h"

ClassImp(PWG::EMCAL::Triggerinfo);

using namespace PWG::EMCAL;

std::string Triggerinfo::ExpandClassName() const {
  std::string result = fTriggerClass + "-" + fBunchCrossing + "-" + fPastFutureProtection + "-" + fTriggerCluster;
  return result;
}

bool Triggerinfo::IsTriggerClass(EMCAL_STRINGVIEW triggerclass) const {
  return fTriggerClass.substr(1) == triggerclass; // remove C from trigger class part
}

std::vector<Triggerinfo> Triggerinfo::DecodeTriggerString(EMCAL_STRINGVIEW triggerstring) {
  std::vector<Triggerinfo> result;
  std::stringstream triggerparser(triggerstring.data());
  std::string currenttrigger;
  while(std::getline(triggerparser, currenttrigger, ' ')){
    if(!currenttrigger.length()) continue;
    std::vector<std::string> tokens;
    std::stringstream triggerdecoder(currenttrigger);
    std::string token;
    while(std::getline(triggerdecoder, token, '-')) tokens.emplace_back(token);
    result.emplace_back(Triggerinfo{tokens[0], tokens[1], tokens[2], tokens[3]});
  }
  return result;
}
