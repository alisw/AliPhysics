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
#include "AliEmcalTriggerPartBadChannelContainer.h"

ClassImp(PWG::EMCAL::TriggerPart::AliEmcalTriggerPartBadChannelContainer);

using namespace PWG::EMCAL::TriggerPart;

void AliEmcalTriggerPartBadChannelContainer::AddChannel(int col, int row){
  if(HasChannel(col, row)) return;
  fChannels.push_back(TriggerChannelPosition(col, row));
}


bool AliEmcalTriggerPartBadChannelContainer::HasChannel(int col, int row){
  TriggerChannelPosition test(col, row);
  bool found(false);
  for(std::vector<TriggerChannelPosition>::iterator channeliter = fChannels.begin(); channeliter != fChannels.end(); ++channeliter){
	  if(*channeliter == test){
		  found = true;
		  break;
	  }
  }
  return found;
}

bool AliEmcalTriggerPartBadChannelContainer::TriggerChannelPosition::operator==(const TriggerChannelPosition &other) const {
  return fCol == other.fCol && fRow == other.fRow;
}

bool AliEmcalTriggerPartBadChannelContainer::TriggerChannelPosition::operator<(const AliEmcalTriggerPartBadChannelContainer::TriggerChannelPosition &other) const {
  if(fCol == other.fCol){
    if(fRow == other.fRow) return 0;
    else if(fRow < other.fRow) return -1;
    else return 1;
  }
  else if(fCol < other.fCol) return -1;
  else return 1;
}
