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
#include "AliEmcalTrackSelResultCombined.h"
#include <TObjArray.h>

/// \cond CLASSIMP
ClassImp(PWG::EMCAL::AliEmcalTrackSelResultCombined)
/// \endcond

using namespace PWG::EMCAL;

AliEmcalTrackSelResultCombined::AliEmcalTrackSelResultCombined():
  TObject(),
  fData()
{

}

AliEmcalTrackSelResultCombined::AliEmcalTrackSelResultCombined(const std::vector<PWG::EMCAL::AliEmcalTrackSelResultPtr>& singleSelPointers):
  TObject(),
  fData(singleSelPointers)
 {

 }


AliEmcalTrackSelResultCombined::AliEmcalTrackSelResultCombined(const TObjArray *singleSelPointers):
  TObject(),
  fData(singleSelPointers->GetEntries())
{
  int index(0);
  for(auto o  : *singleSelPointers) fData[index] = *static_cast<AliEmcalTrackSelResultPtr *>(o);
}

const AliEmcalTrackSelResultPtr &AliEmcalTrackSelResultCombined::operator[](int index) const {
  if(index < 0 || index >= static_cast<Int_t>(fData.size())) throw IndexException(index);
  return fData[index];
}

AliEmcalTrackSelResultPtr &AliEmcalTrackSelResultCombined::operator[](int index) {
  if(index < 0 || index >= static_cast<Int_t>(fData.size())) throw IndexException(index);
  return fData[index];
}

Int_t AliEmcalTrackSelResultCombined::GetNumberOfSelectionResults() const {
  return fData.size();
}