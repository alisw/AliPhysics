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
#include <iostream>
#include <memory>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include "AliEmcalTriggerAlias.h"

ClassImp(PWG::EMCAL::AliEmcalTriggerAlias)

namespace PWG{

namespace EMCAL {

AliEmcalTriggerAlias::AliEmcalTriggerAlias():
  TObject(),
  fTriggerClasses()
{
}

AliEmcalTriggerAlias::AliEmcalTriggerAlias(const char *triggernames):
  TObject(),
  fTriggerClasses()
{
  fTriggerClasses.SetOwner(true);
  DecodeTriggerClasses(triggernames);
}

AliEmcalTriggerAlias::AliEmcalTriggerAlias(const TList &triggernames):
  TObject(),
  fTriggerClasses()
{
  TIter listiter(&triggernames);
  fTriggerClasses.SetOwner(true);
  TObjString *triggerstring(nullptr);
  while((triggerstring = static_cast<TObjString *>(listiter()))){
    fTriggerClasses.Add(new TObjString(*triggerstring));
  }
}

void AliEmcalTriggerAlias::SetTriggerClasses(const char *triggernames) {
  fTriggerClasses.Clear();
  DecodeTriggerClasses(triggernames);
}

void AliEmcalTriggerAlias::SetTriggerClasses(const TList &triggernames) {
  fTriggerClasses.Clear();
  TIter listiter(&triggernames);
  TObjString *triggerstring(nullptr);
  while((triggerstring = static_cast<TObjString *>(listiter()))){
    fTriggerClasses.Add(new TObjString(*triggerstring));
  }
}

void AliEmcalTriggerAlias::DecodeTriggerClasses(const char *triggernames) {
  TString tmpstring(triggernames);
  std::unique_ptr<TObjArray> separated(tmpstring.Tokenize(";"));
  TIter listiter(separated.get());
  TObjString *triggerstring(nullptr);
  while((triggerstring = static_cast<TObjString *>(listiter()))){
    fTriggerClasses.Add(new TObjString(*triggerstring));
  }
}

bool AliEmcalTriggerAlias::HasTriggerClass(const char *triggername) const {
  auto found = fTriggerClasses.FindObject(triggername);
  return found != nullptr;
}

void AliEmcalTriggerAlias::PrintStream(std::ostream &stream) const {
  stream << "Alias for: ";
  TIter listiter(&fTriggerClasses);
  TObjString *tmpstring(nullptr);
  bool first = true;
  while((tmpstring = static_cast<TObjString *>(listiter()))) {
    if(!first) stream << ",";
    stream << " " << tmpstring->String();
    if(first) first = false;
  }
}
    
}

}

std::ostream & operator<< (std::ostream &in, const PWG::EMCAL::AliEmcalTriggerAlias &alias) {
  alias.PrintStream(in);
  return in;
}