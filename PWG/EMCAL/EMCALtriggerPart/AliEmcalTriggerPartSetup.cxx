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
#include "AliEmcalTriggerPartSetup.h"

ClassImp(PWG::EMCAL::TriggerPart::AliEmcalTriggerPartSetup);

using namespace PWG::EMCAL::TriggerPart;

/**
 * Default constructor
 */
AliEmcalTriggerPartSetup::AliEmcalTriggerPartSetup() :
  TObject(),
  fThresholds(),
  fTriggerBitConfig()
{
  for( int i = 0; i < 4; i++) fThresholds[i] = -1;
}

/**
 * Copy constructor
 * \param p Reference for the copy
 */
AliEmcalTriggerPartSetup::AliEmcalTriggerPartSetup(const AliEmcalTriggerPartSetup &p) :
  fTriggerBitConfig()
{
  // Copy constructor.
  for( int i = 0; i < 4; i++ ) fThresholds[i] = p.fThresholds[i];
  fTriggerBitConfig.Initialise(p.fTriggerBitConfig);
}

/**
 * Assignment operator
 * @param p Reference for the assignment
 * @return This object
 */
AliEmcalTriggerPartSetup &AliEmcalTriggerPartSetup::operator=(const AliEmcalTriggerPartSetup &p){
  if (this != &p) {
    for( int i = 0; i < 4; i++ )
      fThresholds[i] = p.fThresholds[i];
    fTriggerBitConfig.Initialise(p.fTriggerBitConfig);
  }

  return *this;
}

/**
 * Cleaning function, resets all arrays to 0
 */
void AliEmcalTriggerPartSetup::Clean(){
  for( int i = 0; i < 4; i++ )
    fThresholds[i] = -1;
}
