/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#include "AliLHCClockPhase.h"
#include "AliDCSValue.h"
#include "AliLog.h"

ClassImp(AliLHCClockPhase)

//______________________________________________________________________________
AliLHCClockPhase::AliLHCClockPhase():
  TObject(),
  fPhaseB1(),
  fPhaseB2()
{
  // default constructor
  // ...
  fPhaseB1.SetOwner();
  fPhaseB2.SetOwner();
}

//______________________________________________________________________________
void AliLHCClockPhase::AddPhaseB1DP(UInt_t timestamp, Float_t phase)
{
  // Add a phase beam1 measurement
  // to the array of data-points

  fPhaseB1.AddLast(new AliDCSValue(phase,timestamp));
}

//______________________________________________________________________________
void AliLHCClockPhase::AddPhaseB2DP(UInt_t timestamp, Float_t phase)
{
  // Add a phase beam2 measurement
  // to the array of data-points

  fPhaseB2.AddLast(new AliDCSValue(phase,timestamp));
}

//______________________________________________________________________________
const AliDCSValue* AliLHCClockPhase::GetPhaseB1DP(Int_t index) const
{
  // Get the value of the phase
  // The argument is the DP index
  AliDCSValue *value = (AliDCSValue*)fPhaseB1.At(index);
  if (!value) {
    AliFatal(Form("Invalid index of the beam1 data point: %d (0 -> %d)",
		  index,fPhaseB1.GetEntries()));
    return NULL;
  }
  return value;
}

//______________________________________________________________________________
const AliDCSValue* AliLHCClockPhase::GetPhaseB2DP(Int_t index) const
{
  // Get the value of the phase
  // The argument is the DP index
  AliDCSValue *value = (AliDCSValue*)fPhaseB2.At(index);
  if (!value) {
    AliFatal(Form("Invalid index of the beam2 data point: %d (0 -> %d)",
		  index,fPhaseB2.GetEntries()));
    return NULL;
  }
  return value;
}
