// clang-format off
/**************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                   *
 * All rights reserved.                                                               *
 *                                                                                    *
 * Redistribution and use in source and binary forms, with or without                 *
 * modification, are permitted provided that the following conditions are met:        *
 *     * Redistributions of source code must retain the above copyright               *
 *       notice, this list of conditions and the following disclaimer.                *
 *     * Redistributions in binary form must reproduce the above copyright            *
 *       notice, this list of conditions and the following disclaimer in the          *
 *       documentation and/or other materials provided with the distribution.         *
 *     * Neither the name of the <organization> nor the                               *
 *       names of its contributors may be used to endorse or promote products         *
 *       derived from this software without specific prior written permission.        *
 *                                                                                    *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND    *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED      *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE             *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY                *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES         *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;       *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND        *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT         *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 **************************************************************************************/
// clang-format on
#include <TH2.h>
#include <TObjArray.h>
#include <cstring>

#include "AliHLTCaloDigitDataStruct.h"
#include "AliHLTEMCALDigitsMonitor.h"
#include "AliHLTEMCALGeometry.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliHLTEMCALDigitsMonitor)
/// \endcond CLASSIMP

AliHLTEMCALDigitsMonitor::AliHLTEMCALDigitsMonitor()
  : TObject(), fListOfHistograms(NULL), fGeometry(NULL)
{
  memset(fHIDvsAmp, 0, sizeof(TH2*) * 2);
  memset(fHIDvsTime, 0, sizeof(TH2*) * 2);
}

AliHLTEMCALDigitsMonitor::~AliHLTEMCALDigitsMonitor() { delete fListOfHistograms; }

void AliHLTEMCALDigitsMonitor::Init()
{
  fListOfHistograms = new TObjArray;
  fListOfHistograms->SetName("EMCALDigitMonitorHists");

  TString tag[2] = { "LG", "HG" };
  for (int i = 0; i < 2; i++) {
    fHIDvsAmp[i] = new TH2F(Form("hIDvsAmp%s", tag[i].Data()), Form("Amplitude vs tower ID for %s", tag[i].Data()), 200,
                            0., 20., 20000, -0.5, 19999.5);
    fHIDvsTime[i] = new TH2F(Form("hIDvsTime%s", tag[i].Data()), Form("Time vs. tower ID for %s", tag[i].Data()), 1000.,
                             0., 1000., 20000., -0.5, 19999.5);
    fListOfHistograms->Add(fHIDvsAmp[i]);
    fListOfHistograms->Add(fHIDvsTime[i]);
  }
}

void AliHLTEMCALDigitsMonitor::ProcessDigits(Int_t ndigits, const AliHLTCaloDigitDataStruct* digits)
{
  const int kNsecPerSec = 1e9;
  for (int idig = 0; idig < ndigits; idig++) {
    int gaintype = digits[idig].fHgPresent ? 1 : 0;
    int cellID = fGeometry->GetGeometryPtr()->GetAbsCellIdFromCellIndexes(digits[idig].fModule, digits[idig].fX, digits[idig].fZ);
    AliDebug(1, Form("Digit %d is of type %s\n", idig, gaintype == 1 ? "high gain" : "low gain"));
    AliDebug(1, Form("Digit ID %d in module %d: Absolute cell ID: %d\n", digits[idig].fID, digits[idig].fModule, cellID));
    fHIDvsAmp[gaintype]->Fill(digits[idig].fEnergy, cellID);
    fHIDvsTime[gaintype]->Fill(digits[idig].fTime * kNsecPerSec, cellID);
  }
}