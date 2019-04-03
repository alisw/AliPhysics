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
#include <cstdlib>
#include <cstring>

#include "AliEmcalTriggerPartChannelMap.h"

ClassImp(PWG::EMCAL::TriggerPart::AliEmcalTriggerPartChannelMap);

using namespace PWG::EMCAL::TriggerPart;

/**
 * Constructor, initializing channel map with the dimensions needed
 * @param ncols Number of columns
 * @param nrows Number of rows
 */
AliEmcalTriggerPartChannelMap::AliEmcalTriggerPartChannelMap(int ncols, int nrows):
    TObject(),
    fNADCCols(ncols),
    fNADCRows(nrows),
    fADC(NULL)
{
  fADC = new double[fNADCCols * fNADCRows];
  memset(fADC, 0, sizeof(double) * fNADCCols * fNADCRows);
}

/**
 * Destructor
 */
AliEmcalTriggerPartChannelMap::~AliEmcalTriggerPartChannelMap() {
  delete[] fADC;
}

/**
 * Set ADC value for position (col, row). Checks for boundary.
 * @param col Column of the position
 * @param row Row of the position
 * @param ADC The value to set
 */
void AliEmcalTriggerPartChannelMap::SetADC(int col, int row, double adc) {
  if(row >= fNADCRows || col >= fNADCCols)
	  throw BoundaryException(row, col, fNADCRows, fNADCCols);
  fADC[GetIndexInArray(col, row)] = adc;
}

/**
 * Add ADC value for position (col, row). Checks for boundary.
 * @param col Column of the position
 * @param row Row of the position
 * @param ADC The value to set
 */
void AliEmcalTriggerPartChannelMap::AddADC(int col, int row, double adc) {
  if(row >= fNADCRows || col >= fNADCCols)
	  throw BoundaryException(row, col, fNADCRows, fNADCCols);
  fADC[GetIndexInArray(col, row)] += adc;
}

/**
 * Set the ADC values stored in the 2D map again to 0
 */
void AliEmcalTriggerPartChannelMap::Reset() {
  memset(fADC, 0, sizeof(double) * fNADCCols * fNADCRows);
}

/**
 * Get ADC value at position (col, row). Checks for boundary.
 * @param col Column of the position
 * @param row Row of the position
 * @return ADC value at the given position
 */
double AliEmcalTriggerPartChannelMap::GetADC(int col, int row) const {
  if(row >= fNADCRows || col >= fNADCCols)
	  throw BoundaryException(row, col, fNADCRows, fNADCCols);
  return fADC[GetIndexInArray(col, row)];
}
