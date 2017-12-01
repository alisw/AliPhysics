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
#ifndef AliHLTEMCALDIGITSMONITOR_H
#define AliHLTEMCALDIGITSMONITOR_H

#include <TObject.h>

class TH2;
class TObjArray;
class AliHLTCaloDigitDataStruct;
class AliHLTEMCALGeometry;

/**
 * @class   AliHLTEMCALDigitsMonitor
 * @brief   Monitor class for EMCAL on digits level
 * @author  Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since   June 1, 2017
 *
 * This class handles the montioring on digits level for EMCAL in the HLT reconstruction
 * chain. As such, the class monitors:
 *   - Energy vs cell ID
 *   - Time vs cell ID
 * both for high gain and low gain cells.
 */
class AliHLTEMCALDigitsMonitor : public TObject
{
 public:
  /** Constructor */
  AliHLTEMCALDigitsMonitor();

  /** Destructor */
  virtual ~AliHLTEMCALDigitsMonitor();

  /**
   * @brief Initialize histograms
   */
  void Init();

  /**
   * @brief Set the EMCAL geometry
   *
   * EMCAL geometry is mandatory in order to get the absolute ID of a cell
   *
   * @param[in] geo EMCAL Geometry (HLT wrapped version)
   */
  void SetGeometry(const AliHLTEMCALGeometry* geo) { fGeometry = geo; }

  /**
   * @brief Process digits
   *
   * Fill histograms:
   * - Cell amp vs Cell ID
   * - Cell time vs Cell ID
   *
   * @param[in] ndigits Number of digits in block
   * @param[in] digits Block (array) of digits
   */
  void ProcessDigits(Int_t ndigits, const AliHLTCaloDigitDataStruct* digits);

  /**
   * @brief Get a TObjArray with the histograms handled by this component
   * @return List (TObjArray) of histograms handled by this component
   */
  TObjArray* GetListOfHistograms() const { return fListOfHistograms; }

 private:
  TObjArray* fListOfHistograms;         //!<! Container with histograms
  const AliHLTEMCALGeometry* fGeometry; //!<! EMCAL geometry
  TH2* fHIDvsAmp[2];                    //!<! Amp vs tower ID
  TH2* fHIDvsTime[2];                   //!<! Time vs tower ID

  /// \cond CLASSIMP
  ClassDef(AliHLTEMCALDigitsMonitor, 1);
  /// \endcond
};
#endif
