/************************************************************************************
 * Copyright (C) 2018, Copyright Holders of the ALICE Collaboration                 *
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
#ifndef ALIEMCALCALIBBADCHANNELS_H
#define ALIEMCALCALIBBADCHANNELS_H

#include <TObject.h>

class TH2;
class TObjArray;

/**
 * @class AliEMCALCalibBadChannels
 * @brief EMCAL bad channel map in the format of the OADB container however as OCDB object, for the HLT
 * @ingroup EMCALbase
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since Oct 16, 2018
 */
class AliEMCALCalibBadChannels : public TObject {
public:
    /**
     * @brief Dummy constructor
     */
    AliEMCALCalibBadChannels();
    /**
     * @brief Destructor, removing also histograms
     */
    virtual ~AliEMCALCalibBadChannels();

    /**
     * @brief Get the bad channel map for a given supermodule
     * @param[in] supermoduleID ID of the supermodule (0-19)
     * @return Histogram with bad and warm cells for the supermodule (NULL if outside range or object not initialized)
     */
    TH2 *GetBadChannelHistForSupermodule(Int_t supermoduleID);
    /**
     * @brief Set the histogram bad channel histogram for a certain supermodule
     * @param[in] supermoduleID ID of the supermodule (0-19)
     * @param[in] badchannels Calibration histogram for the given supermodule
     */
    void SetBadChannelHistForSupermodule(Int_t supermodueID, TH2 *badchannels);

private:
    AliEMCALCalibBadChannels(const AliEMCALCalibBadChannels &);
    AliEMCALCalibBadChannels &operator=(const AliEMCALCalibBadChannels &);

    TObjArray           *fHistsCalibSupermodules;                       ///< Calibration histograms for supermodules    

    ClassDef(AliEMCALCalibBadChannels, 1);
};

#endif