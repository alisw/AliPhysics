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
#include <TH2.h>
#include <TObjArray.h>
#include "AliEMCALCalibBadChannels.h"
#include "AliLog.h"

ClassImp(AliEMCALCalibBadChannels);

AliEMCALCalibBadChannels::AliEMCALCalibBadChannels():
    TObject(),
    fHistsCalibSupermodules(NULL)
{

}

AliEMCALCalibBadChannels::~AliEMCALCalibBadChannels(){
    if(fHistsCalibSupermodules) delete fHistsCalibSupermodules;
}

void AliEMCALCalibBadChannels::SetBadChannelHistForSupermodule(Int_t supermoduleID, TH2 *calibhist) {
    if(!fHistsCalibSupermodules) {
        fHistsCalibSupermodules = new TObjArray(20);
        fHistsCalibSupermodules->SetOwner(true);
    }
    if(supermoduleID < 0 || supermoduleID >= 20) {
        AliError(Form("Trying to add histogram for supermodule outside range: %d", supermoduleID));
        return;
    }
    fHistsCalibSupermodules->AddAt(calibhist, supermoduleID);
}

TH2 *AliEMCALCalibBadChannels::GetBadChannelHistForSupermodule(Int_t supermoduleID){
    if(supermoduleID < 0 || supermoduleID >= 20){
        AliError(Form("Trying to access the calibration histogram for a supermodule outside the range: %d", supermoduleID));
        return NULL;
    }
    if(!fHistsCalibSupermodules){
        AliError("No calibration histograms found for the current object");
        return NULL;
    }
    TH2 * histogram = static_cast<TH2 *>(fHistsCalibSupermodules->At(supermoduleID));
    if(histogram) {
        histogram->SetDirectory(NULL);
    }
    return histogram;
}