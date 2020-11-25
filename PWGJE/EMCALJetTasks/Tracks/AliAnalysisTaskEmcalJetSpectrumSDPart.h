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
#ifndef __ALIANALYSISTASKEMCALJETSPECTRUMSDPART_H__
#define __ALIANALYSISTASKEMCALJETSPECTRUMSDPART_H__

#include "AliAnalysisTaskEmcalJet.h"
#include "AliAnalysisEmcalSoftdropHelper.h"

class THistManager;

namespace PWGJE {

namespace EMCALJetTasks {


class AliAnalysisTaskEmcalJetSpectrumSDPart : public AliAnalysisTaskEmcalJet,
                                              public AliAnalysisEmcalSoftdropHelperImpl
{
public:
    AliAnalysisTaskEmcalJetSpectrumSDPart();
    AliAnalysisTaskEmcalJetSpectrumSDPart(const char *name);
    virtual ~AliAnalysisTaskEmcalJetSpectrumSDPart();

    void SetDoSoftDrop(bool doSoftDrop) { fDoSoftDrop = doSoftDrop; }
    void SetSDBeta(double beta) { fBeta = beta; }
    void SetSDZCut(double zcut) { fZcut = zcut; }
    void SetDropMass0Jets(bool doDrop) { fDropMass0Jets = doDrop; }
    void SetSDUseChargedConstituents(Bool_t doUse) { fUseChargedConstituents = doUse; }
    void SetSDUseNeutralConstituents(Bool_t doUse) { fUseNeutralConstituents = doUse; }

    static AliAnalysisTaskEmcalJetSpectrumSDPart *AddTaskEmcalJetSpectrumSDPart(AliJetContainer::EJetType_t jettype, double R, const char *nameparticles, const char *tag = "");

protected:

    virtual void UserCreateOutputObjects();
    virtual bool Run();
    virtual bool IsEventSelected();
    virtual bool CheckMCOutliers();

private:
    THistManager                            *fHistos;                       //!<! Histogram

    // Softdrop settings
    Bool_t                                  fDoSoftDrop;                    ///< Fill SoftDrop histograms
    Bool_t                                  fDropMass0Jets;                 ///< Drop jets with mass 0
    Double_t                                fBeta;                          ///< SoftDrop Beta
    Double_t                                fZcut;                          ///< SoftDrop Zcut
    Bool_t                                  fUseChargedConstituents;        ///< SoftDrop use charged constituents
    Bool_t                                  fUseNeutralConstituents;        ///< SoftDrop use neutral constituents
    Bool_t                                  fUseStandardOutlierRejection;   ///< Fall back to standard outlier rejection

    ClassDef(AliAnalysisTaskEmcalJetSpectrumSDPart, 1);
};

}

}

#endif
