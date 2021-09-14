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
    enum EWeightType_t {
        kCrossSectionWeightType,
        kEventWeightType,
        kNoWeightType 
    };
    AliAnalysisTaskEmcalJetSpectrumSDPart();
    AliAnalysisTaskEmcalJetSpectrumSDPart(const char *name);
    virtual ~AliAnalysisTaskEmcalJetSpectrumSDPart();

    void SetDoSoftDrop(bool doSoftDrop) { fDoSoftDrop = doSoftDrop; }
    void SetSDBeta(double beta) { fBeta = beta; }
    void SetSDZCut(double zcut) { fZcut = zcut; }
    void SetDropMass0Jets(bool doDrop) { fDropMass0Jets = doDrop; }
    void SetSDUseChargedConstituents(Bool_t doUse) { fUseChargedConstituents = doUse; }
    void SetSDUseNeutralConstituents(Bool_t doUse) { fUseNeutralConstituents = doUse; }
    void SetCutHardestPartonPt(double ptmin, double ptmax) { fCutHardPartonPt = true; fMinPtHardParton = ptmin; fMaxPtHardParton = ptmax; }
    void SetUsePtHardPartonInOutlierCut(bool doUse) { fOutlierMode = kOutlierPtParton; }
    void SetUseOutlierPtMaxBin(double ptmaxBin) { fMaxPtHardValBin = ptmaxBin; fOutlierMode = kOutlierPtMax; }
    void SetFillHistosWeighted(EWeightType_t weighttype) { fFillHistosWeighted = weighttype; }

    static AliAnalysisTaskEmcalJetSpectrumSDPart *AddTaskEmcalJetSpectrumSDPart(AliJetContainer::EJetType_t jettype, double R, const char *nameparticles, const char *tag = "");

protected:
    enum OutlierMode_t {
        kOutlierPtHard,
        kOutlierPtParton,
        kOutlierPtMax
    };
    virtual void UserCreateOutputObjects();
    virtual bool Run();
    virtual bool IsEventSelected();
    virtual bool CheckMCOutliers();
    virtual void UserRetrieveEventObjects();

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
    Bool_t                                  fCutHardPartonPt;               ///< Apply cut on the pt of the hardest parton
    EWeightType_t                           fFillHistosWeighted;            ///< Use event weight in order to fill the histograms
    OutlierMode_t                           fOutlierMode;                   ///< Mode to determine the outlier (event pt-hard, hardest part, max of the pt-hard bin)
    Double_t                                fPtHardParton;                  ///< Pt of the hardest
    Int_t                                   fPdgHardParton;                 ///< Pdg of the hardest parton
    Double_t                                fMinPtHardParton;               ///< Min. pt of the hardest parton used to select events
    Double_t                                fMaxPtHardParton;               ///< Max. pt of the hardest parton used to select events
    Double_t                                fMaxPtHardValBin;               ///< Max. value of the pt-hard bin used in event generation

    ClassDef(AliAnalysisTaskEmcalJetSpectrumSDPart, 1);
};

}

}

#endif
