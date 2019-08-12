/************************************************************************************
 * Copyright (C) 2019, Copyright Holders of the ALICE Collaboration                 *
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
#ifndef __ALIANALYSISTASKEMCALHIGHECLUSTERTREE_H__
#define __ALIANALYSISTASKEMCALHIGHECLUSTERTREE_H__

#include "AliAnalysisTaskEmcal.h"
#include "TString.h"

class TNtuple;

namespace PWGJE {

namespace EMCALJetTasks {

class AliAnalysisTaskEmcalHighEClusterTree : public AliAnalysisTaskEmcal {
public:
    AliAnalysisTaskEmcalHighEClusterTree();
    AliAnalysisTaskEmcalHighEClusterTree(const char *name);
    virtual ~AliAnalysisTaskEmcalHighEClusterTree();

    void SetMinE(Double_t minE) { fMinClusterE = minE; }

    static AliAnalysisTaskEmcalHighEClusterTree * AddTaskEmcalHighClusterE(const char *name);

protected:
    virtual void UserCreateOutputObjects();
    virtual bool Run();
    virtual bool IsTriggerSelected();

private:
    TNtuple             *fOutputTree;               //!<! Output n-tuple
    Double_t            fMinClusterE;               ///< Minimum cluster e of the cluster tree
    
    AliAnalysisTaskEmcalHighEClusterTree(const AliAnalysisTaskEmcalHighEClusterTree &);
    AliAnalysisTaskEmcalHighEClusterTree &operator=(const AliAnalysisTaskEmcalHighEClusterTree &);

    ClassDef(AliAnalysisTaskEmcalHighEClusterTree, 1);
};

}

}

#endif