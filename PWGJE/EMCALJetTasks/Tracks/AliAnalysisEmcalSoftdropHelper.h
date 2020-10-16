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
#ifndef ALIANALYSISEMCALSOFTDROPHELPER_H
#define ALIANALYSISEMCALSOFTDROPHELPER_H

#include <vector>
#include <TObject.h>
#include "AliVCluster.h"

class TBinning;
class AliEmcalJet;

namespace PWGJE {

namespace EMCALJetTasks {

class AliAnalysisEmcalSoftdropHelperImpl {
public:
  enum EBinningMode_t {
    kSDModeINT7,
    kSDModeEJ1,
    kSDModeEJ2,
  };
  enum EReclusterizer_t {
    kCAAlgo = 0,
    kKTAlgo = 1,
    kAKTAlgo = 2
  };

  struct SoftdropResults {
    double fZg;
    double fMg;
    double fRg;
    double fPtg;
    double fMug;
    int fNsd;
  };

  struct SoftdropParams {
    EReclusterizer_t fReclusterizer;
    double fBeta;
    double fZcut;
    bool fUseChargedConstituents;
    bool fUseNeutralConstituents;
  };

  AliAnalysisEmcalSoftdropHelperImpl() {}
  virtual ~AliAnalysisEmcalSoftdropHelperImpl() {}

  TBinning *GetDefaultPartLevelPtBinning(EBinningMode_t binmode) const;
  TBinning *GetDefaultDetLevelPtBinning(EBinningMode_t binmode) const;
  TBinning *GetZgBinning(double zcut) const;
  TBinning *GetRgBinning(double R) const;

  SoftdropResults MakeSoftdrop(const AliEmcalJet &jet, double jetradius, bool isPartLevel, SoftdropParams sdparams, AliVCluster::VCluUserDefEnergy_t energydef, double *vertex);

  ClassDef(AliAnalysisEmcalSoftdropHelperImpl, 1);
};

class AliAnalysisEmcalSoftdropHelper : public TObject, public AliAnalysisEmcalSoftdropHelperImpl {
public:
  AliAnalysisEmcalSoftdropHelper() : TObject(), AliAnalysisEmcalSoftdropHelperImpl() {}
  virtual ~AliAnalysisEmcalSoftdropHelper() {}

  ClassDef(AliAnalysisEmcalSoftdropHelper, 1);
};

}

}

#endif
