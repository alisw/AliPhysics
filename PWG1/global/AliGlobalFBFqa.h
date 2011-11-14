/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

//////////////////////////////////////////////////////
// Basic QA task to monitor the observables related to
// flow and balance function analysis. It is intended
// for both pp and PbPb data.
//
// Mainteiners:
//   Carlos Perez (cperez@cern.ch)
//   Alis Rodriguez (alisrm@nikhef.nl)
//////////////////////////////////////////////////////
#ifndef AliGlobalFBFqa_H
#define AliGlobalFBFqa_H

#include "AliAnalysisTaskSE.h"

class TList;
class TH1D;
class TH2D;
class TH3D;
class AliESDtrackCuts;

class AliGlobalFBFqa : public AliAnalysisTaskSE {
  private:
    AliGlobalFBFqa(const AliGlobalFBFqa& analysisTask);
    AliGlobalFBFqa& operator=(const AliGlobalFBFqa& analysisTask);

    AliESDtrackCuts* CutsITSSAGeneric();
    AliESDtrackCuts* CutsTPC50Generic();

    Bool_t fDebugger;
    TList *fOutputList;
    TH2D *fEvents;
    TH2D *fPhiRinVZERO, *fPhiEtaTPC50, *fPhiEtaITSSA;
    TH2D *fPsiCenVZERO, *fPsiCenTPC50, *fPsiCenITSSA;

  public:
    AliGlobalFBFqa();
    AliGlobalFBFqa(const char *name);
    void SetDebugON() {fDebugger = true;}
    virtual ~AliGlobalFBFqa();
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *);

  ClassDef(AliGlobalFBFqa, 1);
};

#endif
