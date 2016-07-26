#ifndef ALIANALYSISFILTERTREE_H
#define ALIANALYSISFILTERTREE_H

#include <Rtypes.h>
#include <TAxis.h>
#include <THnSparse.h>

#include "AliAnalysisTaskSE.h"
#include "AliESDtrackCuts.h"

class AliPIDResponse;
class TList;
class TH1F;

class AliAnalysisTaskAbsorptionStudies : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskAbsorptionStudies(const char *name = "FilterTree");
    virtual ~AliAnalysisTaskAbsorptionStudies();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(const Option_t*) {}

    bool            mTOFtimeAnalysis;
    bool            mMCtrue;
    int             mMagField;  // positive: positive B, negative: negative B,  0: both B.
    AliESDtrackCuts mCuts;
    TAxis           mCentrality;
    TAxis           mMassAxis;
		TAxis           mPtAxis;
		TAxis           mPhiAxis;
		TAxis           mTimeAxis;
		TAxis           mDCAxyAxis;
  private:
    // Private methods
    AliAnalysisTaskAbsorptionStudies(const AliAnalysisTaskAbsorptionStudies&);            //! Not implemented
    AliAnalysisTaskAbsorptionStudies& operator=(const AliAnalysisTaskAbsorptionStudies&); //! Not implemented

    TList*   mOutput; //!

    AliPIDResponse* mPIDresponse; //!

    TH1F      * mEventCounter;
    THnSparseF* mPhiPtMass;

    THnSparseF* mPhiPtTime;
    THnSparseF* mPhiPtTimeParticles;

    THnSparseF* mPhiPtDCAxy;
    THnSparseF* mPhiPtDCAxyWeak;
    THnSparseF* mPhiPtDCAxyMaterial;

    THnSparseF* mPhiPtGen;
    THnSparseF* mPhiPtRec;

    int GetParticleId(int pdg);

    ClassDef(AliAnalysisTaskAbsorptionStudies,1)
};

#endif

