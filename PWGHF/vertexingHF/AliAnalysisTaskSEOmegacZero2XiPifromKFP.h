#ifndef ALIANALYSISTASKSEOMEGACZERO2XIPIFROMKFP_H
#define ALIANALYSISTASKSEOMEGACZERO2XIPIFROMKFP_H

/// \class AliAnalysisTaskSEOmegacZero2XiPifromKFP
/// \brief This is a brief description of my class
///
/// This is a longer description of the class. This longer description is
/// formatted using the Markdown syntax (see below) and can span on multiple
/// lines.
///
/// \author Federica Zanone <federica.zanone@cern.ch>, Ruprecht-Karls-Universität Heidelberg
/// \date March 16, 2022

/* $Id$ */

#ifndef HomogeneousField
#define HomogeneousField
#endif

#include "AliAnalysisTaskSE.h"
#include "AliAODMCParticle.h"
#include "AliRDHFCutsKFP.h"
#include "AliNormalizationCounter.h"
#include "THnSparse.h"
#include "AliPIDResponse.h"
#include "AliAODInputHandler.h"
#include "AliVertexingHFUtils.h"

// includes added to play with KFParticle
#include <vector>
#include "KFParticleBase.h"
#include "KFParticle.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFVertex.h"

class AliAnalysisTaskSEOmegacZero2XiPifromKFP : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskSEOmegacZero2XiPifromKFP();
                                AliAnalysisTaskSEOmegacZero2XiPifromKFP(const char *name, AliRDHFCutsKFP* cuts);
        virtual                 ~AliAnalysisTaskSEOmegacZero2XiPifromKFP();

        virtual void            UserCreateOutputObjects();
        virtual void            Init();
        virtual void            LocalInit() {Init();}
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

        void                   MakeAnaOmegacZero(AliAODEvent *AODEvent, KFParticle PV);
        void                   FillTreeOmegac0(Int_t flagUSorLS, KFParticle kfpOmegac0, AliAODTrack *trackPiFromOmegac0, KFParticle kfpBP, KFParticle kfpXiMinus, KFParticle kfpXiMinus_m, KFParticle kfpPionFromXi, AliAODTrack *trackPionFromXi, AliAODcascade *casc, KFParticle kfpLambda, KFParticle kfpLambda_m, AliAODTrack *trkProton, AliAODTrack *trkPion, KFParticle PV);

    private:
        void                    DefineTreeOmegac0();

        Bool_t                   fIsMC;               // flag of MC analysis
        AliPIDResponse*          fPID;                // PID info
        AliRDHFCutsKFP*          fAnaCuts;            // cuts described in CutObject
        AliAODVertex*            fpVtx;               // primary vertex
        Double_t                 fBzkG;               // magnetic field value [kG]
        TList*                   fOutputList;         // output list (contains only one histo concerning number of events with different requirements)
        TList*                   fListCuts;           // list of applied cuts (contained in cutobject), included in the output
        TTree*                   fTree_Omegac0;       // tree of the candidate variables
        Float_t*                 fVar_Omegac0;        // variables of Omegac0 to be written to the tree
        AliNormalizationCounter* fCounter;            // Counter for normalization
        TH1F*                    fHistEvents;         // Histogram of selected events (added to output list)

        AliAnalysisTaskSEOmegacZero2XiPifromKFP(const AliAnalysisTaskSEOmegacZero2XiPifromKFP &source); // not implemented
        AliAnalysisTaskSEOmegacZero2XiPifromKFP& operator=(const AliAnalysisTaskSEOmegacZero2XiPifromKFP& source); // not implemented

        ClassDef(AliAnalysisTaskSEOmegacZero2XiPifromKFP, 1);
};

#endif
