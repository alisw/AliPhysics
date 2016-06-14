/* Copyright(c) 1998-kPtDim14, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

 // Short comment describing what this class does needed!


#ifndef AliJPartLifetime_H
#define AliJPartLifetime_H

#include "TString.h"
#include "TH1D.h"
#include "AliTriggerAnalysis.h"
#include "AliAnalysisTaskSE.h"
#include "AliESDtrackCuts.h"
#include "AliJConst.h"
#include "AliJCard.h"

typedef unsigned int uint;

class AliJPartLifetime : public AliAnalysisTaskSE {
    public:
        enum {kK0s, kLamb, kALamb, kAll};

        AliJPartLifetime(AliJCard *);
        ~AliJPartLifetime();


        virtual void    UserCreateOutputObjects();
        virtual void    UserExec(Option_t *);
        virtual void    FinishTaskOutput();
        virtual void    Terminate(Option_t *);
        void            StrangenessMeasure();

        //void SetCard(AliJCard *pcard) {fcard = pcard;}
        void SetOption(char * option) {fOption = option;}

    private:

        TString                         fOption;
        TList*                          fOutput; //!
        AliJCard*                       fcard; //!

        AliESDtrackCuts*                fTrackCuts; //!
        AliESDEvent*                    fEsd; //!
        TH1D*                           fEventNumbers; //!
        TH1D*                           fV0Mass[kAll][kPtDim]; //!
        TH1D*                           fhTriggPtBin[kAll][kPtDim]; //!
        TH1D*                           fDecayLength[kAll][kPtDim]; //!
        TH1D*                           fDCA[kAll][kPtDim]; //!
        TH1D*                           fProperTime[kAll][kPtDim]; //!
        TH1D*                           fpTspectra[kAll]; //!

    ClassDef(AliJPartLifetime, 1);
};

#endif
