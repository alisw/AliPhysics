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
#include "AliPIDResponse.h"

typedef unsigned int uint;

class AliJPartLifetime : public AliAnalysisTaskSE {
    public:
	    enum {kK0s, kLamb, kALamb, kAll};
        enum {kPion, kKaon, kProton, kElectron, kUnknown};

	    AliJPartLifetime();
	    AliJPartLifetime(const char *name,  TString inputformat);
	    AliJPartLifetime(const AliJPartLifetime& ap);
	    AliJPartLifetime& operator = (const AliJPartLifetime& ap);
	    virtual ~AliJPartLifetime();


	    virtual void    UserCreateOutputObjects();
	    virtual void    UserExec(Option_t *);
	    virtual void    FinishTaskOutput();
	    virtual void    Terminate(Option_t *);
	    void            StrangenessMeasure();
        Int_t GetPID(AliPIDResponse *pid, const AliVTrack *trk);

	    void SetCard(AliJCard *pcard) {fcard = pcard;}
	    void SetOption(char * option) {fOption = option;}

    private:

	    TString                         fOption;
	    TList*                          fOutput; //!
	    AliJCard*                       fcard;

	    AliESDtrackCuts*                fTrackCuts; //!
	    AliESDEvent*                    fEsd; //!
        AliPIDResponse                 *fPIDResponse; //!
	    TH1D*                           fEventNumbers; //!
	    TH1D*                           fV0Mass[kAll][kPtDim]; //!
	    TH1D*                           fhTriggPtBin[kAll][kPtDim]; //!
	    TH1D*                           fDecayLength[kAll][kPtDim]; //!
	    TH1D*                           fProperTime[kAll][kPtDim]; //!
        //Particle-identified histograms
        TH1D*                           fV0MassPID[kAll][kPtDim]; //!
        TH1D*                           fProperTimePID[kAll][kPtDim]; //!

	    TH1D*                           fpTspectra[kAll]; //!

	    ClassDef(AliJPartLifetime, 1);
};
#endif // AliJPartLifetime_H
