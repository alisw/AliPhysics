// $Id$

#ifndef ALIANALYSISTASKHLTPHOS_H
#define ALIANALYSISTASKHLTPHOS_H

//* This file is property of and copyright by the ALICE HLT Project *
//* ALICE Experiment at CERN, All rights reserved.                  *
//* See cxx source for full Copyright notice                        *

/** @file AliAnalysisTaskHLTTPC.h
    @author Zhongbao Yin, Kalliopi Kanaki
    @date
    @brief An analysis task to compare the offline and HLT esd trees
*/


// forward declarations
class TH1F;
class TH2F;
class TList;
class AliESDEvent;
class AliESDtrack;
class AliESDRun;
class TObjArray;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskHLTPHOS : public AliAnalysisTaskSE {
 
  public:  
    AliAnalysisTaskHLTPHOS(const char *name);
    virtual ~AliAnalysisTaskHLTPHOS() {}
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    //virtual Bool_t Notify();
    virtual void NotifyRun();

 private:

    AliESDRun *fESDRun;  //!
    TList *fOutputList;

    TH2F *fHistOnlTrk2PHOS; //! track to PHOS 2,3,4 modules in (eta, phi)
    TH2F *fHistOfflTrk2PHOS; //! 
    TH2F *fHistOfflTrk2PHOSTrig; //!
    TH2F *fHistOfflTrk2PHOSNoTrig; //!

    Int_t fNevt;
    TObjArray *fTrgClsArray;

    static const Float_t fgkPhiMin[5];
    static const Float_t fgkPhiMax[5];
    static const Float_t fgkEtaMin;
    static const Float_t fgkEtaMax;
    static const Float_t fgkNormX[5];
    static const Float_t fgkNormY[5];
    static const Float_t fgkInitPosX[5];
    static const Float_t fgkInitPosY[5];

    /** copy constructor */
    AliAnalysisTaskHLTPHOS(const AliAnalysisTaskHLTPHOS&); 
    /** assignment operator */
    AliAnalysisTaskHLTPHOS& operator=(const AliAnalysisTaskHLTPHOS&); 

    Bool_t IsInPHOS(Int_t iMod, AliESDtrack * trk, Float_t b, TVector3& v);

    ClassDef(AliAnalysisTaskHLTPHOS, 0);
};

#endif
