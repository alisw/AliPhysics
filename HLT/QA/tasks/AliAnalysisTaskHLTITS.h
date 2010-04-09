// $Id$

#ifndef ALIANALYSISTASKHLTITS_H
#define ALIANALYSISTASKHLTITS_H

//* This file is property of and copyright by the ALICE HLT Project *
//* ALICE Experiment at CERN, All rights reserved.                  *
//* See cxx source for full Copyright notice                        *

/** @file AliAnalysisTaskHLTITS.h
    @author Kalliopi Kanaki
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

class AliAnalysisTaskHLTITS : public AliAnalysisTaskSE {
 
  public:  
    AliAnalysisTaskHLTITS(const char *name);
    virtual ~AliAnalysisTaskHLTITS() {}
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    //virtual Bool_t Notify();
    virtual void NotifyRun();

 private:

    AliESDRun *fESDRun;  //!
    TList *fOutputList;

    TH1F *fHistOnlITSsignal; //!
    TH1F *fHistOfflITSsignal; //!
    TH1F *fHistOfflITSsignalTrig; //!
    TH1F *fHistOfflITSsignalNoTrig; //!

    TH1F *fHistOnlITSncls; //!
    TH1F *fHistOfflITSncls; //!
    TH1F *fHistOfflITSnclsTrig; //!
    TH1F *fHistOfflITSnclsNoTrig; //!

    /** copy constructor */
    AliAnalysisTaskHLTITS(const AliAnalysisTaskHLTITS&); 
    /** assignment operator */
    AliAnalysisTaskHLTITS& operator=(const AliAnalysisTaskHLTITS&); 

    ClassDef(AliAnalysisTaskHLTITS, 0);
};

#endif
