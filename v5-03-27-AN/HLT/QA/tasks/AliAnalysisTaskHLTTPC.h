// $Id$

#ifndef ALIANALYSISTASKHLTTPC_H
#define ALIANALYSISTASKHLTTPC_H

//* This file is property of and copyright by the ALICE HLT Project *
//* ALICE Experiment at CERN, All rights reserved.                  *
//* See cxx source for full Copyright notice                        *

/** @file AliAnalysisTaskHLTTPC.h
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

class AliAnalysisTaskHLTTPC : public AliAnalysisTaskSE {
 
  public:  
    AliAnalysisTaskHLTTPC(const char *name);
    virtual ~AliAnalysisTaskHLTTPC() {}
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    //virtual Bool_t Notify();
    virtual void NotifyRun();

 private:

    /** copy constructor */
    AliAnalysisTaskHLTTPC(const AliAnalysisTaskHLTTPC&); 
    /** assignment operator */
    AliAnalysisTaskHLTTPC& operator=(const AliAnalysisTaskHLTTPC&); 

    AliESDRun *fESDRun;  //!
    TList *fOutputList;

    TH1F *fHistTrigger; //! trigger counters 
    TH1F *fHistHLTTrigger;  //! HLT online trigger counters

    TH1F *fHistOfflTrkDCA; //! Offline track DCA to beam line
    TH1F *fHistOfflTrkDCATrig; //! Offline track DCA for triggered events
    TH1F *fHistOfflTrkDCANoTrig; //! Offline track DCA for not triggered evts
    
    TH1F *fHistOnlTrkDCA;  //! Online track DCA
    TH1F *fHistOnlTrkDCATrig; //!
    
    TH1F *fHistOfflTrkNcls; //! #clusters for offline tracks
    TH1F *fHistOfflTrkNclsTrig; //! #clusters for offline tracks, triggered evts
    TH1F *fHistOfflTrkNclsNoTrig; //! #clusters not triggered
    
    TH1F *fHistOnlTrkNcls; //! #clusters for online tracks
    TH1F *fHistOnlTrkNclsTrig; //!

    TH1F *fHistOfflTrkDCANoTrigNclsCut1; //! with cut on #clusters>=60
    TH1F *fHistOfflTrkDCANoTrigNclsCut2; //! with cut on #clusters<60
    
    TH1F *fHistOfflTrkP; //! momentum for offline tracks
    TH1F *fHistOfflTrkPTrig; //! momentum for triggered events
    TH1F *fHistOfflTrkPNoTrig; //! momentum for not triggered events
    TH1F *fHistOnlTrkP; //! momentum for online tracks

    TH1F *fHistOfflResPtInv; //! resolution on 1/pt for offline tracks
    TH1F *fHistOnlResPtInv; //! resoltion on 1/pt for online tracks

    TH1F *fHistOffldEdx; //!
    TH1F *fHistOnldEdx; //!

    TH2F *fHistOffldEdxVsP;  //! dE/dx versus momentum for offline TPC tracks
    TH2F *fHistOnldEdxVsP; //! dE/dx versus momentum for online TPC tracks

    TH1F *fHistOffldZ;  //! resolution on z 
    TH1F *fHistOnldZ;   //! resolution on z 

    TH1F *fHistOffldX; //! resolution on r 
    TH1F *fHistOnldX;  //! resolution on r 
    
    TH1F *fHistOfflPhi;  //! resolution on azimuthal angle 
    TH1F *fHistOnlPhi;  //! resolution on azimuthal angle 

    TH1F *fHistOfflTheta; //! resolution on polar angle 
    TH1F *fHistOnlTheta; //! resolution on polar angle 

    TH2F *fHistOnlDZ;  //! online trigger tracks distance to beam and Z to IP
    TH2F *fHistOfflDZ; //! offline tracks distance to beam and Z to IP
    TH2F *fHistOfflDZTrig; //!
    TH2F *fHistOfflDZNoTrig; //!

    Int_t fNevt;
    TObjArray *fTrgClsArray;

    /*
    Int_t fNtracksThruZ0;   //#tracks thru central electrode 
    Int_t fNtracksThruZ0Trig; //#tracks thru central electrode being triggered 
    */

    static const Float_t fgkPhiMin[5];
    static const Float_t fgkPhiMax[5];
    static const Float_t fgkEtaMin;
    static const Float_t fgkEtaMax;
    static const Float_t fgkNormX[5];
    static const Float_t fgkNormY[5];
    static const Float_t fgkInitPosX[5];
    static const Float_t fgkInitPosY[5];

    ClassDef(AliAnalysisTaskHLTTPC, 0); // example of analysis
};

#endif
