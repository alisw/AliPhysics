// $Id$

#ifndef ALIANALYSISTASKHLT_H
#define ALIANALYSISTASKHLT_H

//* This file is property of and copyright by the ALICE HLT Project *
//* ALICE Experiment at CERN, All rights reserved.                  *
//* See cxx source for full Copyright notice                        *

/** @file AliAnalysisTaskHLT.h
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

class AliAnalysisTaskHLT : public AliAnalysisTaskSE {
 
  public:  
    AliAnalysisTaskHLT(const char *name);
    virtual ~AliAnalysisTaskHLT() {}
    virtual void  UserCreateOutputObjects();
    virtual void  UserExec(Option_t *option);
    virtual void  Terminate(Option_t *);
    virtual void  NotifyRun();

 private:

    /** copy constructor */
    AliAnalysisTaskHLT(const AliAnalysisTaskHLT&); 
    /** assignment operator */
    AliAnalysisTaskHLT& operator=(const AliAnalysisTaskHLT&); 

    AliESDRun *fESDRun;
    TList *fOutputList;

    TH1F *fHistTrigger, *fHistHLTTrigger; //! trigger counters 
       
    TH1F *fMomentum_off;   //! momentum     
    TH1F *fDCA_off;        //! track DCA to beam line	
    TH1F *fNcluster_off;   //! #clusters per track
    TH1F *fdEdx_off;       //! TPC signal (offline)
    TH2F *fdEdx_vs_P_off;  //! dE/dx versus momentum for offline TPC tracks
    TH1F *fPhi_off;        //! azimuthal angle distribution
    TH1F *fTheta_off;      //! polar angle distribution
    TH1F *fMult_off;       //! track multiplicity of the event
    TH2F *fXYvertex_off;   //! XY primary vertex distribution
    TH1F *fXvertex_off;    //! X primary vertex distribution
    TH1F *fYvertex_off;    //! Y primary vertex distribution
    TH1F *fZvertex_off;    //! Z primary vertex distribution
    
    TH1F  *fMomentum_hlt; 
    TH1F  *fDCA_hlt;	  
    TH1F  *fNcluster_hlt; 
    TH1F  *fdEdx_hlt;	  
    TH2F  *fdEdx_vs_P_hlt;
    TH1F  *fPhi_hlt;	  
    TH1F  *fTheta_hlt;    
    TH1F  *fMult_hlt;	  
    TH2F  *fXYvertex_hlt; 
    TH1F  *fXvertex_hlt;  
    TH1F  *fYvertex_hlt;  
    TH1F  *fZvertex_hlt;  
   
//     TH1F *fDCA_off_trig;      //! track DCA to beam line for triggered events
//     TH1F *fNcluster_off_trig; //! #clusters per track for triggered events
//     
//     TH1F *fDCA_hlt_trig;     
//     TH1F *fNcluster_hlt_trig;
   
   
//     TH1F *fHistOfflTrkDCANoTrigNclsCut1; //! with cut on #clusters>=60
//     TH1F *fHistOfflTrkDCANoTrigNclsCut2; //! with cut on #clusters<60
//     
//     TH1F *fHistOfflResPtInv; //! resolution on 1/pt for offline tracks
//     TH1F *fHistOnlResPtInv; //! resoltion on 1/pt for online tracks
// 
//     TH1F *fHistOffldZ;  //! resolution on z 
//     TH1F *fHistOnldZ;   //! resolution on z 
// 
//     TH1F *fHistOffldX; //! resolution on r 
//     TH1F *fHistOnldX;  //! resolution on r 
//     
//     TH1F *fHistOfflPhi;  //! resolution on azimuthal angle 
//     TH1F *fHistOnlPhi;  //! resolution on azimuthal angle 
// 
//     TH1F *fHistOfflTheta; //! resolution on polar angle 
//     TH1F *fHistOnlTheta; //! resolution on polar angle 
// 
//     TH2F *fHistOnlDZ;  //! online trigger tracks distance to beam and Z to IP
//     TH2F *fHistOfflDZ; //! offline tracks distance to beam and Z to IP
//     TH2F *fHistOfflDZTrig; //!
//     TH2F *fHistOfflDZNoTrig; //!

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

    ClassDef(AliAnalysisTaskHLT, 0);
};

#endif
