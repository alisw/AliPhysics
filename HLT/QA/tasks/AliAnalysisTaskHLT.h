// $Id$
//-*- Mode: C++ -*-
//* This file is property of and copyright by the ALICE HLT Project *
//* ALICE Experiment at CERN, All rights reserved.                  *
//* See cxx source for full Copyright notice                        *

#ifndef ALIANALYSISTASKHLT_H
#define ALIANALYSISTASKHLT_H


/** @file AliAnalysisTaskHLT.h
    @author Kalliopi Kanaki, Hege Erdal
    @date   
    @brief An analysis task to compare the offline and HLT esd trees
*/

// forward declarations
class TH1F;
class TH2F;
class TList;
//class AliESDEvent;
//class AliESDtrack;
//class AliESDRun;
class TObjArray;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskHLT : public AliAnalysisTaskSE {
 
  public: 


  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

    AliAnalysisTaskHLT();
    AliAnalysisTaskHLT(const char *name);
    virtual ~AliAnalysisTaskHLT() {}

  /*
   * ---------------------------------------------------------------------------------
   *                                    Methods
   * ---------------------------------------------------------------------------------
   */

    virtual void  UserCreateOutputObjects();
    virtual void  UserExec(Option_t *option);
    virtual void  Terminate(Option_t *);
    virtual void  NotifyRun();

    //Use only triggered events
    void SetUseHLTTriggerDecision(Bool_t useHLT = kFALSE) {fUseHLTTrigger = useHLT;}
 

private:

    /** copy constructor */
    AliAnalysisTaskHLT(const AliAnalysisTaskHLT&); 
    /** assignment operator */
    AliAnalysisTaskHLT& operator=(const AliAnalysisTaskHLT&); 


    /*
     * ---------------------------------------------------------------------------------
     *                            Setup Methods - private
     * ---------------------------------------------------------------------------------
     */
    
    void SetupESDtrackCuts();  


 
  
    /*
     * ---------------------------------------------------------------------------------
     *                             Members - private
     * ---------------------------------------------------------------------------------
     */


    Bool_t fUseHLTTrigger;             // use HLT Trigger Decision

 
    //----------------------------------------------------------------------------------



    AliESDtrackCuts *fESDOfftrackCuts;   //! ESD cuts  
    AliESDtrackCuts *fESDHLTtrackCuts;   //! ESD cuts - on HLT object 
    
    TList *fOutputList; // list of output histograms

    TH1F *fHistTrigger, *fHistHLTTrigger; //! trigger counters 

    TH1F  *fChargeOff;         //! Charge distribution      
    TH1F  *fMomentumOff;       //! momentum	
    TH1F  *fMomentumOffTpc;         //! pseudorapidity for kTPCin
    TH1F  *fMomentumOffTpcIts;      //! pseudorapidity for kTPCin && kITSin
    TH1F  *fDCArOff;           //! track DCAr to beam line	 
    TH1F  *fDCAzOff;           //! track DCAz to beam line	 
    TH1F  *fNclusterOff;       //! #clusters per track
    TH1F  *fNclusterOffwCut;   //! #clusters per track w cuts
    TH1F  *fdEdxOff;           //! TPC signal (offline)
    TH2F  *fdEdxVSPOff;        //! dE/dx vs. momentum 
    TH1F  *fPhiOff;            //! azimuthal angle distribution
    TH1F  *fThetaOff;          //! polar angle distribution
    TH1F  *fMultOff;           //! track multiplicity of the event
    TH2F  *fXYvertexOff;       //! XY primary vertex distribution
    TH1F  *fXvertexOff;        //! X primary vertex distribution
    TH1F  *fYvertexOff;        //! Y primary vertex distribution
    TH1F  *fZvertexOff;        //! Z primary vertex distribution
    TH1F  *fEtaOff;            //! pseudorapidity
    TH1F  *fEtaMomentumcutOff;      //! pseudorapidity with DCA cut
    TH2F  *fNclusVSphiOff;     //! clusters per track vs. azimuthal angle 
    TH2F  *fNclusVSthetaOff;   //! clusters per track vs. polar angle 
    TH1F  *fStatusOff;         //! Status counters 
    TH1F  *fEventSpecieOff;    //! Event Specie Offline
    
    TH1F  *fChargeHLT;         //! Charge distribution 
    TH1F  *fMomentumHLT;       //! momentum	
    TH1F  *fMomentumHLTTpc;    //! pseudorapidity for kTPCin
    TH1F  *fMomentumHLTTpcIts; //! pseudorapidity for kTPCin && kITSin 
    TH1F  *fDCArHLT;	       //! track DCAr to beam line	 
    TH1F  *fDCAzHLT;	       //! track DCAz to beam line	 
    TH1F  *fDCArHLTSG;	       //! track DCAr to beam line as calculated in the HLT reco	 
    TH1F  *fDCAzHLTSG;	       //! track DCAz to beam line as calculated in the HLT reco	 
    TH1F  *fNclusterHLT;       //! #clusters per track
    TH1F  *fNclusterHLTwCut;   //! #clusters per track with cuts
    TH1F  *fdEdxHLT;	       //! TPC signal (offline)
    TH2F  *fdEdxVSPHLT;        //! dE/dx vs. momentum 
    TH1F  *fPhiHLT;	       //! azimuthal angle distribution
    TH1F  *fThetaHLT;          //! polar angle distribution
    TH1F  *fMultHLT;	       //! track multiplicity of the event   
    TH2F  *fXYvertexHLT;       //! XY primary vertex distribution
    TH1F  *fXvertexHLT;        //! X primary vertex distribution
    TH1F  *fYvertexHLT;        //! Y primary vertex distribution
    TH1F  *fZvertexHLT;        //! Z primary vertex distribution
    TH1F  *fEtaHLT;	       //! pseudorapidity
    TH1F  *fEtaMomentumcutHLT;      //! pseudorapidity with DCA cut
    TH2F  *fNclusVSphiHLT;     //! clusters per track vs. azimuthal angle 
    TH2F  *fNclusVSthetaHLT;   //! clusters per track vs. polar angle 
    TH1F  *fStatusHLT;         //! Status counters 
    TH1F  *fEventSpecieHLT;    //! Event Specie HLT
    
    TObjArray *fTrgClsArray; //! array of trigger classes
    
    ClassDef(AliAnalysisTaskHLT, 0);
};

#endif
