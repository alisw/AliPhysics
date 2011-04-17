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
class TText;
class AliCentrality;
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskHLT : public AliAnalysisTaskSE {
 
  public: 
    AliAnalysisTaskHLT();
    AliAnalysisTaskHLT(const char *name);
    virtual ~AliAnalysisTaskHLT() {}

    virtual void  UserCreateOutputObjects();
    virtual void  UserExec(Option_t *option);
    virtual void  Terminate(Option_t *);

    //Use only triggered events
    void SetUseHLTTriggerDecision(Bool_t useHLT = kFALSE) { fUseHLTTrigger = useHLT; }
    // function to set the beam type
    void SetBeamType(TString beamType) { fBeamType = beamType; }    
 
private:

    /** copy constructor */
    AliAnalysisTaskHLT(const AliAnalysisTaskHLT&); 
    /** assignment operator */
    AliAnalysisTaskHLT& operator=(const AliAnalysisTaskHLT&); 

    //void SetupESDtrackCuts();  

    Bool_t fUseHLTTrigger;             // use HLT Trigger Decision

    //AliESDtrackCuts *fESDOfftrackCuts;   //! ESD cuts  
    //AliESDtrackCuts *fESDHLTtrackCuts;   //! ESD cuts - on HLT object 
    
    TList *fOutputList; // list of output histograms

    TH1F *fHistTrigger; //! trigger counters 

    TH1F  *fChargeOff;         //! Charge distribution      
    TH1F  *fMomentumOff;       //! momentum	
    TH1F  *fDCArOff;           //! track DCAr to beam line	 
    TH1F  *fDCAzOff;           //! track DCAz to beam line	 
    TH1F  *fNclusterOff;       //! #clusters per track
    TH1F  *fNITSclusterOff;    //! # ITS clusters per track
    TH1F  *fNclusterOffwCut;   //! #clusters per track w cuts
    TH1F  *fPhiOff;            //! azimuthal angle distribution
    TH1F  *fMultOff;           //! track multiplicity of the event
    TH2F  *fXYvertexOff;       //! XY primary vertex distribution
    TH1F  *fXvertexOff;        //! X primary vertex distribution
    TH1F  *fYvertexOff;        //! Y primary vertex distribution
    TH1F  *fZvertexOff;        //! Z primary vertex distribution
    TH1F  *fSPDXvertexOff;     //! SPD X primary vertex distribution
    TH1F  *fSPDYvertexOff;     //! SPD Y primary vertex distribution
    TH1F  *fSPDZvertexOff;     //! SPD Z primary vertex distribution
    TH1F  *fEtaOff;            //! pseudorapidity
    TH1F  *fEtaMomentumcutOff;      //! pseudorapidity with DCA cut
    TH2F  *fNclusVSphiOff;     //! clusters per track vs. azimuthal angle 
    TH2F  *fNclusVSthetaOff;   //! clusters per track vs. polar angle 
    TH1F  *fEventSpecieOff;    //! Event Specie Offline
    TH1F  *fV0cent;            //! V0 centrality information
    TH1F  *fNcontOff;          //! # of contributors to the vertex estimate
    
    TH1F  *fChargeHLT;         //! Charge distribution 
    TH1F  *fMomentumHLT;       //! momentum	
    TH1F  *fDCArHLT;	       //! track DCAr to beam line	 
    TH1F  *fDCAzHLT;	       //! track DCAz to beam line	 
    TH1F  *fNclusterHLT;       //! #clusters per track
    TH1F  *fNITSclusterHLT;    //! # ITS clusters per track
    TH1F  *fNclusterHLTwCut;   //! #clusters per track with cuts
    TH1F  *fPhiHLT;	       //! azimuthal angle distribution
    TH1F  *fMultHLT;	       //! track multiplicity of the event   
    TH2F  *fXYvertexHLT;       //! XY primary vertex distribution
    TH1F  *fXvertexHLT;        //! X primary vertex distribution
    TH1F  *fYvertexHLT;        //! Y primary vertex distribution
    TH1F  *fZvertexHLT;        //! Z primary vertex distribution
    TH1F  *fSPDXvertexHLT;     //! SPD X primary vertex distribution
    TH1F  *fSPDYvertexHLT;     //! SPD Y primary vertex distribution
    TH1F  *fSPDZvertexHLT;     //! SPD Z primary vertex distribution
    TH1F  *fEtaHLT;	       //! pseudorapidity
    TH1F  *fEtaMomentumcutHLT;      //! pseudorapidity with DCA cut
    TH2F  *fNclusVSphiHLT;     //! clusters per track vs. azimuthal angle 
    TH2F  *fNclusVSthetaHLT;   //! clusters per track vs. polar angle 
    TH1F  *fEventSpecieHLT;    //! Event Specie HLT
    TH1F  *fNcontHLT;          //! # of contributors to the vertex estimate
        
    TString fBeamType;         //! beam type: p-p, Pb-Pb, No beam
    TText *fTextBox;           //! TText box containing run number info and date
    Bool_t fSwitch;            //! boolean used to execute parts of the code in the UserExec only once, although
                               // the function is called once per event
    AliCentrality *fCentrality;  //! Centrality holder
   
    ClassDef(AliAnalysisTaskHLT, 0);
};

#endif
