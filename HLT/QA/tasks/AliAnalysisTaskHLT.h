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
    AliAnalysisTaskHLT(const char *name, float eta=2, float pt=0, float DCAr=20, float DCAz=20);
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
    //TH1F *fHistTrigger; //! trigger counters 

    TH1F  *fChargeOff;         //! Charge distribution      
    TH1F  *fMomentumOff;       //! momentum	
    TH1F  *fDCArOff;           //! track DCAr to beam line	 
    TH1F  *fDCAzOff;           //! track DCAz to beam line	 
    TH1F  *fNclusterOff;       //! #clusters per track
    TH1F  *fNITSclusterOff;    //! # ITS clusters per track
    TH1F  *fPhiOff;            //! azimuthal angle distribution
    TH1F  *fEtaOff;            //! pseudorapidity  
    TH1F  *fMultOff;           //! track multiplicity of the event
    TH1F  *fXvertexOff;        //! X primary vertex distribution
    TH1F  *fYvertexOff;        //! Y primary vertex distribution
    TH1F  *fZvertexOff;        //! Z primary vertex distribution
    TH1F  *fSPDXvertexOff;     //! SPD X primary vertex distribution
    TH1F  *fSPDYvertexOff;     //! SPD Y primary vertex distribution
    TH1F  *fSPDZvertexOff;     //! SPD Z primary vertex distribution
    TH1F  *fEventSpecieOff;    //! Event Specie Offline
    TH1F  *fV0cent;            //! V0 centrality information
    TH1F  *fNcontOff;          //! # of contributors to the vertex estimate

   //---------------------------------//
    
    TH1F  *fChargeHLT;         //! Charge distribution 
    TH1F  *fMomentumHLT;       //! momentum	
    TH1F  *fDCArHLT;	       //! track DCAr to beam line	 
    TH1F  *fDCAzHLT;	       //! track DCAz to beam line	 
    TH1F  *fNclusterHLT;       //! #clusters per track
    TH1F  *fNITSclusterHLT;    //! # ITS clusters per track
    TH1F  *fPhiHLT;	       //! azimuthal angle distribution
    TH1F  *fEtaHLT;	       //! pseudorapidity
    
    TH1F  *fChargeHLTcut;      //! Charge distribution with cuts for selecting primary tracks
    TH1F  *fMomentumHLTcut;    //! momentum	
    TH1F  *fNclusterHLTcut;    //! #clusters per track
    TH1F  *fNITSclusterHLTcut; //! # ITS clusters per track
    TH1F  *fPhiHLTcut;	       //! azimuthal angle distribution
    TH1F  *fEtaHLTcut;	       //! pseudorapidity
    
    TH1F  *fMultHLT;	       //! track multiplicity of the event   
    TH1F  *fXvertexHLT;        //! X primary vertex distribution
    TH1F  *fYvertexHLT;        //! Y primary vertex distribution
    TH1F  *fZvertexHLT;        //! Z primary vertex distribution
    TH1F  *fSPDXvertexHLT;     //! SPD X primary vertex distribution
    TH1F  *fSPDYvertexHLT;     //! SPD Y primary vertex distribution
    TH1F  *fSPDZvertexHLT;     //! SPD Z primary vertex distribution
    TH1F  *fEventSpecieHLT;    //! Event Specie HLT
    TH1F  *fNcontHLT;          //! # of contributors to the vertex estimate
        
    TString fBeamType;         //! beam type: p-p, Pb-Pb, No beam
    TText *fTextBox;           //! TText box containing run number info and date
    TText *fCuts;              //! TText box containing the cuts
    Bool_t fSwitch;            //! boolean used to execute parts of the code in the UserExec only once, although
                               // the function is called once per event
    AliCentrality *fCentrality;//! Centrality holder
    Float_t fEta;              //! cut value
    Float_t fPt;               //! cut value
    Float_t fDCAr;             //! cut value
    Float_t fDCAz;             //! cut value
   
    ClassDef(AliAnalysisTaskHLT, 0);
};

#endif
