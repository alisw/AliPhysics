// $Id$
//-*- Mode: C++ -*-
//* This file is property of and copyright by the ALICE HLT Project *
//* ALICE Experiment at CERN, All rights reserved.                  *
//* See cxx source for full Copyright notice                        *

#ifndef ALIANALYSISTASKHLTCENTRALBARREL_H
#define ALIANALYSISTASKHLTCENTRALBARREL_H

/** @file AliAnalysisTaskHLTCentralBarrel.h
    @author Per Ivar Lønne, Hege Erdal, Kalliopi Kanaki
    @date   
    @brief An analysis task to compare the offline and HLT esd trees
*/

// forward declarations
class TList;
class TText;
class TString;
class AliESDEvent;
class AliCentrality;

#include "THnSparse.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskHLTCentralBarrel : public AliAnalysisTaskSE {
 
  public: 

    AliAnalysisTaskHLTCentralBarrel();
    AliAnalysisTaskHLTCentralBarrel(const char *name);
    virtual ~AliAnalysisTaskHLTCentralBarrel();
    
    virtual void  UserCreateOutputObjects();
    virtual void  UserExec(Option_t *option);
    virtual void  Terminate(Option_t *);

    // function to select only HLT triggered events
    //void SetUseHLTTriggerDecision(Bool_t useHLT = kFALSE) { fUseHLTTrigger = useHLT;        }
    // function to set the beam type
    void SetBeamType(TString beamType) {  fBeamType = beamType; }    
    // function to create the THnSparse and name the axis
    THnSparseF* CreateEventTHnSparse(const char* name, Int_t size, const int* bins, double* min, double* max);
    // function to create the THnSparse and name the axis
    THnSparseF* CreateTrackTHnSparse(const char* name, Int_t size, const int* bins, double* min, double* max);
    // options for filling HLT or OFF properties, or event and track properties
    void SetOptions(TString options) { fOptions = options; }
    //function to fill the THnSparse
    //void Fill(AliESDevent *esd, THnSparseF* thn);
    
 private:

    /** copy constructor */
    AliAnalysisTaskHLTCentralBarrel(const AliAnalysisTaskHLTCentralBarrel&); 
    /** assignment operator */
    AliAnalysisTaskHLTCentralBarrel& operator=(const AliAnalysisTaskHLTCentralBarrel&); 
                
    Bool_t fUseHLTTrigger;       // Use HLT Trigger Decision
    AliCentrality *fCentrality;  // Centrality holder
    TString fBeamType;           // beam type: p-p, Pb-Pb, No beam
    
    TList *fOutputList;  // list of output THnSparse objects
    
    THnSparse *fEventOFF; //! offline event properties
    THnSparse *fEventHLT; //! HLT event properties

    THnSparse *fTrackOFF; //! offline track properties
    THnSparse *fTrackHLT; //! HLT track properties
    
    TString fOptions; //! options for filling event and/or track properties for hlt and/or offline
    TText *fTextBox;  //! TText box containing run number info and date
    Bool_t fSwitch;   //! boolean used to execute parts of the code in the UserExec only once, although
                      // the function is called once per event
    
    ClassDef(AliAnalysisTaskHLTCentralBarrel, 0);
};
#endif
