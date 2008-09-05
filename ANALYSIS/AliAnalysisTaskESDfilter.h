#ifndef ALIANALYSISTASKESDFILTER_H
#define ALIANALYSISTASKESDFILTER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliAnalysisTaskESDfilter.h 24429 2008-03-12 10:27:50Z jgrosseo $ */

#include <TList.h> 
#include <TF1.h> 
#include "AliAnalysisTaskSE.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODPid.h"
class AliAnalysisFilter;
class AliRunTag;
class AliAODTagCreator;


class AliAnalysisTaskESDfilter : public AliAnalysisTaskSE
{
 public:
    AliAnalysisTaskESDfilter();
    AliAnalysisTaskESDfilter(const char* name, Bool_t tags = kFALSE);
    virtual ~AliAnalysisTaskESDfilter() {;}
    // Implementation of interface methods
    virtual void   UserCreateOutputObjects();
    virtual void   Init();
    virtual void   LocalInit() {Init();}
    virtual void   UserExec(Option_t *option);
    virtual Bool_t Notify();
    virtual void   Terminate(Option_t *option);
    virtual void   FinishTaskOutput();

    virtual void ConvertESDtoAOD();
    virtual void CreateTags();
    virtual void SetCreateTags() {fCreateTags = kTRUE;}
    // Setters
    virtual void SetTrackFilter(AliAnalysisFilter* trackF) {fTrackFilter = trackF;}
    virtual void SetKinkFilter (AliAnalysisFilter*  KinkF) {fKinkFilter  =  KinkF;}
    virtual void SetV0Filter   (AliAnalysisFilter*    V0F) {fV0Filter    =    V0F;}
    virtual void SetPthreshold (Double_t p)                {fHighPthreshold =  p;}
    virtual void SetPshape     (TF1 *func)                 {fPtshape        = func;}

    virtual void SetAODPID(AliESDtrack *esdtrack, AliAODTrack *aodtrack, AliAODPid *detpid, Double_t timezero);
    void SetDetectorRawSignals(AliAODPid *aodpid, AliESDtrack *track, Double_t timezero);

 private:
    AliAnalysisTaskESDfilter(const AliAnalysisTaskESDfilter&);
    AliAnalysisTaskESDfilter& operator=(const AliAnalysisTaskESDfilter&);
    // Filtering
    AliAnalysisFilter* fTrackFilter; //  Track Filter
    AliAnalysisFilter* fKinkFilter;  //  Kink  Filter
    AliAnalysisFilter* fV0Filter;    //  V0    Filter
    // PID
    Double_t     fHighPthreshold;    //  Pt threshold for detector signal setting
    TF1 *        fPtshape;           //  Pt spectrum distribution
    // Tags
    Bool_t                   fCreateTags;             //  Flag for tag creation
    Bool_t                   fFirstFile;              //! To flag the first file   
    AliRunTag               *fRunTag;                 //! Pointer to run tag
    TTree                   *fTreeT;                  //! tree for  aod tags
    AliAODTagCreator        *fTagCreator;             //! The tag creator
    
    ClassDef(AliAnalysisTaskESDfilter, 1); // Analysis task for standard ESD filtering
};
 
#endif
