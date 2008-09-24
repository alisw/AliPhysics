#ifndef ALIANALYSISTASKTAGCREATOR_H
#define ALIANALYSISTASKTAGCREATOR_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliAnalysisTaskSE.h"
class AliRunTag;
class AliAODTagCreator;
class TTree;


class AliAnalysisTaskTagCreator : public AliAnalysisTaskSE
{
 public:
    AliAnalysisTaskTagCreator();
    AliAnalysisTaskTagCreator(const char* name);
    virtual ~AliAnalysisTaskTagCreator() {;}
    // Implementation of interface methods
    virtual void   UserCreateOutputObjects();
    virtual void   Init();
    virtual void   LocalInit() {Init();}
    virtual void   UserExec(Option_t *option);
    virtual Bool_t Notify();
    virtual void   Terminate(Option_t *option);
    virtual void   FinishTaskOutput();
 private:
    AliAnalysisTaskTagCreator(const AliAnalysisTaskTagCreator&);
    AliAnalysisTaskTagCreator& operator=(const AliAnalysisTaskTagCreator&);
    Bool_t                   fCreateTags;             //  Flag for tag creation
    Bool_t                   fFirstFile;              //! To flag the first file   
    AliRunTag               *fRunTag;                 //! Pointer to run tag
    TTree                   *fTreeT;                  //! tree for  aod tags
    AliAODTagCreator        *fTagCreator;             //! The tag creator
    
    ClassDef(AliAnalysisTaskTagCreator, 1); // Analysis task for standard ESD filtering
};
 
#endif
