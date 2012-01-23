#ifndef ALIANALYSISTASKESDMCLABELADDITION_H
#define ALIANALYSISTASKESDMCLABELADDITION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

#include <TList.h> 
#include "AliAnalysisTaskSE.h"

class AliAnalysisFilter;
class AliStack;
class AliESDMuonTrack;
class AliMUONTrack;
class AliMUONVTrackStore;

class AliAnalysisTaskESDMCLabelAddition : public AliAnalysisTaskSE
{
  
  public:
    AliAnalysisTaskESDMCLabelAddition();
    AliAnalysisTaskESDMCLabelAddition(const char* name);
    virtual ~AliAnalysisTaskESDMCLabelAddition() {;}
    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() {Init();}
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);
    
    
  private:
    
    static const Double_t fgkSigmaCut; // sigma cut applied to match a reconstructed cluster with a trackref
    
    AliAnalysisTaskESDMCLabelAddition(const AliAnalysisTaskESDMCLabelAddition&);
    AliAnalysisTaskESDMCLabelAddition& operator=(const AliAnalysisTaskESDMCLabelAddition&);
    
    void AddMCLabel();
    AliMUONTrack* ESDToMUON(AliESDMuonTrack &esdTrack);
    AliMUONTrack* MatchWithTrackRef(AliESDMuonTrack &esdTrack, AliMUONVTrackStore &trackRefStore);
    Bool_t TrackMatched(AliMUONTrack &track, AliMUONTrack &trackRef);
    
    ClassDef(AliAnalysisTaskESDMCLabelAddition, 1); // Analysis task for standard ESD filtering
    
};

#endif
