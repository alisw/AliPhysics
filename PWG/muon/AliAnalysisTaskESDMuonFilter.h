#ifndef ALIANALYSISTASKESDMUONFILTER_H
#define ALIANALYSISTASKESDMUONFILTER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//
// Add the muon tracks to the generic AOD track branch during the 
// filtering of the ESD. 
//
// Authors: R. Arnaldi 5/5/08 and L. Aphecetche January 2011

#ifndef ALIANALYSISTASKSE_H
#  include "AliAnalysisTaskSE.h"
#endif
#ifndef ALIANALYSISCUTS_H
#  include "AliAnalysisCuts.h"
#endif

class AliAnalysisFilter;

class AliAnalysisTaskESDMuonFilter : public AliAnalysisTaskSE
{
 public:
    AliAnalysisTaskESDMuonFilter(Bool_t onlyMuon=kTRUE, Bool_t keepAllEvents=kTRUE, Int_t mcMode=0);
    AliAnalysisTaskESDMuonFilter(const char* name, Bool_t onlyMuon=kTRUE, Bool_t keepAllEvents=kTRUE, Int_t mcMode=0);
    virtual ~AliAnalysisTaskESDMuonFilter() {;}
    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() {Init();}
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);

    virtual void ConvertESDtoAOD();

    // Setters
    virtual void SetTrackFilter(AliAnalysisFilter* trackF) {fTrackFilter = trackF;}
    void SetWriteMuonAOD(Bool_t enableMuonAOD){fEnableMuonAOD = enableMuonAOD;}
    void SetWriteDimuonAOD(Bool_t enableDimuonAOD){fEnableDimuonAOD = enableDimuonAOD;}

  void PrintTask(Option_t *option="", Int_t indent=0) const;

 private:
    AliAnalysisTaskESDMuonFilter(const AliAnalysisTaskESDMuonFilter&);
    AliAnalysisTaskESDMuonFilter& operator=(const AliAnalysisTaskESDMuonFilter&);
  void AddFilteredAOD(const char* aodfilename, const char* title);

  AliAnalysisFilter* fTrackFilter; //  Track Filter
  Bool_t fEnableMuonAOD; // flag for enabling Muon AOD production
  Bool_t fEnableDimuonAOD; // flag for enabling Dimuon AOD production
  Bool_t fOnlyMuon; // flag for disabling branches irrelevant for (most) muon analyses
  Bool_t fKeepAllEvents; // keep even events where there's no muons (to get e.g. unbiased vertex distribution)
  Int_t  fMCMode; // whether and how we're filtering MC data
  
  ClassDef(AliAnalysisTaskESDMuonFilter, 5); // Analysis task for standard ESD filtering
};
 
class AliAnalysisNonMuonTrackCuts : public AliAnalysisCuts
{
public:
  AliAnalysisNonMuonTrackCuts();
  virtual ~AliAnalysisNonMuonTrackCuts() {}
  virtual Bool_t IsSelected(TObject* obj);
  virtual Bool_t IsSelected(TList*   /* list */ ) { return kTRUE; }

  ClassDef(AliAnalysisNonMuonTrackCuts,1); // Select muon spectrometer tracks
};

class AliAnalysisNonPrimaryVertices : public AliAnalysisCuts
{
public:
  AliAnalysisNonPrimaryVertices();
  virtual ~AliAnalysisNonPrimaryVertices() {}
  virtual Bool_t IsSelected(TObject* obj);
  virtual Bool_t IsSelected(TList*   /* list */ ) { return kTRUE; }
  
  ClassDef(AliAnalysisNonPrimaryVertices,1); // Select primary vertices
};

#endif
