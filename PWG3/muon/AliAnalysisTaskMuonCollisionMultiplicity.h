#ifndef ALIANALYSISTASKMUONCOLLISIONMULTIPLICITY_H
#define ALIANALYSISTASKMUONCOLLISIONMULTIPLICITY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \ingroup base
/// \class AliAnalysisTaskMultiplicity
/// \compute the number of Muons tracks as a function of the SPD tracklets multiplicity
/// Author Matthieu LENHARDT - SUBATECH, Nantes

#include "AliAnalysisTaskSE.h"

class AliAODEvent;
class AliESDEvent;
class AliAODTrack;
class AliESDMuonTrack;
class AliAODDimuon;
class TList;
class TArrayD;

class AliAnalysisTaskMuonCollisionMultiplicity : public AliAnalysisTaskSE
{
 public:
  AliAnalysisTaskMuonCollisionMultiplicity();
  AliAnalysisTaskMuonCollisionMultiplicity(const AliAnalysisTaskMuonCollisionMultiplicity& rhs);
  AliAnalysisTaskMuonCollisionMultiplicity& operator=(const AliAnalysisTaskMuonCollisionMultiplicity&);
  AliAnalysisTaskMuonCollisionMultiplicity(const Char_t* name);
  virtual ~AliAnalysisTaskMuonCollisionMultiplicity();
  
  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  virtual void NotifyRun();
  virtual void FinishTaskOutput();
  
  
  Double_t GetZCut()              {return fZCut;};
  Double_t GetEtaCut()            {return fEtaCut;};
  void SetZCut(Double_t zCut)     {fZCut = zCut;};
  void SetEtaCut(Double_t etaCut) {fEtaCut = etaCut;};
  
 private:
  
  Bool_t fIsInit;              //< Is the class initialized?
  
  AliAODEvent *fAOD;                  //!< AOD Event
  AliESDEvent *fESD;                  //!< ESD Event

  Double_t fZCut;                    //< Cut on the |z| of the primary vertex
  Double_t fEtaCut;                  //< Cut on the eta cut of the SPD tracklets

  Int_t fTrackletMultiplicity;       //< SPD tracklets multiplicity in the current event

  TList *fTriggerList;               //< list of all trigger histos
  TList *fSingleMuonList;            //< List of all single muons histos 
  TList *fDimuonList;                //< List of all dimuons histos
  TList *fMonteCarloList;            //< List of all histos containing MC info


  void Init();
  Bool_t CheckEventAOD();
  Bool_t CheckEventESD();
  void ComputeMultiplicity();
  Bool_t IsUsableMuon(AliAODTrack *track);
  Bool_t IsUsableMuon(AliESDMuonTrack *track);
  void FillHistosAOD(Int_t triggerClass);
  void FillHistosESD(Int_t triggerClass);
  void FillHistosMC();

  ClassDef(AliAnalysisTaskMuonCollisionMultiplicity, 1);
};
    
#endif
    
