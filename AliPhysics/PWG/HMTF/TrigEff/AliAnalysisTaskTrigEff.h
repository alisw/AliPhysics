/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
/* AliAnalysisTaskTrigEff.h
 *
 */
#ifndef ALIANALYSISTASKTRIGEFF_H
#define ALIANALYSISTASKTRIGEFF_H

class TH1F;
class TList;
class AliESDtrackCuts;
class AliTriggerAnalysis;

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

enum {kHistoINT11, kHistoINT7, kHistoINT7Offline, kHistoINT11EvSel, kHistoINT7EvSel, kHistoINT7OfflineEvSel, kNHist};
enum {kEvAll,kEvINT11,kEvFOOnline, kEvInelGT0, kNEvSel};

class AliAnalysisTaskTrigEff : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskTrigEff();
  AliAnalysisTaskTrigEff(const char *name);
  virtual ~AliAnalysisTaskTrigEff();
  
  virtual void     UserCreateOutputObjects();
  virtual void     UserExec(Option_t *option);
  virtual void     Terminate(Option_t *);
  Bool_t GetIsMC() { return fIsMC; }
  void  SetIsMC (Bool_t var) { fIsMC = var;}
  Int_t GetNtrkCut() { return fNtrkCut; }
  void  SetNtrkCut (Int_t var) { fNtrkCut = var;}
  
    
private:
  TList *fOutput;       // Output list
  TH1F  *fHistV0M      [kNHist]; // V0M centrality distribution
  TH1F  *fHistNtrk70100[kNHist]; // Tracklets distribution in the 70-100 V0 mult
  TH1I  *fHistEvCount;  // event counter
  AliTriggerAnalysis * fTriggerAnalysis; //! Transient
  Int_t fNtrkCut; // Number of tracklets in |eta| < 1 used to define the "offline" trigger

  // NEW HISTO to be declared here
  
  Bool_t fIsMC; // true if processing monte carlo


  
  AliAnalysisTaskTrigEff(const AliAnalysisTaskTrigEff&); // not implemented
  AliAnalysisTaskTrigEff& operator=(const AliAnalysisTaskTrigEff&); // not implemented
    
  ClassDef(AliAnalysisTaskTrigEff, 1); // example of analysis
};

#endif

