#ifndef ALIANALYSISTASKQGSEP_H
#define ALIANALYSISTASKQGSEP_H

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

class TList;
class TProfile;
class AliAODEvent;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskQGSep : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskQGSep(const char *name="<default name>");
  virtual ~AliAnalysisTaskQGSep() {}
 
  virtual void   UserCreateOutputObjects(); 
  virtual Bool_t   Notify();
  virtual void   UserExec(Option_t* option);
  virtual void   Terminate(Option_t *);

  void LoopAOD(); //AOD loop
  void LoopAODMC(); //loop containing MC information
  
  void UseMC(Bool_t useMC=kFALSE) { fUseMC = useMC;} //sets use of MC
  void UseAOD(Bool_t useAOD=kFALSE) {fUseAOD = useAOD;} //sets use of AOD inoput
  
 private:
  TString       fBranchRec;  // AOD branch name for reconstructe
  Bool_t       fUseMC; //switch to use MC info
  Bool_t        fUseAOD; //swicth between using AOD input
  Double_t     fXsection; // cross-section from pyxsec.root
  Double_t     fWeight; //fXsection/fAvgTrials; weighting factor for different pT hard bins
  AliAODEvent *fMyAODEvent; // aod event

  TList	      *fOutputList; // output list
  TProfile        *fpHistPtAvEQ;   //Quark Pt_av vs. Energy
  TProfile        *fpHistPtAvEG;   //Gluon Pt_av vs Energy
  TProfile        *fpHistDrEQ;   //Quark Dr vs Energy
  TProfile        *fpHistDrEG;   //Gluon Dr vs Energy
  TProfile        *fpHistDrE;   //Dr vs E for all jets
  TProfile        *fpHistPtAvE;   //Pt_av vs E for all jets
  TProfile        *fpHistDrE3;   //Dr vs E for multi jets
  TProfile        *fpHistPtAvE3;   //Pt_av vs E for multi jets

  AliAnalysisTaskQGSep(const AliAnalysisTaskQGSep&); // not implemented
  AliAnalysisTaskQGSep& operator=(const AliAnalysisTaskQGSep&); // not implemented
  
  ClassDef(AliAnalysisTaskQGSep, 1); // example of analysis
};

#endif
