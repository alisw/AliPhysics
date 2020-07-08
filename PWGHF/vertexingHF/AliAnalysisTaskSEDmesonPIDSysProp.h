#ifndef ALIANALYSISTASKSEDMESONPIDSYSPROP_H
#define ALIANALYSISTASKSEDMESONPIDSYSPROP_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. */

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// \class AliAnalysisTaskSEDmesonPIDSysProp                                                              //
// \brief analysis task for PID Systematic uncertainty propagation from the single track to the D mesons //
// \author: A. M. Barbano, anastasia.maria.barbano@cern.ch                                               //
// \author: F. Grosa, fabrizio.grosa@cern.ch                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTaskSE.h"
#include "AliRDHFCuts.h"
#include "AliPIDResponse.h"
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
 
class AliAnalysisTaskSEDmesonPIDSysProp : public AliAnalysisTaskSE {
public:
  enum DecChannel {kDplustoKpipi,kD0toKpi,kDstartoKpipi,kDstoKKpi};
  enum PIDstrategy {kConservativePID,kStrongPID,knSigmaPID};
  
  enum KaonTOFhisto {kKaonTPCtag,kSamePionV0tag};
  enum KaonTPChisto {kKaonTOFtag,kKaonKinkstag};

  enum VarForProp {kPt, kP};

  AliAnalysisTaskSEDmesonPIDSysProp();
  AliAnalysisTaskSEDmesonPIDSysProp(int ch, AliRDHFCuts* cuts);
  virtual ~AliAnalysisTaskSEDmesonPIDSysProp();
  
  virtual void   UserCreateOutputObjects();
  virtual void   Init();
  virtual void   UserExec(Option_t *);
  
  void SetAODMismatchProtection(int opt=1) {fAODProtection=opt;}
  void SetPIDStrategy(int PIDst=kStrongPID) {fPIDstrategy=PIDst;}
  void SetKaonHistoOptions(int tpcopt, int tofopt) {fKaonTPCHistoOpt=tpcopt; fKaonTOFHistoOpt=tofopt;}
  void SetVariableForUncProp(int var=kPt) {fVarForProp=var;}

  int GetDecayChannel()const {return fDecayChannel;}

  bool LoadEffSystFile(TString systFileName);

private:
  double GetDmesonPIDuncertainty(AliAODTrack *track[], const int nDau, TClonesArray* arrayMC, double ptD);
  void GetSingleTrackSystAndProb(TH1F* hSingleTrackSyst, TH1F* hSingleTrackEff, int bin, double &syst, double &prob);

  TList *fOutput;                 //!<! tlist with output
  TH1F *fHistNEvents;             //!<! histo with number of events
  TH2F *fHistPtDauVsD;            //!<! histo with pT daughters vs pT candidate
  TH2F *fHistSystPIDEffD;         //!<! histo with PID systematic uncertainty on the D candidate

  TH1F *fHistEffPionTPC[2];       //-> histo for Pion TPC Nsigma syst
  TH1F *fHistEffPionTOF;          //-> histo for Pion TOF Nsigma syst
  TH1F *fHistEffKaonTPC[2];       //-> histo for Kaon TPC Nsigma syst
  TH1F *fHistEffKaonTOF;          //-> histo for Kaon TOF Nsigma syst

  TH1F *fHistSystPionTPC[2];      //-> histo for Pion TPC Nsigma syst 
  TH1F *fHistSystPionTOF;         //-> histo for Pion TOF Nsigma syst
  TH1F *fHistSystKaonTPC[2];      //-> histo for Kaon TPC Nsigma syst
  TH1F *fHistSystKaonTOF;         //-> histo for Kaon TOF Nsigma syst

  TString fPartName;              /// string for particle name
  AliPIDResponse *fPIDresp;       /// basic pid object

  int fPIDstrategy;               /// PID strategy (conservative, strong, nsigma..)
  double fnSigma;                 /// number of sigma PID if nsigma strategy enabled
    
  int fDecayChannel;              ///identify the decay channel
  int fKaonTPCHistoOpt;           ///option for syst on kaon TPC PID efficiency
  int fKaonTOFHistoOpt;           ///option for syst on kaon TOF PID efficiency

  int fAODProtection;             /// flag to activate protection against AOD-dAOD mismatch.
  
  int fNPtBins;                   /// number of pT bins
  double *fPtLimits;              //! limits of pT bins
  
  AliRDHFCuts* fAnalysisCuts;     /// cuts
  
  int fVarForProp;                /// variable used for propagation (p or pT)

  ClassDef(AliAnalysisTaskSEDmesonPIDSysProp, 3);
};

#endif
