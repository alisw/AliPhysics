#ifndef ALIANALYSISTASKSEDVSRT_H
#define ALIANALYSISTASKSEDVSRT_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//*************************************************************************
/// \class Class AliAnalysisTaskSEDvsRT
/// \brief AliAnalysisTaskSE for the charmed hadron vs. RT analysis
/// \author Authors: Jeremy Wilkinson,
//*************************************************************************
#include <TROOT.h>
#include <TSystem.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TArrayD.h>
#include <TFile.h>
#include "AliRDHFCuts.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliRDHFCutsLctoV0.h"
#include "AliRDHFCutsLctopKpi.h"

class AliAnalysisTaskSEDvsRT : public AliAnalysisTaskSE
{
 public:

   AliAnalysisTaskSEDvsRT();
   AliAnalysisTaskSEDvsRT(   );
   virtual ~AliAnalysisTaskSEDvsMultiplicity();

   void SetMassLimits(Double_t lowlimit, Double_t uplimit);
   void SetMassLimits(Int_t pdg, Double_t range);
   Double_t GetUpperMassLimit() const {return fUpmasslimit;}
   Double_t GetLowerMassLimit() const {return fLowmasslimit;}
   void SetNMassBins(Int_t nbins) {fNMassBins = nbins;}
   Int_t GetNMassBins() const {return fNMassBins;}
   
   void SetReadMC(Bool_t readMC=kTRUE) {fReadMC = readMC;}
   void SetMCOption(Int_t option=0)    {fMCOption = option;}
   
   /// Implementation of interface methods
   virtual void UserCreateOutputObjects();
   virtual void Init();
   virtual void LocalInit() {Init();}
   virtual void UserExec(Option_t *option);
   virtual void Terminate(Option_t *option);

 private:
 
   AliAnalysisTaskSEDvsRT(const AliAnalysisTaskSEDvsRT &source);
   AliAnalysisTaskSEDvsRT& operator=(const AliAnalysisTaskSEDvsRT& source);



   TList *fOutput;         //!<! list to send on output slot 1
   TList *fListCuts;       ///list of cuts
   TList *fOutputCounters; //!<! list to send on output slot 3
   
   Double_t fUpmasslimit;  /// upper inv. mass limit for histos
   Double_t fLowmasslimit; /// lower inv. mass limit for histos
   Int_t fNMassBins;       /// nbins for invariant mass
   
   AliRDHFCuts *fRDCutsAnalysis; /// cuts for analysis
   
   
   TH1F *fHistNEvents;     //!<! hist. for number of events
   AliNormalizationCounter *fCounter;  //!<! Counter for normalisation
   
   Double_t fUpmasslimit;  /// upper inv mass limit
   Double_t fLowmasslimit; /// lower inv mass limit
   Int_t  fNMassBins;      /// nbins for invariant mass
   Bool_t fReadMC;         /// flag for reading MC
   Int_t  fMCOption;       /// 0 = keep all cand, 1 = keep only signal, 2 = keep only bkg
   Bool_t fUsetBit;        /// flag to use bitmask
   Int_t fAODProtection;   /// flag to activate protection against AOD-dAOD mismatch
   /// -1: no protection, 0: check nEvents only, 1: check nEvents + TProcessID names
   
   
   Int_t fPdgSpecies; /// pdg code of analysed particle
   
   /// \cond CLASSIMP
   ClassDef(AliAnalysisTaskSEDvsRT,1); /// charmed hadrons vs. RT task
   /// \endcond
};

#endif