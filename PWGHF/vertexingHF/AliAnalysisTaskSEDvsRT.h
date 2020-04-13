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
#include "AliNormalizationCounter.h"

class AliAnalysisTaskSEDvsRT : public AliAnalysisTaskSE
{
 public:

   AliAnalysisTaskSEDvsRT();
   AliAnalysisTaskSEDvsRT(const char *name, Int_t pdgSpecies, AliRDHFCuts *cuts  );
   virtual ~AliAnalysisTaskSEDvsRT();

   void SetMassLimits(Double_t lowlimit, Double_t uplimit);
   void SetMassLimits(Int_t pdg, Double_t range);
   Double_t GetUpperMassLimit() const {return fUpmasslimit;}
   Double_t GetLowerMassLimit() const {return fLowmasslimit;}
   void SetNMassBins(Int_t nbins) {fNMassBins = nbins;}
   Int_t GetNMassBins() const {return fNMassBins;}
   
   void SetReadMC(Bool_t readMC=kTRUE) {fReadMC = readMC;}
   void SetMCOption(Int_t option=0)    {fMCOption = option;}
   void SetLctoV0(Bool_t option=kTRUE) {fLctoV0 = option;}
   void SetUseBit(Bool_t option=kTRUE) {fUseBit = option;}
   void SetIsPPb(Bool_t option=kTRUE)  {fisPPbData = option;}
   void SetUseNsparse(Bool_t option=kTRUE) {fUseNsparse = option;}
   
   void SetEtaCut(Double_t etacut)        {fEtaCut = etacut;}
   void SetPtLeadMin(Double_t pt)         {fLeadMin = pt;}
   void SetAveMultiInTrans(Double_t mult) {fAveMultiInTrans = mult;}
   
   
   /// Implementation of interface methods
   virtual void UserCreateOutputObjects();
   virtual void Init();
   virtual void LocalInit() {Init();}
   virtual void UserExec(Option_t *option);
   virtual void Terminate(Option_t *option);

 private:
 
   AliAnalysisTaskSEDvsRT(const AliAnalysisTaskSEDvsRT &source);
   AliAnalysisTaskSEDvsRT& operator=(const AliAnalysisTaskSEDvsRT& source);

   virtual Double_t CalculateRTVal(AliAODEvent* esdEvent);
   virtual ULong64_t GetEventIdAsLong(AliVHeader* header);
   virtual TObjArray* FindLeading(TObjArray* array);
   virtual void QSortTracks(TObjArray &a, Int_t first, Int_t last); 
   virtual TObjArray* SortRegions(const AliVParticle* leading, TObjArray *array);
   virtual TObjArray* GetMinMaxRegion(TList *transv1, TList *transv2);
   TList *fOutput;         //!<! list to send on output slot 1
   TList *fListCuts;       ///list of cuts
   TList *fOutputCounters; //!<! list to send on output slot 3
   TList *fListQAhists;    //!<! list of QA plots on output slot 4
   
   
   Double_t fUpmasslimit;  /// upper inv. mass limit for histos
   Double_t fLowmasslimit; /// lower inv. mass limit for histos
   Int_t fNMassBins;       /// nbins for invariant mass
   
   AliRDHFCuts *fRDCutsAnalysis; /// cuts for analysis
   
   
   TH1F *fHistNEvents;     //!<! hist. for number of events
   TH1F *fGlobalRT;        //!<! hist for RT distribution without D-meson selection
   TH1F *fHistPtLead;      //!<! hist for pT distribution of leading track
   TH3F *fRTvsZvtxvsMult;        //!<! distribution of RT as function of z of primary vertex and tracklet multiplicity
   
   
   AliNormalizationCounter *fCounter;  //!<! Counter for normalisation
   Bool_t fReadMC;         /// flag for reading MC
   Int_t  fMCOption;       /// 0 = keep all cand, 1 = keep only signal, 2 = keep only bkg
   Bool_t fUseBit;         /// flag to use bitmask
   Int_t fAODProtection;   /// flag to activate protection against AOD-dAOD mismatch
   /// -1: no protection, 0: check nEvents only, 1: check nEvents + TProcessID names
   
   
   Int_t fPdgSpecies; /// pdg code of analysed particle
   Bool_t fLctoV0;    /// flag for Lc->pK0 analysis
   Bool_t fisPPbData; /// flag for pPb running

  
   Double_t fEtaCut;  /// cut on eta for RT track counting
   Double_t fLeadMin; /// minimum pt for leading
   
   Double_t fAveMultiInTrans; /// average multiplicity in transverse region   

   Double_t fPhiLeading;      ///  phi of leading particle
   TH3D *fPtvsMassvsRTToward; //!<! output hist for mass vs RT (for candidates toward)
   TH3D *fPtvsMassvsRTAway;   //!<! output hist for mass vs RT (for candidates away)
   TH3D *fPtvsMassvsRTTrans;  //!<! output hist for mass vs RT (for transverse candidates)
   AliAnalysisFilter* fTrackFilter[18]; //! track filter 
   Bool_t fUseNsparse;        /// switch to give nsparse in output
   THnSparse *fOutNsparse;    //!<! output THnSparse for RT analysis
   
   /// \cond CLASSIMP
   ClassDef(AliAnalysisTaskSEDvsRT,2); /// charmed hadrons vs. RT task
   /// \endcond
};

#endif
