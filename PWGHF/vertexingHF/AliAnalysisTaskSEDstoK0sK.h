#ifndef ALIANALYSISTASKSEDSTOK0SK_H
#define ALIANALYSISTASKSEDSTOK0SK_H

/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////////
///
/// \class AliAnalysisTaskSEDstoK0sK
/// \brief AliAnalysisTaskSE to produce Ds->K0S+K invariant mass spectra
///        and THnSparse for cut optimisations.
///
/// \author J.Hamon, julien.hamon@cern.ch (IPHC)
///
/////////////////////////////////////////////////////////////


#include <TROOT.h>
#include <TSystem.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TNtuple.h>

#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"
#include "AliAODRecoCascadeHF.h"
#include "AliNormalizationCounter.h"
#include "AliRDHFCutsDstoK0sK.h"



class AliAnalysisTaskSEDstoK0sK : public AliAnalysisTaskSE
{
public:


   AliAnalysisTaskSEDstoK0sK();
   AliAnalysisTaskSEDstoK0sK(const char *name,
                             AliRDHFCutsDstoK0sK *cuts,
                             Bool_t readMC,
                             Bool_t fillNtuple,
                             Int_t nCutsTuple,
                             Float_t *minCutsTuple,
                             Float_t *maxCutsTuple);
   virtual ~AliAnalysisTaskSEDstoK0sK();


   // Implementation of interface methods
   virtual void Init();
   virtual void LocalInit() { Init(); }
   virtual void UserCreateOutputObjects();
   virtual void UserExec(Option_t* /*option*/);
   virtual void Terminate(Option_t* /*option*/);


   // Setters
   void SetUseSelectionBit(Bool_t flag)      { fUseSelectionBit = flag; }
   void SetAODMismatchProtection(Int_t opt)  { fAODProtection = opt;    }
   void SetPtBins(Int_t nBins, Float_t* limitsPt);
   void SetMassRangeAndBinSize(Float_t range, Float_t bin) { fMassRange = range; fMassBinSize = bin; }
   void SetCutsTupleVariables(Int_t nCuts, Float_t *minCuts, Float_t *maxCuts);



protected:


   enum { kMaxPtBins=20, kNTupleVars=34, kNTupleVarsMC=35 };


   AliAnalysisTaskSEDstoK0sK(const AliAnalysisTaskSEDstoK0sK &source);
   AliAnalysisTaskSEDstoK0sK& operator=(const AliAnalysisTaskSEDstoK0sK& source);


   void FillHistogramsVar(AliAODRecoCascadeHF* dCan, AliRDHFCuts::ESelLevel selFlag, TClonesArray* mcArray = 0);
   void FillHistogramsPID(AliAODRecoCascadeHF* dCan, AliRDHFCuts::ESelLevel selFlag, TClonesArray* mcArray = 0);
   void FillTheTree(AliAODRecoCascadeHF* dCan, AliAODEvent* aod, TClonesArray* mcArray);
   Float_t ComputeSigmaVert(const AliAODEvent* aod, AliAODRecoCascadeHF* dCan) const;
   Float_t CosThetaK0sBachRFrame(AliAODRecoCascadeHF* dCan) const;
   Int_t   MatchToMCDstoK0sKSignal(AliAODRecoCascadeHF* dCan, TClonesArray* mcArray);
   Int_t   MatchToMCDplustoK0spiSignal(AliAODRecoCascadeHF* dCan, TClonesArray* mcArray);


   // Outputs
   AliRDHFCutsDstoK0sK*      fAnalysisCuts;  ///    Cut object for Analysis on output slot #1
   AliNormalizationCounter*  fCounter;       //!<!  AliNormalizationCounter sent on output slot #2
   TList*                    fOutputSele;    //!<!  Various histograms of selected events: TList sent to output slot #3
   TList*                    fOutputCand;    //!<!  Candidate level histograms: TList sent to output slot #4
   TList*                    fOutputPID;     //!<!  PID level histograms: TList sent to output slot #5
   TList*                    fOutputMC;      //!<!  PID level histograms: TList sent to output slot #5
   TNtuple*                  fOutputNtuple;  //!<!  TNtuple for candidates on data sent to output slot #4


   // Histograms of output slot #3 (fOutputSele)
   TH1F*  fHisNEvents;                       //!<!  Counter: events and candidates
   TH1F*  fHisCentrality[3];                 //!<!  Centrality: all, selected and rejected
   TH2F*  fHisCentralityVsMult[3];           //!<!  Centrality VS Multiplicity: all, selected and rejected
   TH2F*  fHisRapidity;                      //!<!  Rapidity all candidates
   TH2F*  fHisRapiditySel;                   //!<!  Rapidity selected (kCandidate) candidates
   TH2F*  fHisPidTPCKaonVsPt;                //!<!
   TH2F*  fHisPidTOFKaonVsPt;                //!<!
   TH2F*  fHisPidTPCTOFKaon;                 //!<!
   TH2F*  fHisPidTPCTOFKaonSel;              //!<!
   TH2F*  fHisPidTPCKaonVsPtSel;             //!<!
   TH2F*  fHisPidTOFKaonVsPtSel;             //!<!


   // Histograms of output slot #4 (fOutputCand) and #5 (fOutputPID)
   TH2F*  fHisInvMassDs[5];                  //!<! (kCandidate, kPID, mcSignal, mcBackground, mcReflection)
   TH2F*  fHisInvMassDplus[5];               //!<! (kCandidate, kPID, mcSignal, mcBackground, mcReflection)
   TH2F*  fHisInvMassK0s[5];                 //!<! (kCandidate, kPID, mcSignal, mcBackground, mcReflection)
   TH2F*  fHisPtK0s[5];                      //!<! (kCandidate, kPID, mcSignal, mcBackground, mcReflection)
   TH2F*  fHisPtBachelor[5];                 //!<! (kCandidate, kPID, mcSignal, mcBackground, mcReflection)
   TH2F*  fHisImpParK0s[5];                  //!<! (kCandidate, kPID, mcSignal, mcBackground, mcReflection)
   TH2F*  fHisImpParBach[5];                 //!<! (kCandidate, kPID, mcSignal, mcBackground, mcReflection)
   TH2F*  fHisCTauK0s[5];                    //!<! (kCandidate, kPID, mcSignal, mcBackground, mcReflection)
   TH2F*  fHisCosPointingDs[5];              //!<! (kCandidate, kPID, mcSignal, mcBackground, mcReflection)
   TH2F*  fHisCosPointingXYDs[5];            //!<! (kCandidate, kPID, mcSignal, mcBackground, mcReflection)
   TH2F*  fHisCosThetaStarK0s[5];            //!<! (kCandidate, kPID, mcSignal, mcBackground, mcReflection)
   TH2F*  fHisCosThetaStarBach[5];           //!<! (kCandidate, kPID, mcSignal, mcBackground, mcReflection)
   TH2F*  fHisDCAK0sBach[5];                 //!<! (kCandidate, kPID, mcSignal, mcBackground, mcReflection)
   TH2F*  fHisDecayLxyDs[5];                 //!<! (kCandidate, kPID, mcSignal, mcBackground, mcReflection)
   TH2F*  fHisNormDecayLxyDs[5];             //!<! (kCandidate, kPID, mcSignal, mcBackground, mcReflection)


   // Data members for analysis
   Int_t    fNPtBins;                        ///  Number of Pt bins
   Float_t  fPtLimits[kMaxPtBins+1];         ///  Limits of Pt bins
   Float_t  fMassRange;                      ///  Size range of invariant mass histograms
   Float_t  fMassBinSize;                    ///  Bin size of invariant mass histograms (GeV)
   Bool_t   fReadMC;                         ///  Flag for accessing MC
   Bool_t   fUseSelectionBit;                ///  Flag for using selection bit (to select Ds flags)
   Bool_t   fFillNtuple;                     ///  Flag for using THnSparse
   Float_t  fCutsMinTupleVars[kNTupleVars];  ///  Minimum cut values for tuple variables
   Float_t  fCutsMaxTupleVars[kNTupleVars];  ///  Maximum cut values for tuple variables
   Int_t    fAODProtection;                  ///  Protection against AOD-deltaAOD mismatch
                                             /// -1: no protection,  0: check AOD/dAOD nEvents only,  1: check AOD/dAOD nEvents + TProcessID names


   /// \cond CLASSIMP
   ClassDef(AliAnalysisTaskSEDstoK0sK, 1);   /// class for Ds->K0S+K invariant mass spectra
   /// \endcond
};



inline Int_t AliAnalysisTaskSEDstoK0sK::MatchToMCDstoK0sKSignal(AliAODRecoCascadeHF* dCan, TClonesArray* mcArray)
{
   Int_t pdgDgDstoK0sK[2]   = {321, 310};
   Int_t pdgDgK0stoPions[2] = {211, 211};
   return (Int_t) dCan->MatchToMC(431, pdgDgDstoK0sK[1], pdgDgDstoK0sK, pdgDgK0stoPions, mcArray, kTRUE);
}


inline Int_t AliAnalysisTaskSEDstoK0sK::MatchToMCDplustoK0spiSignal(AliAODRecoCascadeHF* dCan, TClonesArray* mcArray)
{
   Int_t pdgDgDplustoK0spi[2] = {211, 310};
   Int_t pdgDgK0stoPions[2]   = {211, 211};
   return (Int_t) dCan->MatchToMC(411, pdgDgDplustoK0spi[1], pdgDgDplustoK0spi, pdgDgK0stoPions, mcArray, kTRUE);
}


#endif
