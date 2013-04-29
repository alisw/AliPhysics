#ifndef ALIANALYSISTASKSERECOJETCORRELATIONS_H
#define ALIANALYSISTASKSERECOJetCORRELATIONS_H

// $Id$

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskFlavourJetCorrelations
// AliAnalysisTaskSE for Dmesons - jet correlations analysis
// Author: Xiaoming Zhang, xmzhang@lbl.gov
//*************************************************************************

#include "AliAnalysisTaskEmcalJet.h"

class TList;
class TClonesArray;

class AliAnalysisTaskFlavourJetCorrelations : public AliAnalysisTaskEmcalJet {

 public :

  enum {
    kMatchConeCandi,
    kMatchAreaCandi,
    kMatchConeProng,
    kMatchAreaProng,
    kDzeroMatchType
  };

  AliAnalysisTaskFlavourJetCorrelations();
  AliAnalysisTaskFlavourJetCorrelations(const char *name, Bool_t bIsHisto=kTRUE);
  virtual ~AliAnalysisTaskFlavourJetCorrelations();

  virtual void UserCreateOutputObjects();

 private :

  AliAnalysisTaskFlavourJetCorrelations(const AliAnalysisTaskFlavourJetCorrelations &);
  AliAnalysisTaskFlavourJetCorrelations& operator=(const AliAnalysisTaskFlavourJetCorrelations &);

  virtual void   ExecOnce();
  virtual Bool_t FillGeneralHistograms();
  virtual Bool_t FillHistograms();
  virtual Bool_t IsEventSelected();
  virtual Bool_t RetrieveEventObjects();
  virtual Bool_t Run();

  void   MakeControlHistograms();
  Bool_t FillControlHistograms();

  void CreateDzeroHistograms();
  void CreateDstarHistograms();

  void RunDzeroJet(AliEmcalJet const *pJet, const Int_t iJetPtBin, const Bool_t bIsD0);
  void RunDstarJet(AliEmcalJet const *pJet, const Int_t iJetPtBin);

  TClonesArray *fUsedDzeros;  //! input Dzero candidates array
  TClonesArray *fUsedD0bars;  //! input D0bar candidates array
  TClonesArray *fUsedDstars;  //! input Dstar candidates array

  TList *fListControlHistos;  //! list of output contral histograms
  TList *fListAnDzeroHistos;  //! list of output Dzero - jet correlation histograms
  TList *fListAnDstarHistos;  //! list of output Dstar - jet correlation histograms

  ClassDef(AliAnalysisTaskFlavourJetCorrelations, 1);
};

#endif
