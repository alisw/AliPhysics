#ifndef ALIANALYSISTASKSEDMESONSFILTERCJ_H
#define ALIANALYSISTASKSEDMESONSFILTERCJ_H
// $Id$
/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//-----------------------------------------------------------------------
// Author : A. Grelli,  Utrecht University
//          C. Bianchin, Utrecht University
//          X. Zhang, LBNL
//-----------------------------------------------------------------------


#include <TH2F.h>
#include "AliAODEvent.h"
#include "AliPicoTrack.h"
#include "AliAnalysisTaskSE.h"

class TH3F;
class TString;
class TParticle ;
class TClonesArray ;
class AliMCParticle;
class AliAODMCParticle;
class AliRDHFCuts;
class AliAODRecoCascadeHF;

class AliAnalysisTaskSEDmesonsFilterCJ : public AliAnalysisTaskSE 
{

 public :

  enum ECandidateType{ kD0toKpi, kDstartoKpipi };
  
  AliAnalysisTaskSEDmesonsFilterCJ();
  AliAnalysisTaskSEDmesonsFilterCJ(const Char_t* name,AliRDHFCuts* cuts,ECandidateType candtype);
  virtual ~AliAnalysisTaskSEDmesonsFilterCJ();

  virtual void     UserCreateOutputObjects();
  virtual void     UserExec(Option_t *option);
  virtual void     Terminate(Option_t *);
  virtual void     Init();
  virtual void     LocalInit() { Init(); }

  // inizializations
  Bool_t DefineHistoForAnalysis();

  // set MC usage
  void   SetMC(Bool_t theMCon) { fUseMCInfo = theMCon; }
  Bool_t GetMC() const { return fUseMCInfo; }
  
  void SetMassLimits(Double_t range, Int_t pdg);
  void SetMassLimits(Double_t lowlimit, Double_t uplimit);

  // Array of D0 width for the Dstar
  Bool_t SetD0WidthForDStar(Int_t nptbins, Float_t *width);

 private :
  
  AliAnalysisTaskSEDmesonsFilterCJ(const AliAnalysisTaskSEDmesonsFilterCJ &source);
  AliAnalysisTaskSEDmesonsFilterCJ& operator=(const AliAnalysisTaskSEDmesonsFilterCJ& source); 

  Bool_t fUseMCInfo;               //  Use MC info

  UInt_t  fCandidateType;          //  Dstar or D0
  TString fCandidateName;          //  Dstar or D0

  Int_t fPDGmother;                //  PDG code of D meson
  Int_t fNProngs;                  //  number of prong of the decay channel  
  Int_t fPDGdaughters[4];          //  PDG codes of daughters
  Float_t fSigmaD0[30];            //  D0 sigma for Dstar

  TString fBranchName;             //  AOD branch name
  TList  *fOutput;                 //! user output
//TList *fOutputCandidates;        //! output of array of candidates (kExchange)

  AliRDHFCuts *fCuts;              //! Cuts 
  Double_t fMinMass;               //  mass lower limit histogram
  Double_t fMaxMass;               //  mass upper limit histogram

  TClonesArray *fCandidateArray;   //! contains candidates selected by AliRDHFCuts
//TClonesArray *fIsSelectedArray;  //! contains result of IsSelected for candidates which pass the cuts (needed for D0)

  ClassDef(AliAnalysisTaskSEDmesonsFilterCJ,1); // class for charm-jet correlations
};

#endif
