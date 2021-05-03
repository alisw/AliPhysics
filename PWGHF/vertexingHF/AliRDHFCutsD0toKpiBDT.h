#ifndef ALIRDHFCUTSD0TOKPIBDT_H
#define ALIRDHFCUTSD0TOKPIBDT_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//***********************************************************
/// \class Class AliRDHFCutsD0toKpiBDT
/// \brief class for inherited from AliRDHFBDTD0toKpi with additional selection based on AliRDHFBDT class objects 
/// \author Author: M.Cai, cai.mengke@cern.ch
//***********************************************************

#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFBDT.h"

class AliAODEvent;
class AliAODRecoDecayHF;
class AliAODRecoDecayHF2Prong;
class AliRDHFCutsD0toKpi;

class AliRDHFCutsD0toKpiBDT : public AliRDHFCutsD0toKpi 
{
 public:


  AliRDHFCutsD0toKpiBDT(const char* name="CutsD0toKpiBDT");
  virtual ~AliRDHFCutsD0toKpiBDT();
  AliRDHFCutsD0toKpiBDT(const AliRDHFCutsD0toKpi& source);
  AliRDHFCutsD0toKpiBDT(const AliRDHFCutsD0toKpiBDT& source);
  AliRDHFCutsD0toKpiBDT& operator=(const AliRDHFCutsD0toKpiBDT& source); 

  using AliRDHFCutsD0toKpi::AliRDHFCutsD0toKpi;
  
  void SetBDTNames(Int_t nOpt,TString *Names,Bool_t *isUpperCut);
  void SetPtBinsBDT(Int_t nPtBinLimits,Float_t *ptBinLimits);
  void SetNPtBinsBDT(Int_t nptBins){fnPtBinsBDT=nptBins;}
  void SetBDTCutGlobalIndex(){fBDTCutGlobalIndex=fNBDTOpt*fnPtBinsBDT;}
  void SetBDTCutGlobalIndex(Int_t n,Int_t nptBins){fNBDTOpt=n; fnPtBinsBDT=nptBins; SetBDTCutGlobalIndex();} // Set this first to initialize cut
  void SetBDTCuts(Int_t nOpt,Int_t nPtBins,Float_t** cuts);
  void SetBDTCuts(Int_t glIndex, Float_t* cutsRDGlob);
  
  void SetListOfBDT(TList *l) const;
  
  Float_t *GetPtBinLimitsBDT() const {return fPtBinLimitsBDT;}
  Int_t   GetNPtBinsBDT() const {return fnPtBinsBDT;}
  Int_t   GetNBDTOpt() const {return fNBDTOpt;}
  TString *GetBDTNames() const {return fBDTNames;}
  Bool_t  *GetIsUpperCutBDT() const {return fIsUpperCutBDT;}
  Int_t	GetBDTCutGlobalIndex(Int_t iOpt,Int_t iPtBin) const;
  TList *GetListOfBDT() const {return fListOfBDT;}
  
  Int_t IsSelectedBDT(AliAODRecoDecayHF2Prong *d, AliAODEvent *aod) const; // using BDTResp cut to determine whether it is selected
  Double_t GetBDTResponse(AliAODRecoDecayHF2Prong *d, AliAODEvent *aod, Int_t iOpt, Int_t isD0bar) const; // Get BDT output
  Double_t GetRDHFVarsForSel(AliAODRecoDecayHF2Prong *d, AliAODEvent *aod, TString VarName, Int_t isD0bar) const; // Get RDHF2prong variable by name
  
  Int_t PtBinBDT(Float_t pt) const {
    if(pt<fPtBinLimitsBDT[0]) return -1;
    for (Int_t i=0;i<fnPtBinsBDT;i++) if(pt<fPtBinLimitsBDT[i+1]) return i;
    return -1;
  }
  virtual void PrintAll()const;
  

 protected:
  
  TList *fListOfBDT;
  Int_t fNBDTOpt;	/// Number of BDT
  
  Int_t fnPtBinsBDT;  /// allow to use different number of pt bins from cuts
  Int_t fnPtBinLimitsBDT; /// fnPtBinsBDT+1
  Float_t* fPtBinLimitsBDT; //[fnPtBinLimits]
  TString *fBDTNames; //[fNBDTOpt] names of the BDT for each pT bin
  Int_t fBDTCutGlobalIndex; /// fnVars*fnPtBins
  Float_t *fBDTCuts; //[fBDTCutGlobalIndex] the cuts values
  Bool_t  *fIsUpperCutBDT; //[fNBDTOpt] use > or < to select

  /// \cond CLASSIMP    
  ClassDef(AliRDHFCutsD0toKpiBDT,1);  /// class for cuts on AOD reconstructed D0->Kpi
  /// \endcond
};

#endif
