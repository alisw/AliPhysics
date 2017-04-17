#ifndef AlIRDHFCUTSB0TODSTARPI_H
#define AlIRDHFCUTSB0TODSTARPI_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//***********************************************************
/// \class Class AliRDHFCutsB0toDStarPi
/// \brief class for cuts on AOD reconstructed B0->DStarPi->D0PiPi->KPiPiPi
/// \author Author: Lennart van Doremalen, l.v.r.vandoremalen@uu.nl  
//
// Based on work by A.Grelli, alessandro.grelli@uu.nl
// PID method implemented by   Y.Wang, yifei@physi.uni-heidelberg.de
//***********************************************************

#include "AliRDHFCuts.h"

class AliAODEvent;
class AliAODRecoCascadeHF;
class AliAODRecoDecayHF;

class AliRDHFCutsB0toDStarPi : public AliRDHFCuts 
{
 public:

  enum EUpperCut {kCutBelowValue             = 0,
                  kCutAboveValue             = 1
                 };

  AliRDHFCutsB0toDStarPi(const char* name="CutsB0toDStarPi");
  
  virtual ~AliRDHFCutsB0toDStarPi();

  AliRDHFCutsB0toDStarPi(const AliRDHFCutsB0toDStarPi& source);
  AliRDHFCutsB0toDStarPi& operator=(const AliRDHFCutsB0toDStarPi& source); 
 
  using AliRDHFCuts::GetCutVarsForOpt;
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters);

  using AliRDHFCuts::IsSelected;
  virtual Int_t IsSelected(TObject* obj, Int_t selectionLevel, AliAODEvent* aod, Bool_t* bCutArray);
  virtual Int_t IsSelected(TObject* obj, Int_t selectionLevel) {::Error("AliAnalysisTaskB0toDStarPi", "Old selection function"); return 0;}

  Int_t IsDStarFromB0Selected(Double_t pt, TObject* obj,Int_t selectionLevel, AliAODEvent* aod, Bool_t* bCutArray);
  Int_t IsD0FromDStarSelected(Double_t pt, TObject* obj,Int_t selectionLevel, AliAODEvent* aod, Bool_t* bCutArray);
  Int_t IsD0forD0ptbinSelected(TObject* obj,Int_t selectionLevel, AliAODEvent* aod, Bool_t* bCutArray);
  Int_t IsD0forDStarptbinSelected(TObject* obj,Int_t selectionLevel, AliAODEvent* aod, Bool_t* bCutArray);
  Int_t IsDStarforDStarptbinSelected(TObject* obj,Int_t selectionLevel, AliAODEvent* aod, Bool_t* bCutArray);

  virtual Int_t IsSelectedPID(AliAODRecoDecayHF *rd);
  virtual Int_t SelectPID(AliAODTrack *track, Int_t type);
  virtual Bool_t IsInFiducialAcceptance(Double_t pt,Double_t y) const;
  
  void SetMaxPtPid(Float_t maxPt){fMaxPtPid = maxPt;}

  void SetOffHighPtPIDinTPC(Float_t TPCrem =999.){fTPCflag = TPCrem;}

  void SetGetCutInfo(Bool_t value){fGetCutInfo = value;}
  void InitializeCuts();

  Double_t GetCircRadius() { return fCircRadius; }
  void SetCircRadius(Double_t radius) { fCircRadius = radius; }

  Double_t DeltaInvMassDStarKpipi(AliAODRecoCascadeHF * DStar) const;
  Double_t DeltaInvMassB0Kpipipi(AliAODRecoCascadeHF * B0) const;

  void SetCutsD0forD0ptbin(Int_t nVars,Int_t nPtBins,Float_t **cutsRDD0forD0ptbin);
  Int_t PtBinD0forD0ptbin(Double_t pt) const;
  void SetPtBinsD0forD0ptbin(Int_t nPtBinLimits,Float_t *ptBinLimits);

  void SetCutsD0forDStarptbin(Int_t nVars,Int_t nPtBins,Float_t **cutsRDD0forDStarptbin);
  Int_t PtBinD0forDStarptbin(Double_t pt) const;
  void SetPtBinsD0forDStarptbin(Int_t nPtBinLimits,Float_t *ptBinLimits);

  void SetCutsDStarforDStarptbin(Int_t nVars,Int_t nPtBins,Float_t **cutsRDDStarforDStarptbin);
  Int_t PtBinDStarforDStarptbin(Double_t pt) const;
  void SetPtBinsDStarforDStarptbin(Int_t nPtBinLimits,Float_t *ptBinLimits);

  Float_t *GetPtBinLimitsD0forD0ptbin() const {return fPtBinLimitsD0forD0ptbin;}
  Int_t GetNPtBinsD0forD0ptbin() const {return fnPtBinsD0forD0ptbin;}
  Int_t GetNVarsD0forD0ptbin() const {return fnVarsD0forD0ptbin;}
  Int_t GetGlobalIndexD0forD0ptbin(Int_t iVar,Int_t iPtBin) const{return iPtBin*fnVarsD0forD0ptbin+iVar;}
  void  SetGlobalIndexD0forD0ptbin(){fGlobalIndexD0forD0ptbin=fnVarsD0forD0ptbin*fnPtBinsD0forD0ptbin;}
  void  SetNPtBinsD0forD0ptbin(Int_t nptBins){fnPtBinsD0forD0ptbin=nptBins;}
  void  SetNVarsD0forD0ptbin(Int_t nVars){fnVarsD0forD0ptbin=nVars;} 

  Float_t *GetPtBinLimitsD0forDStarptbin() const {return fPtBinLimitsD0forDStarptbin;}
  Int_t GetNPtBinsD0forDStarptbin() const {return fnPtBinsD0forDStarptbin;}
  Int_t GetNVarsD0forDStarptbin() const {return fnVarsD0forDStarptbin;}
  Int_t GetGlobalIndexD0forDStarptbin(Int_t iVar,Int_t iPtBin) const{return iPtBin*fnVarsD0forDStarptbin+iVar;}
  void  SetGlobalIndexD0forDStarptbin(){fGlobalIndexD0forDStarptbin=fnVarsD0forDStarptbin*fnPtBinsD0forDStarptbin;}
  void  SetNPtBinsD0forDStarptbin(Int_t nptBins){fnPtBinsD0forDStarptbin=nptBins;}
  void  SetNVarsD0forDStarptbin(Int_t nVars){fnVarsD0forDStarptbin=nVars;} 

  Float_t *GetPtBinLimitsDStarforDStarptbin() const {return fPtBinLimitsDStarforDStarptbin;}
  Int_t GetNPtBinsDStarforDStarptbin() const {return fnPtBinsDStarforDStarptbin;}
  Int_t GetNVarsDStarforDStarptbin() const {return fnVarsDStarforDStarptbin;}
  Int_t GetGlobalIndexDStarforDStarptbin(Int_t iVar,Int_t iPtBin) const{return iPtBin*fnVarsDStarforDStarptbin+iVar;}
  void  SetGlobalIndexDStarforDStarptbin(){fGlobalIndexDStarforDStarptbin=fnVarsDStarforDStarptbin*fnPtBinsDStarforDStarptbin;}
  void  SetNPtBinsDStarforDStarptbin(Int_t nptBins){fnPtBinsDStarforDStarptbin=nptBins;}
  void  SetNVarsDStarforDStarptbin(Int_t nVars){fnVarsDStarforDStarptbin=nVars;} 

  void SetIsUpperCut(Int_t nCutIndex, Bool_t isUpperCut){fIsUpperCut[nCutIndex] = isUpperCut; return;}
  void SetIsUpperCutD0forD0ptbin(Int_t nCutIndex, Bool_t isUpperCut){fIsUpperCutD0forD0ptbin[nCutIndex] = isUpperCut; return;}
  void SetIsUpperCutD0forDStarptbin(Int_t nCutIndex, Bool_t isUpperCut){fIsUpperCutD0forDStarptbin[nCutIndex] = isUpperCut; return;}
  void SetIsUpperCutDStarforDStarptbin(Int_t nCutIndex, Bool_t isUpperCut){fIsUpperCutDStarforDStarptbin[nCutIndex] = isUpperCut; return;}

  Bool_t GetIsUpperCut(Int_t nCutIndex){return fIsUpperCut[nCutIndex];}
  Bool_t GetIsUpperCutD0forD0ptbin(Int_t nCutIndex){return fIsUpperCutD0forD0ptbin[nCutIndex];}
  Bool_t GetIsUpperCutD0forDStarptbin(Int_t nCutIndex){return fIsUpperCutD0forDStarptbin[nCutIndex];}
  Bool_t GetIsUpperCutDStarforDStarptbin(Int_t nCutIndex){return fIsUpperCutDStarforDStarptbin[nCutIndex];}

  void SetIsCutUsed(Int_t nCutIndex, Int_t ptbin, Bool_t isCutUsed){fIsCutUsed[GetGlobalIndex(nCutIndex,ptbin)] = isCutUsed; return;}
  void SetIsCutUsedD0forD0ptbin(Int_t nCutIndex, Int_t ptbin, Bool_t isCutUsed){fIsCutUsedD0forD0ptbin[GetGlobalIndexD0forD0ptbin(nCutIndex,ptbin)] = isCutUsed; return;}
  void SetIsCutUsedD0forDStarptbin(Int_t nCutIndex, Int_t ptbin, Bool_t isCutUsed){fIsCutUsedD0forDStarptbin[GetGlobalIndexD0forDStarptbin(nCutIndex,ptbin)] = isCutUsed; return;}
  void SetIsCutUsedDStarforDStarptbin(Int_t nCutIndex, Int_t ptbin, Bool_t isCutUsed){fIsCutUsedDStarforDStarptbin[GetGlobalIndexDStarforDStarptbin(nCutIndex,ptbin)] = isCutUsed; return;}

  Bool_t GetIsCutUsed(Int_t nCutIndex, Int_t ptbin){return fIsCutUsed[GetGlobalIndex(nCutIndex,ptbin)];}
  Bool_t GetIsCutUsedD0forD0ptbin(Int_t nCutIndex, Int_t ptbin){return fIsCutUsedD0forD0ptbin[GetGlobalIndexD0forD0ptbin(nCutIndex,ptbin)];}
  Bool_t GetIsCutUsedD0forDStarptbin(Int_t nCutIndex, Int_t ptbin){return fIsCutUsedD0forDStarptbin[GetGlobalIndexD0forDStarptbin(nCutIndex,ptbin)];}
  Bool_t GetIsCutUsedDStarforDStarptbin(Int_t nCutIndex, Int_t ptbin){return fIsCutUsedDStarforDStarptbin[GetGlobalIndexDStarforDStarptbin(nCutIndex,ptbin)];}

  Int_t ApplyCutOnVariable(Int_t nCutIndex, Int_t ptbin, Float_t cutVariableValue, Bool_t bCutArray[85]);
  Int_t ApplyCutOnVariableD0forD0ptbin(Int_t nCutIndex, Int_t ptbin, Float_t cutVariableValue, Bool_t bCutArray[25]);
  Int_t ApplyCutOnVariableD0forDStarptbin(Int_t nCutIndex, Int_t ptbin, Float_t cutVariableValue, Bool_t bCutArray[35]);
  Int_t ApplyCutOnVariableDStarforDStarptbin(Int_t nCutIndex, Int_t ptbin, Float_t cutVariableValue, Bool_t bCutArray[25]);

  void SetVarNamesD0forD0ptbin(Int_t nVars,TString *varNames,Bool_t *isUpperCut);
  void SetVarNamesD0forDStarptbin(Int_t nVars,TString *varNames,Bool_t *isUpperCut);
  void SetVarNamesDStarforDStarptbin(Int_t nVars,TString *varNames,Bool_t *isUpperCut);


  Int_t GetMinITSNclsD0Pion(){return fMinITSNclsD0Pion;}
  Int_t GetMinTPCNclsD0Pion(){return fMinTPCNclsD0Pion;}
  Bool_t UseITSRefitD0Pion(){return fUseITSRefitD0Pion;}
  Bool_t UseTPCRefitD0Pion(){return fUseTPCRefitD0Pion;}
  Bool_t UseFilterBitD0Pion(){return fUseFilterBitD0Pion;}
  Int_t GetFilterBitD0Pion(){return fFilterBitD0Pion;}
  Double_t GetMinPtD0Pion(){return fMinPtD0Pion;}

  void SetMinITSNclsD0Pion(Int_t value){fMinITSNclsD0Pion = value;}
  void SetMinTPCNclsD0Pion(Int_t value){fMinTPCNclsD0Pion = value;}
  void SetUseITSRefitD0Pion(Bool_t option){fUseITSRefitD0Pion = option;}
  void SetUseTPCRefitD0Pion(Bool_t option){fUseTPCRefitD0Pion = option;}
  void SetUseFilterBitD0Pion(Bool_t option){fUseFilterBitD0Pion = option;}
  void SetFilterBitD0Pion(Int_t value){fFilterBitD0Pion = value;}
  void SetMinPtD0Pion(Double_t value){fMinPtD0Pion = value;}

  Int_t GetMinITSNclsD0Kaon(){return fMinITSNclsD0Kaon;}
  Int_t GetMinTPCNclsD0Kaon(){return fMinTPCNclsD0Kaon;}
  Bool_t UseITSRefitD0Kaon(){return fUseITSRefitD0Kaon;}
  Bool_t UseTPCRefitD0Kaon(){return fUseTPCRefitD0Kaon;}
  Bool_t UseFilterBitD0Kaon(){return fUseFilterBitD0Kaon;}
  Int_t GetFilterBitD0Kaon(){return fFilterBitD0Kaon;}
  Double_t GetMinPtD0Kaon(){return fMinPtD0Kaon;}

  void SetMinITSNclsD0Kaon(Int_t value){fMinITSNclsD0Kaon = value;}
  void SetMinTPCNclsD0Kaon(Int_t value){fMinTPCNclsD0Kaon = value;}
  void SetUseITSRefitD0Kaon(Bool_t option){fUseITSRefitD0Kaon = option;}
  void SetUseTPCRefitD0Kaon(Bool_t option){fUseTPCRefitD0Kaon = option;}
  void SetUseFilterBitD0Kaon(Bool_t option){fUseFilterBitD0Kaon = option;}
  void SetFilterBitD0Kaon(Int_t value){fFilterBitD0Kaon = value;}
  void SetMinPtD0Kaon(Double_t value){fMinPtD0Kaon = value;}

  Int_t GetMinITSNclsDStarPion(){return fMinITSNclsDStarPion;}
  Int_t GetMinTPCNclsDStarPion(){return fMinTPCNclsDStarPion;}
  Bool_t UseITSRefitDStarPion(){return fUseITSRefitDStarPion;}
  Bool_t UseTPCRefitDStarPion(){return fUseTPCRefitDStarPion;}
  Bool_t UseFilterBitDStarPion(){return fUseFilterBitDStarPion;}
  Int_t GetFilterBitDStarPion(){return fFilterBitDStarPion;}
  Double_t GetMinPtDStarPion(){return fMinPtDStarPion;}

  void SetMinITSNclsDStarPion(Int_t value){fMinITSNclsDStarPion = value;}
  void SetMinTPCNclsDStarPion(Int_t value){fMinTPCNclsDStarPion = value;}
  void SetUseITSRefitDStarPion(Bool_t option){fUseITSRefitDStarPion = option;}
  void SetUseTPCRefitDStarPion(Bool_t option){fUseTPCRefitDStarPion = option;}
  void SetUseFilterBitDStarPion(Bool_t option){fUseFilterBitDStarPion = option;}
  void SetFilterBitDStarPion(Int_t value){fFilterBitDStarPion = value;}
  void SetMinPtDStarPion(Double_t value){fMinPtDStarPion = value;}

  Int_t GetMinITSNclsB0Pion(){return fMinITSNclsB0Pion;}
  Int_t GetMinTPCNclsB0Pion(){return fMinTPCNclsB0Pion;}
  Bool_t UseITSRefitB0Pion(){return fUseITSRefitB0Pion;}
  Bool_t UseTPCRefitB0Pion(){return fUseTPCRefitB0Pion;}
  Bool_t UseFilterBitB0Pion(){return fUseFilterBitB0Pion;}
  Int_t GetFilterBitB0Pion(){return fFilterBitB0Pion;}
  Double_t GetMinPtB0Pion(){return fMinPtB0Pion;}

  void SetMinITSNclsB0Pion(Int_t value){fMinITSNclsB0Pion = value;}
  void SetMinTPCNclsB0Pion(Int_t value){fMinTPCNclsB0Pion = value;}
  void SetUseITSRefitB0Pion(Bool_t option){fUseITSRefitB0Pion = option;}
  void SetUseTPCRefitB0Pion(Bool_t option){fUseTPCRefitB0Pion = option;}
  void SetUseFilterBitB0Pion(Bool_t option){fUseFilterBitB0Pion = option;}
  void SetFilterBitB0Pion(Int_t value){fFilterBitB0Pion = value;}
  void SetMinPtB0Pion(Double_t value){fMinPtB0Pion = value;}

  void SetCut(Int_t nCutIndex, Int_t ptBin, AliRDHFCutsB0toDStarPi::EUpperCut cutDirection, Float_t cutValue);
  void SetCutD0forD0ptbin(Int_t nCutIndex, Int_t ptBin, AliRDHFCutsB0toDStarPi::EUpperCut cutDirection, Float_t cutValue);
  void SetCutD0forDStarptbin(Int_t nCutIndex, Int_t ptBin, AliRDHFCutsB0toDStarPi::EUpperCut cutDirection, Float_t cutValue);
  void SetCutDStarforDStarptbin(Int_t nCutIndex, Int_t ptBin, AliRDHFCutsB0toDStarPi::EUpperCut cutDirection, Float_t cutValue);

 protected:

  Float_t fMaxPtPid; 
  Float_t fTPCflag;   
  Double_t fCircRadius; /// Radius for circular PID nsigma cut
  Bool_t fGetCutInfo;

  Bool_t * fIsCutUsed;                       //[fGlobalIndex]

  Int_t fnVarsD0forD0ptbin;
  Int_t fnPtBinsD0forD0ptbin;
  Int_t fGlobalIndexD0forD0ptbin;
  Float_t * fCutsRDD0forD0ptbin;             //[fGlobalIndexD0forD0ptbin]
  Int_t fnPtBinLimitsD0forD0ptbin;
  Float_t * fPtBinLimitsD0forD0ptbin;        //[fnPtBinsD0forD0ptbin]
  Bool_t * fIsUpperCutD0forD0ptbin;          //[fnVarsD0forD0ptbin]
  Bool_t * fIsCutUsedD0forD0ptbin;           //[fGlobalIndexD0forD0ptbin]
  TString * fVarNamesD0forD0ptbin;           //[fnVarsD0forD0ptbin]

  Int_t fnVarsD0forDStarptbin;
  Int_t fnPtBinsD0forDStarptbin;
  Int_t fGlobalIndexD0forDStarptbin;
  Float_t * fCutsRDD0forDStarptbin;          //[fGlobalIndexD0forDStarptbin]
  Int_t fnPtBinLimitsD0forDStarptbin;
  Float_t * fPtBinLimitsD0forDStarptbin;     //[fnPtBinsD0forDStarptbin]
  Bool_t * fIsUpperCutD0forDStarptbin;       //[fnVarsD0forDStarptbin]
  Bool_t * fIsCutUsedD0forDStarptbin;        //[fGlobalIndexD0forDStarptbin]
  TString * fVarNamesD0forDStarptbin;        //[fnVarsD0forDStarptbin]

  Int_t fnVarsDStarforDStarptbin;
  Int_t fnPtBinsDStarforDStarptbin;
  Int_t fGlobalIndexDStarforDStarptbin;
  Float_t * fCutsRDDStarforDStarptbin;       //[fGlobalIndexDStarforDStarptbin]
  Int_t fnPtBinLimitsDStarforDStarptbin;
  Float_t * fPtBinLimitsDStarforDStarptbin;  //[fnPtBinsDStarforDStarptbin]
  Bool_t * fIsUpperCutDStarforDStarptbin;    //[fnVarsDStarforDStarptbin]
  Bool_t * fIsCutUsedDStarforDStarptbin;     //[fGlobalIndexDStarforDStarptbin]
  TString * fVarNamesDStarforDStarptbin;     //[fnVarsDStarforDStarptbin]

  Int_t fMinITSNclsD0Pion;
  Int_t fMinTPCNclsD0Pion;
  Bool_t fUseITSRefitD0Pion;
  Bool_t fUseTPCRefitD0Pion;
  Bool_t fUseFilterBitD0Pion;
  Int_t fFilterBitD0Pion;
  Double_t fMinPtD0Pion;

  Int_t fMinITSNclsD0Kaon;
  Int_t fMinTPCNclsD0Kaon;
  Bool_t fUseITSRefitD0Kaon;
  Bool_t fUseTPCRefitD0Kaon;
  Bool_t fUseFilterBitD0Kaon;
  Int_t fFilterBitD0Kaon;
  Double_t fMinPtD0Kaon;

  Int_t fMinITSNclsDStarPion;
  Int_t fMinTPCNclsDStarPion;
  Bool_t fUseITSRefitDStarPion;
  Bool_t fUseTPCRefitDStarPion;
  Bool_t fUseFilterBitDStarPion;
  Int_t fFilterBitDStarPion;
  Double_t fMinPtDStarPion;

  Int_t fMinITSNclsB0Pion;
  Int_t fMinTPCNclsB0Pion;
  Bool_t fUseITSRefitB0Pion;
  Bool_t fUseTPCRefitB0Pion;
  Bool_t fUseFilterBitB0Pion;
  Int_t fFilterBitB0Pion;
  Double_t fMinPtB0Pion;

  Double_t fMagneticField = 0.0;

  /// \cond CLASSIMP    
  ClassDef(AliRDHFCutsB0toDStarPi,1)
  /// \endcond
};

#endif
