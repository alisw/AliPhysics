#ifndef AlIRDHFCUTSB0TODSTARPI_H
#define AlIRDHFCUTSB0TODSTARPI_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///***********************************************************
/// \class Class AliRDHFCutsB0toDStarPi
/// \brief class for cuts on AOD reconstructed B0->DStarPi->D0PiPi->KPiPiPi
/// \author Author: Lennart van Doremalen, l.v.r.vandoremalen@uu.nl  
///
//
//                 Author Lennart van Doremalen
//           Utrecht University - l.v.r.vandoremalen@uu.nl
//
//     Several AliPhysics classes have been used as a basis for this code
//
///***********************************************************

#include "AliRDHFCuts.h"

class AliAODEvent;
class AliAODRecoDecayHF2Prong;
class AliAODRecoDecayHF;

class AliRDHFCutsB0toDStarPi : public AliRDHFCuts 
{
 public:

  enum EUpperCut {kCutBelowValue             = 0,
                  kCutAboveValue             = 1
                 };

  AliRDHFCutsB0toDStarPi(const char* name="B0toDStarPiCuts");
  
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

  Double_t DeltaInvMassDStarKpipi(AliAODRecoDecayHF2Prong *DStar) const;
  Double_t DeltaInvMassB0Kpipipi(AliAODRecoDecayHF2Prong *Bzero) const;

  void SetCutsD0forD0ptbin(Int_t nVars,Int_t nPtBins,Float_t **cutsRDD0forD0ptbin);
  void SetCutsD0forD0ptbin(Int_t glIndex,Float_t *cutsRDD0forD0ptbin);
  Int_t PtBinD0forD0ptbin(Double_t pt) const;
  void SetPtBinsD0forD0ptbin(Int_t nPtBinLimits,Float_t *ptBinLimits);

  void SetCutsD0forDStarptbin(Int_t nVars,Int_t nPtBins,Float_t **cutsRDD0forDStarptbin);
  void SetCutsD0forDStarptbin(Int_t glIndex,Float_t *cutsRDD0forDStarptbin);
  Int_t PtBinD0forDStarptbin(Double_t pt) const;
  void SetPtBinsD0forDStarptbin(Int_t nPtBinLimits,Float_t *ptBinLimits);

  void SetCutsDStarforDStarptbin(Int_t nVars,Int_t nPtBins,Float_t **cutsRDDStarforDStarptbin);
  void SetCutsDStarforDStarptbin(Int_t glIndex,Float_t *cutsRDDStarforDStarptbin);
  Int_t PtBinDStarforDStarptbin(Double_t pt) const;
  void SetPtBinsDStarforDStarptbin(Int_t nPtBinLimits,Float_t *ptBinLimits);

  Float_t *GetPtBinLimitsD0forD0ptbin() const {return fPtBinLimitsD0forD0ptbin;}
  Int_t GetNPtBinsD0forD0ptbin() const {return fnPtBinsD0forD0ptbin;}
  Int_t GetNVarsD0forD0ptbin() const {return fnVarsD0forD0ptbin;}
  Int_t GetGlobalIndexD0forD0ptbin(Int_t iVar,Int_t iPtBin) const{return iPtBin*fnVarsD0forD0ptbin+iVar;}
  void  SetGlobalIndexD0forD0ptbin(){fGlobalIndexD0forD0ptbin=fnVarsD0forD0ptbin*fnPtBinsD0forD0ptbin; return;}
  void  SetNPtBinsD0forD0ptbin(Int_t nptBins){fnPtBinsD0forD0ptbin=nptBins; return;}
  void  SetNVarsD0forD0ptbin(Int_t nVars){fnVarsD0forD0ptbin=nVars; return;} 

  Float_t *GetPtBinLimitsD0forDStarptbin() const {return fPtBinLimitsD0forDStarptbin;}
  Int_t GetNPtBinsD0forDStarptbin() const {return fnPtBinsD0forDStarptbin;}
  Int_t GetNVarsD0forDStarptbin() const {return fnVarsD0forDStarptbin;}
  Int_t GetGlobalIndexD0forDStarptbin(Int_t iVar,Int_t iPtBin) const{return iPtBin*fnVarsD0forDStarptbin+iVar;}
  void  SetGlobalIndexD0forDStarptbin(){fGlobalIndexD0forDStarptbin=fnVarsD0forDStarptbin*fnPtBinsD0forDStarptbin; return;}
  void  SetNPtBinsD0forDStarptbin(Int_t nptBins){fnPtBinsD0forDStarptbin=nptBins; return;}
  void  SetNVarsD0forDStarptbin(Int_t nVars){fnVarsD0forDStarptbin=nVars; return;} 

  Float_t *GetPtBinLimitsDStarforDStarptbin() const {return fPtBinLimitsDStarforDStarptbin;}
  Int_t GetNPtBinsDStarforDStarptbin() const {return fnPtBinsDStarforDStarptbin;}
  Int_t GetNVarsDStarforDStarptbin() const {return fnVarsDStarforDStarptbin;}
  Int_t GetGlobalIndexDStarforDStarptbin(Int_t iVar,Int_t iPtBin) const{return iPtBin*fnVarsDStarforDStarptbin+iVar;}
  void  SetGlobalIndexDStarforDStarptbin(){fGlobalIndexDStarforDStarptbin=fnVarsDStarforDStarptbin*fnPtBinsDStarforDStarptbin; return;}
  void  SetNPtBinsDStarforDStarptbin(Int_t nptBins){fnPtBinsDStarforDStarptbin=nptBins; return;}
  void  SetNVarsDStarforDStarptbin(Int_t nVars){fnVarsDStarforDStarptbin=nVars; return;} 

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

  Bool_t GetIsCutUsed(Int_t nCutIndex, Int_t ptbin) const {return fIsCutUsed[GetGlobalIndex(nCutIndex,ptbin)];}
  Bool_t GetIsCutUsedD0forD0ptbin(Int_t nCutIndex, Int_t ptbin) const {return fIsCutUsedD0forD0ptbin[GetGlobalIndexD0forD0ptbin(nCutIndex,ptbin)];}
  Bool_t GetIsCutUsedD0forDStarptbin(Int_t nCutIndex, Int_t ptbin) const {return fIsCutUsedD0forDStarptbin[GetGlobalIndexD0forDStarptbin(nCutIndex,ptbin)];}
  Bool_t GetIsCutUsedDStarforDStarptbin(Int_t nCutIndex, Int_t ptbin) const {return fIsCutUsedDStarforDStarptbin[GetGlobalIndexDStarforDStarptbin(nCutIndex,ptbin)];}

  Int_t ApplyCutOnVariable(Int_t nCutIndex, Int_t ptbin, Float_t cutVariableValue, Bool_t bCutArray[97]);
  Int_t ApplyCutOnVariableD0forD0ptbin(Int_t nCutIndex, Int_t ptbin, Float_t cutVariableValue, Bool_t bCutArray[29]);
  Int_t ApplyCutOnVariableD0forDStarptbin(Int_t nCutIndex, Int_t ptbin, Float_t cutVariableValue, Bool_t bCutArray[39]);
  Int_t ApplyCutOnVariableDStarforDStarptbin(Int_t nCutIndex, Int_t ptbin, Float_t cutVariableValue, Bool_t bCutArray[29]);

  void SetVarNamesD0forD0ptbin(Int_t nVars,TString *varNames,Bool_t *isUpperCut);
  void SetVarNamesD0forDStarptbin(Int_t nVars,TString *varNames,Bool_t *isUpperCut);
  void SetVarNamesDStarforDStarptbin(Int_t nVars,TString *varNames,Bool_t *isUpperCut);


  Int_t GetMinITSNclsD0FirstDaughter(){return fMinITSNclsD0FirstDaughter;}
  Int_t GetMinTPCNclsD0FirstDaughter(){return fMinTPCNclsD0FirstDaughter;}
  Bool_t UseITSRefitD0FirstDaughter(){return fUseITSRefitD0FirstDaughter;}
  Bool_t UseTPCRefitD0FirstDaughter(){return fUseTPCRefitD0FirstDaughter;}
  Bool_t UseFilterBitD0FirstDaughter(){return fUseFilterBitD0FirstDaughter;}
  Int_t GetFilterBitD0FirstDaughter(){return fFilterBitD0FirstDaughter;}
  Double_t GetMinPtD0FirstDaughter(){return fMinPtD0FirstDaughter;}
  Double_t GetMaxAbsEtaD0FirstDaughter(){return fMaxAbsEtaD0FirstDaughter;}
  void GetHardSelectionArrayITSD0FirstDaughter(Bool_t array[7] = 0){for(Int_t i=0;i<7;i++){array[i] = fHardSelectionArrayITSD0FirstDaughter[i];} return;}
  void GetSoftSelectionArrayITSD0FirstDaughter(Bool_t array[7] = 0){for(Int_t i=0;i<7;i++){array[i] = fSoftSelectionArrayITSD0FirstDaughter[i];} return;}
  Int_t GetNSoftITSCutD0FirstDaughter(){return fNSoftITSCutD0FirstDaughter;}

  void SetMinITSNclsD0FirstDaughter(Int_t value){fMinITSNclsD0FirstDaughter = value; return;}
  void SetMinTPCNclsD0FirstDaughter(Int_t value){fMinTPCNclsD0FirstDaughter = value; return;}
  void SetUseITSRefitD0FirstDaughter(Bool_t option){fUseITSRefitD0FirstDaughter = option; return;}
  void SetUseTPCRefitD0FirstDaughter(Bool_t option){fUseTPCRefitD0FirstDaughter = option; return;}
  void SetUseFilterBitD0FirstDaughter(Bool_t option){fUseFilterBitD0FirstDaughter = option; return;}
  void SetFilterBitD0FirstDaughter(Int_t value){fFilterBitD0FirstDaughter = value; return;}
  void SetMinPtD0FirstDaughter(Double_t value){fMinPtD0FirstDaughter = value; return;}
  void SetMaxAbsEtaD0FirstDaughter(Double_t value){fMaxAbsEtaD0FirstDaughter = value;}
  void SetHardSelectionArrayITSD0FirstDaughter(const Bool_t array[7] = 0){for(Int_t i=0;i<7;i++){fHardSelectionArrayITSD0FirstDaughter[i] = array[i];} return;}
  void SetSoftSelectionArrayITSD0FirstDaughter(const Bool_t array[7] = 0){for(Int_t i=0;i<7;i++){fSoftSelectionArrayITSD0FirstDaughter[i] = array[i];} return;}
  void SetNSoftITSCutD0FirstDaughter(Int_t value){fNSoftITSCutD0FirstDaughter = value;}

  Int_t GetMinITSNclsD0SecondDaughter(){return fMinITSNclsD0SecondDaughter;}
  Int_t GetMinTPCNclsD0SecondDaughter(){return fMinTPCNclsD0SecondDaughter;}
  Bool_t UseITSRefitD0SecondDaughter(){return fUseITSRefitD0SecondDaughter;}
  Bool_t UseTPCRefitD0SecondDaughter(){return fUseTPCRefitD0SecondDaughter;}
  Bool_t UseFilterBitD0SecondDaughter(){return fUseFilterBitD0SecondDaughter;}
  Int_t GetFilterBitD0SecondDaughter(){return fFilterBitD0SecondDaughter;}
  Double_t GetMinPtD0SecondDaughter(){return fMinPtD0SecondDaughter;}
  Double_t GetMaxAbsEtaD0SecondDaughter(){return fMaxAbsEtaD0SecondDaughter;}
  void GetHardSelectionArrayITSD0SecondDaughter(Bool_t array[7] = 0){for(Int_t i=0;i<7;i++){array[i] = fHardSelectionArrayITSD0SecondDaughter[i];} return;}
  void GetSoftSelectionArrayITSD0SecondDaughter(Bool_t array[7] = 0){for(Int_t i=0;i<7;i++){array[i] = fSoftSelectionArrayITSD0SecondDaughter[i];} return;}
  Int_t GetNSoftITSCutD0SecondDaughter(){return fNSoftITSCutD0SecondDaughter;}

  void SetMinITSNclsD0SecondDaughter(Int_t value){fMinITSNclsD0SecondDaughter = value; return;}
  void SetMinTPCNclsD0SecondDaughter(Int_t value){fMinTPCNclsD0SecondDaughter = value; return;}
  void SetUseITSRefitD0SecondDaughter(Bool_t option){fUseITSRefitD0SecondDaughter = option; return;}
  void SetUseTPCRefitD0SecondDaughter(Bool_t option){fUseTPCRefitD0SecondDaughter = option; return;}
  void SetUseFilterBitD0SecondDaughter(Bool_t option){fUseFilterBitD0SecondDaughter = option; return;}
  void SetFilterBitD0SecondDaughter(Int_t value){fFilterBitD0SecondDaughter = value; return;}
  void SetMinPtD0SecondDaughter(Double_t value){fMinPtD0SecondDaughter = value; return;}
  void SetMaxAbsEtaD0SecondDaughter(Double_t value){fMaxAbsEtaD0SecondDaughter = value; return;}
  void SetHardSelectionArrayITSD0SecondDaughter(const Bool_t array[7] = 0){for(Int_t i=0;i<7;i++){fHardSelectionArrayITSD0SecondDaughter[i] = array[i];} return;}
  void SetSoftSelectionArrayITSD0SecondDaughter(const Bool_t array[7] = 0){for(Int_t i=0;i<7;i++){fSoftSelectionArrayITSD0SecondDaughter[i] = array[i];} return;}
  void SetNSoftITSCutD0SecondDaughter(Int_t value){fNSoftITSCutD0SecondDaughter = value;}

  Int_t GetMinITSNclsDStarPion(){return fMinITSNclsDStarPion;}
  Int_t GetMinTPCNclsDStarPion(){return fMinTPCNclsDStarPion;}
  Bool_t UseITSRefitDStarPion(){return fUseITSRefitDStarPion;}
  Bool_t UseTPCRefitDStarPion(){return fUseTPCRefitDStarPion;}
  Bool_t UseFilterBitDStarPion(){return fUseFilterBitDStarPion;}
  Int_t GetFilterBitDStarPion(){return fFilterBitDStarPion;}
  Double_t GetMinPtDStarPion(){return fMinPtDStarPion;}
  Double_t GetMaxAbsEtaDStarPion(){return fMaxAbsEtaDStarPion;}
  void GetHardSelectionArrayITSDStarPion(Bool_t array[7] = 0){for(Int_t i=0;i<7;i++){array[i] = fHardSelectionArrayITSDStarPion[i];} return;}
  void GetSoftSelectionArrayITSDStarPion(Bool_t array[7] = 0){for(Int_t i=0;i<7;i++){array[i] = fSoftSelectionArrayITSDStarPion[i];} return;}
  Int_t GetNSoftITSCutDStarPion(){return fNSoftITSCutDStarPion;}

  void SetMinITSNclsDStarPion(Int_t value){fMinITSNclsDStarPion = value; return;}
  void SetMinTPCNclsDStarPion(Int_t value){fMinTPCNclsDStarPion = value; return;}
  void SetUseITSRefitDStarPion(Bool_t option){fUseITSRefitDStarPion = option; return;}
  void SetUseTPCRefitDStarPion(Bool_t option){fUseTPCRefitDStarPion = option; return;}
  void SetUseFilterBitDStarPion(Bool_t option){fUseFilterBitDStarPion = option; return;}
  void SetFilterBitDStarPion(Int_t value){fFilterBitDStarPion = value; return;}
  void SetMinPtDStarPion(Double_t value){fMinPtDStarPion = value; return;}
  void SetMaxAbsEtaDStarPion(Double_t value){fMaxAbsEtaDStarPion = value; return;}
  void SetHardSelectionArrayITSDStarPion(const Bool_t array[7] = 0){for(Int_t i=0;i<7;i++){fHardSelectionArrayITSDStarPion[i] = array[i];} return;}
  void SetSoftSelectionArrayITSDStarPion(const Bool_t array[7] = 0){for(Int_t i=0;i<7;i++){fSoftSelectionArrayITSDStarPion[i] = array[i];} return;}
  void SetNSoftITSCutDStarPion(Int_t value){fNSoftITSCutDStarPion = value; return;}

  Int_t GetMinITSNclsB0Pion(){return fMinITSNclsB0Pion;}
  Int_t GetMinTPCNclsB0Pion(){return fMinTPCNclsB0Pion;}
  Bool_t UseITSRefitB0Pion(){return fUseITSRefitB0Pion;}
  Bool_t UseTPCRefitB0Pion(){return fUseTPCRefitB0Pion;}
  Bool_t UseFilterBitB0Pion(){return fUseFilterBitB0Pion;}
  Int_t GetFilterBitB0Pion(){return fFilterBitB0Pion;}
  Double_t GetMinPtB0Pion(){return fMinPtB0Pion;}
  Double_t GetMaxAbsEtaB0Pion(){return fMaxAbsEtaB0Pion;}
  void GetHardSelectionArrayITSB0Pion(Bool_t array[7] = 0){for(Int_t i=0;i<7;i++){array[i] = fHardSelectionArrayITSB0Pion[i];} return;}
  void GetSoftSelectionArrayITSB0Pion(Bool_t array[7] = 0){for(Int_t i=0;i<7;i++){array[i] = fSoftSelectionArrayITSB0Pion[i];} return;}
  Int_t GetNSoftITSCutB0Pion(){return fNSoftITSCutB0Pion;}

  void SetMinITSNclsB0Pion(Int_t value){fMinITSNclsB0Pion = value; return;}
  void SetMinTPCNclsB0Pion(Int_t value){fMinTPCNclsB0Pion = value; return;}
  void SetUseITSRefitB0Pion(Bool_t option){fUseITSRefitB0Pion = option; return;}
  void SetUseTPCRefitB0Pion(Bool_t option){fUseTPCRefitB0Pion = option; return;}
  void SetUseFilterBitB0Pion(Bool_t option){fUseFilterBitB0Pion = option; return;}
  void SetFilterBitB0Pion(Int_t value){fFilterBitB0Pion = value; return;}
  void SetMinPtB0Pion(Double_t value){fMinPtB0Pion = value; return;}
  void SetMaxAbsEtaB0Pion(Double_t value){fMaxAbsEtaB0Pion = value; return;}
  void SetHardSelectionArrayITSB0Pion(const Bool_t array[7] = 0){for(Int_t i=0;i<7;i++){fHardSelectionArrayITSB0Pion[i] = array[i];}}
  void SetSoftSelectionArrayITSB0Pion(const Bool_t array[7] = 0){for(Int_t i=0;i<7;i++){fSoftSelectionArrayITSB0Pion[i] = array[i];}}
  void SetNSoftITSCutB0Pion(Int_t value){fNSoftITSCutB0Pion = value; return;}

  void SetCut(Int_t nCutIndex, Int_t ptBin, AliRDHFCutsB0toDStarPi::EUpperCut cutDirection, Float_t cutValue);
  void SetCutD0forD0ptbin(Int_t nCutIndex, Int_t ptBin, AliRDHFCutsB0toDStarPi::EUpperCut cutDirection, Float_t cutValue);
  void SetCutD0forDStarptbin(Int_t nCutIndex, Int_t ptBin, AliRDHFCutsB0toDStarPi::EUpperCut cutDirection, Float_t cutValue);
  void SetCutDStarforDStarptbin(Int_t nCutIndex, Int_t ptBin, AliRDHFCutsB0toDStarPi::EUpperCut cutDirection, Float_t cutValue);

  Double_t GetMaxDCADStarPionD0(){return fMaxDCADStarPionD0;}
  Double_t GetMaxDCADStarPionB0Pion(){return fMaxDCADStarPionB0Pion;}
  Double_t GetMaxDCAB0PionD0(){return fMaxDCAB0PionD0;}
  Double_t GetMaxDCACombined(){return fMaxDCACombined;}

  void SetMaxDCADStarPionD0(Double_t value){fMaxDCADStarPionD0 = value; return;}
  void SetMaxDCADStarPionB0Pion(Double_t value){fMaxDCADStarPionB0Pion = value; return;}
  void SetMaxDCAB0PionD0(Double_t value){fMaxDCAB0PionD0 = value; return;}
  void SetMaxDCACombined(Double_t value){fMaxDCACombined = value; return;}

  Double_t GetMind0D0FirstDaughter(){return fMind0D0FirstDaughter;}
  Double_t GetMind0D0SecondDaughter(){return fMind0D0SecondDaughter;}
  Double_t GetMind0DStarPion(){return fMind0DStarPion;}
  Double_t GetMind0B0Pion(){return fMind0B0Pion;}

  void SetMind0D0FirstDaughter(Double_t value){fMind0D0FirstDaughter = value; return;}
  void SetMind0D0SecondDaughter(Double_t value){fMind0D0SecondDaughter = value; return;}
  void SetMind0DStarPion(Double_t value){fMind0DStarPion = value; return;}
  void SetMind0B0Pion(Double_t value){fMind0B0Pion = value; return;}

  Double_t GetMaxPtDStarPion(){return fMaxPtDStarPion;}
  void SetMaxPtDStarPion(Double_t value){fMaxPtDStarPion = value; return;}

 protected:

  Float_t fMaxPtPid;                                  ///
  Float_t fTPCflag;                                   ///
  Double_t fCircRadius;                               /// Radius for circular PID nsigma cut
  Bool_t fGetCutInfo;                                 ///

  Bool_t * fIsCutUsed;                                //[fGlobalIndex]

  Int_t fnVarsD0forD0ptbin;                           ///
  Int_t fnPtBinsD0forD0ptbin;                         ///
  Int_t fGlobalIndexD0forD0ptbin;                     ///
  Float_t * fCutsRDD0forD0ptbin;                      //[fGlobalIndexD0forD0ptbin]
  Int_t fnPtBinLimitsD0forD0ptbin;                    ///
  Float_t * fPtBinLimitsD0forD0ptbin;                 //[fnPtBinLimitsD0forD0ptbin]
  Bool_t * fIsUpperCutD0forD0ptbin;                   //[fnVarsD0forD0ptbin]
  Bool_t * fIsCutUsedD0forD0ptbin;                    //[fGlobalIndexD0forD0ptbin]
  TString * fVarNamesD0forD0ptbin;                    //[fnVarsD0forD0ptbin]

  Int_t fnVarsD0forDStarptbin;                        ///
  Int_t fnPtBinsD0forDStarptbin;                      ///
  Int_t fGlobalIndexD0forDStarptbin;                  ///
  Float_t * fCutsRDD0forDStarptbin;                   //[fGlobalIndexD0forDStarptbin]
  Int_t fnPtBinLimitsD0forDStarptbin;                 ///
  Float_t * fPtBinLimitsD0forDStarptbin;              //[fnPtBinLimitsD0forDStarptbin]
  Bool_t * fIsUpperCutD0forDStarptbin;                //[fnVarsD0forDStarptbin]
  Bool_t * fIsCutUsedD0forDStarptbin;                 //[fGlobalIndexD0forDStarptbin]
  TString * fVarNamesD0forDStarptbin;                 //[fnVarsD0forDStarptbin]

  Int_t fnVarsDStarforDStarptbin;                     ///
  Int_t fnPtBinsDStarforDStarptbin;                   ///
  Int_t fGlobalIndexDStarforDStarptbin;               ///
  Float_t * fCutsRDDStarforDStarptbin;                //[fGlobalIndexDStarforDStarptbin]
  Int_t fnPtBinLimitsDStarforDStarptbin;              ///
  Float_t * fPtBinLimitsDStarforDStarptbin;           //[fnPtBinLimitsDStarforDStarptbin]
  Bool_t * fIsUpperCutDStarforDStarptbin;             //[fnVarsDStarforDStarptbin]
  Bool_t * fIsCutUsedDStarforDStarptbin;              //[fGlobalIndexDStarforDStarptbin]
  TString * fVarNamesDStarforDStarptbin;              //[fnVarsDStarforDStarptbin]

  Int_t fMinITSNclsD0FirstDaughter;                   ///
  Int_t fMinTPCNclsD0FirstDaughter;                   ///
  Bool_t fUseITSRefitD0FirstDaughter;                 ///
  Bool_t fUseTPCRefitD0FirstDaughter;                 ///
  Bool_t fUseFilterBitD0FirstDaughter;                ///
  Int_t fFilterBitD0FirstDaughter;                    ///
  Double_t fMinPtD0FirstDaughter;                     ///
  Double_t fMaxAbsEtaD0FirstDaughter;                 ///
  Bool_t fHardSelectionArrayITSD0FirstDaughter[7];    ///
  Bool_t fSoftSelectionArrayITSD0FirstDaughter[7];    ///
  Int_t fNSoftITSCutD0FirstDaughter;                  ///

  Int_t fMinITSNclsD0SecondDaughter;                  ///
  Int_t fMinTPCNclsD0SecondDaughter;                  ///
  Bool_t fUseITSRefitD0SecondDaughter;                ///
  Bool_t fUseTPCRefitD0SecondDaughter;                ///
  Bool_t fUseFilterBitD0SecondDaughter;               ///
  Int_t fFilterBitD0SecondDaughter;                   ///
  Double_t fMinPtD0SecondDaughter;                    ///
  Double_t fMaxAbsEtaD0SecondDaughter;                ///
  Bool_t fHardSelectionArrayITSD0SecondDaughter[7];   ///
  Bool_t fSoftSelectionArrayITSD0SecondDaughter[7];   ///
  Int_t fNSoftITSCutD0SecondDaughter;                 ///

  Int_t fMinITSNclsDStarPion;                         ///
  Int_t fMinTPCNclsDStarPion;                         ///
  Bool_t fUseITSRefitDStarPion;                       ///
  Bool_t fUseTPCRefitDStarPion;                       ///
  Bool_t fUseFilterBitDStarPion;                      ///
  Int_t fFilterBitDStarPion;                          ///
  Double_t fMinPtDStarPion;                           ///
  Double_t fMaxAbsEtaDStarPion;                       ///
  Bool_t fHardSelectionArrayITSDStarPion[7];          ///
  Bool_t fSoftSelectionArrayITSDStarPion[7];          ///
  Int_t fNSoftITSCutDStarPion;                        ///

  Int_t fMinITSNclsB0Pion;                            ///
  Int_t fMinTPCNclsB0Pion;                            ///
  Bool_t fUseITSRefitB0Pion;                          ///
  Bool_t fUseTPCRefitB0Pion;                          ///
  Bool_t fUseFilterBitB0Pion;                         ///
  Int_t fFilterBitB0Pion;                             ///
  Double_t fMinPtB0Pion;                              ///
  Double_t fMaxAbsEtaB0Pion;                          ///
  Bool_t fHardSelectionArrayITSB0Pion[7];             ///
  Bool_t fSoftSelectionArrayITSB0Pion[7];             ///
  Int_t fNSoftITSCutB0Pion;                           ///

  Double_t fMaxDCADStarPionD0;                        ///
  Double_t fMaxDCADStarPionB0Pion;                    ///
  Double_t fMaxDCAB0PionD0;                           ///
  Double_t fMaxDCACombined;                           ///

  Double_t fMind0D0FirstDaughter;                     ///
  Double_t fMind0D0SecondDaughter;                    ///
  Double_t fMind0DStarPion;                           ///
  Double_t fMind0B0Pion;                              ///

  Double_t fMaxPtDStarPion;                           ///

  /// \cond CLASSIMP    
  ClassDef(AliRDHFCutsB0toDStarPi,8) ///
  /// \endcond
};

#endif


