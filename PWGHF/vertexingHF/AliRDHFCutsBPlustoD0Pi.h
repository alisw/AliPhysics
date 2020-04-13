#ifndef AlIRDHFCUTSBPlustoD0Pi_H
#define AlIRDHFCUTSBPlustoD0Pi_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///***********************************************************
/// \class Class AliRDHFCutsBPlustoD0Pi
/// \brief class for cuts on AOD reconstructed BPlus->D0Pi->KPiPi
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

class AliRDHFCutsBPlustoD0Pi : public AliRDHFCuts 
{
 public:

  enum EUpperCut {kCutBelowValue             = 0,
                  kCutAboveValue             = 1
                 };

  AliRDHFCutsBPlustoD0Pi(const char* name="BPlustoD0PiCuts");
  
  virtual ~AliRDHFCutsBPlustoD0Pi();

  AliRDHFCutsBPlustoD0Pi(const AliRDHFCutsBPlustoD0Pi& source);
  AliRDHFCutsBPlustoD0Pi& operator=(const AliRDHFCutsBPlustoD0Pi& source); 
 
  using AliRDHFCuts::GetCutVarsForOpt;
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters);

  //First two for MVA, last one for standard analysis (TBD: Merge them?)
  using AliRDHFCuts::IsSelected;
  virtual Int_t IsSelected(TObject* obj, Int_t selectionLevel, AliAODEvent* aod);
  virtual Int_t IsSelected(TObject* obj, Int_t selectionLevel){ return IsSelected(obj,selectionLevel,0); }
  virtual Int_t IsSelected(TObject* obj, Int_t selectionLevel, AliAODEvent* aod, Bool_t* bCutArray);

  Bool_t IsThisDaughterSelected(AliAODTrack *track, AliAODVertex *primary, const AliAODEvent* aod);

  //Last three for MVA, first two for standard analysis (TBD: Merge them?)
  Int_t IsD0FromBPlusSelected(Double_t ptBPlus, TObject* obj,Int_t selectionLevel, AliAODEvent* aod, Bool_t bCutArray[75]);
  Int_t IsD0forD0ptbinSelected(TObject* obj,Int_t selectionLevel, AliAODEvent* aod, Bool_t* bCutArray);
  Int_t IsBplusPionSelectedMVA(TObject* obj,Int_t selectionLevel, AliAODEvent* aod, AliAODVertex *primaryVertex, Double_t bz);
  Int_t IsD0FromBPlusSelectedMVA(Double_t ptBPlus, TObject* obj,Int_t selectionLevel, AliAODEvent* aod, AliAODVertex *primaryVertex, Double_t bz);
  Int_t IsD0forD0ptbinSelectedMVA(TObject* obj,Int_t selectionLevel, AliAODEvent* aod, AliAODVertex *primaryVertex, Double_t bz);
  Int_t IsD0SelectedPreRecVtxMVA(AliAODRecoDecayHF2Prong* d, AliAODTrack* pion, AliAODVertex *primaryVertex, Double_t bz, Int_t selLevel);

  virtual Int_t IsSelectedPID(AliAODRecoDecayHF *rd);
  virtual Int_t SelectPID(AliAODTrack *track, Int_t type);
  virtual Bool_t IsInFiducialAcceptance(Double_t pt,Double_t y) const;
  
  void SetMaxPtPid(Float_t maxPt){fMaxPtPid = maxPt;}

  void SetOffHighPtPIDinTPC(Float_t TPCrem =999.){fTPCflag = TPCrem;}

  //Pion from Bplus not really soft, but re-using code from Dstar (hence the name)
  void AddTrackCutsSoftPi(const AliESDtrackCuts *cuts){ fTrackCutsSoftPi = new AliESDtrackCuts(*cuts); return;}
  virtual AliESDtrackCuts *GetTrackCutsSoftPi() const {return fTrackCutsSoftPi;}

  void InitializeCuts();
  void InitializeCutsForCutOptimization(Int_t nCutsForOptimization, Int_t nVariables);
  void SetCutsForCutOptimization(Int_t glIndex,Float_t *cutsRDForCutOptimization);

  Double_t GetCircRadius() { return fCircRadius; }
  void SetCircRadius(Double_t radius) { fCircRadius = radius; }

  Double_t DeltaInvMassBPlusKpipi(AliAODRecoDecayHF2Prong * BPlus) const;

  void SetCutsD0forD0ptbin(Int_t nVars,Int_t nPtBins,Float_t **cutsRDD0forD0ptbin);
  void SetCutsD0forD0ptbin(Int_t glIndex,Float_t *cutsRDD0forD0ptbin);
  Int_t PtBinD0forD0ptbin(Double_t pt) const;
  void SetPtBinsD0forD0ptbin(Int_t nPtBinLimits,Float_t *ptBinLimits);

  Float_t *GetPtBinLimitsD0forD0ptbin() const {return fPtBinLimitsD0forD0ptbin;}
  Int_t GetNPtBinsD0forD0ptbin() const {return fnPtBinsD0forD0ptbin;}
  Int_t GetNVarsD0forD0ptbin() const {return fnVarsD0forD0ptbin;}
  Int_t GetGlobalIndexD0forD0ptbin(Int_t iVar,Int_t iPtBin) const{return iPtBin*fnVarsD0forD0ptbin+iVar;}
  void  SetGlobalIndexD0forD0ptbin(){fGlobalIndexD0forD0ptbin=fnVarsD0forD0ptbin*fnPtBinsD0forD0ptbin; return;}
  void  SetNPtBinsD0forD0ptbin(Int_t nptBins){fnPtBinsD0forD0ptbin=nptBins; return;}
  void  SetNVarsD0forD0ptbin(Int_t nVars){fnVarsD0forD0ptbin=nVars; return;} 

  void SetIsUpperCut(Int_t nCutIndex, Bool_t isUpperCut){fIsUpperCut[nCutIndex] = isUpperCut; return;}
  void SetIsUpperCutD0forD0ptbin(Int_t nCutIndex, Bool_t isUpperCut){fIsUpperCutD0forD0ptbin[nCutIndex] = isUpperCut; return;}

  Bool_t GetIsUpperCut(Int_t nCutIndex){return fIsUpperCut[nCutIndex];}
  Bool_t GetIsUpperCutD0forD0ptbin(Int_t nCutIndex){return fIsUpperCutD0forD0ptbin[nCutIndex];}

  void SetIsCutUsed(Int_t nCutIndex, Int_t ptbin, Bool_t isCutUsed){fIsCutUsed[GetGlobalIndex(nCutIndex,ptbin)] = isCutUsed; return;}
  void SetIsCutUsedD0forD0ptbin(Int_t nCutIndex, Int_t ptbin, Bool_t isCutUsed){fIsCutUsedD0forD0ptbin[GetGlobalIndexD0forD0ptbin(nCutIndex,ptbin)] = isCutUsed; return;}

  Bool_t GetIsCutUsed(Int_t nCutIndex, Int_t ptbin) const {return fIsCutUsed[GetGlobalIndex(nCutIndex,ptbin)];}
  Bool_t GetIsCutUsedD0forD0ptbin(Int_t nCutIndex, Int_t ptbin) const {return fIsCutUsedD0forD0ptbin[GetGlobalIndexD0forD0ptbin(nCutIndex,ptbin)];}

  //Last two for MVA, first two for standard analysis (TBD: Merge them?)
  Int_t ApplyCutOnVariable(Int_t nCutIndex, Int_t ptbin, Float_t cutVariableValue, Bool_t bCutArray[75]);
  Int_t ApplyCutOnVariableD0forD0ptbin(Int_t nCutIndex, Int_t ptbin, Float_t cutVariableValue, Bool_t bCutArray[29]);
  Int_t ApplyCutOnVariableMVA(Int_t nCutIndex, Int_t ptbin, Float_t cutVariableValue);
  Int_t ApplyCutOnVariableD0forD0ptbinMVA(Int_t nCutIndex, Int_t ptbin, Float_t cutVariableValue);

  void SetVarNamesD0forD0ptbin(Int_t nVars,TString *varNames,Bool_t *isUpperCut);

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

  Int_t GetMinITSNclsBPlusPion(){return fMinITSNclsBPlusPion;}
  Int_t GetMinTPCNclsBPlusPion(){return fMinTPCNclsBPlusPion;}
  Bool_t UseITSRefitBPlusPion(){return fUseITSRefitBPlusPion;}
  Bool_t UseTPCRefitBPlusPion(){return fUseTPCRefitBPlusPion;}
  Bool_t UseFilterBitBPlusPion(){return fUseFilterBitBPlusPion;}
  Int_t GetFilterBitBPlusPion(){return fFilterBitBPlusPion;}
  Double_t GetMinPtBPlusPion(){return fMinPtBPlusPion;}
  Double_t GetMaxAbsEtaBPlusPion(){return fMaxAbsEtaBPlusPion;}
  void GetHardSelectionArrayITSBPlusPion(Bool_t array[7] = 0){for(Int_t i=0;i<7;i++){array[i] = fHardSelectionArrayITSBPlusPion[i];} return;}
  void GetSoftSelectionArrayITSBPlusPion(Bool_t array[7] = 0){for(Int_t i=0;i<7;i++){array[i] = fSoftSelectionArrayITSBPlusPion[i];} return;}
  Int_t GetNSoftITSCutBPlusPion(){return fNSoftITSCutBPlusPion;}

  void SetMinITSNclsBPlusPion(Int_t value){fMinITSNclsBPlusPion = value; return;}
  void SetMinTPCNclsBPlusPion(Int_t value){fMinTPCNclsBPlusPion = value; return;}
  void SetUseITSRefitBPlusPion(Bool_t option){fUseITSRefitBPlusPion = option; return;}
  void SetUseTPCRefitBPlusPion(Bool_t option){fUseTPCRefitBPlusPion = option; return;}
  void SetUseFilterBitBPlusPion(Bool_t option){fUseFilterBitBPlusPion = option; return;}
  void SetFilterBitBPlusPion(Int_t value){fFilterBitBPlusPion = value; return;}
  void SetMinPtBPlusPion(Double_t value){fMinPtBPlusPion = value; return;}
  void SetMaxAbsEtaBPlusPion(Double_t value){fMaxAbsEtaBPlusPion = value; return;}
  void SetHardSelectionArrayITSBPlusPion(const Bool_t array[7] = 0){for(Int_t i=0;i<7;i++){fHardSelectionArrayITSBPlusPion[i] = array[i];}}
  void SetSoftSelectionArrayITSBPlusPion(const Bool_t array[7] = 0){for(Int_t i=0;i<7;i++){fSoftSelectionArrayITSBPlusPion[i] = array[i];}}
  void SetNSoftITSCutBPlusPion(Int_t value){fNSoftITSCutBPlusPion = value; return;}

  void SetCut(Int_t nCutIndex, Int_t ptBin, AliRDHFCutsBPlustoD0Pi::EUpperCut cutDirection, Float_t cutValue);
  void SetCutD0forD0ptbin(Int_t nCutIndex, Int_t ptBin, AliRDHFCutsBPlustoD0Pi::EUpperCut cutDirection, Float_t cutValue);
  void SetCutForCutOptimization(Int_t nCutIndex, Int_t nVariable, Int_t ptBin, AliRDHFCutsBPlustoD0Pi::EUpperCut cutDirection, Float_t * cutValues);
  Float_t GetCutForCutOptimization(Int_t nCutIndex, Int_t nVariable, Int_t ptBin){return fCutsRDForCutOptimization[GetGlobalIndexForCutOptimization(nCutIndex,nVariable,ptBin)];}

  Double_t GetMind0D0FirstDaughter(){return fMind0D0FirstDaughter;}
  Double_t GetMind0D0SecondDaughter(){return fMind0D0SecondDaughter;}
  Double_t GetMind0BPlusPion(){return fMind0BPlusPion;}
  Double_t GetFiducialYCut() const {return fFiducialYCut;}

  void SetMind0D0FirstDaughter(Double_t value){fMind0D0FirstDaughter = value; return;}
  void SetMind0D0SecondDaughter(Double_t value){fMind0D0SecondDaughter = value; return;}
  void SetMind0BPlusPion(Double_t value){fMind0BPlusPion = value; return;}
  void SetFiducialYCut(Double_t value){fFiducialYCut = value; return;}

  void SetnVariablesForCutOptimization(Double_t value){fnVariablesForCutOptimization = value; return;}
  Int_t GetnVariablesForCutOptimization(){return fnVariablesForCutOptimization;}

  void SetnCutsForOptimization(Double_t value){fnCutsForOptimization = value; return;}
  Int_t GetnCutsForOptimization(){return fnCutsForOptimization;}

  void SetGlobalIndexForCutOptimization(){fGlobalIndexCutOptimization = fnVariablesForCutOptimization*fnCutsForOptimization*fnPtBins; return;}
  Int_t GetGlobalIndexForCutOptimization(Int_t iCut, Int_t iVar,Int_t iPtBin) {return iCut + iVar * fnCutsForOptimization + iPtBin * fnCutsForOptimization * fnVariablesForCutOptimization;}

  void SetIsUpperCutForCutOptimization(Int_t nVariable, Bool_t isUpperCut){fIsUpperCutForCutOptimization[nVariable] = isUpperCut; return;}
  Bool_t GetIsUpperCutForCutOptimization(Int_t nVariable) const {return fIsUpperCutForCutOptimization[nVariable];}

  void SetCutIndexForCutOptimization(Int_t nVariable, Int_t nCutIndex){fCutIndexForCutOptimization[nVariable] = nCutIndex; return;}
  Int_t GetCutIndexForCutOptimization(Int_t nVariable) const {return fCutIndexForCutOptimization[nVariable];}

  void SetSigmaForCutOptimization(Double_t value, Int_t iPtBin){fSigmaForCutOptimization[iPtBin] = value; return;}
  Double_t GetSigmaForCutOptimization(Int_t iPtBin) const {return fSigmaForCutOptimization[iPtBin];}

  void SetNumberOfSigmaBinsForCutOptimization(Int_t nSigma){fNumberOfSigmaBinsForCutOptimization = nSigma; return;}
  Int_t GetNumberOfSigmaBinsForCutOptimization() const {return fNumberOfSigmaBinsForCutOptimization;}

 protected:

  //Pion from Bplus not really soft, but re-using code from Dstar
  AliESDtrackCuts * fTrackCutsSoftPi;                   /// cuts for pion from Bplus (AOD converted to ESD on the flight!)

  Float_t fMaxPtPid;                                  ///
  Float_t fTPCflag;                                   ///
  Double_t fCircRadius;                               /// Radius for circular PID nsigma cut

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

  Int_t fMinITSNclsBPlusPion;                         ///
  Int_t fMinTPCNclsBPlusPion;                         ///
  Bool_t fUseITSRefitBPlusPion;                       ///
  Bool_t fUseTPCRefitBPlusPion;                       ///
  Bool_t fUseFilterBitBPlusPion;                      ///
  Int_t fFilterBitBPlusPion;                          ///
  Double_t fMinPtBPlusPion;                           ///
  Double_t fMaxAbsEtaBPlusPion;                       ///
  Bool_t fHardSelectionArrayITSBPlusPion[7];          ///
  Bool_t fSoftSelectionArrayITSBPlusPion[7];          ///
  Int_t fNSoftITSCutBPlusPion;                        ///

  Double_t fMind0D0FirstDaughter;                     ///
  Double_t fMind0D0SecondDaughter;                    ///
  Double_t fMind0BPlusPion;                           ///
  Double_t fFiducialYCut;                             ///

  Int_t fnVariablesForCutOptimization;                ///
  Int_t fnCutsForOptimization;                        ///
  Int_t fGlobalIndexCutOptimization;                  ///
  Float_t * fCutsRDForCutOptimization;                //[fGlobalIndexCutOptimization]
  Bool_t * fIsUpperCutForCutOptimization;             //[fnVariablesForCutOptimization]
  Int_t * fCutIndexForCutOptimization;                //[fnVariablesForCutOptimization]
  Float_t * fSigmaForCutOptimization;                 //[fnPtBins]
  Int_t fNumberOfSigmaBinsForCutOptimization;         ///

  /// \cond CLASSIMP    
  ClassDef(AliRDHFCutsBPlustoD0Pi,5) ///
  /// \endcond
};

#endif


