#ifndef AlIRDHFCUTSB0toDPi_H
#define AlIRDHFCUTSB0toDPi_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///***********************************************************
/// \class Class AliRDHFCutsB0toDPi
/// \brief class for cuts on AOD reconstructed B0->DPlusPi->KPiPi
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
class AliAODRecoDecayHF3Prong;
class AliAODRecoDecayHF;

class AliRDHFCutsB0toDPi : public AliRDHFCuts 
{
 public:

  enum EUpperCut {kCutBelowValue             = 0,
                  kCutAboveValue             = 1
                 };

  AliRDHFCutsB0toDPi(const char* name="B0toDPiCuts");
  
  virtual ~AliRDHFCutsB0toDPi();

  AliRDHFCutsB0toDPi(const AliRDHFCutsB0toDPi& source);
  AliRDHFCutsB0toDPi& operator=(const AliRDHFCutsB0toDPi& source); 
 
  using AliRDHFCuts::GetCutVarsForOpt;
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters);

  using AliRDHFCuts::IsSelected;
  virtual Int_t IsSelected(TObject* obj, Int_t selectionLevel) {return 0;}
  virtual Int_t IsSelected(TObject* obj, Int_t selectionLevel, AliAODEvent* aod, Bool_t* bCutArray);

  Int_t IsDPlusFromB0Selected(Double_t ptB0, TObject* obj,Int_t selectionLevel, AliAODEvent* aod, Bool_t bCutArray[78]);
  Int_t IsDPlusforDPlusptbinSelected(TObject* obj,Int_t selectionLevel, AliAODEvent* aod, Bool_t* bCutArray);

  virtual Int_t IsSelectedPID(AliAODRecoDecayHF *rd);
  virtual Int_t SelectPID(AliAODTrack *track, Int_t type);
  virtual Bool_t IsInFiducialAcceptance(Double_t pt,Double_t y) const;
  
  void SetMaxPtPid(Float_t maxPt){fMaxPtPid = maxPt;}

  void SetOffHighPtPIDinTPC(Float_t TPCrem =999.){fTPCflag = TPCrem;}

  //Pion from B0 not really soft, but re-using code from Dstar (hence the name)
  // void AddTrackCutsSoftPi(const AliESDtrackCuts *cuts){ fTrackCutsSoftPi = new AliESDtrackCuts(*cuts); return;}
  // virtual AliESDtrackCuts *GetTrackCutsSoftPi() const {return fTrackCutsSoftPi;}

  void InitializeCuts();
  void InitializeCutsForCutOptimization(Int_t nCutsForOptimization, Int_t nVariables);
  void SetCutsForCutOptimization(Int_t glIndex,Float_t *cutsRDForCutOptimization);

  Double_t GetCircRadius() { return fCircRadius; }
  void SetCircRadius(Double_t radius) { fCircRadius = radius; }

  Double_t DeltaInvMassB0Kpipipi(AliAODRecoDecayHF2Prong *Bzero) const;

  void SetCutsDPlusforDPlusptbin(Int_t nVars,Int_t nPtBins,Float_t **cutsRDDPlusforDPlusptbin);
  void SetCutsDPlusforDPlusptbin(Int_t glIndex,Float_t *cutsRDDPlusforDPlusptbin);
  Int_t PtBinDPlusforDPlusptbin(Double_t pt) const;
  void SetPtBinsDPlusforDPlusptbin(Int_t nPtBinLimits,Float_t *ptBinLimits);

  Float_t *GetPtBinLimitsDPlusforDPlusptbin() const {return fPtBinLimitsDPlusforDPlusptbin;}
  Int_t GetNPtBinsDPlusforDPlusptbin() const {return fnPtBinsDPlusforDPlusptbin;}
  Int_t GetNVarsDPlusforDPlusptbin() const {return fnVarsDPlusforDPlusptbin;}
  Int_t GetGlobalIndexDPlusforDPlusptbin(Int_t iVar,Int_t iPtBin) const{return iPtBin*fnVarsDPlusforDPlusptbin+iVar;}
  void  SetGlobalIndexDPlusforDPlusptbin(){fGlobalIndexDPlusforDPlusptbin=fnVarsDPlusforDPlusptbin*fnPtBinsDPlusforDPlusptbin; return;}
  void  SetNPtBinsDPlusforDPlusptbin(Int_t nptBins){fnPtBinsDPlusforDPlusptbin=nptBins; return;}
  void  SetNVarsDPlusforDPlusptbin(Int_t nVars){fnVarsDPlusforDPlusptbin=nVars; return;} 

  void SetIsUpperCut(Int_t nCutIndex, Bool_t isUpperCut){fIsUpperCut[nCutIndex] = isUpperCut; return;}
  void SetIsUpperCutDPlusforDPlusptbin(Int_t nCutIndex, Bool_t isUpperCut){fIsUpperCutDPlusforDPlusptbin[nCutIndex] = isUpperCut; return;}

  Bool_t GetIsUpperCut(Int_t nCutIndex){return fIsUpperCut[nCutIndex];}
  Bool_t GetIsUpperCutDPlusforDPlusptbin(Int_t nCutIndex){return fIsUpperCutDPlusforDPlusptbin[nCutIndex];}

  void SetIsCutUsed(Int_t nCutIndex, Int_t ptbin, Bool_t isCutUsed){fIsCutUsed[GetGlobalIndex(nCutIndex,ptbin)] = isCutUsed; return;}
  void SetIsCutUsedDPlusforDPlusptbin(Int_t nCutIndex, Int_t ptbin, Bool_t isCutUsed){fIsCutUsedDPlusforDPlusptbin[GetGlobalIndexDPlusforDPlusptbin(nCutIndex,ptbin)] = isCutUsed; return;}

  Bool_t GetIsCutUsed(Int_t nCutIndex, Int_t ptbin) const {return fIsCutUsed[GetGlobalIndex(nCutIndex,ptbin)];}
  Bool_t GetIsCutUsedDPlusforDPlusptbin(Int_t nCutIndex, Int_t ptbin) const {return fIsCutUsedDPlusforDPlusptbin[GetGlobalIndexDPlusforDPlusptbin(nCutIndex,ptbin)];}

  Int_t ApplyCutOnVariable(Int_t nCutIndex, Int_t ptbin, Float_t cutVariableValue, Bool_t bCutArray[78]);
  Int_t ApplyCutOnVariableDPlusforDPlusptbin(Int_t nCutIndex, Int_t ptbin, Float_t cutVariableValue, Bool_t bCutArray[39]);

  void SetVarNamesDPlusforDPlusptbin(Int_t nVars,TString *varNames,Bool_t *isUpperCut);

  Int_t GetMinITSNclsDaughterType(Int_t nDaughterType) const {return fMinITSNclsDaughterType[nDaughterType];}
  Int_t GetMinTPCNclsDaughterType(Int_t nDaughterType) const {return fMinTPCNclsDaughterType[nDaughterType];}
  Bool_t UseITSRefitDaughterType(Int_t nDaughterType) const {return fUseITSRefitDaughterType[nDaughterType];}
  Bool_t UseTPCRefitDaughterType(Int_t nDaughterType) const {return fUseTPCRefitDaughterType[nDaughterType];}
  Bool_t UseFilterBitDaughterType(Int_t nDaughterType) const {return fUseFilterBitDaughterType[nDaughterType];}
  Int_t GetFilterBitDaughterType(Int_t nDaughterType) const {return fFilterBitDaughterType[nDaughterType];}
  Double_t GetMinPtDaughterType(Int_t nDaughterType) const {return fMinPtDaughterType[nDaughterType];}
  Double_t GetMaxAbsEtaDaughterType(Int_t nDaughterType) const {return fMaxAbsEtaDaughterType[nDaughterType];}
  void GetHardSelectionArrayITSDaughterType(Int_t nDaughterType, Bool_t array[7] = 0) const {for(Int_t i=0;i<7;i++){array[i] = fHardSelectionArrayITSDaughterType[nDaughterType][i];} return;}
  void GetSoftSelectionArrayITSDaughterType(Int_t nDaughterType, Bool_t array[7] = 0) const {for(Int_t i=0;i<7;i++){array[i] = fSoftSelectionArrayITSDaughterType[nDaughterType][i];} return;}
  Int_t GetNSoftITSCutDaughterType(Int_t nDaughterType) const {return fNSoftITSCutDaughterType[nDaughterType];}

  void SetMinITSNclsDaughterType(Int_t nDaughterType, Int_t value){fMinITSNclsDaughterType[nDaughterType] = value; return;}
  void SetMinTPCNclsDaughterType(Int_t nDaughterType, Int_t value){fMinTPCNclsDaughterType[nDaughterType] = value; return;}
  void SetUseITSRefitDaughterType(Int_t nDaughterType, Bool_t option){fUseITSRefitDaughterType[nDaughterType] = option; return;}
  void SetUseTPCRefitDaughterType(Int_t nDaughterType, Bool_t option){fUseTPCRefitDaughterType[nDaughterType] = option; return;}
  void SetUseFilterBitDaughterType(Int_t nDaughterType, Bool_t option){fUseFilterBitDaughterType[nDaughterType] = option; return;}
  void SetFilterBitDaughterType(Int_t nDaughterType, Int_t value){fFilterBitDaughterType[nDaughterType] = value; return;}
  void SetMinPtDaughterType(Int_t nDaughterType, Double_t value){fMinPtDaughterType[nDaughterType] = value; return;}
  void SetMaxAbsEtaDaughterType(Int_t nDaughterType, Double_t value){fMaxAbsEtaDaughterType[nDaughterType] = value;}
  void SetHardSelectionArrayITSDaughterType(Int_t nDaughterType, const Bool_t array[7] = 0){for(Int_t i=0;i<7;i++){fHardSelectionArrayITSDaughterType[nDaughterType][i] = array[i];} return;}
  void SetSoftSelectionArrayITSDaughterType(Int_t nDaughterType, const Bool_t array[7] = 0){for(Int_t i=0;i<7;i++){fSoftSelectionArrayITSDaughterType[nDaughterType][i] = array[i];} return;}
  void SetNSoftITSCutDaughterType(Int_t nDaughterType, Int_t value){fNSoftITSCutDaughterType[nDaughterType] = value;}

  void SetCut(Int_t nCutIndex, Int_t ptBin, AliRDHFCutsB0toDPi::EUpperCut cutDirection, Float_t cutValue);
  void SetCutDPlusforDPlusptbin(Int_t nCutIndex, Int_t ptBin, AliRDHFCutsB0toDPi::EUpperCut cutDirection, Float_t cutValue);
  void SetCutForCutOptimization(Int_t nCutIndex, Int_t nVariable, Int_t ptBin, AliRDHFCutsB0toDPi::EUpperCut cutDirection, Float_t * cutValues);
  Float_t GetCutForCutOptimization(Int_t nCutIndex, Int_t nVariable, Int_t ptBin){return fCutsRDForCutOptimization[GetGlobalIndexForCutOptimization(nCutIndex,nVariable,ptBin)];}

  Double_t GetMind0DaughterType(Int_t nDaughterType) const {return fMind0DaughterType[nDaughterType];}
  Double_t GetMinNormd0DaughterType(Int_t nDaughterType) const {return fMinNormd0DaughterType[nDaughterType];}
  Double_t GetFiducialYCut() const {return fFiducialYCut;}

  void SetMind0DaughterType(Int_t nDaughterType, Double_t value){fMind0DaughterType[nDaughterType] = value; return;}
  void SetMinNormd0DaughterType(Int_t nDaughterType, Double_t value){fMinNormd0DaughterType[nDaughterType] = value; return;}  
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

  Float_t fMaxPtPid;                                  ///
  Float_t fTPCflag;                                   ///
  Double_t fCircRadius;                               /// Radius for circular PID nsigma cut

  Bool_t * fIsCutUsed;                                //[fGlobalIndex]

  Int_t fnVarsDPlusforDPlusptbin;                           ///
  Int_t fnPtBinsDPlusforDPlusptbin;                         ///
  Int_t fGlobalIndexDPlusforDPlusptbin;                     ///
  Float_t * fCutsRDDPlusforDPlusptbin;                      //[fGlobalIndexDPlusforDPlusptbin]
  Int_t fnPtBinLimitsDPlusforDPlusptbin;                    ///
  Float_t * fPtBinLimitsDPlusforDPlusptbin;                 //[fnPtBinLimitsDPlusforDPlusptbin]
  Bool_t * fIsUpperCutDPlusforDPlusptbin;                   //[fnVarsDPlusforDPlusptbin]
  Bool_t * fIsCutUsedDPlusforDPlusptbin;                    //[fGlobalIndexDPlusforDPlusptbin]
  TString * fVarNamesDPlusforDPlusptbin;                    //[fnVarsDPlusforDPlusptbin]

  Int_t fMinITSNclsDaughterType[3];                   ///
  Int_t fMinTPCNclsDaughterType[3];                   ///
  Bool_t fUseITSRefitDaughterType[3];                 ///
  Bool_t fUseTPCRefitDaughterType[3];                 ///
  Bool_t fUseFilterBitDaughterType[3];                ///
  Int_t fFilterBitDaughterType[3];                    ///
  Double_t fMinPtDaughterType[3];                     ///
  Double_t fMaxAbsEtaDaughterType[3];                 ///
  Bool_t fHardSelectionArrayITSDaughterType[3][7];    ///
  Bool_t fSoftSelectionArrayITSDaughterType[3][7];    ///
  Int_t fNSoftITSCutDaughterType[3];                  ///

  Double_t fMind0DaughterType[3];                     ///
  Double_t fMinNormd0DaughterType[3];                 ///
  Double_t fFiducialYCut;                             ///

  Int_t fnVariablesForCutOptimization;                ///
  Int_t fnCutsForOptimization;                        ///
  Int_t fGlobalIndexCutOptimization;                  ///
  Float_t * fCutsRDForCutOptimization;                //[fGlobalIndexCutOptimization]
  Bool_t * fIsUpperCutForCutOptimization;             //[fnVariablesForCutOptimization]
  Int_t * fCutIndexForCutOptimization;                //[fnVariablesForCutOptimization]
  Float_t * fSigmaForCutOptimization;                 //[fnPtBins]
  Int_t fNumberOfSigmaBinsForCutOptimization;             ///

  /// \cond CLASSIMP    
  ClassDef(AliRDHFCutsB0toDPi,1) ///
  /// \endcond
};

#endif


