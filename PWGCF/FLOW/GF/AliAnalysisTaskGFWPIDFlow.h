/*
Author: Vytautas Vislavicius
Extention of Generic Flow (https://arxiv.org/abs/1312.3572)
*/
#ifndef AliAnalysisTaskGFWPIDFlow__H
#define AliAnalysisTaskGFWPIDFlow__H
#include "AliAnalysisTaskSE.h"
#include "TComplex.h"
#include "AliEventCuts.h"
#include "AliVEvent.h"
#include "AliGFW.h"
#include "AliPID.h"

class TList;
class TH1D;
class TH2D;
class TH3D;
class TProfile;
class TProfile2D;
class TComplex;
class AliAODEvent;
class AliVTrack;
class AliVVertex;
class AliInputEventHandler;
class AliAODTrack;
class TClonesArray;
class AliAODVertex;
class AliAnalysisUtils;
class TProfile;
class AliGFWWeights;
class AliVParticle;
class AliGFWCuts;
class AliGFWFlowContainer;
class AliPIDResponse;
class AliPIDCombined;
class TRandom;

class AliAnalysisTaskGFWPIDFlow : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskGFWPIDFlow();
  AliAnalysisTaskGFWPIDFlow(const char *name, Bool_t IsMC=kTRUE, TString StageSwitch="");
  virtual ~AliAnalysisTaskGFWPIDFlow();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  Bool_t CheckTrigger(Double_t);
  Bool_t AcceptAOD(AliAODEvent*, Double_t lvtxXYZ[3]);
  Bool_t AcceptParticle(AliVParticle*);
  void SetTriggerType(UInt_t newval) {fTriggerType = newval; };
  void FillWeights(AliAODEvent*, Double_t vz, Double_t l_Cent);
  void FillMeanPtCounter(Double_t l_pt, Double_t &l_sum, Double_t &l_count, AliGFWWeights *inWeight); //passing by ref., considering how ofter this is called
  void FillMeanPt(AliAODEvent*, Double_t vz, Double_t l_Cent);
  void FillCK(AliAODEvent *fAOD, Double_t vz, Double_t l_Cent);
  Int_t GetStageSwitch(TString instr);
  AliGFW::CorrConfig GetConf(TString head, TString desc, Bool_t ptdif) { return fGFW->GetCorrelatorConfig(desc,head,ptdif);};
  void CreateCorrConfigs();
  void LoadWeightAndMPT(AliAODEvent*);
  void GetSingleWeightFromList(AliGFWWeights **inWeights, Int_t runno, TString pf="");
  Bool_t WithinSigma(Double_t SigmaCut, AliAODTrack *inTrack, AliPID::EParticleType partType);
  //In development
  Double_t GetZMWeight(Double_t eta, Double_t phi, Int_t PIDIndex);
  void DevFunction(AliAODEvent *fAOD, Double_t vz, Double_t l_Cent);
  void FillCustomWeights(AliAODEvent *fAOD, Double_t vz, Double_t l_Cent);
  void LoadMyWeights(AliAODEvent*);
  void SetUseRunAvgWeights(Bool_t newval) { fUseRunAveragedWeights = newval; };
 protected:
  AliEventCuts fEventCuts;
 private:
  AliAnalysisTaskGFWPIDFlow(const AliAnalysisTaskGFWPIDFlow&);
  AliAnalysisTaskGFWPIDFlow& operator=(const AliAnalysisTaskGFWPIDFlow&);
  Int_t fStageSwitch;
  Bool_t fIsMC;
  AliPIDResponse *fPIDResponse; //!
  TList *fMPTList; //!
  TProfile *fmPT; //!
  TProfile *fmPT_pi; //!
  TProfile *fmPT_ka; //!
  TProfile *fmPT_pr; //!
  TH1D *fMultiDist;
  TProfile *fptvar;
  TProfile *fCovariance;
  Bool_t fmptSet;
  UInt_t fTriggerType; //! No need to store
  Bool_t fUseRunAveragedWeights; //!
  TList *fWeightList; //!
  AliGFWWeights *fWeights;//! This should be stored in TList
  AliGFWWeights *fWeights_pi;//! This should be stored in TList
  AliGFWWeights *fWeights_ka;//! This should be stored in TList
  AliGFWWeights *fWeights_pr;//! This should be stored in TList
  Int_t fRunNo; //!
  AliGFWCuts *fMidSelection; //!
  AliGFWCuts *fFWSelection; //!
  AliGFWFlowContainer *fFC;
  AliGFW *fGFW; //! not stored
  vector<AliGFW::CorrConfig> corrconfigs; //! do not store
  Bool_t FillFCs(AliGFW::CorrConfig corconf, Double_t cent, Double_t rndmn, Bool_t EnableDebug=kFALSE); //Pending implementation: possibility to pass pre-calculated values (e.g. for ref flow)
  Bool_t FillCovariance(AliGFW::CorrConfig corconf, Double_t cent, Double_t d_mpt, Double_t dw_mpt);
  Bool_t AcceptAODTrack(AliAODTrack *lTr, Double_t*);
  //In development
  TH2D **fZMWeights; //
  AliPIDCombined *fBayesPID;
  Bool_t HasTPCPID(AliAODTrack* l_track);
  Bool_t HasTOFPID(AliAODTrack* l_track);
  Int_t GetBayesPIDIndex(AliAODTrack *l_track);
  Bool_t devAcceptAODTrack(AliAODTrack *lTr, Double_t*);
  void AddToOBA(TObjArray *oba, TString l_name, Int_t nPT=0);
  TAxis *fPtAxis; //!
  Bool_t GetIntValAndDNX(AliGFW::CorrConfig corconf, Double_t &l_val, Double_t &l_dnx);
  TRandom *fRndm; //!
  ClassDef(AliAnalysisTaskGFWPIDFlow,1);
};

#endif
