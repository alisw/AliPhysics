/*
Author: Vytautas Vislavicius
Extention of Generic Flow (https://arxiv.org/abs/1312.3572)
*/
#ifndef ALIANALYSISTASKGFWFLOW__H
#define ALIANALYSISTASKGFWFLOW__H
#include "AliAnalysisTaskSE.h"
// #include "TComplex.h"
#include "AliVParticle.h"
#include "TAxis.h"
#include "TStopwatch.h"
#include "AliGFW.h"
#include "AliVEvent.h"
#include "GFWFlags.h"
#include "AliGFWFilter.h"
// #include "TGrid.h"

class TList;
class TH1D;
class TH2D;
class TH3D;
class TProfile;
class TProfile2D;
// class TComplex;
class AliVEvent;
class AliAODEvent;
class AliVTrack;
class AliVVertex;
class AliInputEventHandler;
class AliAODTrack;
class TTree;
class TClonesArray;
class AliMCEvent;
class AliGFWWeights;
class AliGFWFlowContainer;
class TObjArray;
class TNamed;
class AliAODVertex;
class AliAnalysisUtils;
using namespace GFWFlags;
class AliAnalysisTaskGFWFlow : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskGFWFlow();
  AliAnalysisTaskGFWFlow(const char *name, Bool_t ProduceWeights=kTRUE, Bool_t IsMC=kTRUE, Bool_t IsTrain=kFALSE);
  virtual ~AliAnalysisTaskGFWFlow();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void NotifyRun();
  virtual void Terminate(Option_t *);
  void SetPtBins(Int_t nBins, Double_t *bins, Double_t RFpTMin=-1, Double_t RFpTMax=-1); //Also set the RF pT acceptance
  //In case we want custom nominal flags (defaults are the first ones)
  void SetNominalFlags(Int_t lEvFlagIndex, Int_t lTrFlagIndex) { fEvNomFlag=(1<<lEvFlagIndex); fTrNomFlag=(1<<lTrFlagIndex); };
  void SetupFlagsByIndex(Int_t ind); //Local envelope for the function below
  static void SetupFlagsByIndex(const Int_t &ind, UInt_t &l_EvFlag, UInt_t &l_TrFlag); //Function to setup flags. Static, so one is able to call from the outside
  void SetCustomNoFlags(Int_t nEvFlags, Int_t nTrFlags) {fTotTrackFlags=nTrFlags; fTotEvFlags=nEvFlags; };
  vector<AliGFW::CorrConfig> corrconfigs; //! do not store
  AliGFW::CorrConfig GetConf(TString head, TString desc, Bool_t ptdif) { return fGFW->GetCorrelatorConfig(desc.Data(),head.Data(),ptdif);};
  void CreateCorrConfigs();
  void SetTriggerType(UInt_t newval) { fTriggerType = newval; };
  Bool_t CheckTriggerVsCentrality(Double_t l_cent); //Hard cuts on centrality for special triggers
  void SetBypassCalculations(Bool_t newval) { fBypassCalculations = newval; };
  void SetCollisionSystem(Int_t newval) { fCollisionsSystem = newval; };
 private:
  AliAnalysisTaskGFWFlow(const AliAnalysisTaskGFWFlow&);
  AliAnalysisTaskGFWFlow& operator=(const AliAnalysisTaskGFWFlow&);
  UInt_t fTriggerType; //Need to store this for it to be able to work on trains
  Bool_t fProduceWeights;
  TList *fWeightList; //! Stored via PostData
  TH1D *fCentMap; //! centrality map for on-fly trains
  AliGFWWeights *fWeights; //! these are stored in a list now
  AliGFWFlowContainer *fFC; // Flow container
  AliGFW *fGFW; //! no need to store this
  TTree *fOutputTree; //! Not stored and not needed
  AliMCEvent *fMCEvent; //! Not stored
  Bool_t fIsMC;
  Bool_t fIsTrain;
  UInt_t fEvNomFlag; //Nominal event selection flag
  UInt_t fTrNomFlag; //Nominal track selection flag
  TAxis *fPtAxis; // No need to store this
  Double_t fPOIpTMin; //pT min for POI
  Double_t fPOIpTMax; //pT max for POI
  Double_t fRFpTMin; //pT min for RF
  Double_t fRFpTMax; //pT max for RF
  //Double_t fPtBins; //! Not stored
  Int_t fTotFlags; //1 for normal, plus 1 per each flag
  Int_t fTotTrackFlags; //Total number of track flags
  Int_t fTotEvFlags; //Total number of event flags
  Int_t fRunNo;
  Bool_t fBypassCalculations; //Flag to bypass all the calculations, so only event selection is performed (for QA)
  TH1D *fMultiDist;
  Int_t fCollisionsSystem; //0 for pp, 1 for pPb, 2 for PbPb
  Double_t fOverrideCentrality; //Relevant for when running on trains with fixed centrality classes.
  Bool_t LoadWeights(Int_t runno);
  Bool_t FillFCs(AliGFW::CorrConfig corconf, Double_t cent, Double_t rndm, Bool_t DisableOverlap=kFALSE);
  AliMCEvent *FetchMCEvent(Double_t &impactParameter);
  Double_t GetCentFromIP(Double_t impactParameter) { return fCentMap->GetBinContent(fCentMap->FindBin(impactParameter)); };
  void OverrideCentralityValue(Double_t newval=-1) { fOverrideCentrality = newval; };
 // TStopwatch mywatch;
 // TStopwatch mywatchFill;
 // TStopwatch mywatchStore;
  ClassDef(AliAnalysisTaskGFWFlow,1);
};

#endif
