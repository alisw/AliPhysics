#ifndef ALIANALYSISTASKGFWFLOW__H
#define ALIANALYSISTASKGFWFLOW__H
#include "AliAnalysisTaskSE.h"
#include "TComplex.h"
#include "AliEventCuts.h"
#include "AliVParticle.h"
#include "AliGFWCuts.h"
#include "TAxis.h"
#include "TStopwatch.h"
#include "AliGFW.h"
#include "AliVEvent.h"


class TList;
class TH1D;
class TH2D;
class TH3D;
class TProfile;
class TProfile2D;
class TComplex;
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

class AliAnalysisTaskGFWFlow : public AliAnalysisTaskSE {
 public:
  Int_t debugpar;
  AliAnalysisTaskGFWFlow();
  AliAnalysisTaskGFWFlow(const char *name, Bool_t ProduceWeights=kTRUE, Bool_t IsMC=kTRUE, Bool_t AddQA=kFALSE);
  virtual ~AliAnalysisTaskGFWFlow();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  Bool_t AcceptEvent();
  Bool_t AcceptAODVertex(AliAODEvent*);
  void SetPtBins(Int_t nBins, Double_t *bins, Double_t RFpTMin=-1, Double_t RFpTMax=-1); //Also set the RF pT acceptance
  void SetCurrSystFlag(Int_t newval) { fCurrSystFlag = newval; };
  void SetWeightDir(const char *newval) { fWeightDir.Clear(); fWeightDir.Append(newval); };
  Bool_t SetInputWeightList(TList *inList);
  vector<AliGFW::CorrConfig> corrconfigs; //! do not store
  AliGFW::CorrConfig GetConf(TString head, TString desc, Bool_t ptdif) { return fGFW->GetCorrelatorConfig(desc,head,ptdif);};
  void CreateCorrConfigs();
  void SetTriggerType(AliVEvent::EOfflineTriggerTypes newval) { fTriggerType = newval; };
  Bool_t CheckTriggerVsCentrality(Double_t l_cent); //Hard cuts on centrality for special triggers
  void SetBypassCalculations(Bool_t newval) { fBypassCalculations = newval; };
 protected:
  AliEventCuts fEventCuts, fEventCutsForPU;
 private:
  AliAnalysisTaskGFWFlow(const AliAnalysisTaskGFWFlow&);
  AliAnalysisTaskGFWFlow& operator=(const AliAnalysisTaskGFWFlow&);
  AliVEvent::EOfflineTriggerTypes fTriggerType; //! No need to store
  Bool_t fProduceWeights;
  AliGFWCuts **fSelections; //! Selection array; not store
  TList *fWeightList; //! Stored via PostData
  AliGFWWeights *fWeights; //! these are stored in a list now
  AliGFWWeights *fExtraWeights; //! to fetch ITS weights, if required
  AliGFWFlowContainer *fFC; // Flow container
  AliGFW *fGFW; //! no need to store this
  TTree *fOutputTree; //! Not stored and not needed
  AliMCEvent *fMCEvent; //! Not stored
  Bool_t fIsMC;
  TAxis *fPtAxis; // No need to store this
  Double_t fPOIpTMin; //pT min for POI
  Double_t fPOIpTMax; //pT max for POI
  Double_t fRFpTMin; //pT min for RF
  Double_t fRFpTMax; //pT max for RF
  TString fWeightPath; //! No need to store this
  TString fWeightDir; //Directory where to find weights
  //Double_t fPtBins; //! Not stored
  Int_t fTotFlags; //1 for normal, plus 1 per each flag
  Int_t fTotTrackFlags; //Total number of track flags
  Int_t fRunNo;
  Int_t fCurrSystFlag;
  Bool_t fAddQA; // Add AliEventSelection QA plots
  TList *fQAList;
  Bool_t fBypassCalculations; //Flag to bypass all the calculations, so only event selection is performed (for QA)
  Int_t AcceptedEventCount;
  Int_t GetVtxBit(AliAODEvent *mev);
  Int_t GetParticleBit(AliVParticle *mpa);
  Int_t GetTrackBit(AliAODTrack *mtr, Double_t *lDCA);
  Int_t CombineBits(Int_t VtxBit, Int_t TrkBit);
  Bool_t AcceptParticle(AliVParticle *mPa);
  Bool_t InitRun();
  Bool_t LoadWeights(Int_t runno);
  Bool_t FillFCs(AliGFW::CorrConfig corconf, Double_t cent, Double_t rndm, Bool_t DisableOverlap=kFALSE);
  Bool_t FillFCs(TString head, TString hn, Double_t cent, Bool_t diff, Double_t rndmn);
 // TStopwatch mywatch;
 // TStopwatch mywatchFill;
 // TStopwatch mywatchStore;
  ClassDef(AliAnalysisTaskGFWFlow,1);
};

#endif
