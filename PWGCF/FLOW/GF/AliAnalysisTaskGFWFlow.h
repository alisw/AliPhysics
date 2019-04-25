#ifndef ALIANALYSISTASKGFWFLOW__H
#define ALIANALYSISTASKGFWFLOW__H
#include "AliAnalysisTaskSE.h"
#include "TComplex.h"
#include "AliEventCuts.h"
#include "AliVParticle.h"
#include "AliGFWCuts.h"
#include "TAxis.h"

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
class AliGFW;
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
  void SetPtBins(Int_t nBins, Double_t *bins) { fPtAxis->Set(nBins,bins); };
  void SetCurrSystFlag(Int_t newval) { fCurrSystFlag = newval; };
  void SetWeightDir(const char *newval) { fWeightDir.Clear(); fWeightDir.Append(newval); };
  Bool_t SetInputWeightList(TList *inList);
 protected:
  AliEventCuts fEventCuts, fEventCutsForPU;
 private:
  AliAnalysisTaskGFWFlow(const AliAnalysisTaskGFWFlow&);
  AliAnalysisTaskGFWFlow& operator=(const AliAnalysisTaskGFWFlow&);
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
  TString fWeightPath; //! No need to store this
  TString fWeightDir; //Directory where to find weights
  //Double_t fPtBins; //! Not stored
  Int_t fTotFlags; //1 for normal, plus 1 per each flag
  Int_t fTotTrackFlags; //Total number of track flags
  Int_t fRunNo;
  Int_t fCurrSystFlag;
  Bool_t fAddQA; // Add AliEventSelection QA plots
  TList *fQAList;
  Int_t AcceptedEventCount;
  Int_t GetVtxBit(AliAODEvent *mev);
  Int_t GetParticleBit(AliVParticle *mpa);
  Int_t GetTrackBit(AliAODTrack *mtr, Double_t *lDCA);
  Int_t CombineBits(Int_t VtxBit, Int_t TrkBit);
  Bool_t AcceptParticle(AliVParticle *mPa);
  Bool_t InitRun();
  Bool_t LoadWeights(Int_t runno);
  Bool_t FillFCs(TString head, TString hn, Double_t cent, Bool_t diff, Double_t rndmn);
  ClassDef(AliAnalysisTaskGFWFlow,1);
};

#endif
