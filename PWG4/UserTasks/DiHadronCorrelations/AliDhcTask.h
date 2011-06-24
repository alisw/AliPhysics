// Dihadron correlations task - simple task to read ESD or AOD input,
// calculate same- and mixed-event correlations, and fill THnSparse
// output. -A. Adare, Apr 2011


#ifndef AliDhcTask_cxx
#define AliDhcTask_cxx

class TH1;
class TH2;
class THnSparse;
class TObject;
class TObjArray;
class AliESDEvent;
class AliAODEvent;
class AliESDtrackCuts;
class KiddiePoolManager;

#include "AliAnalysisTaskSE.h"
#include "KiddiePoolClasses.h"

class AliDhcTask : public AliAnalysisTaskSE {
 public:
  AliDhcTask() : 
    AliAnalysisTaskSE(), fVerbosity(0), fESD(0), fAOD(0), fOutputList(0), 
    fHistPt(0), fHEvt(0), fHTrk(0), fHS(0), fHM(0), fPoolMgr(0), fCentrality(99), 
    fZVertex(99), fZVtxMax(10), fPtMin(0), fPtMax(100), fInputHandler(0), 
    fEsdTrackCutsTPCOnly(0) {}
  AliDhcTask(const char *name);
  virtual ~AliDhcTask() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  void SetVerbosity(const Int_t v) { fVerbosity = v; }

 protected:
  enum ePairHistAxes  {kDeta, kPtAssc, kPtTrig, kCent, kDphi,
		       kZvtx, kChargeComb};
  enum eEventHistAxes {kZvtxEvt, kCentV0M, kCentCL1};
  enum eTrackHistAxes {kPhiTrk, kEtaTrk};
  enum ePairingScheme {kSameEvt, kDiffEvt};
  enum eDataType {kESD, kAOD};

  void BookHistos();
  void InitEventMixer();
  MiniEvent* GetESDTrax() const;
  MiniEvent* GetAODTrax() const;
  Bool_t VertexOk(TObject* obj) const;
  Double_t DeltaPhi(Double_t phia, Double_t phib, 
		    Double_t rangeMin = -TMath::Pi()/2, 
		    Double_t rangeMax = 3*TMath::Pi()/2) const;
  Int_t Correlate(const MiniEvent &arr1, const MiniEvent &arr2, 
		  Int_t pairing = kSameEvt, Double_t weight = 1.);

 private:
  Int_t        fVerbosity;       // 0 = silence
  AliESDEvent *fESD;             //! ESD object
  AliAODEvent *fAOD;             //! AOD object
  TList       *fOutputList;      //! Output list
  TH1F        *fHistPt;          //! Pt spectrum
  TH2         *fHEvt;            //! Cent, vtx, etc.
  TH2         *fHTrk;            //! Phi, Eta, etc.
  THnSparse   *fHS;              //! Same-evt correlations
  THnSparse   *fHM;              //! Diff-evt correlations
  KiddiePoolManager* fPoolMgr;   //! Event mixer
  Double_t    fCentrality;       //! V0M for now
  Double_t    fZVertex;          //! Of current event
  Double_t    fZVtxMax;          //! Max |z| cut (cm)
  Double_t    fPtMin;            //! Min pt cut
  Double_t    fPtMax;            //! Max pt cut
  AliInputEventHandler*  fInputHandler;
  AliESDtrackCuts* fEsdTrackCutsTPCOnly;

  AliDhcTask(const AliDhcTask&);            // not implemented
  AliDhcTask &operator=(const AliDhcTask&); // not implemented

  ClassDef(AliDhcTask, 1);
};

#endif
