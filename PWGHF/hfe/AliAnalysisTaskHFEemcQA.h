#ifndef AliAnalysisTaskHFEemcQA_cxx
#define AliAnalysisTaskHFEemcQA_cxx

//QA task for EMCAL electron analysis 

class TH1F;
class THnSparse;
class AliESDEvent;
class AliAODEvent;
class AliHFEcontainer;
class AliHFEcuts;
class AliHFEpid;
class AliHFEpidQAmanager;
class AliCFManager;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskHFEemcQA : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskHFEemcQA();
  AliAnalysisTaskHFEemcQA(const char *name);
  virtual ~AliAnalysisTaskHFEemcQA();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetAODAnalysis() { SetBit(kAODanalysis, kTRUE); };
  void SetESDAnalysis() { SetBit(kAODanalysis, kFALSE); };
  Bool_t IsAODanalysis() const { return TestBit(kAODanalysis); };

  Bool_t GetElecIDsparse() {return fFlagSparse;};
  void SetElecIDsparse(Bool_t flagelecIDsparse){fFlagSparse = flagelecIDsparse;};

 private:
  enum{
    kAODanalysis = BIT(20),
  };

  AliVEvent   *fVevent;  //!event object
  AliESDEvent *fESD;    //!ESD object
  AliAODEvent *fAOD;    //!AOD object
  AliPIDResponse *fpidResponse; //!pid response

  Bool_t      fFlagSparse;// switch to THnspare

  TList       *fOutputList; //!Output list
  TH1F        *fNevents;//! no of events
  TH1F        *fVtxZ;//!Vertex z 
  TH1F        *fVtxX;//!Vertex x 
  TH1F        *fVtxY;//!Vertex y 
  TH2F        *fTrigMulti;//!trigger multiplicity 
  TH1F        *fHistClustE;//!cluster energy
  TH2F        *fEMCClsEtaPhi;//! EMC cluster eta and phi 
  TH1F        *fNegTrkIDPt;//!neg track ID
  TH1F        *fTrkPt;//!track pt
  TH1F        *fTrketa;//!track eta
  TH1F        *fTrkphi;//!track phi 
  TH2F        *fdEdx;//!dedx vs pt
  TH2F        *fTPCNpts;//!TPC Npoints used for dedx
  TH2F        *fTPCnsig;//!TPC Nsigma
  TH1F        *fHistPtMatch;//!tracks matched to EMCAL 
  TH2F        *fEMCTrkMatch;//!Distance of EMC cluster to closest track in phi and z
  TH1F        *fEMCTrkPt;//!tracks with EMCAL cluster
  TH1F        *fEMCTrketa;//!EMC trk eta
  TH1F        *fEMCTrkphi;//!EMC trk phi
  TH2F        *fEMCdEdx;//!EMC trk dedx
  TH2F        *fEMCTPCnsig;//! EMC trk nsig
  TH2F        *fEMCTPCNpts;//!EMC Npoints used for dedx
  TH2F        *fHistdEdxEop;//!E/p vs dedx
  TH2F        *fHistNsigEop;//!E/p vs dedx
  TH2F        *fHistEop;//!pt vs E/p
  TH2F        *fEleCanTPCNpts;//!ele cand TPC Npoints used for dedx
  TH2F        *fEleCanTPCNCls;//!ele cand TPC N clusters
  TH2F        *fEleCanITSNCls;//!ele cand ITS N clusters
  TH1F        *fEleCanITShit;//!ele cand ITS hit map
  TH2F        *fEleCanSPD1;//!ele cand hit SPD layer 1
  TH2F        *fEleCanSPD2;//!ele cand hit SPD layer 2
  TH2F        *fEleCanSPDBoth;//!ele cand SPD both layer
  TH2F        *fEleCanSPDOr;//!ele cand SPD or

  THnSparse  *fSparseElectron;//!Electron info 
  Double_t *fvalueElectron;//!Electron info

  AliAnalysisTaskHFEemcQA(const AliAnalysisTaskHFEemcQA&); // not implemented
  AliAnalysisTaskHFEemcQA& operator=(const AliAnalysisTaskHFEemcQA&); // not implemented

  ClassDef(AliAnalysisTaskHFEemcQA, 1); // example of analysis
};

#endif


