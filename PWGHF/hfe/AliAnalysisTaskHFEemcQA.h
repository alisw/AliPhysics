#ifndef AliAnalysisTaskHFEemcQA_cxx
#define AliAnalysisTaskHFEemcQA_cxx

//QA task for EMCAL electron analysis 

class TH1F;
class AliESDEvent;
class AliAODEvent;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskHFEemcQA : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskHFEemcQA() : AliAnalysisTaskSE(), fVevent(0),fESD(0),fAOD(0),fOutputList(0),fVtxZ(0),fVtxX(0),fVtxY(0),fTrigMulti(0),fHistClustE(0),fEMCClsEtaPhi(0),fNegTrkIDPt(0),fTrkPt(0),fTrketa(0),fTrkphi(0),fdEdx(0),fTPCNpts(0),fHistPtMatch(0),fEMCTrkMatch(0),fEMCTrkPt(0),fEMCTrketa(0),fEMCTrkphi(0),fEMCdEdx(0),fEMCTPCNpts(0),fHistdEdxEop(0),fHistEop(0),fEleCanTPCNpts(0),fEleCanTPCNCls(0),fEleCanITSNCls(0),fEleCanITShit(0),fEleCanSPD1(0),fEleCanSPD2(0),fEleCanSPDBoth(0),fEleCanSPDOr(0) {}
  AliAnalysisTaskHFEemcQA(const char *name);
  virtual ~AliAnalysisTaskHFEemcQA() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:
  AliVEvent   *fVevent;  //!V event object
  AliESDEvent *fESD;    //!ESD object
  AliAODEvent *fAOD;    //!AOD object
  TList       *fOutputList; //!Output list
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
  TH1F        *fHistPtMatch;//!tracks matched to EMCAL 
  TH2F        *fEMCTrkMatch;//!Distance of EMC cluster to closest track in phi and z
  TH1F        *fEMCTrkPt;//!tracks with EMCAL cluster
  TH1F        *fEMCTrketa;//!EMC trk eta
  TH1F        *fEMCTrkphi;//!EMC trk phi
  TH2F        *fEMCdEdx;//!EMC trk dedx
  TH2F        *fEMCTPCNpts;//!EMC Npoints used for dedx
  TH2F        *fHistdEdxEop;//!E/p vs dedx
  TH2F        *fHistEop;//!pt vs E/p
  TH2F        *fEleCanTPCNpts;//!ele cand TPC Npoints used for dedx
  TH2F        *fEleCanTPCNCls;//!ele cand TPC N clusters
  TH2F        *fEleCanITSNCls;//!ele cand ITS N clusters
  TH1F        *fEleCanITShit;//!ele cand ITS hit map
  TH2F        *fEleCanSPD1;//!ele cand hit SPD layer 1
  TH2F        *fEleCanSPD2;//!ele cand hit SPD layer 2
  TH2F        *fEleCanSPDBoth;//!ele cand SPD both layer
  TH2F        *fEleCanSPDOr;//!ele cand SPD or

  AliAnalysisTaskHFEemcQA(const AliAnalysisTaskHFEemcQA&); // not implemented
  AliAnalysisTaskHFEemcQA& operator=(const AliAnalysisTaskHFEemcQA&); // not implemented
  
  ClassDef(AliAnalysisTaskHFEemcQA, 1); // example of analysis
};

#endif


