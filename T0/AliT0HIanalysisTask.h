#ifndef AliT0HIanalysisTask_cxx
#define AliT0HIanalysisTask_cxx

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"
#include "TTree.h"
#include "TString.h"

#define NPMT0 24  //number T0 of photomultipliers

class AliESDEvent;
#include "AliAnalysisTaskSE.h"

class AliT0HIanalysisTask : public AliAnalysisTaskSE {
 public:
  AliT0HIanalysisTask() : AliAnalysisTaskSE(), 
    fESD(0), fOutputList(0), fT0OutTree(0) {};
  AliT0HIanalysisTask(const char *name);
  virtual ~AliT0HIanalysisTask();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  Bool_t UserNotify();  
 private:
  AliESDEvent *fESD;    //! ESD object
  TList     *fOutputList; //! Output list
  TTree     *fT0OutTree;  //output tree
  Int_t fEvent;
  Int_t fOrbit;
  Int_t fBC;
  Int_t fTrackletSPD;
  Int_t fClustersSPD;
  Int_t fNcont;
  Int_t fNcontTPC;
  Float_t fVertex;
  Float_t fVertexPrim;
  Float_t fVertexSPD;
  Float_t fVertexTPC;
  Float_t fMeanAC;
  Float_t fMeanA;
  Float_t fMeanC;
  Float_t fMeanACcalc;
  Float_t fMeanAcalc;
  Float_t fMeanCcalc;
  Float_t fMultV0A;
  Float_t fMultV0C;
  Float_t fTimeV0A;
  Float_t fTimeV0C;
  Float_t fSumampA;
  Float_t fSumampC;
  UInt_t ftimestamp;
  Float_t fSep2;
  Bool_t fZDCcut;
  Int_t fT0Trigger;
  Bool_t fpileup;
  TObjString fTrigger;
  TH1F    **fT0_amplitude; //! Amplitudes
  TH1F    **fT0_time;      //! Time
  Float_t fcentralityV0M;
  Float_t fcentralityZDC;
  Float_t fcentralityTRK;
  Float_t fcentralityCLA;
  Int_t fESDtracks;
  Float_t fOrA[5];
  Float_t fOrC[5];
  Float_t fTVDC[5];

  Float_t famp[24];
  Float_t ftime[24];
  Float_t fRawTime[24][5];
  Bool_t fT0pileup[3];
  Int_t fMultiplicity;
  TObjString fTriggerinput ;  
  Int_t fTOFtracks;
  Float_t  t0tofTrack;
  AliESDpid* fESDpid;  //! esd pid 

  AliT0HIanalysisTask(const AliT0HIanalysisTask&); // not implemented
  AliT0HIanalysisTask& operator=(const AliT0HIanalysisTask&); // not implemented
  
  ClassDef(AliT0HIanalysisTask, 3); // example of analysis
};

#endif
