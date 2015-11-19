#ifndef AliT0CalibAnalysisTask_cxx
#define AliT0CalibAnalysisTask_cxx

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"
#include "TTree.h"
#include "TString.h"

#define NPMT0 24  //number T0 of photomultipliers

class AliESDEvent;
#include  "AliESDpid.h"
#include "AliAnalysisTaskSE.h"

class AliT0CalibAnalysisTask : public AliAnalysisTaskSE {
 public:
  AliT0CalibAnalysisTask() : AliAnalysisTaskSE(), 
    fESD(0), fOutputList(0), fT0OutTree(0), fEvent(-99999),   fOrbit(-99999), fBC(-99999),
  fTrackletSPD(-99999),  fNcont(-99999),
  fVertex(-99999),  fVertexSPD(-99999), 
  fMeanAC(-99999), fMeanA(-99999),fMeanC(-99999), 
  fMultV0A(-99999), fMultV0C(-99999),fTimeV0A(-99999),fTimeV0C(-99999), fSumampA(-99999), fSumampC(-99999),
  ftimestamp(0), fSep2(0),
  fZDCcut(kFALSE), fT0Trigger(-99999), fpileup(kFALSE), fTrigger(0x0),
  fT0_amplitude(0x0), fT0_time(0x0),
  fcentralityV0M(0), fcentralityZDC(0), fcentralityTRK(0), 
  fESDtracks(-99999),
    fMultiplicity(-99999),
  fTriggerinput(0x0), fZDCbg(kFALSE),
  fTOFtracks(0), fT0tofTrack(0),
  fESDpid(new AliESDpid())
    {};

  AliT0CalibAnalysisTask(const char *name);
  virtual ~AliT0CalibAnalysisTask();
  
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
  Int_t fNcont;
  Float_t fVertex;
  Float_t fVertexSPD;
  Float_t fMeanAC;
  Float_t fMeanA;
  Float_t fMeanC;
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
  Int_t fESDtracks;
  Float_t fOrA[5];
  Float_t fOrC[5];
  Float_t fTVDC[5];

  Float_t famp[24];
 Float_t famp_new[24];	
  Float_t ftime[24];
  Float_t fRawTime[24];
  Bool_t fT0pileup[3];
  Int_t fMultiplicity;
  TObjString fTriggerinput ;
  Bool_t fZDCbg;   //ZDC BG flag
  Int_t fTOFtracks;
  Float_t fT0tofTrack;
  AliESDpid* fESDpid;  //! esd pid
  TBits fPFPbit; //PFP bits
  
  AliT0CalibAnalysisTask(const AliT0CalibAnalysisTask&); // not implemented
  AliT0CalibAnalysisTask& operator=(const AliT0CalibAnalysisTask&); // not implemented
  
  ClassDef(AliT0CalibAnalysisTask, 1); // example of analysis
};

#endif
