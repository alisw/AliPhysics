#ifndef ALIANALYSISTASKTOFQA_h
#define ALIANALYSISTASKTOFQA_h

class TString;
class TList;
class AliESDEvent;
class AliAnalysisFilter;
class TDatabasePDG;


#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskTOFqa : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskTOFqa();
  AliAnalysisTaskTOFqa(const char *name);
  AliAnalysisTaskTOFqa(const AliAnalysisTaskTOFqa& copy);
  AliAnalysisTaskTOFqa& operator= (const AliAnalysisTaskTOFqa& copy);
  virtual ~AliAnalysisTaskTOFqa();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

    
  Int_t GetStripIndex(const Int_t * const in);
  void SetTrackFilter(AliAnalysisFilter *filter) {fTrackFilter = filter;};
  void EnableAdvancedCheck(Bool_t enable){fEnableAdvancedCheck=enable;};
  void SetExpTimeHistoRange(Float_t min, Float_t max){fExpTimeRangeMin=min; fExpTimeRangeMax=max;return;};
  void SetExpTimeHistoSmallRange(Float_t min, Float_t max){fExpTimeSmallRangeMin=min; fExpTimeSmallRangeMax=max;return;};
  void SetExpTimeBinWidth(Float_t width){fExpTimeBinWidth=width;return;};
 private: 
  UInt_t fRunNumber; //run number
  AliESDEvent *fESD;    //ESD object
  AliAnalysisFilter *fTrackFilter; //track filter object
  AliESDVertex *fVertex; //pointer to the vertex object
  AliESDpid *fESDpid; //pointer to the PID object
    
  Int_t fNTOFtracks; //number of tracks matching with TOF
  //Int_t fNPrimaryTracks; //number of primary tracks
  Float_t fT0[3]; //event time
  Float_t fSigmaSpecie[5]; //number of TOF PID sigmas, ie.fSigmaPion, fSigmaKaon, fSigmaProton;
  Double_t fTrkExpTimes[5]; //expected times from tracking for 5 mass hypothesis
  Double_t fThExpTimes[5]; //theoretical expected times for 5 mass hypothesis
  Bool_t fEnableAdvancedCheck; //flag to enable advanced checks
  Float_t fExpTimeBinWidth;//bin width for t-texp histos
  Float_t fExpTimeRangeMin, fExpTimeRangeMax; //range of t-texp histogram
  Float_t fExpTimeSmallRangeMin, fExpTimeSmallRangeMax; //reduced range of t-texp histogram

  //output objects
  TList *fHlist;  //list of general histos
  TList *fHlistTimeZero; //list of timeZero related histos
  TList *fHlistPID; //list of PID-related histos
  TList *fHpos;  //list of general histos for positive tracks
  TList *fHneg;  //list of general histos for negative tracks
  

  ClassDef(AliAnalysisTaskTOFqa, 3); // example of analysis
};

#endif
