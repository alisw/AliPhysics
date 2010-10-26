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

 private: 
  UInt_t fRunNumber; //run number
  AliESDEvent *fESD;    //ESD object
  AliAnalysisFilter *fTrackFilter; //track filter object
  AliESDVertex *fVertex; //pointer to the vertex object
  
  Int_t fNTOFtracks; //number of tracks matching with TOF
  Int_t fNPrimaryTracks; //number of primary tracks
  Float_t fT0; //event time
  
  //output objects
  TList *fHlist;  //list of general histos
  TList *fHlistExperts; //list of experts histos

  ClassDef(AliAnalysisTaskTOFqa, 1); // example of analysis
};

#endif
