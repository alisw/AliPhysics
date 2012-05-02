#ifndef AliAnalysisTaskScale_cxx
#define AliAnalysisTaskScale_cxx

// $Id$

class TList;
class TH1F;
class TH1I;
class TH2F;
class AliESDEvent;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskScale : public AliAnalysisTaskSE {
 public:
    AliAnalysisTaskScale() : AliAnalysisTaskSE(), fTracksName(), fClustersName(),
      fOutputList(0), fHistCentrality(0), fHistPtTPCvsCent(0), fHistPtEMCALvsCent(0), fHistEtvsCent(0),  
      fHistScalevsCent(0),  fHistDeltaScalevsCent(0), fHistPtTPCvsNtrack(0), fHistPtEMCALvsNtrack(0), 
      fHistEtvsNtrack(0), fHistScalevsNtrack(0), fHistDeltaScalevsNtrack(0) {}
  AliAnalysisTaskScale(const char *name);
  virtual ~AliAnalysisTaskScale() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  virtual void   SetTracksName(const char *n)    { fTracksName   = n; }
  virtual void   SetClustersName(const char *n)  { fClustersName = n; }
  
 private:
  TString                fTracksName;             // name of track collection
  TString                fClustersName;           // name of clusters collection
  TList                 *fOutputList;             //!output list
  TH1F                  *fHistCentrality;         //!output histogram
  TH2F                  *fHistPtTPCvsCent;        //!output histogram
  TH2F                  *fHistPtEMCALvsCent;      //!output histogram
  TH2F                  *fHistEtvsCent;           //!output histogram
  TH2F                  *fHistScalevsCent;        //!output histogram
  TH2F                  *fHistDeltaScalevsCent;   //!output histogram
  TH2F                  *fHistPtTPCvsNtrack;      //!output histogram
  TH2F                  *fHistPtEMCALvsNtrack;    //!output histogram
  TH2F                  *fHistEtvsNtrack;         //!output histogram
  TH2F                  *fHistScalevsNtrack;      //!output histogram
  TH2F                  *fHistDeltaScalevsNtrack; //!output histogram

  AliAnalysisTaskScale(const AliAnalysisTaskScale&); // not implemented
  AliAnalysisTaskScale& operator=(const AliAnalysisTaskScale&); // not implemented
  
  ClassDef(AliAnalysisTaskScale, 2); // example of analysis
};
#endif
