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
    AliAnalysisTaskScale() : AliAnalysisTaskSE(), fTracksName("tracks"), fClustersName("clusters"), fESD(0),
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
  AliESDEvent           *fESD;                    //!ESD object
  TList                 *fOutputList;             //!Output list
  TH1F                  *fHistCentrality;         //!
  TH2F                  *fHistPtTPCvsCent;        //!
  TH2F                  *fHistPtEMCALvsCent;      //!
  TH2F                  *fHistEtvsCent;           //!
  TH2F                  *fHistScalevsCent;        //!
  TH2F                  *fHistDeltaScalevsCent;   //!
  TH2F                  *fHistPtTPCvsNtrack;      //!
  TH2F                  *fHistPtEMCALvsNtrack;    //!
  TH2F                  *fHistEtvsNtrack;         //!
  TH2F                  *fHistScalevsNtrack;      //!
  TH2F                  *fHistDeltaScalevsNtrack; //!

  AliAnalysisTaskScale(const AliAnalysisTaskScale&); // not implemented
  AliAnalysisTaskScale& operator=(const AliAnalysisTaskScale&); // not implemented
  
  ClassDef(AliAnalysisTaskScale, 1); // example of analysis
};
#endif
