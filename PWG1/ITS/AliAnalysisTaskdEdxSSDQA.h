#ifndef AliAnalysisTaskdEdxSSDQA_cxx
#define AliAnalysisTaskdEdxSSDQA_cxx

#include "AliAnalysisTaskSE.h"

class TH1F;
class TH2F;
class TList;

class AliAnalysisTaskdEdxSSDQA : public AliAnalysisTaskSE {

 public:
  AliAnalysisTaskdEdxSSDQA(const char *name = "AliAnalysisTaskdEdxSSDQA");
  virtual ~AliAnalysisTaskdEdxSSDQA() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  virtual void   LocalInit();
  
  void SetPcut(Float_t pcut){ fPcut=pcut;}
  Float_t GetPcut() const{return fPcut;}
 private:

  TH2F*   fHist1;         // signal differenc between ITSsignal stored in ESD and calcuated after using AliTrackPoint
  TH2F*   fHist2;         // landau distributions for each module	  
  TList*  fListOfHistos;  // output list	
  Float_t fPcut;          // Momentum cut



 AliAnalysisTaskdEdxSSDQA(const AliAnalysisTaskdEdxSSDQA&); // not implemented
 AliAnalysisTaskdEdxSSDQA& operator=(const AliAnalysisTaskdEdxSSDQA&); // not implemented
 ClassDef(AliAnalysisTaskdEdxSSDQA, 1); // example of analysis
};

#endif
