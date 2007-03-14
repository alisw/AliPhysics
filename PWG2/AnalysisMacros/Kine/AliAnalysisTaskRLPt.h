#include "TH1.h"

#include "AliESD.h"

#include "AliAnalysisTaskRL.h"

class AliAnalysisTaskRLPt : public AliAnalysisTaskRL {
 public:
  AliAnalysisTaskRLPt(const char *name);
  virtual ~AliAnalysisTaskRLPt() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:
  AliESD *fESD; //ESD object
  TH1F   *fHistPt; //Pt spectrum
   
  ClassDef(AliAnalysisTaskRLPt, 0); // example of analysis
};

