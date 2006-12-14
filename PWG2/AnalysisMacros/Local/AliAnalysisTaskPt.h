#include "TH1.h"

#include "AliESD.h"

#include "AliAnalysisTask.h"

class AliAnalysisTaskPt : public AliAnalysisTask {
 public:
  AliAnalysisTaskPt(const char *name);
  virtual ~AliAnalysisTaskPt() {}
  
  virtual void   Init(Option_t *);
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:
  AliESD *fESD; //ESD object
  TH1F   *fHistPt; //Pt spectrum
   
  ClassDef(AliAnalysisTaskPt, 0); // example of analysis
};

