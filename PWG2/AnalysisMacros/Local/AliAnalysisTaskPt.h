#include "TH1.h"

#include "AliESD.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

class AliAnalysisTaskPt : public AliAnalysisTask {
 public:
  AliAnalysisTaskPt() : AliAnalysisTask(), fESD(0), fHistPt(0) {}
  AliAnalysisTaskPt(const char *name);
  virtual ~AliAnalysisTaskPt() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:
  AliESD *fESD; //ESD object
  TH1F   *fHistPt; //Pt spectrum
   
  ClassDef(AliAnalysisTaskPt, 1); // example of analysis
};
