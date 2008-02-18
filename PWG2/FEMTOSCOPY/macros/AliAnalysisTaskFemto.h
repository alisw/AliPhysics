#include "TH1.h"

#include "AliESDEvent.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"

#include "AliFemtoEventReaderESDChain.h"
#include "AliFemtoManager.h"

class AliAnalysisTaskFemto : public AliAnalysisTask {
 public:
  AliAnalysisTaskFemto() : AliAnalysisTask(), fESD(0), fHistPt(0) {}
  AliAnalysisTaskFemto(const char *name);
  virtual ~AliAnalysisTaskFemto() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetFemtoReader(AliFemtoEventReaderESDChain *aReader);
  void SetFemtoManager(AliFemtoManager *aManager);
  //  void SetFriendAddress(AliESDfriend **aFriendAddress);

 private:
  AliESDEvent   *fESD;       // ESD object
  TH1F          *fHistPt;    // Pt spectrum
  AliFemtoEventReaderESDChain *fReader;   // AliFemto reader for ESD given by the chain
  AliFemtoManager             *fManager;  // AliFemto top-level manager 

  ClassDef(AliAnalysisTaskFemto, 1); // example of analysis
};
