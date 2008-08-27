#ifndef  ALIANALYSISHELPERJETTASKS_H
#define  ALIANALYSISHELPERJETTASKS_H


#include "TObject.h"
class AliMCEvent;
class AliGenPythiaEventHeader;

// Helper Class that contains a lot of usefull static functions (i.e. for Flavor selection.

class AliAnalysisHelperJetTasks : public TObject {
 public:
  AliAnalysisHelperJetTasks() : TObject() {;}
  virtual ~AliAnalysisHelperJetTasks(){;}
  
  static AliGenPythiaEventHeader*  GetPythiaEventHeader(AliMCEvent *mcEvent);
  

  private:
  
  ClassDef(AliAnalysisHelperJetTasks, 1) // 
};

#endif // ALIANALYSISHELPERJETTASKS_H
