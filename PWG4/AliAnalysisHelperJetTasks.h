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
  static void PrintStack(AliMCEvent *mcEvent,Int_t iFirst = 0,Int_t iLast = 0,Int_t iMaxPrint = 10);
  

  private:
  
  ClassDef(AliAnalysisHelperJetTasks, 1) // 
};

#endif // ALIANALYSISHELPERJETTASKS_H
