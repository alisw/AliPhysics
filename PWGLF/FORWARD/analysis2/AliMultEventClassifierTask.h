#ifndef ALIMULTEVENTCLASSIFIERTASK_H
#define ALIMULTEVENTCLASSIFIERTASK_H
#include <AliAnalysisTaskSE.h>
#include "AliMultEventClassifier.h"
#include "AliAODMultEventClass.h"

class AliMultEventClassifierTask : public AliAnalysisTaskSE
{
public:
  /**
   * Default CTOR - for ROOT I/O only
   */
  AliMultEventClassifierTask();
  /**
   * User CTOR
   *
   * @param name Name of task 
   */
  AliMultEventClassifierTask(const char* name);
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliMultEventClassifierTask(const AliMultEventClassifierTask& o);
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   *
   * @return reference to this object
   */
  AliMultEventClassifierTask& operator=(const AliMultEventClassifierTask& o);
  /** 
   * Connect to train 
   * 
   * @param sumFile Output file 
   * 
   * @return true on usccess 
   */
  Bool_t Connect(const char* sumFile);
  /** 
   * Create output objects 
   * 
   */
  void UserCreateOutputObjects();
  /** 
   * Process an event 
   * 
   * @param option  Note used
   */
  void UserExec(Option_t* option="");
  /** 
   * Called at end. 
   * 
   * @param option Not used
   */
  void Terminate(Option_t* option);
protected:
  AliMultEventClassifier fClassifier;
  AliAODMultEventClass   fData;
  TList*                 fList;

  ClassDef(AliMultEventClassifierTask,1);
};

#endif
// Local Variables:
//  mode: C++
// End:
