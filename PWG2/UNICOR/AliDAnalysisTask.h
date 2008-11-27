// Author: Dariusz Miskowiec 2007

#ifndef ALIDANALYSISTASK_H
#define ALIDANALYSISTASK_H

#include "AliESDEvent.h"
#include "AliAnalysisTask.h"
class AliDEventAliceESD;
class AliDAnalGlobal;
class AliDAnalSingle;
class AliDAnalCorrel;
class AliDAnalPtfluc;

/*****************************************************************************/
class AliDAnalysisTask : public AliAnalysisTask {
   
 public:                                        
  AliDAnalysisTask();                            // constructor
  virtual ~AliDAnalysisTask() {}                 // destructor
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  virtual void   LocalInit() {}
  virtual Bool_t Notify() {return kTRUE;}
  virtual Bool_t NotifyBinChange() {return kTRUE;}
  virtual void   FinishTaskOutput() {}

 protected:
  AliESDEvent    *fESD;                       //! ESD event
  AliDEventAliceESD *fEv0;                       //! data/analysis interface
  AliDEventAliceESD *fEv1;                       //! another for event mixing
  AliDAnalGlobal    *fDag;                       //! global analysis
  AliDAnalSingle    *fAll;                       //! single analysis
  AliDAnalSingle    *fPim;                       //! single analysis
  AliDAnalSingle    *fPip;                       //! single analysis
  AliDAnalCorrel    *fCnn;                       //! correlation analysis pi-pi-
  AliDAnalCorrel    *fCpp;                       //! correlation analysis pi+pi+
  AliDAnalPtfluc    *fPtf;                       //! pt-fluctuation analysis
  TList          *fOutputList;                //  list of output objects
  AliDAnalysisTask(const AliDAnalysisTask&); 
  AliDAnalysisTask& operator=(const AliDAnalysisTask&); 

  ClassDef(AliDAnalysisTask,1)
};
/*****************************************************************************/

#endif
