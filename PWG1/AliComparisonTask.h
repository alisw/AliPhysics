#ifndef ALICOMPARISONRESTASK_H
#define ALICOMPARISONRESTASK_H

//------------------------------------------------------------------------------
// Class to compare properties of reconstructed and MC particle tracks. 
// 
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

class AliComparisonObject;
class AliMagF;
class TList;

#include "AliAnalysisTask.h"

class AliComparisonTask : public AliAnalysisTask {
 public:
  AliComparisonTask(const char *name = "AliComparisonTask");
  virtual ~AliComparisonTask();
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  // Read TTree entry (event by event)
  Bool_t  ReadEntry(Int_t evt);

  // Set comparison objects
  Bool_t AddComparisonObject(AliComparisonObject* comp);

 private:
  TTree* fTree;                   //! input tree
  AliMCInfo *fInfoMC;             //! AliMCInfo object
  AliESDRecInfo *fInfoRC;         //! AliESDRecInfo object

  TList* fOutput;                 //! list send on output slot 0
  static Int_t fEvtNumber;        //! event number
  TIterator *fPitList;            //! iterator over the output objetcs  
  TList *fCompList;               // list of comparison objects

  AliComparisonTask(const AliComparisonTask&); // not implemented
  AliComparisonTask& operator=(const AliComparisonTask&); // not implemented
  
  ClassDef(AliComparisonTask, 1); // example of analysis
};

#endif
