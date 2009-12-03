#ifndef ALIANALYSISTASKGLOBALQA_H
#define ALIANALYSISTASKGLOBALQA_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//                 AliAnalysisTaskGlobalQA class
// This task is for running the GlobalQA over already existing ESDs
//          Origin:  I.Belikov, Iouri.Belikov@cern.ch, June 2009
//-----------------------------------------------------------------

#include <TObjArray.h>
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskGlobalQA : public AliAnalysisTaskSE {
public:
  enum {
    kClr0,kClr1,kClr2,kClr3,
    kTrk0,kTrk1,kTrk2,kTrk3,kTrk4,kTrk5,kTrk6,kTrk7,kTrk8,kTrk9,kTrk10,
    kK0on,kK0off,kL0on,kL0off,
    kPid0,kPid1,kPid2,kPid3,
    kMlt0,kMlt1,
    kLast
  };
  AliAnalysisTaskGlobalQA();
  virtual ~AliAnalysisTaskGlobalQA() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
private:
  void Add2ESDsList(TObject *o,Int_t index) {fArrayQA->AddAt(o,index);}
  TH1 *GetESDsData(Int_t i) {return dynamic_cast<TH1*>(fArrayQA->At(i));}

  TObjArray   *fArrayQA;       //! the array of output histos
   
  AliAnalysisTaskGlobalQA(const AliAnalysisTaskGlobalQA&);
  AliAnalysisTaskGlobalQA& operator=(const AliAnalysisTaskGlobalQA&);
  
  ClassDef(AliAnalysisTaskGlobalQA, 1); // GlobalQA analysis
};

#endif
