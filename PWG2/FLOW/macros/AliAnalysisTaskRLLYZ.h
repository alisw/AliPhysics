#ifndef AliAnalysisTaskRLLYZ_H
#define AliAnalysisTaskRLLYZ_H

 
#include "AliAnalysisTaskRL.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "TObjArray.h"
#include "TList.h"
#include "AliESD.h"
//#include "AliFlowConstants.h"
//#include "AliFlowLYZConstants.h"
#include "AliFlowSelection.h"
#include "AliFlowMaker.h"
#include "AliFlowLeeYangZerosMaker.h"
#include "AliFlowEvent.h"

 

class AliAnalysisTaskRLLYZ : public AliAnalysisTaskRL {
 public:
  AliAnalysisTaskRLLYZ(const char *name, Bool_t firstrun);
  virtual ~AliAnalysisTaskRLLYZ();
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  //lyz flags
  void           SetFirstRunLYZ(Bool_t kt)    { this->fFirstRunLYZ = kt ;  }
  Bool_t         GetFirstRunLYZ() const       { return this->fFirstRunLYZ ; }
  void           SetUseSumLYZ(Bool_t kt)      { this->fUseSumLYZ = kt ;  }
  Bool_t         GetUseSumLYZ() const         { return this->fUseSumLYZ ; }


 private:

  TFile*                    fFirstRunFile ;   //! pointer to file from first run
  AliESD*                   fESD;             //! ESD object
  AliFlowEvent*             fFlowEvent;       //! flowevent object
  AliFlowSelection*         fFlowSelect;      //! flowselection object
  AliFlowMaker*             fFlowMaker;       //! flowmaker object
  AliFlowLeeYangZerosMaker* fFlowLYZ;         //! lyz maker object

  //flags
  Bool_t                    fFirstRunLYZ ;    //! flag for lyz analysis 
  Bool_t                    fUseSumLYZ ;      //! flag for lyz analysis 
  
   
  ClassDef(AliAnalysisTaskRLLYZ, 0);          // lyz analysis
};

 #endif
