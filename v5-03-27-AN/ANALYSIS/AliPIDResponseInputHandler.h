//
// Class AliPIDResponseInputHandler
//
// AliPIDResponseInputHandler
// TODO example
// authors:
//        Jens Wiechula (jens.wiechula@cern.ch)
//        Martin Vala (martin.vala@cern.ch)
//

#ifndef ALIPIDRESPONSEINPUTHANDLER_H
#define ALIPIDRESPONSEINPUTHANDLER_H

#include "AliInputEventHandler.h"
class AliPIDResponse;
class AliMultiInputEventHandler;

class AliPIDResponseInputHandler : public AliInputEventHandler {
  
public:
  AliPIDResponseInputHandler(const char *name = "PIDResoponseIH");
  virtual ~AliPIDResponseInputHandler();
  
   // From the interface
  virtual Bool_t  Init(Option_t *opt);
  virtual Bool_t  Init(TTree *tree, Option_t *opt);
  virtual Bool_t  BeginEvent(Long64_t entry);
  virtual Bool_t  FinishEvent();
  virtual Bool_t  Notify();
  virtual Bool_t  Notify(const char *path);
  virtual Bool_t  GetEntry();
  
  void SetIsMC(Bool_t isMC=kTRUE) { fIsMC=isMC; }
private:
  
  Bool_t fIsMC;                        //  If we run on MC data
  
  AliPIDResponse *fPIDResponse;        //! PID response Handler
  Int_t   fRun;                        //! current run number
  Int_t   fOldRun;                     //! current run number
  Int_t   fRecoPass;                   //! reconstruction pass
  
  AliMultiInputEventHandler *fMCurrentMutliIH;  //! input handler
  
  //
  void SetRecoInfo();
  
  
  AliPIDResponseInputHandler(const AliPIDResponseInputHandler& handler);
  AliPIDResponseInputHandler &operator=(const AliPIDResponseInputHandler &handler);
  
  ClassDef(AliPIDResponseInputHandler, 1)
    
};

#endif
