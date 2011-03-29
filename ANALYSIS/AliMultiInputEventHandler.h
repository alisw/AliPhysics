//
// Class AliMultiInputEventHandler
//
// Multi input event handler
// TODO example
// author:
//        Martin Vala (martin.vala@cern.ch)
//

#ifndef ALIMULTIINPUTEVENTHANDLER_H
#define ALIMULTIINPUTEVENTHANDLER_H

#include <TObjArray.h>

#include "AliInputEventHandler.h"

class AliMCEventHandler;
class AliMultiInputEventHandler : public AliInputEventHandler {

public:
   AliMultiInputEventHandler(const Int_t size = 0, const char *name = "name");
   virtual ~AliMultiInputEventHandler();

   // From the interface
   virtual Bool_t  Init(Option_t *opt);
   virtual Bool_t  Init(TTree *tree, Option_t *opt);
   virtual Bool_t  BeginEvent(Long64_t entry);
   virtual Bool_t  GetEntry();
   virtual Bool_t  FinishEvent();
   virtual Bool_t  Notify();
   virtual Bool_t  Notify(const char *path);
   virtual UInt_t  IsEventSelected();
   // add input handler
   virtual void                AddInputEventHandler(AliVEventHandler*inHandler);
   AliVEventHandler           *InputEventHandler(const Int_t index);
   Int_t                       BufferSize() const { return fBufferSize; }
   TObjArray                  *InputEventHandlers() { return &fInputHandlers; }

   AliInputEventHandler       *GetFirstInputEventHandler();
   AliMCEventHandler          *GetFirstMCEventHandler();
   AliMultiInputEventHandler  *GetFirstMultiInputHandler();

   Option_t                   *GetDataType() const;
   
   //PID response
   virtual AliPIDResponse* GetPIDResponse();
   virtual void CreatePIDResponse(Bool_t isMC);

protected:

   Int_t                   fBufferSize;    // Size of the buffer
   TObjArray               fInputHandlers; // buffer of input handlers
   Option_t               *fAnalysisType;  //! local, proof, grid
private:
   AliMultiInputEventHandler(const AliMultiInputEventHandler& handler);
   AliMultiInputEventHandler &operator=(const AliMultiInputEventHandler &handler);

   ClassDef(AliMultiInputEventHandler, 1)
};

#endif
