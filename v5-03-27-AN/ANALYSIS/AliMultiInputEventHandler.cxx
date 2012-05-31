//
// Class AliMultiInputEventHandler
//
// Multi input event handler
// TODO example
// author:
//        Martin Vala (martin.vala@cern.ch)
//

#include "AliLog.h"
#include "AliMCEventHandler.h"

#include "AliMultiInputEventHandler.h"

ClassImp(AliMultiInputEventHandler)

static Option_t *gCurrentMultiDataType = "ESD";

//_____________________________________________________________________________
AliMultiInputEventHandler::AliMultiInputEventHandler(const Int_t size, const char *name) :
   AliInputEventHandler(name, name),
   fBufferSize(size),
   fInputHandlers(),
   fAnalysisType(0)
{
//
// Default constructor.
//
   AliDebug(AliLog::kDebug + 10, "<-");
   fInputHandlers.SetOwner(kTRUE);
   AliDebug(AliLog::kDebug + 10, "->");
}

//_____________________________________________________________________________
AliMultiInputEventHandler::~AliMultiInputEventHandler()
{
   //
   // Destructor
   //
   AliDebug(AliLog::kDebug + 10, "<-");
   AliDebug(AliLog::kDebug + 10, "->");
}


//_____________________________________________________________________________
AliVEventHandler *AliMultiInputEventHandler::InputEventHandler(const Int_t index)
{
   //
   // Returns input handler
   //
   AliDebug(AliLog::kDebug + 5, Form("<-"));
   if ((index >= 0) && (index < fBufferSize)) {
      return (AliVEventHandler *) fInputHandlers.At(index);
   }
   AliDebug(AliLog::kDebug + 5, Form("->"));
   return 0;
}

//_____________________________________________________________________________
void AliMultiInputEventHandler::AddInputEventHandler(AliVEventHandler*inHandler)
{
   //
   // Create N (fBufferSize) copies of input handler
   //
   if (inHandler->InheritsFrom("AliESDInputHandler")) gCurrentMultiDataType = "ESD";
   if (inHandler->InheritsFrom("AliAODInputHandler")) gCurrentMultiDataType = "AOD";
   AliDebug(AliLog::kDebug + 5, Form("<-"));
   AliDebug(AliLog::kDebug + 5, Form("Creating %d input event handlers ...", fBufferSize));
   AliDebug(AliLog::kDebug + 5, Form("Adding input handler with index %d ...", fBufferSize));
   fInputHandlers.Add(inHandler);
   fBufferSize++;
   AliDebug(AliLog::kDebug + 5, Form("->"));
}

//_____________________________________________________________________________
Bool_t AliMultiInputEventHandler::Init(Option_t *opt)
{
   //
   // Init() is called for all mix input handlers.
   //
   fAnalysisType = opt;
   AliDebug(AliLog::kDebug + 5, Form("<- \"%s\"", opt));
   AliInputEventHandler *eh = 0;
   TObjArrayIter next(&fInputHandlers);
   while ((eh = (AliInputEventHandler *) next())) {
      eh->Init(fAnalysisType);
   }
   AliDebug(AliLog::kDebug + 5, Form("->"));
   return AliInputEventHandler::Init(opt);
}
//_____________________________________________________________________________
Bool_t AliMultiInputEventHandler::Init(TTree *tree, Option_t *opt)
{
   //
   // Init(const char*path) is called for all mix input handlers.
   // Create event pool if needed
   //
   fAnalysisType = opt;
   AliDebug(AliLog::kDebug + 5, Form("<- %p %s", (void *) tree, tree->GetName()));
   if (!tree) {
      AliError(Form("-> tree is null"));
      return kFALSE;
   }
   AliInputEventHandler *eh = 0;
   TObjArrayIter next(&fInputHandlers);
   while ((eh = (AliInputEventHandler *) next())) {
      // using mixing input hadnler from Base class
      // for me fParentHandler would be better name
      eh->SetParentHandler(this);
      eh->Init(tree, fAnalysisType);
   }
   AliDebug(AliLog::kDebug + 5, Form("->"));
   return AliInputEventHandler::Init(tree, opt);
}
//_____________________________________________________________________________
Bool_t AliMultiInputEventHandler::Notify()
{
   //
   // Notify() is called for all mix input handlers
   //
   AliDebug(AliLog::kDebug + 5, Form("<-"));
   AliInputEventHandler *eh = 0;
   TObjArrayIter next(&fInputHandlers);
   while ((eh = (AliInputEventHandler *) next())) {
      eh->Notify();
   }
   AliDebug(AliLog::kDebug + 5, Form("->"));
   return AliInputEventHandler::Notify();
}

//_____________________________________________________________________________
Bool_t AliMultiInputEventHandler::Notify(const char *path)
{
   //
   // Notify(const char*path) is called for all mix input handlers
   //
   AliDebug(AliLog::kDebug + 5, Form("<- %s", path));
   AliInputEventHandler *eh = 0;
   TObjArrayIter next(&fInputHandlers);
   while ((eh = (AliInputEventHandler *) next())) {
      eh->Notify(path);
   }
   AliDebug(AliLog::kDebug + 5, Form("->"));
//   return AliInputEventHandler::Notify(path);
   return AliInputEventHandler::Notify(path);
}
//_____________________________________________________________________________
Bool_t AliMultiInputEventHandler::BeginEvent(Long64_t entry)
{
   //
   // BeginEvent(Long64_t entry) is called for all mix input handlers
   //
   AliDebug(AliLog::kDebug + 5, Form("<- %lld", entry));
   AliInputEventHandler *eh = 0;
   TObjArrayIter next(&fInputHandlers);
   while ((eh = (AliInputEventHandler *) next())) {
      eh->BeginEvent(entry);
   }
   GetEntry();
   AliDebug(AliLog::kDebug + 5, Form("->"));
   return AliInputEventHandler::BeginEvent(entry);
}


//_____________________________________________________________________________
Bool_t AliMultiInputEventHandler::GetEntry()
{
   //
   // Sets correct events to every mix events
   //
   AliDebug(AliLog::kDebug + 5, Form("<-"));
   AliInputEventHandler *eh = 0;
   TObjArrayIter next(&fInputHandlers);
   while ((eh = (AliInputEventHandler *) next())) {
      eh->GetEntry();
   }
   AliDebug(AliLog::kDebug + 5, Form("->"));
   return AliInputEventHandler::GetEntry();
}
//_____________________________________________________________________________
Bool_t AliMultiInputEventHandler::FinishEvent()
{
   //
   // FinishEvent() is called for all mix input handlers
   //
   AliDebug(AliLog::kDebug + 5, Form("<-"));
   AliInputEventHandler *eh = 0;
   TObjArrayIter next(&fInputHandlers);
   while ((eh = (AliInputEventHandler *) next())) {
      eh->FinishEvent();
   }
   AliDebug(AliLog::kDebug + 5, Form("->"));
   return AliInputEventHandler::FinishEvent();
}

AliInputEventHandler *AliMultiInputEventHandler::GetFirstInputEventHandler()
{
   //
   // Return first InputEventHandler
   //
   AliDebug(AliLog::kDebug + 5, Form("<-"));
   AliVEventHandler *eh = 0;
   AliInputEventHandler *handler = 0;
   TObjArrayIter next(&fInputHandlers);
   while ((eh = (AliVEventHandler *) next())) {
      handler = dynamic_cast<AliInputEventHandler *>(eh);
      if (handler) return handler;
   }
   AliDebug(AliLog::kDebug + 5, Form("->"));
   return 0;
}
AliMCEventHandler *AliMultiInputEventHandler::GetFirstMCEventHandler()
{
   //
   // Return first MCEventHandler
   //
   AliDebug(AliLog::kDebug + 5, Form("<-"));
   AliVEventHandler *eh = 0;
   AliMCEventHandler *handler = 0;
   TObjArrayIter next(&fInputHandlers);
   while ((eh = (AliVEventHandler *) next())) {
      handler = dynamic_cast<AliMCEventHandler *>(eh);
      if (handler) return handler;
   }
   AliDebug(AliLog::kDebug + 5, Form("->"));
   return 0;
}

AliMultiInputEventHandler *AliMultiInputEventHandler::GetFirstMultiInputHandler()
{
   //
   // Return first MultiInputHandler
   //
   AliDebug(AliLog::kDebug + 5, Form("<-"));
   AliVEventHandler *eh = 0;
   AliMultiInputEventHandler *handler = 0;
   TObjArrayIter next(&fInputHandlers);
   while ((eh = (AliVEventHandler *) next())) {
      handler = dynamic_cast<AliMultiInputEventHandler *>(eh);
      if (handler) return handler;
   }
   AliDebug(AliLog::kDebug + 5, Form("->"));
   return 0;
}

//______________________________________________________________________________
Option_t *AliMultiInputEventHandler::GetDataType() const
{
   // Returns handled data type.
   return gCurrentMultiDataType;
}

//______________________________________________________________________________
UInt_t  AliMultiInputEventHandler::IsEventSelected() 
{
  // returns if event is selected
  
  AliInputEventHandler *firstIH = dynamic_cast<AliInputEventHandler*> (GetFirstInputEventHandler());
  if (firstIH) {
    return firstIH->IsEventSelected();
  }
  
  return fIsSelectedResult;
}

//______________________________________________________________________________
AliPIDResponse* AliMultiInputEventHandler::GetPIDResponse()
{
   // retrieve PID response
   
   AliInputEventHandler *firstIH = dynamic_cast<AliInputEventHandler*> (GetFirstInputEventHandler());
   if (firstIH) {
      return firstIH->GetPIDResponse();
   }
   
   return 0x0;
}
   
//______________________________________________________________________________
void AliMultiInputEventHandler::CreatePIDResponse(Bool_t isMC)
{
   // create PID response
   AliInputEventHandler *firstIH = dynamic_cast<AliInputEventHandler*> (GetFirstInputEventHandler());
   if (firstIH) {
      firstIH->CreatePIDResponse(isMC);
   }
}
