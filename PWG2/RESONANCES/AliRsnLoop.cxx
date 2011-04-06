//
// Base class to implement any computation within the RSN package.
// It contains only an array of output objects which must derive 
// from AliRsnOutput.
// Its core functions ar Init() and DoLoop() which must be 
// overloaded by any class which inherits from this.
//

#include <TList.h>
#include <TEntryList.h>

#include "AliLog.h"

#include "AliRsnDaughterSelector.h"

#include "AliRsnLoop.h"

ClassImp(AliRsnLoop)

//_____________________________________________________________________________
AliRsnLoop::AliRsnLoop(const char *name, Bool_t isMixed) :
   TNamed(name, ""),
   fIsMixed(isMixed),
   fEventCuts(0x0),
   fOutputs("AliRsnListOutput", 0)
{
//
// Default constructor
//
}

//_____________________________________________________________________________
AliRsnLoop::AliRsnLoop(const AliRsnLoop& copy) :
   TNamed(copy),
   fIsMixed(copy.fIsMixed),
   fEventCuts(copy.fEventCuts),
   fOutputs(copy.fOutputs)
{
//
// Default constructor
//
}

//_____________________________________________________________________________
AliRsnLoop& AliRsnLoop::operator=(const AliRsnLoop& copy)
{
//
// Assignment operator
//

   fIsMixed = copy.fIsMixed;
   fEventCuts = copy.fEventCuts;
   fOutputs = copy.fOutputs;
   return (*this);
}

//_____________________________________________________________________________
AliRsnLoop::~AliRsnLoop()
{
//
// Destructor
//
}

//_____________________________________________________________________________
void AliRsnLoop::AddOutput(TObject *object)
{
//
// Adds an object to any of the collections.
// The target depends on the object type.
// Returns kFALSE if the addition failed.
//

   //fOutputs.AddLast(out);
   AliRsnListOutput *out = (AliRsnListOutput*)object;
   Int_t n = fOutputs.GetEntries();
   new (fOutputs[n]) AliRsnListOutput(*out);
}

//_____________________________________________________________________________
void AliRsnLoop::Print(Option_t*) const
{
//
// Prints info about pair
//
   
   TObjArrayIter next(&fOutputs);
   AliRsnListOutput *out = 0x0;
   
   while ( (out = (AliRsnListOutput*)next()) ) {
      out->Print();
   }
}

//_____________________________________________________________________________
Bool_t AliRsnLoop::OkEvent(AliRsnEvent *rsn)
{
//
// If event cuts are defined, check event against them
//

   if (fEventCuts) 
      return fEventCuts->IsSelected(rsn);
   else
      return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliRsnLoop::Init(const char *prefix, TList *list)
{
//
// Initialization function.
// Loops on all outputs, and initialize each of them.
// Returns kTRUE only if all initializations were successful.
//

   TObjArrayIter next(&fOutputs);
   AliRsnListOutput *out;
   Bool_t globalOK = kTRUE;
   
   while ( (out = (AliRsnListOutput*)next()) ) {
      globalOK = globalOK && out->Init(prefix, list);
   }
   
   AliInfo(Form("[%s] Object initialization: %s", GetName(), (globalOK ? "successful" : "failed")));
   return globalOK;
}

//_____________________________________________________________________________
Int_t AliRsnLoop::DoLoop
(AliRsnEvent *, AliRsnDaughterSelector *, AliRsnEvent *, AliRsnDaughterSelector *)
{
//
// Main loop.
// Performs all the computations, looping on the passed event(s) and using the lists
// of selected daughters which are provided, for allowing the user to choose what to do
// with them.
//

   AliWarning("Implement this method in derived class");
   return 0;
}
