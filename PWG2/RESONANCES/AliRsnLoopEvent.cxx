//
// Computator for single daughters.
// Implements a simple loop on tracks from one of the entry lists
// filled by the task AliRsnInputHandler, adding a check on their
// definition specified in the daughter def.
//

#include "AliLog.h"

#include "AliRsnEvent.h"

#include "AliRsnLoopEvent.h"

ClassImp(AliRsnLoopEvent)

//_____________________________________________________________________________
AliRsnLoopEvent::AliRsnLoopEvent(const char *name) :
   AliRsnLoop(name)
{
//
// Default constructor
//
}

//_____________________________________________________________________________
AliRsnLoopEvent::AliRsnLoopEvent(const AliRsnLoopEvent& copy) :
   AliRsnLoop(copy)
{
//
// Copy constructor
//
}

//_____________________________________________________________________________
AliRsnLoopEvent& AliRsnLoopEvent::operator=(const AliRsnLoopEvent& copy)
{
//
// Assignment operator
//

   AliRsnLoop::operator=(copy);

   return (*this);
}

//_____________________________________________________________________________
AliRsnLoopEvent::~AliRsnLoopEvent()
{
//
// Destructor
//
}

//_____________________________________________________________________________
void AliRsnLoopEvent::Print(Option_t* /*option*/) const
{
//
// Prints info about pair
//

   AliRsnLoop::Print();
}

//_____________________________________________________________________________
Bool_t AliRsnLoopEvent::Init(const char *prefix, TList *list)
{
//
// Initialization function.
// Loops on all functions and eventual the ntuple, to initialize output objects.
//

   return AliRsnLoop::Init(Form("%s_%s", prefix, GetName()), list);
}

//_____________________________________________________________________________
Int_t AliRsnLoopEvent::DoLoop
(AliRsnEvent *evMain, AliRsnDaughterSelector *, AliRsnEvent *, AliRsnDaughterSelector *)
{
//
// Loop function.
// Computes what is needed from passed events.
// Returns the number of pairs successfully processed.
//

   TObjArrayIter next(&fOutputs);
   AliRsnListOutput *out = 0x0;
   
   // check cuts
   if (!OkEvent(evMain)) return 0;
   
   while ( (out = (AliRsnListOutput*)next()) ) {
      out->Fill(evMain);
   }
   
   return 1;
}
