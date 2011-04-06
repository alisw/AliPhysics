//
// *** Class AliRsnPair ***
//
// "Core" method for defining the work on a pari of particles.
// For one analysis, one must setup one of this for each pair he wants to analyze,
// adding to it all analysis which he desires to do.
// Here he defines the cuts, and the particle types and charges, and can add
// functions which do different operations on the same pair, and some binning
// with respect to some kinematic variables (eta, momentum)
//
// authors: A. Pulvirenti (email: alberto.pulvirenti@ct.infn.it)
//          M. Vala (email: martin.vala@cern.ch)
//

#include "AliRsnLoopEff.h"

ClassImp(AliRsnLoopEff)

//_____________________________________________________________________________
AliRsnLoopEff::AliRsnLoopEff(const char *name, Int_t nSteps) :
   AliRsnLoop(name),
   fAddSteps(nSteps),
   fSteps(0)
{
//
// Default constructor
//
}

//_____________________________________________________________________________
AliRsnLoopEff::AliRsnLoopEff(const AliRsnLoopEff& copy) :
   AliRsnLoop(copy),
   fAddSteps(copy.fAddSteps),
   fSteps(copy.fSteps)
{
//
// Copy constructor
//
}

//_____________________________________________________________________________
AliRsnLoopEff& AliRsnLoopEff::operator=(const AliRsnLoopEff& copy)
{
   AliRsnLoop::operator=(copy);

   fAddSteps = copy.fAddSteps;
   fSteps = copy.fSteps;

   return (*this);
}

//_____________________________________________________________________________
AliRsnLoopEff::~AliRsnLoopEff()
{
//
// Destructor
//

   fSteps.Delete();
}

//_____________________________________________________________________________
void AliRsnLoopEff::CreateOutput()
{
//
// Create the unique output object of this loop
//

   if (!fOutputs.IsEmpty()) {
      AliInfo("Clearing container of this efficiency loop.");
      fOutputs.Delete();
   }

   AliRsnListOutput out(Form("%s_out", GetName()), AliRsnListOutput::kCFContainer);
   out.SetSteps(NStepsAll());
   
   AddOutput(&out);
}

//_____________________________________________________________________________
void AliRsnLoopEff::AddStep(TObject *cuts)
{
//
// Add a step on reconstruction
//

   fSteps.AddLast(cuts);
}

//_____________________________________________________________________________
Bool_t AliRsnLoopEff::Init(const char *prefix, TList *list)
{
//
// Initialization function.
// Loops on all functions and eventual the ntuple, to initialize output objects.
//

   CreateOutput();
   return AliRsnLoop::Init(Form("%s_%s", prefix, GetName()), list);
}

//_____________________________________________________________________________
Int_t AliRsnLoopEff::FindTrack(Int_t label, AliVEvent *event)
{
//
// Loops an event and find all tracks which have a label
// equal to that passed as first argument.
//

   Int_t   i = 0;
   Int_t   ntracks = event->GetNumberOfTracks();
   TArrayI array(100);

   for (i = 0; i < ntracks; i++) {
      AliVParticle *track = event->GetTrack(i);
      if (TMath::Abs(track->GetLabel()) != label) continue;
      return i;
   }
   
   return -1;
}
