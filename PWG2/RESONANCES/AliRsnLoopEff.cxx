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

#include "Riostream.h"
#include "AliLog.h"
#include "AliRsnLoopEff.h"

ClassImp(AliRsnLoopEff)

//_____________________________________________________________________________
AliRsnLoopEff::AliRsnLoopEff(const char *name, Int_t nSteps, Double_t maxDist) :
   AliRsnLoop(name),
   fAddSteps(nSteps),
   fSteps(0),
   fOutput(0),
   fMaxDistPV(maxDist)
{
//
// Default constructor
//

   fVertex[0] = fVertex[1] = fVertex[2] = 0.0;
}

//_____________________________________________________________________________
AliRsnLoopEff::AliRsnLoopEff(const AliRsnLoopEff& copy) :
   AliRsnLoop(copy),
   fAddSteps(copy.fAddSteps),
   fSteps(copy.fSteps),
   fOutput(copy.fOutput)
{
//
// Copy constructor
//

   fVertex[0] = fVertex[1] = fVertex[2] = 0.0;
}

//_____________________________________________________________________________
AliRsnLoopEff& AliRsnLoopEff::operator=(const AliRsnLoopEff& copy)
{
   AliRsnLoop::operator=(copy);

   fAddSteps = copy.fAddSteps;
   fSteps = copy.fSteps;
   fOutput = copy.fOutput;

   return (*this);
}

//_____________________________________________________________________________
AliRsnLoopEff::~AliRsnLoopEff()
{
//
// Destructor
//

   fSteps.Delete();
   delete fOutput;
}

//_____________________________________________________________________________
void AliRsnLoopEff::CreateOutput()
{
//
// Create the unique output object of this loop
//

   fOutput = new AliRsnListOutput(Form("%s_out", GetName()), AliRsnListOutput::kCFContainer);
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

   if (!fOutputs.IsEmpty()) {
      AliInfo("Clearing container of this efficiency loop.");
      fOutputs.Delete();
   }
   
   Int_t nSteps = (Int_t)fSteps.GetEntries();
   nSteps += fAddSteps;

   fOutput->SetSteps(nSteps);
   fOutput->SetSkipFailed(kFALSE);
   AliRsnLoop::AddOutput(fOutput);
   
   if (AliRsnLoop::Init(Form("%s_%s", prefix, GetName()), list)) {
      fOutput = (AliRsnListOutput*)fOutputs[0];
      return kTRUE;
   } else {
      fOutput = 0x0;
      return kFALSE;
   }
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

//__________________________________________________________________________________________________
Int_t AliRsnLoopEff::GetMatchedDaughter(Int_t label, AliRsnEvent *event)
{
//
// Searches an object among all possible daughters which matches the corresponding label
// and if it is found, assigns to the daughter and returns it
//

   if (!event) return -1;
   
   AliRsnDaughter out;
   
   Int_t i, imax = event->GetAbsoluteSum();
   for (i = 0; i < imax; i++) {
      event->SetDaughter(out, i);
      if (out.IsOK() && out.GetLabel() == label) return i;
   }
   
   return -1;
}

//__________________________________________________________________________________________________
Double_t AliRsnLoopEff::DistanceFromPV(Double_t x, Double_t y, Double_t z)
{
//
// Compute distance from current primary vertex
//

   AliDebugClass(1, Form("Vertex = %.3f %.3f %.3f -- vprod = %.3f %.3f %.3f", fVertex[0], fVertex[1], fVertex[2], x, y, z));

   x -= fVertex[0];
   y -= fVertex[1];
   z -= fVertex[2];
   
   return TMath::Sqrt(x*x + y*y + z*z);
}

   
