//
// Class AliRsnAnalysisTaskEff
//
// Base class for efficiency computation tasks
// which should be inherited by different efficiency computators
//
// author: Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include "AliVEvent.h"

#include "AliRsnEvent.h"
#include "AliRsnCutSet.h"

#include "AliRsnAnalysisTaskEff.h"

ClassImp(AliRsnAnalysisTaskEff)

//_____________________________________________________________________________
AliRsnAnalysisTaskEff::AliRsnAnalysisTaskEff(const char *name) :
   AliRsnVAnalysisTask(name),
   fDefs(0),
   fStepsMC(0),
   fStepsRec(0),
   fAxes("AliRsnValue", 0),
   fOutList(0x0),
   fEventCuts("eventCuts", AliRsnCut::kEvent),
   fVar(0)
{
//
// Default constructor.
//

   DefineOutput(2, TList::Class());
}

//_____________________________________________________________________________
AliRsnAnalysisTaskEff::AliRsnAnalysisTaskEff(const AliRsnAnalysisTaskEff& copy) :
   AliRsnVAnalysisTask(copy),
   fDefs(copy.fDefs),
   fStepsMC(copy.fStepsMC),
   fStepsRec(copy.fStepsRec),
   fAxes(copy.fAxes),
   fOutList(0x0),
   fEventCuts(copy.fEventCuts),
   fVar(0)
{
//
// Copy constrtuctor
//
}

//_____________________________________________________________________________
AliRsnAnalysisTaskEff& AliRsnAnalysisTaskEff::operator=(const AliRsnAnalysisTaskEff& copy)
{
//
// Assignment operator
//

   AliRsnVAnalysisTask::operator=(copy);
   
   fDefs = copy.fDefs;
   fStepsMC = copy.fStepsMC;
   fStepsRec = copy.fStepsRec;
   fAxes = copy.fAxes;
   fEventCuts = copy.fEventCuts;
   
   return (*this);
}

//_____________________________________________________________________________
void AliRsnAnalysisTaskEff::AddDef(TObject* def)
{
//
//  Adds pair definition
//
   fDefs.AddLast(def);
}

//_____________________________________________________________________________
void AliRsnAnalysisTaskEff::AddAxis(AliRsnValue *axis)
{
//
// Add a new axis
//

   Int_t n = fAxes.GetEntries();
   new (fAxes[n]) AliRsnValue(*axis);
}

//_____________________________________________________________________________
void AliRsnAnalysisTaskEff::AddStepMC(TObject *cuts)
{
//
// Add a step on montecarlo
//

   fStepsMC.AddLast(cuts);
}

//_____________________________________________________________________________
void AliRsnAnalysisTaskEff::AddStepRec(TObject *cuts)
{
//
// Add a step on ESD
//

   fStepsRec.AddLast(cuts);
}

//_____________________________________________________________________________
void AliRsnAnalysisTaskEff::RsnUserCreateOutputObjects()
{
//
// Creation of output objects.
// These are created through the utility methods in the analysis manager,
// which produces a list of histograms for each specified set of pairs.
// Each of these lists is added to the main list of this task.
//

   // get number of steps and axes
   Int_t iaxis  = 0;
   Int_t nAxes  = fAxes.GetEntries();
   Int_t nSteps = (Int_t)fStepsMC.GetEntries() + (Int_t)fStepsRec.GetEntries();

   if (!nSteps) {
      AliError("No steps defined");
      return;
   }
   if (!nAxes) {
      AliError("No axes defined");
      return;
   }

   // initialize variable list
   fVar.Set(nAxes);

   // retrieve number of bins for each axis
   Int_t *nBins = new Int_t[nAxes];
   for (iaxis = 0; iaxis < nAxes; iaxis++) {
      AliRsnValue *fcnAxis = (AliRsnValue*)fAxes.At(iaxis);
      nBins[iaxis] = fcnAxis->GetArray().GetSize() - 1;
   }

   // initialize output list
   OpenFile(2);
   fOutList = new TList();
   fOutList->SetOwner();

   // create the containers
   Int_t i, nDef = (Int_t)fDefs.GetEntries();
   for (i = 0; i < nDef; i++) {
      // instantiate a new container
      AliCFContainer *cont = new AliCFContainer(fDefs[i]->GetName(), "", nSteps, nAxes, nBins);
      // set the bin limits for each axis
      for (iaxis = 0; iaxis < nAxes; iaxis++) {
         AliRsnValue *fcnAxis = (AliRsnValue*)fAxes.At(iaxis);
         cont->SetBinLimits(iaxis, fcnAxis->GetArray().GetArray());
      }
      // add the container to output list
      fOutList->Add(cont);
   }

   PostData(2, fOutList);

   // clear heap
   delete [] nBins;
}

//_____________________________________________________________________________
void AliRsnAnalysisTaskEff::RsnUserExec(Option_t*)
{
//
// Execution of the analysis task.
// Recovers the input event and processes it with all included pair objects.
// In this case, we NEED to have reconstruction and MC, otherwise we cannot do anything.
//

   if (fRsnEvent[0].IsESD()) ProcessEventESD();
   if (fRsnEvent[0].IsAOD()) ProcessEventAOD();

   PostData(2, fOutList);
}

//_____________________________________________________________________________
void AliRsnAnalysisTaskEff::RsnTerminate(Option_t*)
{
//
// Termination.
// Could be added some monitor histograms here.
//
}

//______________________________________________________________________________
Bool_t AliRsnAnalysisTaskEff::RsnEventProcess()
{
//
// Customized event pre-processing.
// First checks if the current event passes all cuts,
// and if it does, updates the informations and then
// call the operations which are already defined in the
// omonyme function in mother class
//

   // initially, an event is expected to be bad
   fTaskInfo.SetEventUsed(kFALSE);

   // check the event cuts and update the info data accordingly
   // events not passing the cuts must be rejected
   if (!fEventCuts.IsSelected(&fRsnEvent[0])) {
      fTaskInfo.SetEventUsed(kFALSE);
      return kFALSE;
   }

   // if we reach this point, cuts were passed;
   // then additional operations can be done

   // find leading particle (without any PID/momentum restriction)
   fRsnEvent[0].SelectLeadingParticle(0);

   // final return value is positive
   // but call the mother class method which updates info object
   fTaskInfo.SetEventUsed(kTRUE);
   return AliRsnVAnalysisTask::RsnEventProcess();
}

//_____________________________________________________________________________
TArrayI AliRsnAnalysisTaskEff::FindTracks(Int_t label, AliVEvent *event)
{
//
// Loops an event and find all tracks which have a label
// equal to that passed as first argument.
//

   Int_t   i = 0, nfound = 0;
   Int_t   ntracks = event->GetNumberOfTracks();
   TArrayI array(100);

   for (i = 0; i < ntracks; i++) {
      AliVParticle *track = event->GetTrack(i);
      if (TMath::Abs(track->GetLabel()) != label) continue;
      array[nfound++] = i;
   }

   array.Set(nfound);
   return array;
}

//_____________________________________________________________________________
Int_t AliRsnAnalysisTaskEff::NGoodSteps()
{
//
// Tells the maximum step where cuts for the checked objects are passed
//

   AliInfo("Must be overloaded in inheriting tasks");
   
   return 0;
}

//_____________________________________________________________________________
void AliRsnAnalysisTaskEff::ProcessEventESD()
{
//
// Process current event with the definitions of the specified step in MC list
// and store results in the container slot defined in second argument.
// It is associated with the AliCFContainer with the name of the pair.
//

   AliInfo("Must be overloaded in inheriting tasks");
}

//_____________________________________________________________________________
void AliRsnAnalysisTaskEff::ProcessEventAOD()
{
//
// Process current event with the definitions of the specified step in MC list
// and store results in the container slot defined in second argument.
// It is associated with the AliCFContainer with the name of the pair.
//

   AliInfo("Must be overloaded in inheriting tasks");
}

//_____________________________________________________________________________
void AliRsnAnalysisTaskEff::FillContainer(Bool_t, TObject*)
{
//
// Fill the containers
//

   AliInfo("Must be overloaded in inheriting tasks");
}
