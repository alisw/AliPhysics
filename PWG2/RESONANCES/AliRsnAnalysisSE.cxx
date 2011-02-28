//
// Class AliRsnAnalysisSE
//
// Virtual Class derivated from AliRsnVAnalysisTaskSE which will be base class
// for all RSN SE tasks
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include <Riostream.h>
#include <TList.h>
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliAODEvent.h"

#include "AliRsnCutSet.h"
#include "AliRsnVATProcessInfo.h"
#include "AliRsnAnalysisSE.h"

ClassImp(AliRsnAnalysisSE)

//_____________________________________________________________________________
AliRsnAnalysisSE::AliRsnAnalysisSE(const char *name, Bool_t useKine) :
   AliRsnVAnalysisTaskSE(name, useKine),
   fRsnAnalysisManager(),
   fEventCuts("eventCuts", AliRsnCut::kEvent),
   fOutList(0x0),
   fZeroEventPercentWarning(100),
   fUseZeroEventWarning(kTRUE)
{
//
// Default constructor.
// Defines another output slot for histograms/ntuples
//

   DefineOutput(2, TList::Class());
}

//_____________________________________________________________________________
AliRsnAnalysisSE::AliRsnAnalysisSE(const AliRsnAnalysisSE& copy) :
   AliRsnVAnalysisTaskSE(copy),
   fRsnAnalysisManager(copy.fRsnAnalysisManager),
   fEventCuts(copy.fEventCuts),
   fOutList(0x0),
   fZeroEventPercentWarning(copy.fZeroEventPercentWarning),
   fUseZeroEventWarning(copy.fUseZeroEventWarning)
{
//
// Copy constructor.
//
}

//_____________________________________________________________________________
AliRsnAnalysisSE& AliRsnAnalysisSE::operator=(const AliRsnAnalysisSE& copy)
{
//
// Assigment operator.
//

   AliRsnVAnalysisTaskSE::operator=(copy);

   fRsnAnalysisManager = copy.fRsnAnalysisManager;
   fEventCuts = copy.fEventCuts;
   if (fOutList) fOutList->Clear();
   fZeroEventPercentWarning = copy.fZeroEventPercentWarning;
   fUseZeroEventWarning = copy.fUseZeroEventWarning;

   return (*this);
}

//_____________________________________________________________________________
void AliRsnAnalysisSE::RsnUserCreateOutputObjects()
{
//
// Creation of output objects.
// These are created through the utility methods in the analysis manager,
// which asks all the AliRsnPair objects to initialize their output which
// is then linked to the TList data member of this, which will contain all the output.
//

   if (!fOutList) fOutList = new TList;
   fOutList->Clear();

   fRsnAnalysisManager.InitAllPairs(fOutList);

   PostData(2, fOutList);
}

//_____________________________________________________________________________
void AliRsnAnalysisSE::RsnUserExec(Option_t*)
{
//
// Execution of the analysis task.
// Recovers the input event and processes it with all included pair objects,
// using 'reconstructed' or 'MonteCarlo' functions depending on MC-only flag.
//

   fRsnAnalysisManager.ProcessAll(fMCOnly);

   PostData(2, fOutList);
}

//_____________________________________________________________________________
void AliRsnAnalysisSE::RsnTerminate(Option_t*)
{
//
// Termination.
// Could be added some monitor histograms here.
//
}

//______________________________________________________________________________
Bool_t AliRsnAnalysisSE::EventProcess()
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

   // check #1: number of tracks in event (reject empty events)
   Int_t    ntracks = fRsnEvent.GetMultiplicity();
   Double_t zeroEventPercent = 0.0;
   if (ntracks < 1) {
      // if using the checker for amount of empty events, update it
      if (fUseZeroEventWarning) {
         TH1I *hist = (TH1I*)fInfoList->FindObject(fTaskInfo.GetEventHistogramName());
         if (hist) {
            if (hist->Integral() > 1) zeroEventPercent = (Double_t)hist->GetBinContent(1) / hist->Integral() * 100;
            if ((zeroEventPercent > fZeroEventPercentWarning) && (fEntry > 100))
               AliWarning(Form("%3.2f%% Events are with zero tracks (CurrentEvent=%d)!!!", zeroEventPercent, fEntry));
         }
      }

      // empty events are rejected by default
      fTaskInfo.SetEventUsed(kFALSE);
      AliDebug(AliLog::kDebug, "Empty event. Skipping...");
      return kFALSE;
   }

   // check the event cuts and update the info data accordingly
   // events not passing the cuts must be rejected
   if (!fEventCuts.IsSelected(&fRsnEvent)) {
      fTaskInfo.SetEventUsed(kFALSE);
      return kFALSE;
   }

   // if we reach this point, cuts were passed;
   // then additional operations can be done

   // find leading particle (without any PID/momentum restriction)
   fRsnEvent.SelectLeadingParticle(0);

   // final return value is positive
   // but call the mother class method which updates info object
   fTaskInfo.SetEventUsed(kTRUE);
   return AliRsnVAnalysisTaskSE::EventProcess();
}
