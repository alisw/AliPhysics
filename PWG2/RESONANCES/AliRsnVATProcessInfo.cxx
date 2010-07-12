//
// *** Class AliRsnVATProcessInfo ***
//
//  TODO
//  TODO
//
// authors: A. Pulvirenti (email: alberto.pulvirenti@ct.infn.it)
//          M. Vala (email: martin.vala@cern.ch)
//

#include <TList.h>
#include <TH1.h>

#include "AliLog.h"

#include "AliRsnVATProcessInfo.h"

ClassImp(AliRsnVATProcessInfo)

//_____________________________________________________________________________
AliRsnVATProcessInfo::AliRsnVATProcessInfo(const char *name) :
    TNamed(name,name),
    fHistUsedEvents(0x0),
    fEventUsed(kFALSE),
    fPrintInfoNumber(1000)
{
//
// Constructor.
//

  AliDebug(AliLog::kDebug+2,"<-");
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
AliRsnVATProcessInfo::AliRsnVATProcessInfo(const AliRsnVATProcessInfo& copy) :
    TNamed(copy),
    fHistUsedEvents(copy.fHistUsedEvents),
    fEventUsed(copy.fEventUsed),
    fPrintInfoNumber(copy.fPrintInfoNumber)
{
//
// Copy constructor.
//

  AliDebug(AliLog::kDebug+2,"<-");
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
AliRsnVATProcessInfo::~AliRsnVATProcessInfo()
{
//
// Destructor.
//

  AliDebug(AliLog::kDebug+2,"<-");
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnVATProcessInfo::GenerateInfoList(TList *list)
{
//
// Creates a TList containing all objects created in this class.
// This list will be passe to the AnalysisTask object and will generate
// a separate directory in the output file, with these informations.
// This TList will have the same name of this object and will own its content.
//

  AliDebug(AliLog::kDebug+2,"<-");

  // create list
  AliDebug(AliLog::kWarning,"Doing new TList(), so make sure you delete this list ... ");
//   TList* list = new TList();
//   list->SetName(GetName());
//   list->SetOwner();

  // initialize contained objects
  fHistUsedEvents = new TH1I(GetEventHistogramName(), "Skipped/used events", 2, 0, 2);

  // ad objects to list
  list->Add(fHistUsedEvents);

  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
void AliRsnVATProcessInfo::FillInfo()
{
//
// This method defines how the information histograms must be filled.
// The structure of this class is auto-consistend, but in case of inheritance
// this method must be modified accordingly.
// Current implementation uses the 'fEventUsed' flag to choose if the event
// has been used or not, and increments the corrisponding bin in the related
// histogram (bin '0' = skipped, bin '1' = used).
//

  if (fEventUsed)
    fHistUsedEvents->Fill(1);
  else
    fHistUsedEvents->Fill(0);
}

//_____________________________________________________________________________
void AliRsnVATProcessInfo::PrintInfo(const Long64_t &num)
{
//
// This method is used in some cases
// to inform about number of events processed
//

  if ((num+1) % fPrintInfoNumber == 0) AliInfo(Form("Events processed %l",num+1));
}


Long64_t AliRsnVATProcessInfo::GetNumerOfEventsProcessed()
{
//
// returns number of events from histogram
//
  if (fHistUsedEvents)
    return fHistUsedEvents->Integral();
  return 0;
}

