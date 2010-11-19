//
// *** Class AliRsnVATProcessInfo ***
//
// Virtual class which makes computations at the event level,
// in order to return a list of histograms useful to have a look
// of the characteristics of used events.
// If can be inherited and customized for the needs of the analysis.
//
// authors: A. Pulvirenti (email: alberto.pulvirenti@ct.infn.it)
//          M. Vala (email: martin.vala@cern.ch)
//

#include <TList.h>
#include <TH1.h>

#include "AliLog.h"

#include "AliRsnVATProcessInfo.h"

ClassImp(AliRsnVATProcessInfo)

//______________________________________________________________________________
AliRsnVATProcessInfo::AliRsnVATProcessInfo(const char *name) :
  TNamed(name,name),
  fHistUsedEvents(0x0),
  fEventUsed(kFALSE),
  fPrintInfoNumber(1000)
{
//
// Constructor.
// Does nothing more than initialization of data members.
//

  AliDebug(AliLog::kDebug+2, "Entering");
  AliDebug(AliLog::kDebug+2, "Exiting");
}

//______________________________________________________________________________
AliRsnVATProcessInfo::AliRsnVATProcessInfo(const AliRsnVATProcessInfo& copy) :
  TNamed(copy),
  fHistUsedEvents(0x0),
  fEventUsed(copy.fEventUsed),
  fPrintInfoNumber(copy.fPrintInfoNumber)
{
//
// Copy constructor.
// Clones the histogram and copies the values of other data members
//

  AliDebug(AliLog::kDebug+2, "Entering");
  
  fHistUsedEvents  = (TH1I*)copy.fHistUsedEvents->Clone();
  
  AliDebug(AliLog::kDebug+2, "Exiting");
}

//______________________________________________________________________________
AliRsnVATProcessInfo& AliRsnVATProcessInfo::operator= 
(const AliRsnVATProcessInfo& copy)
{
//
// Assignment operator.
// Clones the histogram and copies the values of other data members.
//

  AliDebug(AliLog::kDebug+2, "Entering");

  fHistUsedEvents  = (TH1I*)copy.fHistUsedEvents->Clone();
  fEventUsed       = copy.fEventUsed;
  fPrintInfoNumber = copy.fPrintInfoNumber;
  
  AliDebug(AliLog::kDebug+2, "Exiting");
  
  return (*this);
}

//______________________________________________________________________________
AliRsnVATProcessInfo::~AliRsnVATProcessInfo()
{
//
// Destructor.
// Does nothing, since the histogram it creates is usually owned
// by another object (TList output of AnalysisTask's), but sets
// the data member pointers to NULL.
//

  AliDebug(AliLog::kDebug+2, "Entering");
  
  fHistUsedEvents  = 0x0;
  fEventUsed       = 0;
  fPrintInfoNumber = 0;
  
  AliDebug(AliLog::kDebug+2, "Exiting");
}

//______________________________________________________________________________
void AliRsnVATProcessInfo::GenerateInfoList(TList *list)
{
//
// Allocate in memory the histograms created in this class and store them
// inside the TList object passed as argument, which usually belongs to 
// an AnalysisTask object, and represents one of its outputs.
// If the histogram was already initialized, it is deleted and recreated.
//

  AliDebug(AliLog::kDebug+2, "Entering");

  // delete existing allocation of stored objects
  if (fHistUsedEvents) delete fHistUsedEvents;
  
  // create stored objects
  fHistUsedEvents = new TH1I(GetEventHistogramName(), "Skipped/used events", 2, 0, 2);

  // ad objects to list
  list->Add(fHistUsedEvents);

  AliDebug(AliLog::kDebug+2, "Exiting");
}

//______________________________________________________________________________
void AliRsnVATProcessInfo::FillInfo()
{
//
// This method defines how the information histograms must be filled.
// The structure of this class is auto-consistent, but in case of inheritance
// this method must be modified accordingly.
// Current implementation uses the 'fEventUsed' flag to choose if the event
// has been used or not, and increments the corresponding bin in the related
// histogram (bin '0' = skipped, bin '1' = used).
//

  if (fEventUsed)
    fHistUsedEvents->Fill(1);
  else
    fHistUsedEvents->Fill(0);
}

//______________________________________________________________________________
void AliRsnVATProcessInfo::PrintInfo(const Long64_t &num)
{
//
// This method is used in some cases
// to inform about number of events processed
//

  if ((num+1) % fPrintInfoNumber == 0) AliInfo(Form("Events processed %lld", num+1));
}

//______________________________________________________________________________
Long64_t AliRsnVATProcessInfo::GetNumerOfEventsProcessed()
{
//
// returns number of events from histogram
//

  if (fHistUsedEvents) return fHistUsedEvents->Integral();
  return 0;
}
