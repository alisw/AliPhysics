//
// *** Class AliRsnVATProcessInfo ***
//
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
AliRsnVATProcessInfo::AliRsnVATProcessInfo(const char *name) : TNamed(name,name),
        fHistUsedEvents(0x0),
        fPrintInfoNumber(1000)
{
    AliDebug(AliLog::kDebug+2,"<-");
    AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
AliRsnVATProcessInfo::AliRsnVATProcessInfo(const AliRsnVATProcessInfo& copy) : TNamed(copy),
        fHistUsedEvents(copy.fHistUsedEvents),
        fPrintInfoNumber(copy.fPrintInfoNumber)

{
    AliDebug(AliLog::kDebug+2,"<-");
    AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
AliRsnVATProcessInfo::~AliRsnVATProcessInfo()
{
    AliDebug(AliLog::kDebug+2,"<-");
    AliDebug(AliLog::kDebug+2,"->");
}

TList* AliRsnVATProcessInfo::GenerateInfoList() {
    AliDebug(AliLog::kDebug+2,"<-");

    AliDebug(AliLog::kWarning,"Doing new TList(), so make sure you delete this list ... ");

    TList* list = new TList();
    list->SetName(GetName());
    list->SetOwner();

    fHistUsedEvents = new TH1I(GetEventHistogramName(), "skipped and used events in this analysis", 2, 0, 2);
    list->Add(fHistUsedEvents);

    AliDebug(AliLog::kDebug+2,"->");
    return list;
}

void AliRsnVATProcessInfo::FillInfo() {
    if (!fNumOfTracks)
        fHistUsedEvents->Fill(fNumOfTracks);
    else
        fHistUsedEvents->Fill(1);
}

void AliRsnVATProcessInfo::PrintInfo(const Long64_t &num) {
  if ((num+1)%fPrintInfoNumber == 0){
    AliInfo(Form("Events processed %d",num+1));
  }
    
}

