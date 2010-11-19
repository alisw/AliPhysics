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

#ifndef ALIRSNVATPROCESSINFO_H
#define ALIRSNVATPROCESSINFO_H

#include <TNamed.h>

class TH1I;

class AliRsnVATProcessInfo : public TNamed
{
  public:
    AliRsnVATProcessInfo(const char *name = "RSNInfo");
    AliRsnVATProcessInfo(const AliRsnVATProcessInfo& copy);
    AliRsnVATProcessInfo& operator= (const AliRsnVATProcessInfo& copy);
    ~AliRsnVATProcessInfo();

    void          GenerateInfoList(TList* list);
    virtual void  FillInfo();
    virtual void  PrintInfo(const Long64_t &num);

    const char*   GetEventHistogramName() { return Form("hEventsUsed_%s",GetName()); };
    Long64_t      GetNumerOfEventsProcessed();
    void          SetEventUsed(Int_t flag) { fEventUsed = flag; }
    Int_t         IsEventUsed() const { return fEventUsed; };

    void          SetPrintInfoNumber(const Long64_t &num=1) { fPrintInfoNumber = num; }

  private:

    TH1I         *fHistUsedEvents;      // hist of used events
    Int_t         fEventUsed;           // number of used events

    Long64_t      fPrintInfoNumber;     // print info number 

    ClassDef(AliRsnVATProcessInfo, 1)
};

#endif
