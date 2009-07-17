//
// *** Class AliRsnVATProcessInfo ***
//
// TODO
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
    AliRsnVATProcessInfo(const char*name = "AT_RSNInfo");
    AliRsnVATProcessInfo(const AliRsnVATProcessInfo& copy);
    AliRsnVATProcessInfo& operator= (const AliRsnVATProcessInfo& /*copy*/) {return *this;}
    ~AliRsnVATProcessInfo();

    void         GenerateInfoList(TList* list);
    virtual void FillInfo();
    virtual void PrintInfo(const Long64_t &num);

    const char* GetEventHistogramName() { return Form("hEventsUsed_%s",GetName()); };
    Long64_t    GetNumerOfEventsProcessed();
    void        SetEventUsed(Bool_t isUsed = kTRUE) { fEventUsed = isUsed; }
    Bool_t      IsEventUsed() const { return fEventUsed; };

    void        SetPrintInfoNumber(const Long64_t &num=1) { fPrintInfoNumber = num; }

  private:

    TH1I         *fHistUsedEvents;
    Int_t         fEventUsed;

    Long64_t      fPrintInfoNumber;

    ClassDef(AliRsnVATProcessInfo, 1)
};

#endif
