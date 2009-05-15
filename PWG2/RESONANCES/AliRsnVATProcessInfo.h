//
// *** Class AliRsnVATProcessInfo ***
//
// TODO
//
// authors: A. Pulvirenti (email: alberto.pulvirenti@ct.infn.it)
//          M. Vala (email: martin.vala@cern.ch)
//

#ifndef AliRsnVATProcessInfo_H
#define AliRsnVATProcessInfo_H

#include <TNamed.h>
#include <TH1.h>

class AliRsnVATProcessInfo : public TNamed
{
public:
    AliRsnVATProcessInfo(const char*name="AT_RSNInfo");
    AliRsnVATProcessInfo(const AliRsnVATProcessInfo& copy);
    AliRsnVATProcessInfo& operator= (const AliRsnVATProcessInfo& /*copy*/) {
        return *this;
    }
    ~AliRsnVATProcessInfo();

    TList* GenerateInfoList();
    virtual void FillInfo();
    virtual void PrintInfo(const Long64_t &num);

    void SetNumberOfTracks(const Int_t & num) { fNumOfTracks = num; }
    Int_t GetNumberOfTracks() { return fNumOfTracks; };

    const char* GetEventHistogramName() { return Form("hEventsUsed_%s",GetName()); };

    void SetPrintInfoNumber(const Long64_t &num=1) { fPrintInfoNumber = num; }

private:

    TH1I         *fHistUsedEvents;
    Int_t         fNumOfTracks;

    Long64_t      fPrintInfoNumber;

    ClassDef(AliRsnVATProcessInfo, 1)
};

#endif
