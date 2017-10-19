#ifndef AliPHOSEventCuts_cxx
#define AliPHOSEventCuts_cxx

//Author: Daiki Sekihata (Hiroshima University)
//this analsyis class provides event selection.
//e.g., reject pileup event, Incomplete event, |zvtx|>10cm.

#include "AliAnalysisCuts.h"

class AliPHOSEventCuts : public AliAnalysisCuts {

  public:
    AliPHOSEventCuts(const char *name = "AliPHOSEventCuts");
    virtual ~AliPHOSEventCuts(); 

    virtual Bool_t IsSelected(TObject* obj) {return AcceptEvent((AliVEvent*)obj);}
    virtual Bool_t IsSelected(TList* /*list*/) {return kTRUE;}
    Bool_t AcceptEvent(AliVEvent *event);

    void SetMCFlag(Bool_t mc) {fIsMC = mc;}
    void SetMaxAbsZvtx(Double_t maxZ) {fMaxAbsZvtx = maxZ;}
    void SetRejectPileup(Bool_t reject) {fRejectPileup = reject;}
    void SetRejectDAQIncompleteEvent(Bool_t reject) {fRejectDAQIncomplete = reject;}

    Bool_t IsMC() {return fIsMC;}

  private:
    Bool_t fIsMC;
    Double_t fMaxAbsZvtx;
    Bool_t fRejectPileup;
    Bool_t fRejectDAQIncomplete;

  private:
    AliPHOSEventCuts(const AliPHOSEventCuts&);
    AliPHOSEventCuts& operator=(const AliPHOSEventCuts&);

    ClassDef(AliPHOSEventCuts, 8);
};

#endif

