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

    enum PileupFinder {
      kNone = -1,
      kSPD = 0,
      kSPDInMultBins = 1,
      kMultiVertexer = 2
    };

    virtual Bool_t IsSelected(TObject* obj) {return AcceptEvent((AliVEvent*)obj);}
    virtual Bool_t IsSelected(TList* /*list*/) {return kTRUE;}
    Bool_t AcceptEvent(AliVEvent *event);

    void SetMCFlag(Bool_t mc) {fIsMC = mc;}
    void SetMaxAbsZvtx(Double_t maxZ) {fMaxAbsZvtx = maxZ;}
    void SetRejectPileup(Bool_t reject) {fRejectPileup = reject;}
    void SetRejectDAQIncompleteEvent(Bool_t reject) {fRejectDAQIncomplete = reject;}
    void SetPileupFinder(AliPHOSEventCuts::PileupFinder pf) {fPF = pf;}

    Bool_t IsMC() {return fIsMC;}

  private:
    Bool_t fIsMC;
    Double_t fMaxAbsZvtx;
    Bool_t fRejectPileup;
    Bool_t fRejectDAQIncomplete;
    PileupFinder fPF;

  private:
    AliPHOSEventCuts(const AliPHOSEventCuts&);
    AliPHOSEventCuts& operator=(const AliPHOSEventCuts&);

    ClassDef(AliPHOSEventCuts, 9);
};

#endif

