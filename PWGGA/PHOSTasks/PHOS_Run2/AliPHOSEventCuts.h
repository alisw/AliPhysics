#ifndef AliPHOSEventCuts_cxx
#define AliPHOSEventCuts_cxx

// Author: Daiki Sekihata (Hiroshima University)
//this analsyis class provides event selection.
//e.g., reject pileup event, Incomplete event, |zvtx|>10cm.

class TH2;
class AliPHOSGeometry;
class AliPHOSTriggerHelper;

#include "AliPHOSTriggerHelper.h"
#include "AliAnalysisCuts.h"

class AliPHOSEventCuts : public AliAnalysisCuts {

  public:
    AliPHOSEventCuts(const char *name = "AliPHOSEventCuts");
    virtual ~AliPHOSEventCuts(); 

    virtual Bool_t IsSelected(TObject* obj) {return AcceptEvent((AliVEvent*)obj);}
    virtual Bool_t IsSelected(TList* /*list*/) {return kTRUE;}
    Bool_t AcceptEvent(AliVEvent *event);
    Bool_t IsPHOSTriggerAnalysis() {return fIsPHOSTriggerAnalysis;}

    void SetMCFlag(Bool_t mc) {fIsMC = mc;}
    void SetMaxAbsZvtx(Double_t maxZ) {fMaxAbsZvtx = maxZ;}
    void SetRejectPileup(Bool_t reject) {fRejectPileup = reject;}
    void SetRejectDAQIncompleteEvent(Bool_t reject) {fRejectDAQIncomplete = reject;}

    void DoPHOSTriggerAnalysis(Bool_t flag, TObject *obj=0x0){
      fIsPHOSTriggerAnalysis = flag;
      if(flag) fTriggerHelper = (AliPHOSTriggerHelper*)obj;
    }

    Bool_t IsMC() {return fIsMC;}

    AliPHOSTriggerHelper *GetPHOSTriggerHelper() {return fTriggerHelper;}

  private:
    Bool_t fIsMC;
    Bool_t fUsePHOSTender;
    Double_t fMaxAbsZvtx;
    Bool_t fRejectPileup;
    Bool_t fRejectDAQIncomplete;
    Bool_t fIsPHOSTriggerAnalysis;
    AliPHOSTriggerHelper *fTriggerHelper;
    AliPHOSGeometry *fPHOSGeo;
    TH2I* fPHOSTRUBadMap[6];

  private:
    AliPHOSEventCuts(const AliPHOSEventCuts&);
    AliPHOSEventCuts& operator=(const AliPHOSEventCuts&);

    ClassDef(AliPHOSEventCuts, 5); // example of analysis
};

#endif

