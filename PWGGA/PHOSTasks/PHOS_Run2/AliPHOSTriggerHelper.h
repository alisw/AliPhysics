#ifndef AliPHOSTriggerHelper_cxx
#define AliPHOSTriggerHelper_cxx

//Author: Daiki Sekihata (Hiroshima University)

#include "TObject.h"
#include "TVector3.h"
#include "AliPHOSGeometry.h"
#include "AliVEvent.h"
#include "AliVCaloTrigger.h"
#include "AliPHOSClusterCuts.h"

class AliPHOSGeometry;
class AliVEvent;

class AliPHOSTriggerHelper : public TObject {

  public:
    AliPHOSTriggerHelper();
    AliPHOSTriggerHelper(TString trigger, Bool_t isMC);
    AliPHOSTriggerHelper(Int_t L1triggerinput, Int_t L0triggerinput, Bool_t isMC);
    virtual ~AliPHOSTriggerHelper();

    void SetPHOSTRUBadMap(Int_t mod, TH2I *h)
    {
      fIsUserTRUBadMap = kTRUE;
      if(fPHOSTRUBadMap[mod]){
        delete fPHOSTRUBadMap[mod];
        fPHOSTRUBadMap[mod] = 0x0;
      }
      fPHOSTRUBadMap[mod] = new TH2I(*h);
    }

    TH2I* GetPHOSTRUBadMap(Int_t mod) {return fPHOSTRUBadMap[mod];}
 
    Bool_t IsOnActiveTRUChannel(AliCaloPhoton *ph);
    Bool_t IsGoodTRUChannel(const char * det, Int_t mod,Int_t ix, Int_t iz);
    Int_t WhichTRU(Int_t cellx, Int_t cellz);
    Int_t WhichTRUChannel(Int_t cellx, Int_t cellz, Int_t &chX, Int_t &chZ);
    Int_t FindHighestAmplitudeCellAbsId(AliVCluster *clu, AliVCaloCells *cells);
    Bool_t IsMatched(Int_t *trgrelid, Int_t *clurelid);
    Bool_t IsMatchedDeltaR(Int_t *trgrelid, TVector3 position);

    Int_t GetL1TriggerInput(){return fTriggerInputL1;}
    Int_t GetL0TriggerInput(){return fTriggerInputL0;}

    void SetMatchingDistance(Int_t xmin, Int_t zmin, Int_t xmax, Int_t zmax){
      fUseDeltaRMatching = kFALSE;
      fXmin = xmin;
      fXmax = xmax;
      fZmin = zmin;
      fZmax = zmax;
    }

    void SetMatchingDeltaR(Double_t DeltaR){
      fUseDeltaRMatching = kTRUE;
      fMatchingDeltaR = DeltaR;
    }

    Double_t GetMatchingDeltaR(){return fMatchingDeltaR;}
    Bool_t IsDeltaRUsed() {return fUseDeltaRMatching;}
    Bool_t IsPHI7(AliVEvent *event, AliPHOSClusterCuts *cuts);
    Double_t GetDistanceToClosestTRUChannel(AliCaloPhoton *ph);

  private:
    AliPHOSGeometry *fPHOSGeo;
    TH2I* fPHOSTRUBadMap[6];
    Int_t fXmin;
    Int_t fXmax;
    Int_t fZmin;
    Int_t fZmax;
    Double_t fMatchingDeltaR;
    AliVEvent *fEvent;
    AliESDEvent* fESDEvent;
    AliAODEvent* fAODEvent;
    Int_t fTriggerInputL1;//1PHL:5,1PHM:6,1PHH:7
    Int_t fTriggerInputL0;//0PH0 should be 9
    Bool_t fIsMC;
    AliVCaloTrigger* fCaloTrigger;
    Bool_t fIsUserTRUBadMap;
    Int_t fRunNumber;
    Bool_t fUseDeltaRMatching;

  private:
    AliPHOSTriggerHelper(const AliPHOSTriggerHelper&);
    AliPHOSTriggerHelper& operator=(const AliPHOSTriggerHelper&);

    ClassDef(AliPHOSTriggerHelper, 14);

};

#endif

