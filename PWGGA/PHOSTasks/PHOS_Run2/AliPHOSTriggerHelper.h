#ifndef AliPHOSTriggerHelper_cxx
#define AliPHOSTriggerHelper_cxx

//Author: Daiki Sekihata (Hiroshima University)

#include "TObject.h"
#include "AliPHOSGeometry.h"
#include "AliVEvent.h"
#include "AliVCaloTrigger.h"

class AliPHOSGeometry;
class AliVEvent;

class AliPHOSTriggerHelper : public TObject {

  public:
    AliPHOSTriggerHelper();
    AliPHOSTriggerHelper(Int_t inputL1, Int_t inputL0);
    //AliPHOSTriggerHelper(Int_t xmin, Int_t zmin, Int_t xmax, Int_t zmax);
    virtual ~AliPHOSTriggerHelper();

    void SetPHOSTRUBadMap(Int_t mod, TH2I *h)
    {
      if(fPHOSTRUBadMap[mod]) delete fPHOSTRUBadMap[mod];
      fPHOSTRUBadMap[mod] = new TH2I(*h);
      AliInfo(Form("Setting Bad Map Histogram  %s",fPHOSTRUBadMap[mod]->GetName()));
    }

    TH2I* GetPHOSTRUBadMap(Int_t mod) {return fPHOSTRUBadMap[mod];}
 
    Bool_t IsGoodTRUChannel(const char * det, Int_t mod,Int_t ix, Int_t iz);
    Int_t WhichTRU(Int_t cellx, Int_t cellz);
    Int_t WhichTRUChannel(Int_t cellx, Int_t cellz, Int_t &chX, Int_t &chZ);
    Int_t FindHighestAmplitudeCellAbsId(AliVCluster *clu, AliVCaloCells *cells);
    Bool_t IsMatched(Int_t *trgrelid, Int_t *clurelid);

    Int_t GetL1TriggerInput(){return fTriggerInputL1;}
    Int_t GetL0TriggerInput(){return fTriggerInputL0;}

    void SetMatchingDistance(Int_t xmin, Int_t zmin, Int_t xmax, Int_t zmax){
      fXmin = xmin;
      fXmax = xmax;
      fZmin = zmin;
      fZmax = zmax;
    }

    Bool_t IsPHI7(AliVEvent *event);

  private:
    AliPHOSGeometry *fPHOSGeo;
    TH2I* fPHOSTRUBadMap[6];
    Int_t fXmin;
    Int_t fXmax;
    Int_t fZmin;
    Int_t fZmax;
    AliVEvent *fEvent;
    AliESDEvent* fESDEvent;
    AliAODEvent* fAODEvent;
    Int_t fTriggerInputL1;//1PHL:5,1PHM:6,1PHH:7
    Int_t fTriggerInputL0;//0PH0 should be 9
    Bool_t fIsMC;
    AliVCaloTrigger* fCaloTrigger;

  private:
    AliPHOSTriggerHelper(const AliPHOSTriggerHelper&);
    AliPHOSTriggerHelper& operator=(const AliPHOSTriggerHelper&);

    ClassDef(AliPHOSTriggerHelper, 4);

};

#endif

