#ifndef ALIREDUCEDQNFILLEVENT_H
#define ALIREDUCEDQNFILLEVENT_H

#include <TNamed.h>

class AliQnCorrectionsManager;
//class AliQnCorrectionsHistos;
class AliReducedEventInfo;
class AliReducedBaseEvent;
class AliHistogramManager;
class AliReducedEventInfo;

class AliReducedQnFillEvent : public TNamed{
public:

  AliReducedQnFillEvent();
  ~AliReducedQnFillEvent();

void Process(AliReducedEventInfo* event, Float_t* values);
void FillDetectors(Float_t* values);
void FillTPC(Float_t* values);
void FillVZERO();
void FillTZERO();
void FillZDC();
void FillFMD();
//void FillEventInfo(Float_t* values); 
//void FillTrackInfo(AliReducedTrackInfo* p, Float_t* values);

//static void SetEvent() {fEvent=InputEvent();}
void SetEvent(AliReducedEventInfo* ev) {fEvent=ev;}
void SetQnCorrectionsManager(AliQnCorrectionsManager* EPmanager) {fEventPlaneManager=EPmanager; SetDetectors();}
void SetHistogramManager(AliHistogramManager* EPmanager) {fEventPlaneHistos=EPmanager;}
void SetDetectors();

private:

  AliReducedQnFillEvent(const AliReducedQnFillEvent &c);
  AliReducedQnFillEvent& operator= (const AliReducedQnFillEvent &c);

  AliReducedEventInfo* fEvent;
  AliQnCorrectionsManager* fEventPlaneManager;
  AliHistogramManager* fEventPlaneHistos;

  Bool_t fFillVZERO;
  Bool_t fFillTPC;
  Bool_t fFillZDC;
  Bool_t fFillTZERO;
  Bool_t fFillFMD;

ClassDef(AliReducedQnFillEvent, 1);

};





#endif
