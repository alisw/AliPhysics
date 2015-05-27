#ifndef ALIEVENTPLANEESDVARMANAGER_H
#define ALIEVENTPLANEESDVARMANAGER_H



#include "AliEventPlaneManager.h"
#include <AliVEvent.h>
#include <AliESDEvent.h>


class AliEventPlaneESDVarManager : public TNamed{
public:


  AliEventPlaneESDVarManager();
  ~AliEventPlaneESDVarManager();


  
static void FillTPC(AliEventPlaneManager* EPmanager,AliVEvent* event, Float_t* values);
static void FillVZERO(AliEventPlaneManager* EPmanager,AliVEvent* event);
static void FillTZERO(AliEventPlaneManager* EPmanager,AliVEvent* event);
static void FillZDC(AliEventPlaneManager* EPmanager,AliVEvent* event);
static void FillEventInfo(AliVEvent* ev, Float_t* values); 
static void FillTrackInfo(AliESDtrack* p, Float_t* values);

ClassDef(AliEventPlaneESDVarManager, 1);

};





#endif
