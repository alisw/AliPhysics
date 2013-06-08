#ifndef AliToyMCEventGenerator_H
#define AliToyMCEventGenerator_H


#include "AliToyMCEvent.h"
#include "AliToyMCTrack.h"
#include <AliTPCSpaceCharge3D.h>
#include <AliTPCParam.h>
class AliToyMCEventGenerator : public TObject {
 public:
  AliToyMCEventGenerator();
  AliToyMCEventGenerator(const AliToyMCEventGenerator &gen);
  virtual ~AliToyMCEventGenerator();

  virtual AliToyMCEvent* Generate(Double_t time) = 0;

  Bool_t DistortTrack(AliToyMCTrack &trackIn, Double_t t0);

 protected:
  AliTPCParam *fTPCParam;
  
 private:
  AliTPCSpaceCharge3D *fSpaceCharge;
  
  ClassDef(AliToyMCEventGenerator, 1)
     
};

#endif

