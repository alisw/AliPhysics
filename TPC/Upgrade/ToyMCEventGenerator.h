#ifndef ToyMCEventGenerator_H
#define ToyMCEventGenerator_H


#include "ToyMCEvent.h"
#include "ToyMCTrack.h"
#include <AliTPCSpaceCharge3D.h>
#include <AliTPCParam.h>
class ToyMCEventGenerator : public TObject {
 public:
  ToyMCEventGenerator();
  ToyMCEventGenerator(const ToyMCEventGenerator &gen);
  virtual ~ToyMCEventGenerator();

  virtual ToyMCEvent* Generate(Double_t time) = 0;

  Bool_t DistortTrack(ToyMCTrack &trackIn, Double_t t0);

 protected:
  AliTPCParam *fTPCParam;
  
 private:
  AliTPCSpaceCharge3D *fSpaceCharge;
  
  ClassDef(ToyMCEventGenerator, 1)
     
};

#endif

