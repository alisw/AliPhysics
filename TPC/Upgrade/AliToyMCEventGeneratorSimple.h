#ifndef AliToyMCEventGeneratorSimple_H
#define AliToyMCEventGeneratorSimple_H


#include "AliToyMCEvent.h"
#include "AliToyMCTrack.h"
#include "AliToyMCEventGenerator.h"
class AliToyMCEventGeneratorSimple : public AliToyMCEventGenerator {
 public:
  AliToyMCEventGeneratorSimple();
  AliToyMCEventGeneratorSimple(const AliToyMCEventGeneratorSimple &gen);
  virtual ~AliToyMCEventGeneratorSimple();
  AliToyMCEventGeneratorSimple & operator = (const AliToyMCEventGeneratorSimple &gen);

  AliToyMCEvent* Generate(Double_t time);
  void SetParameters(Double_t vertexMean, Double_t vertexSigma);
 private:
  
  Double_t fVertexMean;
  Double_t fVertexSigma;

  ClassDef(AliToyMCEventGeneratorSimple, 1)

};






#endif

