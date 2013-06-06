#ifndef ToyMCEventGeneratorSimple_H
#define ToyMCEventGeneratorSimple_H


#include "ToyMCEvent.h"
#include "ToyMCTrack.h"
#include "ToyMCEventGenerator.h"
class ToyMCEventGeneratorSimple : public ToyMCEventGenerator {
 public:
  ToyMCEventGeneratorSimple();
  ToyMCEventGeneratorSimple(const ToyMCEventGeneratorSimple &gen);
  virtual ~ToyMCEventGeneratorSimple();
  ToyMCEventGeneratorSimple & operator = (const ToyMCEventGeneratorSimple &gen);

  ToyMCEvent* Generate(Double_t time);
  void SetParameters(Double_t vertexMean, Double_t vertexSigma);
 private:
  
  Double_t fVertexMean;
  Double_t fVertexSigma;

  ClassDef(ToyMCEventGeneratorSimple, 1)

};






#endif

