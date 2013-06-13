#ifndef AliToyMCEventGeneratorSimple_H
#define AliToyMCEventGeneratorSimple_H


#include "AliToyMCEventGenerator.h"

class AliToyMCEvent;

class AliToyMCEventGeneratorSimple : public AliToyMCEventGenerator {
 public:
  AliToyMCEventGeneratorSimple();
  AliToyMCEventGeneratorSimple(const AliToyMCEventGeneratorSimple &gen);
  virtual ~AliToyMCEventGeneratorSimple();
  AliToyMCEventGeneratorSimple & operator = (const AliToyMCEventGeneratorSimple &gen);

  AliToyMCEvent* Generate(Double_t time);
  void SetParameters(Double_t vertexMean, Double_t vertexSigma);

  void RunSimulation(const Int_t nevents=10, const Int_t ntracks=20);
 private:
  
  Double_t fVertexMean;
  Double_t fVertexSigma;

  Int_t fNtracks;

  ClassDef(AliToyMCEventGeneratorSimple, 1)

};






#endif

