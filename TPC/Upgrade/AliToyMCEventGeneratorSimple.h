#ifndef AliToyMCEventGeneratorSimple_H
#define AliToyMCEventGeneratorSimple_H

#include <AliESDEvent.h>

#include "AliToyMCEventGenerator.h"
#include <TString.h>
class AliToyMCEvent;
class AliESDtrackCuts;
class TTree;


class AliToyMCEventGeneratorSimple : public AliToyMCEventGenerator {
 public:
  AliToyMCEventGeneratorSimple();
  AliToyMCEventGeneratorSimple(const AliToyMCEventGeneratorSimple &gen);
  virtual ~AliToyMCEventGeneratorSimple();
  AliToyMCEventGeneratorSimple & operator = (const AliToyMCEventGeneratorSimple &gen);

  AliToyMCEvent* Generate(Double_t time);
  AliToyMCEvent* GenerateESD(AliESDEvent& esdEvent, Double_t time);
  void SetParametersSimple(Double_t vertexMean, Double_t vertexSigma);
  
  void RunSimulation(const Int_t nevents=10, const Int_t ntracks=20);
  void RunSimulationBunchTrain(const Int_t nevents=10, const Int_t ntracks=20);
  void RunSimulationESD(const Int_t nevents=10, const Int_t ntracks=20);
  void SetInputESD(const Char_t* filename) {fInputFileNameESD = filename;}

 private:
  
  Double_t fVertexMean;
  Double_t fVertexSigma;

  Int_t fNtracks;
  TString fInputFileNameESD;

  AliESDtrackCuts *fESDCuts;
  //AliESDEvent* fESDEvent;
  //TTree* fESDTree;
 
  ClassDef(AliToyMCEventGeneratorSimple, 1)

};






#endif

