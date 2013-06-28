#ifndef AliToyMCEventGeneratorSimple_H
#define AliToyMCEventGeneratorSimple_H

#include <AliESDEvent.h>

#include "AliToyMCEventGenerator.h"
#include <TString.h>
class AliToyMCEvent;
class AliESDtrackCuts;
class TTree;
class TFile;

class AliToyMCEventGeneratorSimple : public AliToyMCEventGenerator {
 public:
  AliToyMCEventGeneratorSimple();
  AliToyMCEventGeneratorSimple(const AliToyMCEventGeneratorSimple &gen);
  virtual ~AliToyMCEventGeneratorSimple();
  AliToyMCEventGeneratorSimple & operator = (const AliToyMCEventGeneratorSimple &gen);

  AliToyMCEvent* Generate(Double_t time);
  AliToyMCEvent* GenerateESD(AliESDEvent& esdEvent, Double_t time);
  AliToyMCEvent* GenerateESD2(Double_t time);
  void SetParametersSimple(Double_t vertexMean, Double_t vertexSigma);
  
  void RunSimulation(const Int_t nevents=10, const Int_t ntracks=20);
  void RunSimulationBunchTrain(const Int_t nevents=10, const Int_t ntracks=20);
  void RunSimulationESD(const Int_t nevents=10, const Int_t ntracks=20);
  void SetInputESD(const Char_t* filename) {fInputFileNameESD = filename;}
  Int_t OpenInputAndGetMaxEvents(const Int_t type, const Int_t nevents);
  void RunSimulation2(const Bool_t equalspacing, const Int_t type, const Int_t nevents, const Int_t ntracks);
  void GetNGeneratedEventsAndSpacing(const Bool_t equalSpacing, Int_t &ngen, Double_t &spacing);
  Bool_t CloseInputFile();

 private:
  
  Double_t fVertexMean;
  Double_t fVertexSigma;

  Int_t fNtracks;
  TString fInputFileNameESD;

  AliESDtrackCuts *fESDCuts;
  Int_t fInputIndex;
  AliESDEvent* fESDEvent;
  TTree* fESDTree;
  TFile* fInputFile;



  ClassDef(AliToyMCEventGeneratorSimple, 1)

};






#endif

