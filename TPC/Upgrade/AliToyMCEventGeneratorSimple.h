#ifndef AliToyMCEventGeneratorSimple_H
#define AliToyMCEventGeneratorSimple_H

#include <AliESDEvent.h>

#include "AliToyMCEventGenerator.h"
#include <TString.h>
class AliToyMCEvent;
class AliESDtrackCuts;
class TTree;
class TFile;
class TH1F;
class AliToyMCEventGeneratorSimple : public AliToyMCEventGenerator {
 public:
  AliToyMCEventGeneratorSimple();
  AliToyMCEventGeneratorSimple(const AliToyMCEventGeneratorSimple &gen);
  virtual ~AliToyMCEventGeneratorSimple();
  AliToyMCEventGeneratorSimple & operator = (const AliToyMCEventGeneratorSimple &gen);

  AliToyMCEvent* Generate(Double_t time);
  AliToyMCEvent* GenerateESD(AliESDEvent& esdEvent, Double_t time);
  AliToyMCEvent* GenerateESD2(Double_t time);
  AliToyMCEvent* GenerateLaser(Double_t time);
  
  void SetParametersToyGen(const Char_t* parfilename="$ALICE_ROOT/TPC/Upgrade/files/params.root", Double_t vertexMean = 0., Double_t vertexSigma = 7.);
  void RunSimulation(Int_t nevents=10, Int_t ntracks=20, Int_t rate=50);
  void RunSimulationBunchTrain(Int_t nevents=10, Int_t ntracks=20);
  void RunSimulationESD(Int_t nevents=10, Int_t ntracks=20);
  void RunSimulationLaser(Int_t nevents=1);
  
  void SetInputESD(const Char_t* filename) {fInputFileNameESD = filename;}
  Int_t OpenInputAndGetMaxEvents(Int_t type, Int_t nevents);
  void RunSimulation2(Bool_t equalspacing, Int_t type, Int_t nevents, Int_t ntracks);
  void GetNGeneratedEventsAndSpacing(Bool_t equalSpacing, Int_t &ngen, Double_t &spacing);
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
  TH1F* fHPt;
  TH1F* fHEta;
  TH1I* fHMult;
  Bool_t fHistosSet;
  TFile* fParamFile;

  ClassDef(AliToyMCEventGeneratorSimple, 1)

};






#endif

