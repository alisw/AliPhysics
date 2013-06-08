#ifndef AliToyMCEventGenerator_H
#define AliToyMCEventGenerator_H

#include <TString.h>

class TFile;
class TTree;

class AliTPCParam;
class AliTPCSpaceCharge3D;

class AliToyMCTrack;
class AliToyMCEvent;

class AliToyMCEventGenerator : public TObject {
 public:
  AliToyMCEventGenerator();
  AliToyMCEventGenerator(const AliToyMCEventGenerator &gen);
  virtual ~AliToyMCEventGenerator();

  virtual AliToyMCEvent* Generate(Double_t time) = 0;

  Bool_t DistortTrack(AliToyMCTrack &trackIn, Double_t t0);

  void SetOutputFileName(const char* file) { fOutputFileName=file; }
  const char* GetOutputFileName()    const { return fOutputFileName.Data(); }
 protected:
  AliTPCParam *fTPCParam;
  AliToyMCEvent *fEvent;
  
  Bool_t ConnectOutputFile();
  Bool_t CloseOutputFile();
  void FillTree();
  
 private:
  AliToyMCEventGenerator& operator= (const AliToyMCEventGenerator& );
   
  AliTPCSpaceCharge3D *fSpaceCharge;

  TString fOutputFileName;
  TFile   *fOutFile;
  TTree   *fOutTree;
  
  ClassDef(AliToyMCEventGenerator, 1)
     
};

#endif

