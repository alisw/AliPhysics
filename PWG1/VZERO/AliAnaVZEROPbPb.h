#ifndef AliAnaVZEROPbPb_cxx
#define AliAnaVZEROPbPb_cxx

class TH1F;
class TH2F;
class AliESDEvent;

#include "AliAnalysisTaskSE.h"

class AliAnaVZEROPbPb : public AliAnalysisTaskSE {

public:
  AliAnaVZEROPbPb();
  AliAnaVZEROPbPb(const char *name);
  virtual ~AliAnaVZEROPbPb() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  virtual void Init();

  Bool_t IsGoodEvent();

  void SetClassesNames(const Char_t * names);

  void CreateQAHistos();
  void CreateHistosPerL2Trigger();

  void FillQAHistos();
  void FillPerL2TriggerHistos();
  
  TH1F* CreateHisto1D(const char* name, const char* title, Int_t nBins, Double_t xMin, Double_t xMax,
		      const char* xLabel = NULL, const char* yLabel = NULL);
  TH2F* CreateHisto2D(const char* name, const char* title, Int_t nBinsX, Double_t xMin, Double_t xMax,
		      Int_t nBinsY, Double_t yMin, Double_t yMax,
		      const char* xLabel = NULL, const char* yLabel = NULL);

 private:
  AliESDEvent *fESD;    //! ESD object
  AliESDVZERO* fEsdV0;
  TList       *fOutputList; //! Output list
  Int_t fNClasses;
  TObjArray *fClassesNames;

  TH2F *fNFlags; //!

  TH1F *fhAdcNoTime[2];     // ADC spectra (no time measurement)
  TH1F *fhAdcWithTime[2];   // ADC spectra (with time measurement)

  TH2F *fhAdcPMTNoTime;   // ADC spectra per PMT (no time measurement)
  TH2F *fhAdcPMTWithTime; // ADC spectra per PMT (with time measurement)
 
  TH1F *fhTime[2];          // Time spectra per side

  TH1F *fhWidth[2];         // Signal width per side

  TH2F *fhTimePMT;        // Time spectra per PMT
  TH2F *fhWidthPMT;       // Signal width per PMT

  TH2F *fhAdcWidth[2];      // ADC vs Signal width per side

  TH2F *fhTimeCorr;       // Corrected mean time V0C vs V0A

  TH2F *fhAdcTime[2];       // ADC vs Time per side

  TH2F *fhPmtMult;             // Number of fired PMTs in V0C vs V0A
  TH1F *fhV0ampl;          // ADC spectra for both rings

  TH2F *fhEvents;         // Event statistics histogram

  TH2F *fhVtxXYBB;        // XY vertex for beam-beam events
  TH1F *fhVtxZBB;         // Z vertex for beam-beam events
  TH2F *fhVtxXYBGA;       // XY vertex for beam-gas (A side) events
  TH1F *fhVtxZBGA;        // Z vertex for beam-gas (A side) events
  TH2F *fhVtxXYBGC;       // XY vertex for beam-gas (C side) events
  TH1F *fhVtxZBGC;        // Z vertex for beam-gas (C side) events

  TH1F *fhL2Triggers;
  TH2F **fhOnlineCharge;
  TH2F **fhRecoMult;
  TH2F **fhV0vsSPDCentrality;
  TH1F **fhTriggerBits;
  TH1F **fhTotRecoMult;
  TH1F **fhCentrality;

  AliAnaVZEROPbPb(const AliAnaVZEROPbPb&); // not implemented
  AliAnaVZEROPbPb& operator=(const AliAnaVZEROPbPb&); // not implemented
  
  ClassDef(AliAnaVZEROPbPb, 1); // example of analysis
};

#endif
