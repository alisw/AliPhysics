#ifndef ALIHIGHPTDEDXDATA_H
#define ALIHIGHPTDEDXDATA_H

#include "AliHighPtDeDxBase.h"

class AliHighPtDeDxData : public AliHighPtDeDxBase {
 public:
  AliHighPtDeDxData(); // default constructor  
  AliHighPtDeDxData(const char* name, const char* title); // named constructor  
  virtual ~AliHighPtDeDxData(); // default destructor

  //  virtual void Init(Int_t nPtBins, Double_t* ptBins);
  //  virtual void FillTrackInfo(Float_t weight=1);
  
  TH2D* GetHistDeltaPiVsPt(Int_t pid, Int_t charge);
  TH2D* GetHistDeltaPiVsPtMc(Int_t pid, Int_t charge);

  void SetPionDeDxFunction(TF1* piFunc)     { fDeDxPi    = piFunc; }
  void SetKaonDeDxFunction(TF1* kFunc)      { fDeDxK     = kFunc; }
  void SetProtonDeDxFunction(TF1* pFunc)    { fDeDxP     = pFunc; }
  void SetElectronDeDxFunction(TF1* eFunc)  { fDeDxE     = eFunc; }
  void SetSigmaDeDxFunction(TF1* sigmaFunc) { fSigmaDeDx = sigmaFunc; }

  void Init(Int_t nPtBins, Double_t* ptBins);
  void FillTrackInfo(Float_t weight);

 private:

  TF1*   fDeDxPi;          //! dE/dx vs p for pions
  TF1*   fDeDxK;           //! dE/dx vs p for kaons
  TF1*   fDeDxP;           //! dE/dx vs p for protons
  TF1*   fDeDxE;           //! dE/dx vs p for electrons
  TF1*   fSigmaDeDx;       //! sigma dE/dx vs ncl

  // histograms
  TH2D* hDeltaPiVsPt;      // Delta pi vs pt (both q)
  TH2D* hDeltaPiVsPtNeg;   // Delta pi vs pt (q < 0)
  TH2D* hDeltaPiVsPtPos;   // Delta pi vs pt (q > 0)

  TH2D* hDeltaPiVsPtPiGen;      // Delta pi vs pt (both q) - generated pions
  TH2D* hDeltaPiVsPtPiGenNeg;   // Delta pi vs pt (q < 0)  - generated pions
  TH2D* hDeltaPiVsPtPiGenPos;   // Delta pi vs pt (q > 0)  - generated pions

  TH2D* hDeltaPiVsPtKGen;       // Delta pi vs pt (both q) - generated kaons
  TH2D* hDeltaPiVsPtKGenNeg;    // Delta pi vs pt (q < 0)  - generated kaons
  TH2D* hDeltaPiVsPtKGenPos;    // Delta pi vs pt (q > 0)  - generated kaons

  TH2D* hDeltaPiVsPtPGen;       // Delta pi vs pt (both q) - generated protons
  TH2D* hDeltaPiVsPtPGenNeg;    // Delta pi vs pt (q < 0)  - generated protons
  TH2D* hDeltaPiVsPtPGenPos;    // Delta pi vs pt (q > 0)  - generated protons

  TH2D* hDeltaPiVsPtEGen;       // Delta pi vs pt (both q) - generated electrons
  TH2D* hDeltaPiVsPtEGenNeg;    // Delta pi vs pt (q < 0)  - generated electrons
  TH2D* hDeltaPiVsPtEGenPos;    // Delta pi vs pt (q > 0)  - generated electrons

  TH2D* hDeltaPiVsPtPiMc;     // Delta pi vs pt for MC pions
  TH2D* hDeltaPiVsPtPiMcNeg;  // Delta pi vs pt for MC pions
  TH2D* hDeltaPiVsPtPiMcPos;  // Delta pi vs pt for MC pions
  TH2D* hDeltaPiVsPtKMc;      // Delta pi vs pt for MC Kaons
  TH2D* hDeltaPiVsPtKMcNeg;   // Delta pi vs pt for MC Kaons
  TH2D* hDeltaPiVsPtKMcPos;   // Delta pi vs pt for MC Kaons
  TH2D* hDeltaPiVsPtPMc;      // Delta pi vs pt for MC protons
  TH2D* hDeltaPiVsPtPMcNeg;   // Delta pi vs pt for MC protons
  TH2D* hDeltaPiVsPtPMcPos;   // Delta pi vs pt for MC protons
  TH1D* hPtPi;             // pt distribution for pions (from fits)
  TH1D* hPtK;              // pt distribution for kaons (from fits)
  TH1D* hPtP;              // pt distribution for protons (from fits)
  TH3D* hPrimaryVsPidVsPt; // pt disrtibutions for efficiency
  /* TH1D* hPtPiMc; */
  /* TH1D* hPtKMc; */
  /* TH1D* hPtPMc; */
  
  //  void Print(Option_t* option) const;
  
  ClassDef(AliHighPtDeDxData, 4)  // AliHighPtDeDxData information
    };

#endif
	
