#ifndef ALIITSSIMULATIONSPDBARI_H
#define ALIITSSIMULATIONSPDBARI_H

#include "AliITSsimulation.h"

class AliITSMapA2;
class AliITSsegmentation;
class AliITSresponse;
class AliITSmodule;

//-------------------------------------------------------------------

class AliITSsimulationSPDbari : public AliITSsimulation {

public:
        
  AliITSsimulationSPDbari();
  AliITSsimulationSPDbari(AliITSsegmentation *seg, AliITSresponse *res);
  ~AliITSsimulationSPDbari();
  AliITSsimulationSPDbari(const AliITSsimulationSPDbari &source); // copy constructor
  AliITSsimulationSPDbari& operator=(const AliITSsimulationSPDbari &source); // ass. operator

  void DigitiseModule(AliITSmodule *mod,Int_t module, Int_t dummy);
  void SetFluctuations(Float_t **pList);
  void HitToDigit(AliITSmodule *mod, Int_t hitpos, Int_t module, 
              Int_t *frowpixel, Int_t *fcolpixel, Double_t *fenepixel,
	      Float_t **pList);
	     
  void UpdateMap( Int_t row, Int_t col, Double_t ene); 
  void ChargeSharing(Float_t x1l,Float_t z1l,Float_t x2l,
                    Float_t z2l,Int_t c1,Int_t r1,Int_t c2,
				    Int_t r2,Float_t etot,
				    Int_t &npixel,Int_t *frowpixel,
				    Int_t *fcolpixel,Double_t *fenepixel);
  
  void SetCoupling(Int_t row, Int_t col, Int_t ntrack, Float_t **pList);
  void CreateDigit(Int_t nhits, Int_t module, Float_t **pList);
  void GetList(Int_t track, Float_t **pList, Int_t row, Int_t col);

  void CreateHistograms();
  void ResetHistograms();
  TObjArray*  GetHistArray() {
    // get hist array
    return fHis;
  }

private:

  AliITSMapA2  *fMapA2;        // MapA2
  Float_t      fThresh;        // Threshold
  Float_t      fSigma;         // Noise 
  Float_t      fCouplCol;      // Coupling along columns
  Float_t      fCouplRow;      // Coupling along rows
  Int_t        fNPixelsX;      // NPixelsX
  Int_t        fNPixelsZ;      // NPixelsZ

  TObjArray *fHis;             // just in case for histogramming
    
  ClassDef(AliITSsimulationSPDbari,1)  // Simulation of SPD clusters

};

#endif 



