#ifndef ALIITSSIMULATIONSSD_H
#define ALIITSSIMULATIONSSD_H

#include <TArrayF.h>

#include "AliITSdcsSSD.h"
#include "AliITSsimulation.h"
#include "AliITSdictSSD.h"


class AliITSdictSSD;
class AliITSdcsSSD;

class AliITSsimulationSSD: public AliITSsimulation {

public:

  AliITSsimulationSSD() {}
  AliITSsimulationSSD(AliITSsimulationSSD &source); // copy constructor
  AliITSsimulationSSD& operator=(AliITSsimulationSSD &source); // operator =
  AliITSsimulationSSD(AliITSsegmentation *seg, AliITSresponse *resp);
  virtual ~AliITSsimulationSSD();
    
  void DigitiseModule(AliITSmodule *mod, Int_t module, Int_t dummy);  
  void HitToDigit(Int_t &hit,Int_t idtrack,Int_t nhits,TObjArray *hits);            

  TArrayF* GetSignalP() {
                         // return the signal of P strip
                         return fP;
								}
  TArrayF* GetSignalN() {
                         // return the signal of N strip
                         return fN;
								}

protected:
  
  Int_t IntegrateGaussian(Double_t par, Double_t av, Double_t sigma, 
			  Double_t inf, Double_t sup);
  void  NormalizeCharge(Int_t k, Double_t pair);
  Int_t NumOfSteps(Double_t x, Double_t y, Double_t z,
		   Double_t &dex,Double_t &dey,Double_t &dez);
  void  ApplyNoise();
  void  ApplyCoupling();
  void  ApplyThreshold();
  void  ApplyDAQ();
  
  Float_t  F(Float_t x, Float_t s); 
  Float_t  Get2Strip(Int_t flag, Int_t istrip, Float_t x, Float_t z);
  
  // Data members 
  
protected:
  
  AliITSdcsSSD    *fDCS;
  Int_t   fNstrips;
  Float_t fPitch;
  TArrayF *fN;         // for signal N side
  TArrayF *fP;         // for signal P side
  AliITSdictSSD *fTracksP;  //!
  AliITSdictSSD *fTracksN;  //! 
  //______________________________________________________________
  //  
  // Parameters for simulation
  //______________________________________________________________
  Int_t    fSteps;                //Number of steps 
  
  ClassDef(AliITSsimulationSSD,1)
    
    
    };


#endif
