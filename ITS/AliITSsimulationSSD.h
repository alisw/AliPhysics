#ifndef ALIITSSIMULATIONSSD_H
#define ALIITSSIMULATIONSSD_H

#include <TArrayF.h>

#include "AliITSdcsSSD.h"
//#include "AliITSdictSSD.h"
#include "AliITSsimulation.h"
#include "AliITSsegmentationSSD.h"


class AliITSMapA2;
//class AliITSdictSSD;
class AliITSdcsSSD;

class AliITSsimulationSSD: public AliITSsimulation {

public:

  AliITSsimulationSSD() {}
  AliITSsimulationSSD(AliITSsimulationSSD &source); // copy constructor
  AliITSsimulationSSD& operator=(AliITSsimulationSSD &source); // operator =
  AliITSsimulationSSD(AliITSsegmentation *seg, AliITSresponse *resp);
  virtual ~AliITSsimulationSSD();
    
  void DigitiseModule(AliITSmodule *mod, Int_t imod, Int_t dummy);  
  //void HitToDigit(Double_t x0, Double_t y0, Double_t z0, 
  void HitToDigit(Int_t module, Double_t x0, Double_t y0, Double_t z0, //b.b. 
		  Double_t x, Double_t y, Double_t z, Double_t de,
		  Int_t *indexRange, Bool_t first);            
  Int_t NumOfSteps(Double_t x, Double_t y, Double_t z,
		   Double_t  & dex,Double_t & dey,Double_t & dez);
  void GetList(Int_t track, Float_t **pList, Int_t *IndexRange);
  void ChargeToSignal(Float_t **pList);
  
  //TArrayF* GetSignalP() {return fP;}
  //TArrayF* GetSignalN() {return fN;}
  
  AliITSsegmentationSSD *GetSegmentation() {return (AliITSsegmentationSSD*)fSegmentation;}

 protected:
  
  void  IntegrateGaussian(Int_t k, Double_t par, Double_t av, Double_t sigma, 
			  Double_t inf, Double_t sup,
			  Int_t *indexRange, Bool_t first);
  void  ApplyNoise();
  void  ApplyCoupling();
  //  void  ApplyThreshold();
  //void  ApplyDAQ();
  
  Float_t  F(Float_t av, Float_t x, Float_t s); 
  //  Float_t  Get2Strip(Int_t flag, Int_t istrip, Float_t x, Float_t z);
  
  // Data members 
  
 protected:
  
  AliITSdcsSSD    *fDCS;
  Int_t   fNstrips;
  Float_t fPitch;
  //TArrayF *fN;         // for signal N side
  //TArrayF *fP;         // for signal P side
  //AliITSdictSSD *fTracksP;  //!
  //AliITSdictSSD *fTracksN;  //! 
  
  //______________________________________________________________
  //  
  // Parameters for simulation
  //______________________________________________________________
  Int_t    fSteps;                //Number of steps 
  
 private:

  AliITSMapA2  *fMapA2;        // MapA2
  
  ClassDef(AliITSsimulationSSD,1)
    
    };


#endif
