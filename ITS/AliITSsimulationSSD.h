#ifndef ALIITSSIMULATIONSSD_H
#define ALIITSSIMULATIONSSD_H

#include <TArrayF.h>

#include "AliITSsimulation.h"


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
  void HitToDigit(Int_t &hit,Int_t nhits,TObjArray *hits);            

  // return the pointer to the signal array of P strip
  TArrayF* GetSignalP() {return fP;}
  // return the pointer to the signal array of N strip
  TArrayF* GetSignalN() {return fN;}

protected:
  void  ApplyNoise();
  void  ApplyCoupling();
  void  ApplyThreshold();
  void  ApplyDAQ();
  
  Float_t  F(Float_t x, Float_t s); 
  Float_t  Get2Strip(Int_t flag, Int_t istrip, Float_t x, Float_t z);
  
  // Data members 
  
protected:
  
  AliITSdcsSSD    *fDCS; //!
  Int_t   fNstrips;      //!
  Float_t fPitch;        //!
  TArrayF *fN;         //! for signal N side
  TArrayF *fP;         //! for signal P side
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
