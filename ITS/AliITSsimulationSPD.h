#ifndef ALIITSSIMULATIONSPD_H
#define ALIITSSIMULATIONSPD_H

#include "AliITSsimulation.h"

class AliITSMapA2;
class AliITSsegmentation;
class AliITSresponse;
class AliITSmodule;

//-------------------------------------------------------------------

class AliITSsimulationSPD : public AliITSsimulation {

public:
        
  AliITSsimulationSPD();
  AliITSsimulationSPD(AliITSsegmentation *seg, AliITSresponse *res);
  ~AliITSsimulationSPD();
  AliITSsimulationSPD(const AliITSsimulationSPD &source); // copy constructor
  AliITSsimulationSPD& operator=(const AliITSsimulationSPD &source); // ass. operator

  void DigitiseModule(AliITSmodule *mod,Int_t module,Int_t dummy);
  void ChargeToSignal(Float_t **pList);
  void GetList(Int_t track, Int_t hit, Float_t **pList, Int_t *IndexRange);

  void CreateHistograms();
  void ResetHistograms();
  TObjArray*  GetHistArray() {
    // get hist array
    return fHis;
  }

private:

  AliITSMapA2  *fMapA2;        // MapA2
  Float_t      fNoise;         // Noise
  Float_t      fBaseline;      // Baseline
  Int_t        fNPixelsX;      // NPixelsX
  Int_t        fNPixelsZ;      // NPixelsZ

  TObjArray *fHis;             // just in case for histogramming
    
  ClassDef(AliITSsimulationSPD,1)  // Simulation of SPD clusters

};

#endif 



