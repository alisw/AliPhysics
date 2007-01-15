#ifndef ALIITSSIMULATIONFASTPOINTS_H
#define ALIITSSIMULATIONFASTPOINTS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//////////////////////////////////////////////////////////
// implements fast simulation 
/////////////////////////////////////////////////////////
//

#include <TMath.h>

#include "AliITSsimulation.h"
class TClonesArray;
class AliITSmodule;
class TRandom;

class AliITSsimulationFastPoints : public AliITSsimulation
{

public:
  AliITSsimulationFastPoints(); // default constructor  
  virtual ~AliITSsimulationFastPoints() {;} 
  virtual AliITSsimulation& operator=(const AliITSsimulation&){return *this;}
  void CreateFastRecPoints(AliITSmodule *mod,Int_t module,TRandom *rndm, 
			   TClonesArray* recp);
  void CreateFastRecPoints(Int_t module,TClonesArray* recp);

  virtual void SetSegmentationModel(Int_t dt, AliITSsegmentation *seg){fDetType->SetSegmentationModel(dt,seg);}
  virtual AliITSsegmentation* GetSegmentationModel(Int_t dt){return fDetType->GetSegmentationModel(dt);}

 private:

    
    virtual void SetSigmaRPhi(Double_t sigmarphi[6]);  
    virtual void SetSigmaZ(Double_t sigmaz[6]);  
    virtual void SetSigmaDe(Double_t sigmade[6]);  
    virtual void SetThrDe(Double_t thrde[6]); 
    Double_t SigmaRPhi(Int_t layer) const {return fSigmaRPhi[layer-1];}  
    Double_t SigmaZ(Int_t layer) const  {return fSigmaZ[layer-1];}  
    Double_t SigmaDe(Int_t layer) const {return fSigmaDe[layer-1];} 
    Double_t ThrDe(Int_t layer) const {return fThrDe[layer-1];} 


    Double_t fSigmaRPhi[6];              // Sigmas in rphi for the 6 layers
    Double_t fSigmaZ[6];                 // Sigmas in Z for the 6 layers
    Double_t fSigmaDe[6];                // Sigmas in energy loss for the 6 layers
    Double_t fThrDe[6];                  // Energy thresholds for the 6 layers


  ClassDef(AliITSsimulationFastPoints,1) // Fast point simulator.

};

#endif
