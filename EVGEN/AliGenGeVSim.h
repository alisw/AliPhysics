#ifndef ALIGENGEVSIM_H
#define ALIGENGEVSIM_H

#include "TObjArray.h"
#include "TF1.h"
#include "TF2.h"
#include "AliGenerator.h"
#include "AliGeVSimParticle.h"


class AliGenGeVSim : public AliGenerator {

  Int_t   fModel;            // Selected model (1-5)
  Float_t fPsi;              // Reaction Plane angle (0-2pi)

  TF1 *fPtFormula;           // Pt formula for model (1)
  TF1 *fYFormula;            // Y formula for model (1)
  TF2 *fPtYFormula[4];       // Pt,Y formulae for model (2)-(4) and custom
  TF1 *fPhiFormula;          // phi formula 

  TObjArray *fPartTypes;     // Registered particles

  void InitFormula();        
  //void PlotDistributions();
  
  Bool_t CheckPtYPhi(Float_t pt, Float_t y, Float_t phi);
  Bool_t CheckP(Float_t p[3]);

  Float_t FindScaler(Int_t paramId, Int_t pdg);

 public:

  AliGenGeVSim();
  AliGenGeVSim(Int_t model, Float_t psi);

  virtual ~AliGenGeVSim();

  /////////////////////////////////////////////////////////////////

  void AddParticleType(AliGeVSimParticle *part);
  
  void Init();
  void Generate();

  /////////////////////////////////////////////////////////////////

  ClassDef(AliGenGeVSim, 1)

};

#endif
