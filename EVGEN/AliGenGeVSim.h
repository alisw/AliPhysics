#ifndef ALIGENGEVSIM_H
#define ALIGENGEVSIM_H

////////////////////////////////////////////////////////////////////////////////
//
// AliGenGeVSim is a class implementing simple Monte-Carlo event generator for 
// testing algorythms and detector performance.
//
// In this event generator particles are generated from thermal distributions 
// without any dynamics and addicional constrains. Distribution parameters like
// multiplicity, particle type yields, inverse slope parameters, flow coeficients 
// and expansion velocities are expleicite defined by the user.
//
// GeVSim contains four thermal distributions the same as
// MevSim event generator developed for STAR experiment.
//
// In addition custom distributions can be used be the mean of TF2 function
// named "gevsimPtY". 
//  
// Azimuthal distribution is deconvoluted from (Pt,Y) distribution
// and is described by two Fourier coefficients representing 
// Directed and Elliptical flow. 
// 
// To apply flow to event ganerated by an arbitraly event generator
// refer to AliGenAfterBurnerFlow class.
// For examples, parameters and testing macros refer to:
// http:/home.cern.ch/radomski
//  
// Author:
// Sylwester Radomski,
// GSI, March 2002
//  
// S.Radomski@gsi.de
//
////////////////////////////////////////////////////////////////////////////////

class TFormula;
class TF1;
class TF2;
class TObjArray;
class AliGeVSimParticle;


class AliGenGeVSim : public AliGenerator {

 public:
  
  AliGenGeVSim();
  AliGenGeVSim(Int_t model, Float_t psi);
  
  virtual ~AliGenGeVSim();
  
  /////////////////////////////////////////////////////////////////
  
  void AddParticleType(AliGeVSimParticle *part);
  
  void Init();
  void Generate();
  
  /////////////////////////////////////////////////////////////////
  
 private:
  
  Int_t   fModel;            // Selected model (1-5)
  Float_t fPsi;              // Reaction Plane angle (0-2pi)
  
  TF1 *fPtFormula;           // Pt formula for model (1)
  TF1 *fYFormula;            // Y formula for model (1)
  TF2 *fPtYFormula[4];       // Pt,Y formulae for model (2)-(4) and custom
  TF1 *fPhiFormula;          // phi formula 
  
  TObjArray *fPartTypes;     // Registered particles
  
  void InitFormula();        
  void DetermineReactionPlane();
 
  TFormula* DetermineModel();
 
  //void PlotDistributions();
  
  Bool_t CheckPtYPhi(Float_t pt, Float_t y, Float_t phi);
  Bool_t CheckP(Float_t p[3]);
  
  Float_t FindScaler(Int_t paramId, Int_t pdg);
  
  /////////////////////////////////////////////////////////////////

 public:

  ClassDef(AliGenGeVSim, 1)

};

#endif
