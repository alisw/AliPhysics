#ifndef AliGenAfterBurnerFlow_h
#define AliGenAfterBurnerFlow_h

////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// AliGenAfterBurnerFlow is a After Burner event generator applying flow.
// The generator changes Phi coordinate of the particle momentum.
// Flow (directed and elliptical) can be defined on particle type level
//
// For examples, parameters and testing macros refer to:
// http:/home.cern.ch/radomski
//
// Author:
// Sylwester Radomski,
// GSI, April 2002
// 
// S.Radomski@gsi.de
//
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "AliGenerator.h"

class TObjArray;

class AliGenAfterBurnerFlow : public AliGenerator {

 public:

  AliGenAfterBurnerFlow();
  AliGenAfterBurnerFlow(Float_t reactionPlane);

  ~AliGenAfterBurnerFlow();

  void SetDirected(Int_t pdg, Float_t v11, Float_t v12 = 0, Float_t v13 = 1, Float_t v14 = 0);
  void SetElliptic(Int_t pdg, Float_t v21, Float_t v22 = 0, Float_t v23 = 0);

  void SetDefDirected(Float_t v11, Float_t v12 = 0, Float_t v13 = 1, Float_t v14 = 0);
  void SetDefElliptic(Float_t v21, Float_t v22 = 0, Float_t v23 = 0);

  void Init();
  void Generate(); 

 private:

  static const Int_t kN = 30;

  Float_t GetCoeff(Int_t pdg, Int_t n, Float_t Pt, Float_t Y);
  void SetFlowParameters(Int_t pdg, Int_t order, Float_t v1, Float_t v2, Float_t v3, Float_t v4);

  Float_t fReactionPlane;      // Reaction plane angle (in rad)
  Float_t fParams[kN][6];      // parameters (0: pdg, 1: order, 2-5: actual parameters) 
  Int_t   fCounter;            // counter

 public:

  ClassDef(AliGenAfterBurnerFlow,1)

};

////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
