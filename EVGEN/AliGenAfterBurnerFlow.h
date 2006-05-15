#ifndef AliGenAfterBurnerFlow_h
#define AliGenAfterBurnerFlow_h

////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// AliGenAfterBurnerFlow is a After Burner event generator applying flow.
// The generator changes Phi coordinate of the particle momentum.
// Flow (directed and elliptic) can be defined on particle type level
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

  void SetDirectedSimple(Int_t pdg, Float_t v1);
  void SetDirectedParam(Int_t pdg, Float_t v11, Float_t v12 = 0, Float_t v13 = 1, Float_t v14 = 0);

  void SetEllipticSimple(Int_t pdg, Float_t v2);
  void SetEllipticParamPion(Int_t pdg, Float_t v21, Float_t pTmax, Float_t v22=0.);
  void SetEllipticParamOld(Int_t pdg, Float_t v21, Float_t v22, Float_t v23);

  void Init();
  void Generate(); 

 private:

  static const Int_t fgkN = 30; // Size of array fParams 

  Float_t GetCoefficient(Int_t pdg, Int_t n, Float_t Pt, Float_t Y);
  void SetFlowParameters(Int_t pdg, Int_t order, Int_t type, Float_t v1, Float_t v2, Float_t v3, Float_t v4);

  Float_t fReactionPlane;      // Reaction plane angle (in rad)
  Float_t fParams[fgkN][7];    // parameters (0: pdg, 1: order, 2: type,  3-6: actual parameters) 
  Int_t   fCounter;            // counter


  ClassDef(AliGenAfterBurnerFlow,3)

};

////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
