#ifndef AliGenAfterBurnerFlow_h
#define AliGenAfterBurnerFlow_h

////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// AliGenAfterBurnerFlow is a After Burner event generator applying flow.
// The generator changes Phi coordinate of the particle momentum.
// Flow (directed and elliptic) can be defined on particle type level
//
////////////////////////////////////////////////////////////////////////////////////////////////////


#include "AliGenerator.h"
#include "AliStack.h"

class TObjArray;

class AliGenAfterBurnerFlow : public AliGenerator {

 public:

  AliGenAfterBurnerFlow();
  AliGenAfterBurnerFlow(Float_t reactionPlane);
  ~AliGenAfterBurnerFlow();

  void      AarbitraryReactionPlaneAngle() {fHow = 3;}
  void      SetDirectedSimple(Int_t pdg, Float_t v1);
  void      SetDirectedParam(Int_t pdg, Float_t v11, Float_t v12 = 0, Float_t v13 = 1, Float_t v14 = 0);
  void      SetEllipticSimple(Int_t pdg, Float_t v2);
  void      SetEllipticParamPion(Int_t pdg, Float_t v21, Float_t pTmax, Float_t v22=0.);
  void      SetEllipticParamOld(Int_t pdg, Float_t v21, Float_t v22, Float_t v23);
  void      SetEllipticParam(Int_t pdg, Float_t v00, Float_t v10, Float_t v11, Float_t v22);
  void      SetNpParams(Int_t order=-1, Float_t p0=-1, Float_t p1=-1, Float_t p2=-1, Float_t p3=-1);
  void      SetNpDefault() { SetNpParams(2,1,-5e-8,-5e-6); } 
  Bool_t    IsPrimary(Int_t pdg) const;
  void      Init();
  void      Generate(); 
  void      NeglectFlow(Int_t pdg) {fIsPrim[pdg]=kFALSE;}

 private:
  AliGenAfterBurnerFlow(const AliGenAfterBurnerFlow&);
  AliGenAfterBurnerFlow& operator=(const AliGenAfterBurnerFlow&);

  static const Int_t fgkN = 300;   // Size of array fParams 
  static const Int_t fgkPDG= 5226; // Size of PDG code array 

  Bool_t    fIsPrim[fgkPDG];   // array of primaries
  Float_t   fReactionPlane;    // Reaction plane angle (in rad)
  Int_t     fHow;              // Choose reaction plane angle
  Float_t   fParams[fgkN][7];  // parameters (0: pdg, 1: order, 2: type,  3-6: actual parameters) 
  Float_t   fNpParams[5];      // parameters (0: order, 1-4: actual parameters) 
  Int_t     fCounter;          // counter
  AliStack *fStack;            //!

  Float_t   GetCoefficient(Int_t pdg, Int_t n, Float_t Pt, Float_t Y) const;
  Float_t   GetNpNorm(Int_t npart) const;
  void      InitPrimaries();
  void      Rotate(Int_t i, Double_t phi, Bool_t IsPrim=kTRUE);
  void      SetFlowParameters(Int_t pdg, Int_t order, Int_t type, Float_t v1, Float_t v2, Float_t v3, Float_t v4);

  ClassDef(AliGenAfterBurnerFlow,4)
};
#endif
