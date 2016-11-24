#ifndef ALIPICOHEADERCJ_H
#define ALIPICOHEADERCJ_H

#include <TNamed.h>
#include <TString.h>

class TString;

class AliVEventHandler;

class AliPicoHeaderCJ : public TNamed {

 public :

  enum {
    kPP = 0,  // pp    collisions
    kPA = 1,  // p-Pb  collisions
    kAP = 2,  // Pb-p  collisions
    kAA = 3   // Pb-Pb collisions
  };

  enum {
    kKshort     = BIT(0),  //
    kLambda     = BIT(1),  //
    kAntiLambda = BIT(2),  //
    kXiPos      = BIT(3),  //
    kXiNeg      = BIT(4),  //
    kDzero      = BIT(5),  //
    kDstar      = BIT(6)   //
  };

  enum {
    kEventAccCheck   = BIT(0),  //
    kEventAccMult    = BIT(1),  //
    kEventAccTrigger = BIT(2),  //
    kEventAccVertex  = BIT(3),  //
    kEventAccPileup  = BIT(4)   //
  };

  enum {
    kPrimary                = BIT(0), //
    kPhysicalPrimary        = BIT(1), //
    kSecondaryFromWeakDecay = BIT(2), //
    kSecondaryFromMaterial  = BIT(3)  //
  };
//=============================================================================

  AliPicoHeaderCJ();
  AliPicoHeaderCJ(const AliPicoHeaderCJ &src);
  AliPicoHeaderCJ& operator=(const AliPicoHeaderCJ &src);
  virtual ~AliPicoHeaderCJ();
//=============================================================================

  void SetEventInfo(AliVEventHandler* const pEH);

  void BackgroundRhoRD02(Double_t d) { fBackgroundRhoRD02 = d; }
  void BackgroundRhoRD03(Double_t d) { fBackgroundRhoRD03 = d; }
  void BackgroundRhoRD04(Double_t d) { fBackgroundRhoRD04 = d; }
  void BackgroundRhoMC02(Double_t d) { fBackgroundRhoMC02 = d; }
  void BackgroundRhoMC03(Double_t d) { fBackgroundRhoMC03 = d; }
  void BackgroundRhoMC04(Double_t d) { fBackgroundRhoMC04 = d; }

  void Reset();
//=============================================================================

  UInt_t  PhysSelMask()       const { return fPhysSelMask; }
  TString FiredTriggerClass() const { return fFiredTriggerClass; }

  void Vertex(Double_t d[3]) { for (Int_t i=3; i--;) d[i] = fVtx[i]; }

  Float_t CentralityV0M() const { return fCentralityV0M; }
  Float_t CentralityV0A() const { return fCentralityV0A; }
  Float_t CentralityCL1() const { return fCentralityCL1; }
  Float_t CentralityZNA() const { return fCentralityZNA; }
  Double32_t EventPlane() const { return fEventPlane;    }
  Double_t  BackgroundRho(const TString sJet) const;
//=============================================================================

 private :

  UInt_t  fPhysSelMask;  //
  TString fFiredTriggerClass;  //

  Double_t fVtx[3];  //!

  Float_t fCentralityV0M;  //
  Float_t fCentralityV0A;  //
  Float_t fCentralityCL1;  //
  Float_t fCentralityZNA;  //

  Double32_t fEventPlane;  //

  Double_t fBackgroundRhoRD02;  //
  Double_t fBackgroundRhoRD03;  //
  Double_t fBackgroundRhoRD04;  //

  Double_t fBackgroundRhoMC02;  //
  Double_t fBackgroundRhoMC03;  //
  Double_t fBackgroundRhoMC04;  //

  ClassDef(AliPicoHeaderCJ, 4);
};


#endif
