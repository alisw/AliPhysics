#ifndef ALIPICOHEADERV0_H
#define ALIPICOHEADERV0_H

#include <TNamed.h>

#include <TMap.h>
#include <TString.h>
#include <TParameter.h>

class AliInputEventHandler;

class AliPicoHeaderV0 : public TNamed {

 public :

  AliPicoHeaderV0(const TString s="");
  AliPicoHeaderV0(const AliPicoHeaderV0 &src);
  AliPicoHeaderV0& operator=(const AliPicoHeaderV0 &src);
  virtual ~AliPicoHeaderV0();
//=============================================================================

  UInt_t  PhysSelMask()       const { return fPS;  }
  TString FiredTriggerClass() const { return fTrg; }

  Float_t MultiplicityPercentile(const TString &s);
  Double32_t EventPlane() const { return fEP; }

  void Vertex(Double_t d[3]) { for (Int_t i=0; i<3; ++i) d[i] = fVtx[i]; }
//=============================================================================

  void AddMultEstimator(const TString &s) {
    if (!fAct) fAct = new TMap();
    fAct->Add(new TObjString(s.Data()),
              new TParameter<Float_t>(s.Data(),-999.));
    return;
  }
//=============================================================================

  void Reset();
  void SetEventInfo(AliVEventHandler* const pH);
//=============================================================================

 private :

  UInt_t  fPS;   //
  TString fTrg;  //

  Double_t fVtx[3];  //
  Double32_t fEP;    //

  TMap  *fAct;  //!
//=============================================================================

  ClassDef(AliPicoHeaderV0, 6);
};

#endif
