#ifndef ALIPICOJETHEADER_H
#define ALIPICOJETHEADER_H

#include <TNamed.h>
#include <TString.h>

class TString;

class AliInputEventHandler;

class AliPicoJetHeader : public TNamed {

 public :

  AliPicoJetHeader();
  AliPicoJetHeader(const AliPicoJetHeader &src);
  AliPicoJetHeader& operator=(const AliPicoJetHeader &src);
  virtual ~AliPicoJetHeader();

  void SetEventInfo(AliInputEventHandler* const pEH);

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


  Double_t   BackgroundRhoRD03() const { return fBackgroundRhoRD03;  }

  Double_t   BackgroundRhoRD04() const { return fBackgroundRhoRD04;  }

  Double_t   BackgroundRhoMC02() const { return fBackgroundRhoMC02;  }

  Double_t   BackgroundRhoMC03() const { return fBackgroundRhoMC03;  }

  Double_t   BackgroundRhoMC04() const { return fBackgroundRhoMC04;  }

 private :

  UInt_t  fPhysSelMask;  //
  TString fFiredTriggerClass;  //

  Double_t fVtx[3];  //

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

  ClassDef(AliPicoJetHeader, 1);
};

#endif
