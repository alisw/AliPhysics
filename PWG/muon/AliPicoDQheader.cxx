#include "AliPicoDQheader.h"

ClassImp(AliPicoDQheader)

//_____________________________________________________________________________
AliPicoDQheader::AliPicoDQheader(const TString s) :
TNamed(s,s),
fPSmask(0),
fTriggerInputs(0),
fFiredTriggerClass(""),
fVz(-9999.),
fMultSPDtrkls(-1.),
fMultV0M(-1.),
fMultV0C(-1.),
fNSPDtrkls(-1)
{
//
//  AliPicoDQheader::AliPicoDQheader(const TString s) :
//
}

//_____________________________________________________________________________
AliPicoDQheader::AliPicoDQheader(const AliPicoDQheader &a) :
TNamed(a),
fPSmask(a.fPSmask),
fTriggerInputs(a.fTriggerInputs),
fFiredTriggerClass(a.fFiredTriggerClass),
fVz(a.fVz),
fMultSPDtrkls(a.fMultSPDtrkls),
fMultV0M(a.fMultV0M),
fMultV0C(a.fMultV0C),
fNSPDtrkls(a.fNSPDtrkls)
{
//
//  AliPicoDQheader::AliPicoDQheader(const AliPicoDQheader &a) :
//
}

//_____________________________________________________________________________
AliPicoDQheader &AliPicoDQheader::operator=(const AliPicoDQheader &a)
{
//
//  AliPicoDQheader &AliPicoDQheader::operator=(const AliPicoDQheader &a)
//

  if(&a==this) return *this;
//=============================================================================

  TNamed::operator=(a);
//=============================================================================

  fPSmask            = a.fPSmask;
  fTriggerInputs     = a.fTriggerInputs;
  fFiredTriggerClass = a.fFiredTriggerClass;

  fVz   = a.fVz;
  fMultSPDtrkls = a.fMultSPDtrkls;
  fMultV0M = a.fMultV0M;
  fMultV0C = a.fMultV0C;
  fNSPDtrkls = a.fNSPDtrkls;
//=============================================================================

  return *this;
}

//_____________________________________________________________________________
AliPicoDQheader::~AliPicoDQheader()
{
//
//  AliPicoDQheader::~AliPicoDQheader()
//
}

//_____________________________________________________________________________
void AliPicoDQheader::SetEventInfo(
  UInt_t   wMask,
  UInt_t   wTrgIn,
  TString  sTrgCls,
  Double_t dVz,
  Float_t  dMultSPDtrkls,
  Float_t  dMultV0M,
  Float_t  dMultV0C,
  Int_t    nTrkls)
{
//
//  void AliPicoDQheader::SetEventInfo( UInt_t, UInt_t, TString, Double_t, Float_t, Int_t)
//

  fPSmask = wMask;
  fTriggerInputs = wTrgIn;
  fFiredTriggerClass = sTrgCls;

  fVz   = dVz;
  fMultSPDtrkls = dMultSPDtrkls;
  fMultV0M = dMultV0M;
  fMultV0C = dMultV0C;
  fNSPDtrkls = nTrkls;
//=============================================================================

  return;
}

//_____________________________________________________________________________
void AliPicoDQheader::Reset()
{
//
//  void AliPicoDQheader:: Reset()
//

  fPSmask = 0;
  fTriggerInputs = 0;
  fFiredTriggerClass = "";

  fVz = -9999.;
  fMultSPDtrkls = -1.;
  fMultV0M = -1.;
  fMultV0C = -1.;
  fNSPDtrkls = -1;
//=============================================================================

  return;
}
