/**********************************************************************
 Created on : 21/11/2013
 Purpose    : pp 8 TeV paraset for JPsi
 Author     : Indranil Das, IPN Orsay
 Email      : indranil.das@cern.ch | indra.ehep@gmail.com
**********************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TRandom.h"
#include "AliGenParam.h"
#include "AliGenMUONlib.h"
#endif 

//           generator functions
//


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


Double_t  V2Zero_M( const Double_t* /*px*/, const Double_t */*dummy*/ )
{

  return 0.0;

}
 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


Double_t YJpsiPPdummy_M(Double_t x, Double_t energy)
{
// J/Psi y
// pp
// from the fit of RHIC + LHC data, see arXiv:1103.2394
//
    x = x/TMath::Log(energy/3.097);
    x = x*x;
    Double_t y = TMath::Exp(-x/0.4/0.4/2);
    if(x > 1) y=0;
    return y;
}

//---------------------------------------------------------------------

Double_t YJpsiPPpoly_M(Double_t x, Double_t energy)
{
// J/Psi y
// pp
// from the fit of RHIC + LHC data, see arXiv:1103.2394
//
    x = x/TMath::Log(energy/3.097);
    x = x*x;
    Double_t y = 1 - 6.9*x*x;
    if(y < 0) y=0;
    return y;
}

//---------------------------------------------------------------------

Double_t YJpsiPP8000_M(const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi y
// pp 7 TeV
//
  return YJpsiPPdummy_M(*px, 8000);
}

//---------------------------------------------------------------------

Double_t YJpsi_M(const Double_t *py, const Double_t */*dummy*/)
{
// J/psi y
  const Double_t ky0 = 4.;
  const Double_t kb=1.;
  Double_t yj;
  Double_t y=TMath::Abs(*py);
  //
  if (y < ky0)
    yj=kb;
  else
    yj=kb*TMath::Exp(-(y-ky0)*(y-ky0)/2);
  return yj;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Int_t IpJpsi_M(TRandom *)
{
// J/Psi composition
    return 443;
}

//---------------------------------------------------------------------

Int_t  IpMuon_M(TRandom* /*ran*/) {

  //muon composition
  // if (ran->Rndm() < 0.5 ) {
  //   return 13;
  // }
  return 443;
  
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Double_t PtJpsiPPdummy_M(Double_t x, Double_t energy)
{
// J/Psi pT
// pp
// from the fit of RHIC, CDF & LHC data, see arXiv:1103.2394
//
  const Double_t kpt0 = 1.04*TMath::Power(energy,0.101);
  const Double_t kxn  = 3.9;
  //
  Double_t pass1 = 1.+0.363*(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}

//---------------------------------------------------------------------

Double_t PtJpsiPP8000_M(const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi pT
// pp 7 TeV
//
  return PtJpsiPPdummy_M(*px,8000);
}

//---------------------------------------------------------------------

Double_t PtJpsi_M( const Double_t *px, const Double_t */*dummy*/)
{
// J/Psi pT
  const Double_t kpt0 = 4.;
  const Double_t kxn  = 3.6;
  Double_t x=*px;
  //
  Double_t pass1 = 1.+(x/kpt0)*(x/kpt0);
  return x/TMath::Power(pass1,kxn);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//           generators
//
AliGenerator* JPsi7TeV()
{
  printf("\nProcessing config setup : JPsi7TeV\n\n");
  AliGenParam *jpsi7TeV = new AliGenParam(1, AliGenMUONlib::kJpsi,"pp 8");
  jpsi7TeV->SetMomentumRange(0,999);
  //jpsi7TeV->SetPtRange(0.,999.);
  //jpsi7TeV->SetYRange(-4.2, -2.3);
  jpsi7TeV->SetPtRange(0,50.);
  jpsi7TeV->SetYRange(-4.5,-2.);
  jpsi7TeV->SetPhiRange(0., 360.);
  jpsi7TeV->SetCutOnChild(1);
  jpsi7TeV->SetChildPhiRange(0.,360.);
  jpsi7TeV->SetChildThetaRange(0.,180.);
  jpsi7TeV->SetForceDecay(kDiMuon);
  jpsi7TeV->SetTrackingFlag(1);
  
  return jpsi7TeV;
}

AliGenerator* GenJPsi8TeVParaSet()
{
  printf("\nProcessing config setup : JPsi8TeVParaSet1\n\n");
  AliGenParam *genJpsi8TeV = new AliGenParam(1,-1, PtJpsiPP8000_M, YJpsiPP8000_M, V2Zero_M, IpJpsi_M);
  genJpsi8TeV->SetMomentumRange(0,999);
  genJpsi8TeV->SetPtRange(0,50.);
  genJpsi8TeV->SetYRange(-4.5,-2.0);
  genJpsi8TeV->SetPhiRange(0., 360.);
  genJpsi8TeV->SetCutOnChild(1);
  genJpsi8TeV->SetChildPhiRange(0.,360.);
  genJpsi8TeV->SetChildThetaRange(0.0,180.0);
  genJpsi8TeV->SetOrigin(0,0,0);
  genJpsi8TeV->SetForceDecay(kDiMuon);
  genJpsi8TeV->SetTrackingFlag(1);
  //genJpsi8TeV->Init();

  return genJpsi8TeV;
}
