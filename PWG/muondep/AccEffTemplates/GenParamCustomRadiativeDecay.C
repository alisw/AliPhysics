#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include "TRandom.h"
#include "AliGenerator.h"
#include "AliGenParam.h"
#include "AliGenCocktail.h"
#include "AliGenEvtGen.h"
#endif

static Int_t IpCustom( TRandom *ran);
static Double_t PtCustom( const Double_t *px, const Double_t */*dummy*/ );
static Double_t YCustom( const Double_t *py, const Double_t */*dummy*/ );
static Double_t V2Custom( const Double_t *pv, const Double_t */*dummy*/ );

//-------------------------------------------------------------------------
AliGenerator* GenParamCustomRadiativeDecay()
{
  
  AliGenCocktail *gener = new AliGenCocktail();
  
  gener->UsePerEventRates();
  AliGenParam *quarkoniaGener = new AliGenParam(1,-1,PtCustom,YCustom,V2Custom,IpCustom);
  quarkoniaGener->SetMomentumRange(0,1e6);
  quarkoniaGener->SetPtRange(0,999.);
  quarkoniaGener->SetYRange(-4.2, -2.3);
  quarkoniaGener->SetPhiRange(0., 360.);
  
  quarkoniaGener->SetTrackingFlag(1);
  quarkoniaGener->SetForceDecay(kNoDecay);  // quarkonia particles decay are switched-off for Pythia
  
  // evtGen (for radiative decays)
  AliGenEvtGen *evtGener = new AliGenEvtGen();
  evtGener->SetForceDecay(kDiMuon);
  evtGener->SetParticleSwitchedOff(AliGenEvtGen::kAllPart); // will be used by both J/psi and Upsilon

  // add quarkonia generator to cocktail
  gener->AddGenerator(quarkoniaGener,"quarkonia",1);
  // add EvtGen generator to cocktail
  gener->AddGenerator(evtGener,"EvtGen",1);
  
  return gener;
}

//-------------------------------------------------------------------------
Int_t IpCustom( TRandom *)
{
  // particle to simulate (e.g. 443 for J/psi, 553 for Upsilon)
  return VAR_GENPARAMCUSTOM_PDGPARTICLECODE;
}

//-------------------------------------------------------------------------
Double_t PtCustom( const Double_t *px, const Double_t */*dummy*/ )
{
  // pT distribution
  Double_t x=*px;
  Float_t p0,p1,p2,p3;
  p0 = VAR_GENPARAMCUSTOM_PT_P0; // 1.13e9;
  p1 = VAR_GENPARAMCUSTOM_PT_P1; // 18.05;
  p2 = VAR_GENPARAMCUSTOM_PT_P2; // 2.05;
  p3 = VAR_GENPARAMCUSTOM_PT_P3; // 3.34;
  return p0 *x / TMath::Power( p1 + TMath::Power(x,p2), p3 );
//  return p0 *x / TMath::Power( 1. + TMath::Power(x/p1,p2), p3 );
}

//-------------------------------------------------------------------------
Double_t YCustom( const Double_t *py, const Double_t */*dummy*/ )
{
  // y distribution
  Double_t y = *py;
  Float_t p0,p1;
  p0 = VAR_GENPARAMCUSTOM_Y_P0; // 4.08e5;
  p1 = VAR_GENPARAMCUSTOM_Y_P1; // 7.1e4;
  return p0 + p1*y;
//  return p0*TMath::Exp(-(1./2.)*TMath::Power(((y)/p1),2));
}

//-------------------------------------------------------------------------
Double_t V2Custom( const Double_t */*dummy*/, const Double_t */*dummy*/ )
{
  // v2
  return 0.;
}

