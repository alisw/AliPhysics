#include <AliGenerator.h>
#include <AliGenCocktail.h>
#include <AliGenEvtGen.h>
#include <AliGenLib.h>
#include <AliGenParam.h>
#include <iostream>

//_____________________________________________________________
Double_t Square( Double_t x )
{ return x*x; }

//_____________________________________________________________
Double_t PtJpsi_13TeV( const Double_t *x, const Double_t* )
{
  // measured from data
  const Double_t pt = x[0];
  const Double_t par[4] = { 1.0, 4.75208, 1.69247, 4.49224 };
  return par[0]*pt / TMath::Power( 1 + TMath::Power( pt/par[1], par[2] ), par[3] );
}

//_____________________________________________________________
Double_t YJpsi_13TeV( const Double_t *x, const Double_t* )
{
  // measured from data
  const Double_t y = x[0];
  const Double_t par[3] = { 1.0, 0.0, 2.98887 };
  return par[0]*TMath::Exp( -0.5*Square((y-par[1])/par[2]) );
}

//_____________________________________________________________
Double_t V2Jpsi( const Double_t*, const Double_t* )
{ return 0; }

//_____________________________________________________________
Int_t IpJpsi( TRandom* )
{ return 443; }

//_____________________________________________________________
AliGenerator* GenJPsi13TeV( const Double_t ptMin = 0 )
{
  std::cout << "GenJPsi13TeV" << std::endl;

  // cocktail generator
  AliGenCocktail *cocktail = new AliGenCocktail();
  cocktail->UsePerEventRates();
  cocktail->SetVertexSmear( kPerEvent );

  // rapidity range
  const Double_t yMin = -4.2;
  const Double_t yMax = -2.3;

  // jpsi generator
  AliGenParam *jpsiGener = new AliGenParam( 1, -1, PtJpsi_13TeV, YJpsi_13TeV, V2Jpsi, IpJpsi );
  jpsiGener->SetMomentumRange( 0, 999 );
  jpsiGener->SetPtRange( ptMin, 999);
  jpsiGener->SetYRange( yMin, yMax );
  jpsiGener->SetPhiRange( 0., 360. );

  // jpsi particles decay are switched-off for Pythia
  jpsiGener->SetForceDecay(kNoDecay);

  // evtGen
  AliGenEvtGen *evtGener = new AliGenEvtGen();
  evtGener->SetForceDecay(kDiMuon);
  evtGener->SetParticleSwitchedOff(AliGenEvtGen::kCharmPart);

  // add Jpsi generator to cocktail
  cocktail->AddGenerator(jpsiGener,"JPsi",1.);

  // add EvtGen generator to cocktail
  cocktail->AddGenerator(evtGener,"EvtGen",1.);

  return cocktail;
}
