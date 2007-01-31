/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//
// AliGenGeVSim is a class implementing GeVSim event generator.
// 
// GeVSim is a simple Monte-Carlo event generator for testing detector and 
// algorythm performance especialy concerning flow and event-by-event studies
//
// In this event generator particles are generated from thermal distributions 
// without any dynamics and addicional constrains. Distribution parameters like
// multiplicity, particle type yields, inverse slope parameters, flow coeficients 
// and expansion velocities are expleicite defined by the user.
//
// GeVSim contains four thermal distributions the same as
// MevSim event generator developed for STAR experiment.
//
// In addition custom distributions can be used be the mean 
// either two dimensional formula (TF2), a two dimensional histogram or
// two one dimensional histograms.
//  
// Azimuthal distribution is deconvoluted from (Pt,Y) distribution
// and is described by two Fourier coefficients representing 
// Directed and Elliptic flow. 
// 
////////////////////////////////////////////////////////////////////////////////
//
// To apply flow to event ganerated by an arbitraly event generator
// refer to AliGenAfterBurnerFlow class.
//
////////////////////////////////////////////////////////////////////////////////
//
// For examples, parameters and testing macros refer to:
// http:/home.cern.ch/radomski
// 
// for more detailed description refer to ALICE NOTE
// "GeVSim Monte-Carlo Event Generator"
// S.Radosmki, P. Foka.
//  
// Author:
// Sylwester Radomski,
// GSI, March 2002
//  
// S.Radomski@gsi.de
//
////////////////////////////////////////////////////////////////////////////////
//
// Updated and revised: September 2002, S. Radomski, GSI
//
////////////////////////////////////////////////////////////////////////////////


#include <Riostream.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1.h>
#include <TH2.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include <TParticle.h>
#include <TROOT.h>


#include "AliGeVSimParticle.h"
#include "AliGenGeVSim.h"
#include "AliGenGeVSimEventHeader.h"
#include "AliGenerator.h"
#include "AliRun.h"


ClassImp(AliGenGeVSim)

//////////////////////////////////////////////////////////////////////////////////

AliGenGeVSim::AliGenGeVSim() : 
  AliGenerator(-1),
  fModel(0),
  fPsi(0),
  fIsMultTotal(kTRUE),
  fPtFormula(0),
  fYFormula(0),
  fPhiFormula(0),
  fCurrentForm(0),
  fPtYHist(0),
  fPartTypes(0) 
{
  //
  //  Default constructor
  // 

  for (Int_t i=0; i<4; i++)  
    fPtYFormula[i] = 0;
  for (Int_t i=0; i<2; i++)
    fHist[i] = 0;
}

//////////////////////////////////////////////////////////////////////////////////

AliGenGeVSim::AliGenGeVSim(Float_t psi, Bool_t isMultTotal) 
    : AliGenerator(-1),
      fModel(0),
      fPsi(psi),
      fIsMultTotal(isMultTotal),
      fPtFormula(0),
      fYFormula(0),
      fPhiFormula(0),
      fCurrentForm(0),
      fPtYHist(0),
      fPartTypes(0) 
 {
  //
  //  Standard Constructor.
  //  
  //  models - thermal model to be used:
  //          1    - deconvoluted pt and Y source
  //          2,3  - thermalized sphericaly symetric sources  
  //          4    - thermalized source with expansion
  //          5    - custom model defined in TF2 object named "gevsimPtY" 
  //          6    - custom model defined by two 1D histograms 
  //          7    - custom model defined by 2D histogram
  //
  //  psi   - reaction plane in degrees
  //  isMultTotal - multiplicity mode
  //                kTRUE - total multiplicity (default)
  //                kFALSE - dN/dY at midrapidity
  // 

  // checking consistancy
  
  if (psi < 0 || psi > 360 ) 
    Error ("AliGenGeVSim", "Reaction plane angle ( %d )out of range [0..360]", psi);

  fPsi = psi * TMath::Pi() / 180. ;
  fIsMultTotal = isMultTotal;

  // Initialization 

  fPartTypes = new TObjArray();
  InitFormula();
}

//////////////////////////////////////////////////////////////////////////////////

AliGenGeVSim::~AliGenGeVSim() {
  //
  //  Default Destructor
  //  
  //  Removes TObjArray keeping list of registered particle types
  //

  if (fPartTypes != NULL) delete fPartTypes;
}


//////////////////////////////////////////////////////////////////////////////////

Bool_t AliGenGeVSim::CheckPtYPhi(Float_t pt, Float_t y, Float_t phi) const {
  //
  // private function used by Generate()
  //
  // Check bounds of Pt, Rapidity and Azimuthal Angle of a track
  // Used only when generating particles from a histogram
  //

  if ( TestBit(kPtRange) && ( pt < fPtMin || pt > fPtMax )) return kFALSE;
  if ( TestBit(kPhiRange) && ( phi < fPhiMin || phi > fPhiMax )) return kFALSE;
  if ( TestBit(kYRange) && ( y < fYMin || y > fYMax )) return kFALSE;

  return kTRUE;
}

//////////////////////////////////////////////////////////////////////////////////

Bool_t AliGenGeVSim::CheckAcceptance(Float_t p[3]) {
  //
  // private function used by Generate()
  //
  // Check bounds of a total momentum and theta of a track
  //

  if ( TestBit(kThetaRange) ) {
    
    Double_t theta = TMath::ATan2( TMath::Sqrt(p[0]*p[0]+p[1]*p[1]), p[2]);
    if ( theta < fThetaMin || theta > fThetaMax ) return kFALSE;
  }


  if ( TestBit(kMomentumRange) ) {
    
    Double_t momentum = TMath::Sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
    if ( momentum < fPMin || momentum > fPMax) return kFALSE;
  }

  return kTRUE;
}

//////////////////////////////////////////////////////////////////////////////////

// Deconvoluted Pt Y formula

static Double_t aPtForm(Double_t * x, Double_t * par) {
  // ptForm: pt -> x[0] ,  mass -> [0] , temperature -> [1]
  // Description as string: " x * exp( -sqrt([0]*[0] + x*x) / [1] )"

    return x[0] * TMath::Exp( -sqrt(par[0]*par[0] + x[0]*x[0]) / par[1]);
  }

static  Double_t aYForm(Double_t * x, Double_t * par) {
  // y Form: y -> x[0] , sigmaY -> [0]
  // Description as string: " exp ( - x*x / (2 * [0]*[0] ) )"

    return TMath::Exp ( - x[0]*x[0] / (2 * par[0]*par[0] ) );
  }

// Models 1-3
// Description as strings:

//  const char *kFormE = " ( sqrt([0]*[0] + x*x) * cosh(y) ) ";
//  const char *kFormG = " ( 1 / sqrt( 1 - [2]*[2] ) ) ";
//  const char *kFormYp = "( [2]*sqrt(([0]*[0]+x*x)*cosh(y)*cosh(y)-[0]*[0])/([1]*sqrt(1-[2]*[2]))) ";

//  const char* kFormula[3] = {
//    " x * %s * exp( -%s / [1]) ", 
//    " (x * %s) / ( exp( %s / [1]) - 1 ) ",
//    " x*%s*exp(-%s*%s/[1])*((sinh(%s)/%s)+([1]/(%s*%s))*(sinh(%s)/%s-cosh(%s)))"
//  };
//   printf(kFormula[0], kFormE, kFormE);
//   printf(kFormula[1], kFormE, kFormE);
//   printf(kFormula[2], kFormE, kFormG, kFormE, kFormYp, kFormYp, kFormG, kFormE, kFormYp, kFormYp, kFormYp);


static Double_t aPtYFormula0(Double_t *x, Double_t * par) {
  // pt -> x , Y -> y
  // mass -> [0] , temperature -> [1] , expansion velocity -> [2]

  Double_t aFormE = TMath::Sqrt(par[0]*par[0] + x[0]*x[0]) * TMath::CosH(x[1]);
  return x[0] * aFormE * TMath::Exp(-aFormE/par[1]);
}

static Double_t aPtYFormula1(Double_t *x, Double_t * par) {
  // pt -> x , Y -> y
  // mass -> [0] , temperature -> [1] , expansion velocity -> [2]

  Double_t aFormE = TMath::Sqrt(par[0]*par[0] + x[0]*x[0]) * TMath::CosH(x[1]);
  return x[0] * aFormE / ( TMath::Exp( aFormE / par[1]) - 1 );
}

static Double_t aPtYFormula2(Double_t *x, Double_t * par) {
  // pt -> x , Y -> y
  // mass -> [0] , temperature -> [1] , expansion velocity -> [2]

  Double_t aFormE = TMath::Sqrt(par[0]*par[0] + x[0]*x[0]) * TMath::CosH(x[1]);
  Double_t aFormG = 1 / TMath::Sqrt( 1 - par[2]*par[2] );
  Double_t aFormYp = par[2]*TMath::Sqrt( (par[0]*par[0] + x[0]*x[0]) 
					 * TMath::CosH(x[1])*TMath::CosH(x[1])
					 - par[0]*par[0] )
    /( par[1]*TMath::Sqrt(1-par[2]*par[2]));

  return x[0] * aFormE * TMath::Exp( - aFormG * aFormE / par[1])
    *( TMath::SinH(aFormYp)/aFormYp 
       + par[1]/(aFormG*aFormE) 
       * ( TMath::SinH(aFormYp)/aFormYp-TMath::CosH(aFormYp) ) );
}

// Phi Flow Formula

static Double_t aPhiForm(Double_t * x, Double_t * par) {
  // phi -> x
  // Psi -> [0] , Direct Flow -> [1] , Elliptical Flow -> [2]
  // Description as string: " 1 + 2*[1]*cos(x-[0]) + 2*[2]*cos(2*(x-[0])) "

  return 1 + 2*par[1]*TMath::Cos(x[0]-par[0]) 
    + 2*par[2]*TMath::Cos(2*(x[0]-par[0]));
}

void AliGenGeVSim::InitFormula() {
  //
  // private function
  //
  // Initalizes formulas used in GeVSim.

  // Deconvoluted Pt Y formula

  fPtFormula  = new TF1("gevsimPt", &aPtForm, 0, 3, 2);
  fYFormula   = new TF1("gevsimRapidity", &aYForm, -3, 3,1);

  fPtFormula->SetParNames("mass", "temperature");
  fPtFormula->SetParameters(1., 1.);
  
  fYFormula->SetParName(0, "sigmaY");
  fYFormula->SetParameter(0, 1.);

  // Grid for Pt and Y
  fPtFormula->SetNpx(100);
  fYFormula->SetNpx(100);
  

  // Models 1-3

  fPtYFormula[0] = new TF2("gevsimPtY_2", &aPtYFormula0, 0, 3, -2, 2, 2);

  fPtYFormula[1] = new TF2("gevsimPtY_3", &aPtYFormula1, 0, 3, -2, 2, 2);

  fPtYFormula[2] = new TF2("gevsimPtY_4", &aPtYFormula2, 0, 3, -2, 2, 3);

  fPtYFormula[3] = 0;


  // setting names & initialisation

  const char* kParamNames[3] = {"mass", "temperature", "expVel"};

  Int_t i, j;
  for (i=0; i<3; i++) {    

    fPtYFormula[i]->SetNpx(100);        // step 30 MeV  
    fPtYFormula[i]->SetNpy(100);        //

    for (j=0; j<3; j++) {

      if ( i != 2 && j == 2 ) continue; // ExpVel
      fPtYFormula[i]->SetParName(j, kParamNames[j]);
      fPtYFormula[i]->SetParameter(j, 0.5);
    }
  }
  
  // Phi Flow Formula

  fPhiFormula = new TF1("gevsimPhi", &aPhiForm, 0, 2*TMath::Pi(), 3);

  fPhiFormula->SetParNames("psi", "directed", "elliptic");
  fPhiFormula->SetParameters(0., 0., 0.);

  fPhiFormula->SetNpx(180);

}

//////////////////////////////////////////////////////////////////////////////////

void AliGenGeVSim::AddParticleType(AliGeVSimParticle *part) {
  //
  // Adds new type of particles.
  // 
  // Parameters are defeined in AliGeVSimParticle object
  // This method has to be called for every particle type
  //

  if (fPartTypes == NULL) 
     fPartTypes = new TObjArray();

  fPartTypes->Add(part);
}

//////////////////////////////////////////////////////////////////////////////////

void AliGenGeVSim::SetMultTotal(Bool_t isTotal) {
  //
  //
  //
  
  fIsMultTotal = isTotal;
}

//////////////////////////////////////////////////////////////////////////////////

Float_t AliGenGeVSim::FindScaler(Int_t paramId, Int_t pdg) {
  //
  // private function
  // Finds Scallar for a given parameter.
  // Function used in event-by-event mode.
  //
  // There are two types of scallars: deterministic and random
  // and they can work on either global or particle type level.
  // For every variable there are four scallars defined.  
  //  
  // Scallars are named as folowa
  // deterministic global level : "gevsimParam"        (eg. "gevsimTemp")
  // deterinistig type level    : "gevsimPdgParam"     (eg. "gevsim211Mult")
  // random global level        : "gevsimParamRndm"    (eg. "gevsimMultRndm")
  // random type level          : "gevsimPdgParamRndm" (eg. "gevsim-211V2Rndm");
  //
  // Pdg - code of a particle type in PDG standard (see: http://pdg.lbl.gov)
  // Param - parameter name. Allowed parameters:
  //
  // "Temp"   - inverse slope parameter
  // "SigmaY" - rapidity width - for model (1) only
  // "ExpVel" - expansion velocity - for model (4) only
  // "V1"     - directed flow
  // "V2"     - elliptic flow
  // "Mult"   - multiplicity
  //
  
  
  static const char* params[] = {"Temp", "SigmaY", "ExpVel", "V1", "V2", "Mult"};
  static const char* ending[] = {"", "Rndm"};

  static const char* patt1 = "gevsim%s%s";
  static const char* patt2 = "gevsim%d%s%s";

  char buffer[80];
  TF1 *form;
  
  Float_t scaler = 1.;

  // Scaler evoluation: i - global/local, j - determ/random

  Int_t i, j;

  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {
      
      form = 0;
      
      if (i == 0) sprintf(buffer, patt1, params[paramId], ending[j]);      
      else sprintf(buffer, patt2, pdg, params[paramId], ending[j]);
      
      form = (TF1 *)gROOT->GetFunction(buffer);

      if (form != 0) {
	if (j == 0) scaler *= form->Eval(gAlice->GetEvNumber()); 
	if (j == 1) {
	  form->SetParameter(0, gAlice->GetEvNumber());
	  scaler *= form->GetRandom();
	}
      }
    }
  }
  
  return scaler;
}

//////////////////////////////////////////////////////////////////////////////////

void AliGenGeVSim::DetermineReactionPlane() {
  //
  // private function used by Generate()
  //
  // Retermines Reaction Plane angle and set this value 
  // as a parameter [0] in fPhiFormula
  //
  // Note: if "gevsimPsiRndm" function is found it override both 
  //       "gevsimPhi" function and initial fPsi value
  //
  
  TF1 *form;
  
  form = 0;
  form = (TF1 *)gROOT->GetFunction("gevsimPsi");
  if (form) fPsi = form->Eval(gAlice->GetEvNumber()) * TMath::Pi() / 180;
  
  form = 0;
  form = (TF1 *)gROOT->GetFunction("gevsimPsiRndm");
  if (form) fPsi = form->GetRandom() * TMath::Pi() / 180;

  
  cout << "Psi = " << fPsi << "\t" << (Int_t)(fPsi*180./TMath::Pi()) << endl;
  
  fPhiFormula->SetParameter(0, fPsi);
}

//////////////////////////////////////////////////////////////////////////////////

Float_t AliGenGeVSim::GetdNdYToTotal() {
  //
  // Private, helper function used by Generate()
  //
  // Returns total multiplicity to dN/dY ration using current distribution.
  // The function have to be called after setting distribution and its 
  // parameters (like temperature).
  //

  Float_t integ, mag;
  const Double_t kMaxPt = 3.0, kMaxY = 2.; 

  if (fModel == 1) {
    
    integ = fYFormula->Integral(-kMaxY, kMaxY);
    mag = fYFormula->Eval(0);
    return integ/mag;
  }

  // 2D formula standard or custom

  if (fModel > 1 && fModel < 6) {
    
    integ =  ((TF2*)fCurrentForm)->Integral(0,kMaxPt, -kMaxY, kMaxY);
    mag = ((TF2*)fCurrentForm)->Integral(0, kMaxPt, -0.1, 0.1) / 0.2;
    return integ/mag;
  }

  // 2 1D histograms

  if (fModel == 6) {

    integ = fHist[1]->Integral(); 
    mag = fHist[0]->GetBinContent(fHist[0]->FindBin(0.));
    mag /= fHist[0]->GetBinWidth(fHist[0]->FindBin(0.));
    return integ/mag;
  }

  // 2D histogram
  
  if (fModel == 7) {

    // Not tested ...
    Int_t yBins = fPtYHist->GetNbinsY();
    Int_t ptBins = fPtYHist->GetNbinsX();

    integ = fPtYHist->Integral(0, ptBins, 0, yBins);
    mag = fPtYHist->Integral(0, ptBins, (yBins/2)-1, (yBins/2)+1 ) / 2;
    return integ/mag;
  }

  return 1;
}

//////////////////////////////////////////////////////////////////////////////////

void AliGenGeVSim::SetFormula(Int_t pdg) {
  //
  // Private function used by Generate() 
  //
  // Configure a formula for a given particle type and model Id (in fModel).
  // If custom formula or histogram was selected the function tries
  // to find it.
  //  
  // The function implements naming conventions for custom distributions names 
  // 

  char buff[40];
  const char* msg[4] = {
    "Custom Formula for Pt Y distribution not found [pdg = %d]",
    "Histogram for Pt distribution not found [pdg = %d]", 
    "Histogram for Y distribution not found [pdg = %d]",
    "HIstogram for Pt Y dostribution not found [pdg = %d]"
  };

  const char* pattern[8] = {
    "gevsimDistPtY", "gevsimDist%dPtY",
    "gevsimHistPt", "gevsimHist%dPt",
    "gevsimHistY", "gevsimHist%dY",
    "gevsimHistPtY", "gevsimHist%dPtY"
  };

  const char *where = "SetFormula";

  
  if (fModel < 1 || fModel > 7)
    Error("SetFormula", "Model Id (%d) out of range [1-7]", fModel);


  // standard models

  if (fModel == 1) fCurrentForm = fPtFormula;
  if (fModel > 1 && fModel < 5) fCurrentForm = fPtYFormula[fModel-2];


  // custom model defined by a formula 

  if (fModel == 5) {
    
    fCurrentForm = 0;
    fCurrentForm = (TF2*)gROOT->GetFunction(pattern[0]);
    
    if (!fCurrentForm) {

      sprintf(buff, pattern[1], pdg);
      fCurrentForm = (TF2*)gROOT->GetFunction(buff);

      if (!fCurrentForm) Error(where, msg[0], pdg);
    }
  }

  // 2 1D histograms

  if (fModel == 6) {
    
    for (Int_t i=0; i<2; i++) {

      fHist[i] = 0;
      fHist[i] = (TH1D*)gROOT->FindObject(pattern[2+2*i]);
      
      if (!fHist[i]) {
	
	sprintf(buff, pattern[3+2*i], pdg);
	fHist[i] = (TH1D*)gROOT->FindObject(buff);
	
	if (!fHist[i]) Error(where, msg[1+i], pdg);
      }
    }
  }
 
  // 2d histogram

  if (fModel == 7) {

    fPtYHist = 0;
    fPtYHist = (TH2D*)gROOT->FindObject(pattern[6]);
    
    if (!fPtYHist) {
      
      sprintf(buff, pattern[7], pdg);
      fPtYHist = (TH2D*)gROOT->FindObject(buff);
    }

    if (!fPtYHist) Error(where, msg[3], pdg);
  }

}

//////////////////////////////////////////////////////////////////////////////////

void AliGenGeVSim:: AdjustFormula() {
  //
  // Private Function
  // Adjust fomula bounds according to acceptance cuts.
  //
  // Since GeVSim is producing "thermal" particles Pt
  // is cut at 3 GeV even when acceptance extends to grater momenta.
  //
  // WARNING !
  // If custom formula was provided function preserves
  // original cuts.
  //

  const Double_t kMaxPt = 3.0;
  const Double_t kMaxY = 2.0;
  Double_t minPt, maxPt, minY, maxY;

  
  if (fModel > 4) return;

  // max Pt 
  if (TestBit(kPtRange) && fPtMax < kMaxPt ) maxPt = fPtMax;
  else maxPt = kMaxPt;

  // min Pt
  if (TestBit(kPtRange)) minPt = fPtMin;
  else minPt = 0;

  if (TestBit(kPtRange) && fPtMin > kMaxPt ) 
    Warning("Acceptance", "Minimum Pt (%3.2f GeV) greater that 3.0 GeV ", fPtMin);

  // Max Pt < Max P
  if (TestBit(kMomentumRange) && fPtMax < maxPt) maxPt = fPtMax;
 
  // max and min rapidity
  if (TestBit(kYRange)) {
    minY = fYMin;
    maxY = fYMax;
  } else {
    minY = -kMaxY;
    maxY = kMaxY;
  }
  
  // adjust formula

  if (fModel == 1) {
    fPtFormula->SetRange(fPtMin, maxPt);
    fYFormula->SetRange(fYMin, fYMax);
  }
  
  if (fModel > 1)
    ((TF2*)fCurrentForm)->SetRange(minPt, minY, maxPt, maxY);

  // azimuthal cut

  if (TestBit(kPhiRange)) 
    fPhiFormula->SetRange(fPhiMin, fPhiMax);

}

//////////////////////////////////////////////////////////////////////////////////

void AliGenGeVSim::GetRandomPtY(Double_t &pt, Double_t &y) {
  //
  // Private function used by Generate()
  //
  // Returns random values of Pt and Y corresponding to selected
  // distribution.
  //
  
  if (fModel == 1) {
    pt = fPtFormula->GetRandom();
    y = fYFormula->GetRandom();
    return;
  }

  if (fModel > 1 && fModel < 6) {
    ((TF2*)fCurrentForm)->GetRandom2(pt, y);
    return;
  }
  
  if (fModel == 6) {
    pt = fHist[0]->GetRandom();
    y = fHist[1]->GetRandom();
  }
  
  if (fModel == 7) {
    fPtYHist->GetRandom2(pt, y);
    return;
  }
}

//////////////////////////////////////////////////////////////////////////////////

void AliGenGeVSim::Init() {
  //
  // Standard AliGenerator initializer.
  // does nothing
  //
}

//////////////////////////////////////////////////////////////////////////////////

void AliGenGeVSim::Generate() {
  //
  // Standard AliGenerator function
  // This function do actual job and puts particles on stack.
  //

  PDG_t pdg;                    // particle type
  Float_t mass;                 // particle mass
  Float_t orgin[3] = {0,0,0};   // particle orgin [cm]
  Float_t polar[3] = {0,0,0};   // polarisation
  Float_t time = 0;             // time of creation

  Float_t multiplicity = 0;           
  Bool_t isMultTotal = kTRUE;

  Float_t paramScaler;
  Float_t directedScaller = 1., ellipticScaller = 1.;

  TLorentzVector *v = new TLorentzVector(0,0,0,0);

  const Int_t kParent = -1;
  Int_t id;  

  // vertexing 
  VertexInternal();
  orgin[0] = fVertex[0];
  orgin[1] = fVertex[1];
  orgin[2] = fVertex[2];


  // Particle params database

  TDatabasePDG *db = TDatabasePDG::Instance(); 
  TParticlePDG *type;
  AliGeVSimParticle *partType;

  Int_t nType, nParticle, nParam;
  const Int_t kNParams = 6;

  // reaction plane determination and model
  DetermineReactionPlane();

  // loop over particle types

  for (nType = 0; nType < fPartTypes->GetEntries(); nType++) {

    partType = (AliGeVSimParticle *)fPartTypes->At(nType);

    pdg = (PDG_t)partType->GetPdgCode();
    type = db->GetParticle(pdg);
    mass = type->Mass();

    fModel = partType->GetModel();
    SetFormula(pdg);
    fCurrentForm->SetParameter("mass", mass);


    // Evaluation of parameters - loop over parameters

    for (nParam = 0; nParam < kNParams; nParam++) {
      
      paramScaler = FindScaler(nParam, pdg);

      if (nParam == 0)
	fCurrentForm->SetParameter("temperature", paramScaler * partType->GetTemperature());

      if (nParam == 1 && fModel == 1) 
	fYFormula->SetParameter("sigmaY", paramScaler * partType->GetSigmaY());
     
      if (nParam == 2 && fModel == 4) {

	Double_t totalExpVal =  paramScaler * partType->GetExpansionVelocity();

	if (totalExpVal == 0.0) totalExpVal = 0.0001;
	if (totalExpVal == 1.0) totalExpVal = 9.9999;

	fCurrentForm->SetParameter("expVel", totalExpVal);
      }

      // flow
     
      if (nParam == 3) directedScaller = paramScaler;
      if (nParam == 4) ellipticScaller = paramScaler;
      
      // multiplicity
      
      if (nParam == 5) {

	if (partType->IsMultForced()) isMultTotal = partType->IsMultTotal();
	else isMultTotal = fIsMultTotal;

	multiplicity = paramScaler * partType->GetMultiplicity();
	multiplicity *= (isMultTotal)? 1 : GetdNdYToTotal();
      } 
    }

    // Flow defined on the particle type level (not parameterised)
    if (partType->IsFlowSimple()) {
      fPhiFormula->SetParameter(1, partType->GetDirectedFlow(0,0) * directedScaller);
      fPhiFormula->SetParameter(2, partType->GetEllipticFlow(0,0) * ellipticScaller);
    }

    AdjustFormula();


    Info("Generate","PDG = %d \t Mult = %d", pdg, (Int_t)multiplicity);

    // loop over particles
    
    nParticle = 0;
    while (nParticle < multiplicity) {

      Double_t pt, y, phi;       // momentum in [pt,y,phi]
      Float_t p[3] = {0,0,0};    // particle momentum

      GetRandomPtY(pt, y);

      // phi distribution configuration when differential flow defined
      // to be optimised in future release 

      if (!partType->IsFlowSimple()) {
	fPhiFormula->SetParameter(1, partType->GetDirectedFlow(pt,y) * directedScaller);
	fPhiFormula->SetParameter(2, partType->GetEllipticFlow(pt,y) * ellipticScaller);
      }

      phi = fPhiFormula->GetRandom(); 

      if (!isMultTotal) nParticle++;
      if (fModel > 4 && !CheckPtYPhi(pt,y,phi) ) continue;
      
      // coordinate transformation
      v->SetPtEtaPhiM(pt, y, phi, mass);

      p[0] = v->Px();
      p[1] = v->Py();
      p[2] = v->Pz();

      // momentum range test
      if ( !CheckAcceptance(p) ) continue;

      // putting particle on the stack

      PushTrack(fTrackIt, kParent, pdg, p, orgin, polar, time, kPPrimary, id, fTrackIt);     
      if (isMultTotal) nParticle++;
    }
  }

  // prepare and store header

  AliGenGeVSimEventHeader *header = new AliGenGeVSimEventHeader("GeVSim header");
  TArrayF eventVertex(3,orgin);

  header->SetPrimaryVertex(eventVertex);
  header->SetEventPlane(fPsi);
  header->SetEllipticFlow(fPhiFormula->GetParameter(2));

  gAlice->SetGenEventHeader(header);

  delete v;
}

//////////////////////////////////////////////////////////////////////////////////
