////////////////////////////////////////////////////////////////////////////////
//
// AliGenGeVSim is a class implementing simple Monte-Carlo event generator for 
// testing algorythms and detector performance.
//
// In this event generator particles are generated from thermal distributions 
// without any dynamics and addicional constrains. Distribution parameters like
// multiplicity, particle type yields, inverse slope parameters, flow coeficients 
// and expansion velocities are expleicite defined by the user.
//
// GeVSim contains four thermal distributions the same as
// MevSim event generator developed for STAR experiment.
//
// In addition custom distributions can be used be the mean of TF2 function
// named "gevsimPtY". 
//  
// Azimuthal distribution is deconvoluted from (Pt,Y) distribution
// and is described by two Fourier coefficients representing 
// Directed and Elliptic flow. Coefficients depend on Pt and Y.
// 
// To apply flow to event ganerated by an arbitraly event generator
// refer to AliGenAfterBurnerFlow class.
//
// For examples, parameters and testing macros refer to:
// http:/home.cern.ch/radomski
//  
// Author:
// Sylwester Radomski,
// GSI, March 2002
//  
// S.Radomski@gsi.de
//
////////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>

#include "TROOT.h"
#include "TCanvas.h"
#include "TParticle.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TF2.h"

#include "AliRun.h"
#include "AliPDG.h"
#include "AliGenerator.h"

#include "AliGenGeVSim.h"
#include "AliGeVSimParticle.h"


ClassImp(AliGenGeVSim);

//////////////////////////////////////////////////////////////////////////////////

AliGenGeVSim::AliGenGeVSim() : AliGenerator(-1) {
  //
  //  Default constructor
  // 

  fModel = 1;
  fPsi = 0;

  //PH  InitFormula();
  for (Int_t i=0; i<4; i++) 
    fPtYFormula[i] = 0;
}

//////////////////////////////////////////////////////////////////////////////////

AliGenGeVSim::AliGenGeVSim(Int_t model, Float_t psi) : AliGenerator(-1) {
  //
  //  Standard Constructor.
  //  
  //  model - thermal model to be used:
  //          1    - deconvoluted pt and Y source
  //          2,3  - thermalized sphericaly symetric sources  
  //          4    - thermalized source with expansion
  //          5    - custom model defined in TF2 object named "gevsimPtY" 
  //  
  //  psi   - reaction plane in degrees
  //

  // checking consistancy

  if (model < 1 || model > 5) 
    Error("AliGenGeVSim","Model Id ( %d ) out of range [1..5]", model);


  if (psi < 0 || psi > 360 ) 
    Error ("AliGenGeVSim", "Reaction plane angle ( %d )out of range [0..360]", psi);

  // setting parameters

  fModel = model;
  fPsi = psi * 2 * TMath::Pi() / 360. ;

  // initialization 

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

Bool_t AliGenGeVSim::CheckPtYPhi(Float_t pt, Float_t y, Float_t phi) {
  //
  // private function used by Generate()
  //
  // Check bounds of Pt, Rapidity and Azimuthal Angle of a track
  //

  if ( TestBit(kPtRange) && ( pt < fPtMin || pt > fPtMax )) return kFALSE;
  if ( TestBit(kPhiRange) && ( phi < fPhiMin || phi > fPhiMax )) return kFALSE;
  if ( TestBit(kYRange) && ( y < fYMin || y > fYMax )) return kFALSE;

  if ( TestBit(kThetaRange) ) {
    Float_t theta = TMath::ACos( TMath::TanH(y) );
    if ( theta < fThetaMin || theta > fThetaMax ) return kFALSE;
  }

  return kTRUE;
}

//////////////////////////////////////////////////////////////////////////////////

Bool_t AliGenGeVSim::CheckP(Float_t p[3]) {
  //
  // private function used by Generate()
  //
  // Check bounds of a total momentum of a track
  //

  if ( !TestBit(kMomentumRange) ) return kTRUE;

  Float_t momentum = TMath::Sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
  if ( momentum < fPMin || momentum > fPMax) return kFALSE;

  return kTRUE;
}

//////////////////////////////////////////////////////////////////////////////////

void AliGenGeVSim::InitFormula() {
  //
  // private function
  //
  // Initalizes formulas used in GeVSim.
  // Manages strings and creates TFormula objects from strings
  // 

  // Deconvoluted Pt Y formula

  // ptForm: pt -> x ,  mass -> [0] , temperature -> [1]
  // y Form: y -> x , sigmaY -> [0]

  const char* ptForm  = " x * exp( -sqrt([0]*[0] + x*x) / [1] )";
  const char* yForm   = " exp ( - x*x / (2 * [0]*[0] ) )";

  fPtFormula  = new TF1("gevsimPt", ptForm, 0, 3);
  fYFormula   = new TF1("gevsimRapidity", yForm, -3, 3);

  fPtFormula->SetParNames("Mass", "Temperature");
  fPtFormula->SetParameters(1., 1.);
  
  fYFormula->SetParName(0, "Sigma Y");
  fYFormula->SetParameter(0, 1.);

  // Grid for Pt and Y
  fPtFormula->SetNpx(100);
  fYFormula->SetNpx(100);
  

  // Models 1-3

  // pt -> x , Y -> y
  // mass -> [0] , temperature -> [1] , expansion velocity -> [2]

  
  const char *formE = " ( sqrt([0]*[0] + x*x) * cosh(y) ) ";
  const char *formG = " ( 1 / sqrt( 1 - [2]*[2] ) ) ";
  const char *formYp = "( [2]*sqrt(([0]*[0]+x*x)*cosh(y)*cosh(y)-[0]*[0])/([1]*sqrt(1-[2]*[2]))) ";

  const char* formula[3] = {
    " x * %s * exp( -%s / [1]) ", 
    " (x * %s) / ( exp( %s / [1]) - 1 ) ",
    " x*%s*exp(-%s*%s/[1])*((sinh(%s)/%s)+([1]/(%s*%s))*(sinh(%s)/%s-cosh(%s)))"
  };

  const char* paramNames[3] = {"Mass", "Temperature", "ExpVel"};

  char buffer[1024];

  sprintf(buffer, formula[0], formE, formE);
  fPtYFormula[0] = new TF2("gevsimPtY_2", buffer, 0, 4, -3, 3);

  sprintf(buffer, formula[1], formE, formE);
  fPtYFormula[1] = new TF2("gevsimPtY_3", buffer, 0, 4, -3, 3);

  sprintf(buffer, formula[2], formE, formG, formE, formYp, formYp, formG, formE, formYp, formYp, formYp);
  fPtYFormula[2] = new TF2("gevsimPtY_4", buffer, 0, 4, -3, 3);

  fPtYFormula[3] = 0;


  // setting names & initialisation

  Int_t i, j;
  for (i=0; i<3; i++) {    

    fPtYFormula[i]->SetNpx(100);        // 40 MeV  
     fPtYFormula[i]->SetNpy(100);        //

    for (j=0; j<3; j++) {

      if ( i != 2 && j == 2 ) continue; // ExpVel
      fPtYFormula[i]->SetParName(j, paramNames[j]);
      fPtYFormula[i]->SetParameter(j, 0.5);
    }
  }
  
  // Phi Flow Formula

  // phi -> x
  // Psi -> [0] , Direct Flow -> [1] , Ellipticla Flow -> [2]

  const char* phiForm = " 1 + 2*[1]*cos(x-[0]) + 2*[2]*cos(2*(x-[0])) ";
  fPhiFormula = new TF1("gevsimPhi", phiForm, 0, 2*TMath::Pi());

  fPhiFormula->SetParNames("Reaction Plane", "Direct Flow", "Elliptical Flow");
  fPhiFormula->SetParameters(0., 0.1, 0.1);

  fPhiFormula->SetNpx(360);

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
  if (form != 0) fPsi = form->Eval(gAlice->GetEvNumber());
  
   form = 0;
  form = (TF1 *)gROOT->GetFunction("gevsimPsiRndm");
  if (form != 0) fPsi = form->GetRandom();
  
  fPhiFormula->SetParameter(0, fPsi);
}

//////////////////////////////////////////////////////////////////////////////////

TFormula* AliGenGeVSim::DetermineModel() {
  //
  // private function used by Generate() 
  //
  // Determines model to be used in generation
  // if "gevsimModel" function is found its overrides initial value
  // of a fModel variable.
  // 
  // Function return TFormula object corresponding to a selected model.
  // 
  
  TF1 *form = 0;
  form = (TF1 *)gROOT->GetFunction("gevsimModel");
  if (form != 0) fModel = (Int_t)form->Eval(gAlice->GetEvNumber());
  
  if (fModel == 1) return fPtFormula;
  if (fModel > 1) return fPtYFormula[fModel-2];
  return 0;
}

//////////////////////////////////////////////////////////////////////////////////

void AliGenGeVSim::Init() {
  //
  // Standard AliGenerator initializer.
  // 
  // The function check if fModel was set to 5 (custom distribution)
  // If fModel==5 try to find function named "gevsimPtY"
  // If fModel==5 and no "gevsimPtY" formula exist Error is thrown.
  //
  // Griding 100x100 is applied for custom function 
  //

  // init custom formula

  Int_t customId = 3;
  
  if (fModel == 5) {
    
    fPtYFormula[customId] = 0;
    fPtYFormula[customId] = (TF2 *)gROOT->GetFunction("gevsimPtY");
    
    if (fPtYFormula[customId] == 0)
      Error("Init", "No custom Formula 'gevsimPtY' found");

    // Grid
    fPtYFormula[customId]->SetNpx(100);
    fPtYFormula[customId]->SetNpy(100);
  }
}

//////////////////////////////////////////////////////////////////////////////////

void AliGenGeVSim::Generate() {
  //
  // Standard AliGenerator function
  // This function do actual job and puts particles on stack.
  //

  Int_t i;

  PDG_t pdg;                    // particle type
  Float_t mass;                 // particle mass
  Float_t orgin[3] = {0,0,0};   // particle orgin [cm]
  Float_t polar[3] = {0,0,0};   // polarisation
  Float_t p[3] = {0,0,0};       // particle momentum
  Float_t time = 0;             // time of creation
  Double_t pt, y, phi;          // particle momentum in {pt, y, phi}  

  Int_t multiplicity;           
  Float_t paramScaler;
  TFormula *model = 0;

  const Int_t kParent = -1;
  Int_t id;  

  TLorentzVector *v = new TLorentzVector(0,0,0,0);

  // vertexing 

  VertexInternal();

  orgin[0] = fVertex[0];
  orgin[1] = fVertex[1];
  orgin[2] = fVertex[2];

  cout << "Vertex ";
  for (i =0; i<3; i++)
    cout << orgin[i] << "\t";
  cout << endl;


  // Particle params database

  TDatabasePDG *db = TDatabasePDG::Instance(); 
  TParticlePDG *type;
  AliGeVSimParticle *partType;

  Int_t nType, nParticle, nParam;
  Int_t nParams = 6;


  // reaction plane determination and model

  DetermineReactionPlane();
  model = DetermineModel();
  
 
  // loop over particle types

  for (nType = 0; nType < fPartTypes->GetEntries(); nType++) {

    partType = (AliGeVSimParticle *)fPartTypes->At(nType);

    pdg = (PDG_t)partType->GetPdgCode();
    type = db->GetParticle(pdg);
    mass = type->Mass();

    cout << "Particle type: " << pdg << "\tMass: " << mass << "\t" << model << endl;

    model->SetParameter("Mass", mass);    
    multiplicity = 0;

    // Evaluation of parameters - loop over parameters

    for (nParam =0; nParam < nParams; nParam++) {
      
      paramScaler = FindScaler(nParam, pdg);

      if (nParam == 0) {
	model->SetParameter("Temperature", paramScaler * partType->GetTemperature());
	cout << "temperature set to: " << (paramScaler * partType->GetTemperature()) << endl;
      }

      if (nParam == 1 && fModel == 1) 
	fYFormula->SetParameter("Sigma Y", paramScaler * partType->GetSigmaY());
     
      if (nParam == 2 && fModel == 4) {

	Double_t totalExpVal;
	//if ( (partType->GetExpansionVelocity() == 0.) && (paramScaler != 1.0)) {
	//  Warning("generate","Scaler = 0.0 setting scaler = 1.0");
	//  partType->SetExpansionVelocity(1.0);
	//}
	
	totalExpVal = paramScaler * partType->GetExpansionVelocity();

	if (totalExpVal == 0.0) totalExpVal = 0.0001;
	if (totalExpVal == 1.0) totalExpVal = 9.9999;

	cout << "Expansion: " << paramScaler << " " << totalExpVal << endl;
	model->SetParameter("ExpVel", totalExpVal);
      }

      // flow

      if (nParam == 3) fPhiFormula->SetParameter(1, paramScaler * partType->GetDirectedFlow());
      if (nParam == 4) fPhiFormula->SetParameter(2, paramScaler * partType->GetEllipticFlow());
      
      // multiplicity

      if (nParam == 5) multiplicity = (Int_t) ( paramScaler * partType->GetMultiplicity() );
    }


    // loop over particles
    
    nParticle = 0;

    while (nParticle < multiplicity) {

      pt = y = phi = 0.;

      // fModel dependent momentum distribution
      
      if (fModel == 1) {
	pt = fPtFormula->GetRandom();
	y = fYFormula->GetRandom();
      }
  
      if (fModel > 1)
	((TF2*)model)->GetRandom2(pt, y);
     

      // phi distribution
      phi = fPhiFormula->GetRandom(); 

      // checking bounds 

      if ( !CheckPtYPhi(pt, y, phi) ) continue;

      // coordinate transformation

      v->SetPtEtaPhiM(pt, y, phi, mass);

      p[0] = v->Px();
      p[1] = v->Py();
      p[2] = v->Pz();

      // momentum range test
      
      if ( !CheckP(p) ) continue;

      // putting particle on the stack

      SetTrack(fTrackIt, kParent, pdg, p, orgin, polar, time, kPPrimary, id, 1);     
      nParticle++;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////
