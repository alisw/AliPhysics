//-------------------------------------------------------------
//  Generate Beam-Gas-Events 
//  Default: underlying event: pO @ 7 TeV (fixed target)
//
//  underlying event can be changed by adding generators 
//  in the same as to AliGenCocktail before calling the Init()
//-------------------------------------------------------------

/* $Id$ */

#include <TList.h>
#include <TObjArray.h>
#include "TParticle.h"

#include "AliGenBeamGasNew.h"
#include "AliGenCocktail.h" 
#include "AliGenCocktailEntry.h"
#include "AliGenCocktailEventHeader.h"
#include "../THijing/AliGenHijing.h" 
#include "AliCollisionGeometry.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliStack.h"
#include "AliLog.h"

ClassImp(AliGenBeamGasNew)

Float_t EtaToTheta(Float_t arg) { return (180./TMath::Pi())*2.*atan(exp(-arg)); }

AliGenBeamGasNew::AliGenBeamGasNew() : 
    AliGenCocktail(), 
    fTwindow(88.e-6)
{
    // Default constructor
}

AliGenBeamGasNew::~AliGenBeamGasNew()
{
    // Destructor
}

AliGenBeamGasNew::AliGenBeamGasNew(const AliGenBeamGasNew& /*rhs*/) : AliGenCocktail()
{
  AliFatal("Copy Constructor not yet implemented!");
}


void AliGenBeamGasNew::Init()
{
  //  Initialisation of the class
  //  if no generators were added before calling Init()
  //  AliGenHijing is added configured for pO

  fVertexSmear = kPerEvent;
  fVertexSource = kInternal;
  fRandom = kTRUE;

  // Adding default underlying event in case none was specified
  // p-O-collision at 7 TeV (fixed target)
  if (!fEntries) {
    AliGenHijing *gen = new AliGenHijing(-1);
    gen->SetEnergyCMS(7000);
    gen->SetReferenceFrame("LAB");
    gen->SetProjectile("P",1,1);
    gen->SetTarget("A",16,8);
    gen->SetPhiRange(0,360);
    Float_t thmin = EtaToTheta(8);
    Float_t thmax = EtaToTheta(-8);
    gen->SetThetaRange(thmin, thmax);
    gen->SetOrigin(0,0,0);
    gen->SetSigma(0,0,0);
    gen->KeepFullEvent();
    gen->SetShadowing(1);
    gen->SetDecaysOff(1);
    gen->SetSpectators(1);
    gen->SetSelectAll(1);
    gen->SetRandomPz(kTRUE);
    gen->SetPileUpTimeWindow(0.);
    AddGenerator(gen,"Hijing pO",1);
  }

  TIter next(fEntries);
  AliGenCocktailEntry *entry;

  // Loop over generators and initialize
  while((entry = (AliGenCocktailEntry*)next())) {
      if (fStack)  entry->Generator()->SetStack(fStack);
      entry->Generator()->Init();
  }  
  
  next.Reset();
  
  if (fRandom) {
    fProb.Set(fNGenerators);
    next.Reset();
    Float_t sum = 0.;
    while((entry = (AliGenCocktailEntry*)next())) {
      sum += entry->Rate();
    } 
    
    next.Reset();
    Int_t i = 0;
    Float_t psum = 0.;
    while((entry = (AliGenCocktailEntry*)next())) {
      psum +=  entry->Rate() / sum;
      fProb[i++] = psum;
    }
  }
  next.Reset();
}

void AliGenBeamGasNew::VertexInternal()
{
  // calculation of the interaction vertex for the beam gas event
  // both spatial and time coordinate are adjusted.

    Float_t random[2];
    Float_t nbunch;
    Int_t bunch;

    Rndm(random,2);
    fVertex[0] = fVertex[1] = 0.;
    fVertex[2] = random[1] * 4000. - 2000.;
    AliInfo(Form("z-Koord.: %13.3f \n", fVertex[2]));
    nbunch = 2 * fTwindow / 25.e-9;
    bunch = (Int_t) (nbunch * (random[0] - 0.5));
    fItime = fVertex[2] / 3.e8 + bunch * 25.e-9;
    AliInfo(Form("fItime: %13.3e \n", fItime));
}

void AliGenBeamGasNew::Generate()
{
  TIter next(fEntries);
  AliGenCocktailEntry *entry = 0;
  AliGenCocktailEntry *preventry = 0;
  AliGenerator* gen = 0;
  
  if (fHeader) delete fHeader;
  fHeader = new AliGenCocktailEventHeader("Beamgas Header");
  
  TObjArray *partArray = gAlice->GetMCApp()->Particles();
  
  TArrayF eventVertex;
  int numbg = 0; // variable for counting how many bg interactions happened
  
  do {
    numbg++;
    Vertex();
    eventVertex.Set(3);
    for (Int_t j=0; j < 3; j++) eventVertex[j] = fVertex[j];
    
    if (!fRandom) {
      Int_t igen=0;
      while((entry = (AliGenCocktailEntry*)next())) {
	if (fUsePerEventRate && (gRandom->Rndm() > entry->Rate())) continue;
	
	igen++;
	if (igen ==1) {
	  entry->SetFirst(0);
	} else {
	  entry->SetFirst((partArray->GetEntriesFast())+1);
	}
	//
	//      Handle case in which current generator needs collision geometry from previous generator
	//
	gen = entry->Generator();
	if (gen->NeedsCollisionGeometry())
	  {
	    if (preventry && preventry->Generator()->ProvidesCollisionGeometry())
	      {
		gen->SetCollisionGeometry(preventry->Generator()->CollisionGeometry());
	      } else {
		Fatal("Generate()", "No Collision Geometry Provided");
	      }
	  }
	entry->Generator()->SetVertex(fVertex.At(0), fVertex.At(1), fVertex.At(2));
	entry->Generator()->Generate();
	entry->SetLast(partArray->GetEntriesFast());
	preventry = entry;
      }  
    } else {
      //
      // Select a generator randomly
      //
      Int_t i;
      Float_t p0 =  gRandom->Rndm();
      
      for (i = 0; i < fNGenerators; i++) {
	if (p0 < fProb[i]) break;
      }
      
      entry = (AliGenCocktailEntry*) fEntries->At(i);
      entry->SetFirst(0);
      gen = entry->Generator();
      entry->Generator()->SetVertex(fVertex.At(0), fVertex.At(1), fVertex.At(2));
      entry->Generator()->Generate();
      entry->SetLast(partArray->GetEntriesFast());
    } 
    
    AliStack *stack = gAlice->GetRunLoader()->Stack();
    TLorentzVector v;
    for (int k = 0; k < stack->GetNprimary(); k++) {
      stack->Particle(k)->ProductionVertex(v);
      v[3] += fItime;
      stack->Particle(k)->SetProductionVertex(v);
    }
    
    next.Reset();
  } while (gRandom->Rndm() < 1. / 5.); // must be adjusted to the rate of bg interactions
  // Event Vertex
  fHeader->SetPrimaryVertex(eventVertex);
  // Total number of produced an stored particles
  fHeader->CalcNProduced();

  gAlice->SetGenEventHeader(fHeader); 
}
