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
// Container class for AliGenerator and AfterBurners 
// (which are AliGenerators as well) through recursion.
// The container is itself an AliGenerator a
// what is stored are not the pointers to the generators directly 
// but to objects of type
// AliGenCocktailAfterBurner entry.   
// The class provides also iterator functionality.  
// Author: andreas.morsch@cern.ch and piotr.skowronski@cern.ch
//
// 24.09.2001  Piotr Skowronski
//             debug -> gDebug,
//             fNEvents replaced with gAlice->GetEventsPerRun()
//


#include <Riostream.h>

#include <TList.h>
#include <TObjArray.h>
#include <TParticle.h>

#include "AliGenCocktailAfterBurner.h"
#include "AliGenCocktailEntry.h"
#include "AliCollisionGeometry.h"
#include "AliStack.h"
#include "AliMC.h"


ClassImp(AliGenCocktailAfterBurner)
/*********************************************************************/ 
/*********************************************************************/ 

    AliGenCocktailAfterBurner::AliGenCocktailAfterBurner():
	fNAfterBurners(0),
	fAfterBurnerEntries(new TList()),
	fGenerationDone(kFALSE),
	fInternalStacks(0),
	fCollisionGeometries(0),
	fCurrentEvent(0),
	fActiveStack(0),
	fActiveEvent(-1),
	fCurrentGenerator(0),
	fNBgEvents(0)
{
// Constructor
    if (gDebug > 0) 
	cout<<"AliGenCocktailAfterBurner::AliGenCocktailAfterBurner()"<<endl;
    SetName("AliGenCocktailAfterBurner");
    SetTitle("AliGenCocktailAfterBurner");
}

/*********************************************************************/ 

AliGenCocktailAfterBurner::~AliGenCocktailAfterBurner()
  {
//destructor

    if (fInternalStacks) //delete stacks
     { 
       fInternalStacks->SetOwner();
       delete fInternalStacks;
    }
    if (fAfterBurnerEntries) delete fAfterBurnerEntries; //delete entries
    delete[] fCollisionGeometries;
  }
/*********************************************************************/ 
/*********************************************************************/ 

void AliGenCocktailAfterBurner::
AddAfterBurner(AliGenerator *AfterBurner, char* Name, Float_t RateExp)
{
//
//  Forward parameters to the new AfterBurner
    
    if (gDebug>0)cout<<"AliGenCocktailAfterBurner::AddAfterBurner  Named "<<Name<<endl;

    if(TestBit(kPtRange)) 
	AfterBurner->SetPtRange(fPtMin,fPtMax);
    if(TestBit(kMomentumRange))
	AfterBurner->SetMomentumRange(fPMin,fPMax);
    
    AfterBurner->SetYRange(fYMin,fYMax);
    AfterBurner->SetPhiRange(fPhiMin*180/TMath::Pi(),fPhiMax*180/TMath::Pi());
    AfterBurner->SetThetaRange(fThetaMin*180/TMath::Pi(),fThetaMax*180/TMath::Pi());
    AfterBurner->SetOrigin(fOrigin[0], fOrigin[1], fOrigin[2]);
    AfterBurner->SetSigma(fOsigma[0], fOsigma[1], fOsigma[2]);
    AfterBurner->SetVertexSmear(fVertexSmear);
    AfterBurner->SetTrackingFlag(fTrackIt);    
//
//  Add AfterBurner to list   
    
    AliGenCocktailEntry *entry = 
	new AliGenCocktailEntry(AfterBurner, Name, RateExp);
    fAfterBurnerEntries->Add(entry);
    fNAfterBurners++;
//
    
}
/*********************************************************************/ 
/*********************************************************************/ 

void AliGenCocktailAfterBurner::Init()
{
// Initialisation
    fGenerationDone = kFALSE;
    if (fInternalStacks) //delete stacks
     { 
       fInternalStacks->SetOwner();
       fInternalStacks->Delete(); //clean after previous generation cycle
     }

// ANDREAS MORSCH ---------------------------------------------------(
    if (fCollisionGeometries) delete[] fCollisionGeometries;
// ANDREAS MORSCH ---------------------------------------------------)
    
    this->AliGenCocktail::Init(); 
    
    if (gDebug>0) cout<<"AliGenCocktailAfterBurner::Init"<<endl;
    TIter next(fAfterBurnerEntries);
    AliGenCocktailEntry *entry;
    //
    // Loop over generators and initialize
    while((entry = (AliGenCocktailEntry*)next())) {
	entry->Generator()->Init();
    }  
}
/*********************************************************************/ 
/*********************************************************************/ 

void AliGenCocktailAfterBurner::Generate()
{
//
// Generate event
//  Firsts runs each generator for all events
//  than after burners ones for each event
//
//  It generates and processes all events during
//  first call only.
//  In next calls it just returns already generated 
//  and processed events to the gAlice

    if (gDebug>0)
      cout<<"#####################################"<<endl
          <<"#AliGenCocktailAfterBurner::Generate#"<<endl
          <<"#####################################"<<endl;
    
    Int_t i; //iterator
    AliStack * stack;
    
    if (fGenerationDone)
    {//if generation is done (in first call) 
     //just copy particles from the stack to the gAlice
      SetTracks(++fCurrentEvent);
      cout<<"Returning event "<<fCurrentEvent<<endl;
      return;  
    }
    else
    { //Here we are in the first call of the method
      fCurrentEvent=0;
      Int_t numberOfEvents = gAlice->GetEventsPerRun();
      //Create stacks
      fInternalStacks      = new TObjArray(numberOfEvents + fNBgEvents); //Create array of internal stacks
      fCollisionGeometries = new AliCollisionGeometry*[numberOfEvents + fNBgEvents]; //Create array of collision geometries
      for(i=0;i<numberOfEvents + fNBgEvents;i++) 
       {	
        stack = new AliStack(10000);
        stack->Reset();
        fInternalStacks->Add(stack);
       }
/*********************************************************************/ 
      TIter next(fEntries);
      AliGenCocktailEntry *entry;
      AliGenCocktailEntry *e1;
      AliGenCocktailEntry *e2;
      TObjArray *partArray;
  //
  // Loop over generators and generate events
      Int_t igen=0;
      while((entry = (AliGenCocktailEntry*)next())) 
      {
        igen++;
        cout<<"Generator "<<igen<<"  : "<<entry->GetName()<<endl;
/***********************************************/
//First generator for all evenets, than second for all events, etc...
        for(i=0;i<numberOfEvents + fNBgEvents;i++) 
          {  
            cout<<"                  EVENT "<<i<<endl;
            stack = GetStack(i);
            partArray = stack->Particles();
            fCurrentGenerator = entry->Generator();
            fCurrentGenerator->SetStack(stack);
            if (igen ==1) 
              {
                entry->SetFirst(0);
              } 
            else 
              {
                entry->SetFirst((partArray->GetEntriesFast())+1);
              }
                fCurrentGenerator->Generate();
                entry->SetLast(partArray->GetEntriesFast());
		
// ANDREAS MORSCH ---------------------------------------------------(
		if (fCurrentGenerator->ProvidesCollisionGeometry())  fCollisionGeometries[i] = fCurrentGenerator->CollisionGeometry();
// ANDREAS MORSCH ---------------------------------------------------)
		
           }
/***********************************************/
      }
      next.Reset();
      while((entry = (AliGenCocktailEntry*)next())) 
        {
          entry->PrintInfo();
        }
      for ( entry=FirstGenerator();entry;entry=NextGenerator() ) 
        {
          entry->PrintInfo();
        }
      for (FirstGeneratorPair(e1,e2); (e1&&e2); NextGeneratorPair(e1,e2) )
        {
          printf("\n -----------------------------");
          e1->PrintInfo();
          e2->PrintInfo();
        }
	
	
      /***********************************************/
      /*******After Burners Processing****************/
      /***********************************************/
      TIter nextAfterBurner(fAfterBurnerEntries);
      AliGenCocktailEntry *afterBurnerEntry;
      Int_t iab =0; //number of current after burner / counter
      
      cout<<"\n\nRunning After Burners"<<endl;
      while((afterBurnerEntry = (AliGenCocktailEntry*)nextAfterBurner()))
        {
          cout<<"After Burner "<<iab++<<"  :"<<afterBurnerEntry->GetName()<<endl;
          fCurrentGenerator = afterBurnerEntry->Generator();
          fCurrentGenerator->Generate();
        }
      cout<<endl<<"Finished. Processed "<<iab<<" After Burners"<<endl;

      /***********************************************/
      /***********************************************/
      /***********************************************/       
        
      fGenerationDone=kTRUE; 
      SetTracks(0); //copy event 0 to gAlice stack
	
/*********************************************************************/
	
    }//else generated
}
/*********************************************************************/
/*********************************************************************/ 

AliStack* AliGenCocktailAfterBurner::GetStack(Int_t n) const
{
//Returns the pointer to the N'th stack (event)
  if( ( n<0 ) || ( n>=GetNumberOfEvents() ) )
    {
      Fatal("AliGenCocktailAfterBurner::GetStack","Asked for non existing stack (%d)",n);
      return 0; 
    }
    return ((AliStack*) fInternalStacks->At(n) );
}

/*********************************************************************/ 
/*********************************************************************/ 

// ANDREAS MORSCH ---------------------------------------------------(

AliCollisionGeometry* AliGenCocktailAfterBurner::GetCollisionGeometry(Int_t n) const
{
//Returns the pointer to the N'th stack (event)
  if( ( n<0 ) || ( n>=GetNumberOfEvents() ) )
    {
      Fatal("AliGenCocktailAfterBurner::GetCollisionGeometry","Asked for non existing stack (%d)",n);
      return 0; 
    }
    return fCollisionGeometries[n];
}

// ANDREAS MORSCH ---------------------------------------------------)

/*********************************************************************/ 
/*********************************************************************/ 

void AliGenCocktailAfterBurner::SetActiveEventNumber(Int_t actev)
{
//Set Active Events Number and Active Stack
//There is only one active event number
//Made fo convinience of work with AfterBurners (HBT processor)

    fActiveEvent = actev;
    fActiveStack = GetStack(actev);
}
/*********************************************************************/ 
/*********************************************************************/ 

void AliGenCocktailAfterBurner::SetTracks(Int_t stackno)
{
//Method which copies tracks from given stack to the
//gAlice's stack
    AliStack* instack = GetStack(stackno);
    Int_t done;
    Int_t parent; 
    Int_t pdg;
    Double_t px, py, pz, e, vx, vy, vz, tof, polx, poly, polz;
    TMCProcess mech;
    Int_t ntr;
    Float_t weight;
    TVector3 pol;
    
    TParticle * p;
    Int_t n = instack->GetNtrack();
    if (gDebug) 
    {
      cout<<"AliGenCocktailAfterBurner::SetTracks("<<stackno<<"). Number of particles is: "<<n<<"\n";
    }
    
    for(Int_t i = 0; i < n; i++)
    {
	
      p = instack->Particle(i);
      done = !p->TestBit(kDoneBit);
      parent = p->GetMother(0);
      pdg = p->GetPdgCode();
      px = p->Px();
      py = p->Py();
      pz = p->Pz();
      e  = p->Energy();
      vx = p->Vx();
      vy = p->Vy();
      vz = p->Vz();
      tof = p->T();
      p->GetPolarisation(pol);
      polx = pol.X();
      poly = pol.Y();
      polz = pol.Z();
      mech = AliGenCocktailAfterBurner::IntToMCProcess(p->GetUniqueID());
      weight = p->GetWeight();

      gAlice->GetMCApp()->PushTrack(done, parent, pdg, px, py, pz, e, vx, vy, vz, tof,polx, poly, polz, mech, ntr, weight);

// ANDREAS MORSCH ---------------------------------------------------(
      SetHighWaterMark(ntr) ; 
// ANDREAS MORSCH ---------------------------------------------------)

    }
}
/*********************************************************************/ 
/*********************************************************************/ 

TMCProcess AliGenCocktailAfterBurner::IntToMCProcess(Int_t no)
{
 //Mothod used to convert uniqueID (integer) to TMCProcess type
    const TMCProcess kMCprocesses[kMaxMCProcess] = 
    {
     kPNoProcess, kPMultipleScattering, kPEnergyLoss, kPMagneticFieldL, 
     kPDecay, kPPair, kPCompton, kPPhotoelectric, kPBrem, kPDeltaRay,
     kPAnnihilation, kPHadronic, kPNoProcess, kPEvaporation, kPNuclearFission,
     kPNuclearAbsorption, kPPbarAnnihilation, kPNCapture, kPHElastic, 
     kPHInhelastic, kPMuonNuclear, kPTOFlimit,kPPhotoFission, kPNoProcess, 
     kPRayleigh, kPNoProcess, kPNoProcess, kPNoProcess, kPNull, kPStop
    };
    
    for (Int_t i = 0;i<kMaxMCProcess;i++)
    {
      if (kMCprocesses[i] == no)
        {
          return kMCprocesses[i];
        }
    } 
    return kPNoProcess;
}

