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

#include "AliGenCocktailAfterBurner.h"
#include "AliGenCocktailEntry.h"
#include "AliRun.h"
#include "AliStack.h"
#include <TObjArray.h>
#include <TList.h>
#include <TParticle.h>
#include <iostream.h>

static const Bool_t debug = kFALSE;
ClassImp(AliGenCocktailAfterBurner)

AliGenCocktailAfterBurner::AliGenCocktailAfterBurner()
{
// Constructor
    if (debug) 
	cout<<"AliGenCocktailAfterBurner::AliGenCocktailAfterBurner()"<<endl;
    SetName("AliGenCocktailAfterBurner");
    SetTitle("AliGenCocktailAfterBurner");
    fInternalStacks =0;
    fCurrentEvent =0;
    fAfterBurnerEntries = new TList();
    fNAfterBurners = 0;
    fGenerationDone = kFALSE;
    
    fActiveEvent = -1;  
}


AliGenCocktailAfterBurner::~AliGenCocktailAfterBurner()
  {
    fInternalStacks->SetOwner();
    delete fInternalStacks;
    delete fAfterBurnerEntries;
  }

void AliGenCocktailAfterBurner::
AddAfterBurner(AliGenerator *AfterBurner, char* Name, Float_t RateExp)
{
//
//  Forward parameters to the new AfterBurner
    
    if (debug)cout<<"AliGenCocktailAfterBurner::AddAfterBurner  Named "<<Name<<endl;

    if(TestBit(kPtRange)) 
	AfterBurner->SetPtRange(fPtMin,fPtMax);
    if(TestBit(kMomentumRange))
	AfterBurner->SetMomentumRange(fPMin,fPMax);
    
    AfterBurner->SetYRange(fYMin,fYMax);
    AfterBurner->
	SetPhiRange(fPhiMin*180/TMath::Pi(),fPhiMax*180/TMath::Pi());
    AfterBurner->
	SetThetaRange(fThetaMin*180/TMath::Pi(),fThetaMax*180/TMath::Pi());
    AfterBurner->
	SetOrigin(fOrigin[0], fOrigin[1], fOrigin[2]);
    AfterBurner->
	SetSigma(fOsigma[0], fOsigma[1], fOsigma[2]);
    AfterBurner->SetVertexSmear(fVertexSmear);
    AfterBurner->SetTrackingFlag(fTrackIt);    
//
//  Add AfterBurner to list   
    
    AliGenCocktailEntry *entry = 
	new AliGenCocktailEntry(AfterBurner, Name, RateExp);
    fAfterBurnerEntries->Add(entry);
    fNAfterBurners++;
//
    SetNumberOfEvents(gAlice->GetEventsPerRun());
    
}

void AliGenCocktailAfterBurner::Init()
{
// Initialisation
    this->AliGenCocktail::Init(); 
    
    if (debug) cout<<"AliGenCocktailAfterBurner::Init"<<endl;
    TIter next(fAfterBurnerEntries);
    AliGenCocktailEntry *entry;
    //
    // Loop over generators and initialize
    while((entry = (AliGenCocktailEntry*)next())) {
	entry->Generator()->Init();
    }  
}

void AliGenCocktailAfterBurner::Generate()
{
//
// Generate event 
    if (debug)cout<<"#####################################"<<endl
                  <<"#AliGenCocktailAfterBurner::Generate#"
		  <<endl<<"#####################################"<<endl;
    
    Int_t i; //iterator
    AliStack * stack;
    if (fGenerationDone)
    { 
	SetTracks(++fCurrentEvent);
	cout<<"Returning event "<<fCurrentEvent<<endl;
	return;  
    }
    else
    { 
	fCurrentEvent=0;
	fInternalStacks = new TObjArray(fNumberOfEvents); //Create array of internal stacks
	for(i=0;i<fNumberOfEvents;i++) 
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
	    cout<<"Generator number "<<igen<<endl;
	    /***********************************************/
//First generator for all evenets, than second for all events, etc...
	    for(i=0;i<fNumberOfEvents;i++) 
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
	Int_t iab =0;
	cout<<"\n\nRunning After Burners"<<endl;
	while((afterBurnerEntry = (AliGenCocktailEntry*)nextAfterBurner()))
	{
	    cout<<"After Burner number "<<iab++<<endl;
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


AliGenCocktailAfterBurner& AliGenCocktailAfterBurner::operator=(const  AliGenCocktailAfterBurner& rhs)
{
// Assignment operator
    return *this;
}


AliStack* AliGenCocktailAfterBurner::GetStack(Int_t n)
{
    if( (n<0) || (n>=fNumberOfEvents) )
    {
	Fatal("AliGenCocktailAfterBurner::GetStack","Asked for non existing stack (%d)",n);
	return 0; 
    }
    //PH    return ((AliStack*) ((*fInternalStacks)[n]) );
    return ((AliStack*) fInternalStacks->At(n) );
}

void AliGenCocktailAfterBurner::SetActiveEventNumber(Int_t actev)
{
    fActiveEvent = actev;
    fActiveStack = GetStack(actev);
}

void AliGenCocktailAfterBurner::SetTracks(Int_t stackno)
{
    AliStack* instack = GetStack(stackno);
    Int_t done;
    Int_t parent; 
    Int_t pdg;
    Double_t px, py, pz, e, vx, vy, vz, tof, polx, poly, polz;
    AliMCProcess mech;
    Int_t ntr;
    Float_t weight;
    TVector3 pol;
    
    TParticle * p;
    Int_t N = instack->GetNtrack();
    if (debug) 
    {
	cout<<"AliGenCocktailAfterBurner::SetTracks("<<stackno<<"). Number of particles is: "<<N<<"\n";
    }
    
    for(Int_t i = 0; i < N; i++)
    {
	
	p = instack->Particle(i);
	done = !p->TestBit(kDoneBit);
	if (debug) {cout<<i<<"  "<<done<<"  "; fflush(0);}
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
	
	gAlice->SetTrack(done, parent, pdg, px, py, pz, e, vx, vy, vz, tof,
			 polx, poly, polz, mech, ntr, weight);
    }
}

AliMCProcess AliGenCocktailAfterBurner::IntToMCProcess(Int_t no)
{
    const AliMCProcess MCprocesses[kMaxMCProcess] = 
    {kPNoProcess, kPMultipleScattering, kPEnergyLoss, kPMagneticFieldL, 
     kPDecay, kPPair, kPCompton, kPPhotoelectric, kPBrem, kPDeltaRay,
     kPAnnihilation, kPHadronic, kPNoProcess, kPEvaporation, kPNuclearFission,
     kPNuclearAbsorption, kPPbarAnnihilation, kPNCapture, kPHElastic, 
     kPHInhelastic, kPMuonNuclear, kPTOFlimit,kPPhotoFission, kPNoProcess, 
     kPRayleigh, kPNoProcess, kPNoProcess, kPNoProcess, kPNull, kPStop};
    
    for (Int_t i = 0;i<kMaxMCProcess;i++)
    {
	if (MCprocesses[i] == no)
	{
	    //if (debug) cout<<"IntToMCProcess("<<no<<") returned AliMCProcess Named \""<<AliMCProcessName[MCprocesses[i]]<<"\""<<endl;
	    return MCprocesses[i];
	}
    } 
    return kPNoProcess;
}
