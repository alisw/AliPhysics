//------------------------------------------------------------------------------
// Impementation of AliMCInfoCuts class. It keeps selection cuts for MC tracks. 
// 
// 
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

#include <iostream>
#include <TArrayI.h>
#include <TList.h>

#include "AliLog.h"
#include "AliMCInfoCuts.h"

using namespace std;

ClassImp(AliMCInfoCuts)

//_____________________________________________________________________________
AliMCInfoCuts::AliMCInfoCuts(const Char_t* name,const Char_t *title) : 
AliAnalysisCuts(name, title)
, fMinRowsWithDigits(0)
, fMaxR(0)
, fMaxVz(0)
, fMinTPCSignal(0)
, fMaxTPCSignal(0)
, aTrackParticles(0)
{
  // default constructor 
  
  // init data members with defaults
  Init();
}

//_____________________________________________________________________________
AliMCInfoCuts::~AliMCInfoCuts()  
{
  // destructor
  if(aTrackParticles != 0) 
  {
    delete aTrackParticles;
	aTrackParticles = 0;
  }
}

//_____________________________________________________________________________
void AliMCInfoCuts::Init()  
{
  // set default values
  SetMinRowsWithDigits();
  SetMaxR();
  SetMaxVz();
  SetRangeTPCSignal();

  // create aTrackParticles array
  aTrackParticles = new TArrayI(kNParticles); // max nb. of particles
  aTrackParticles->Reset(0);

  // create an array of track particles: e, muons, pions, kaons, protons
  if(aTrackParticles != 0)
  {
     // keep order adding a new particles
     AddPdgParticle(0,ep);   // e+
     AddPdgParticle(1,em);   // e-
     AddPdgParticle(2,mup);  // mu+
     AddPdgParticle(3,mum);  // mu-
     AddPdgParticle(4,pip);  // pi+ 
     AddPdgParticle(5,pim);  // pi-
     AddPdgParticle(6,kp);   // K+ 
     AddPdgParticle(7,km);   // K-
     AddPdgParticle(8,prot);    // p 
     AddPdgParticle(9,protbar); // p_bar
  }
}

//_____________________________________________________________________________
void AliMCInfoCuts::AddPdgParticle(Int_t idx, Int_t pdgcode) const
{
  // add particle to the array
  if(aTrackParticles != 0) aTrackParticles->AddAt(pdgcode,idx);
  else AliDebug(AliLog::kError, "ERROR: Cannot add particle to the array");
}

//_____________________________________________________________________________
Bool_t AliMCInfoCuts::IsPdgParticle(Int_t pdgcode) const
{
  // check PDG particle 
  if(aTrackParticles == 0) { 
    AliDebug(AliLog::kError, "ERROR: Cannot get particle array");
    return kFALSE;
  }

  Int_t size = aTrackParticles->GetSize();
  for(int i=0; i<size; ++i) {
    if(pdgcode == aTrackParticles->At(i)) return kTRUE;
  }
  return kFALSE;
}

//_____________________________________________________________________________
Long64_t AliMCInfoCuts::Merge(TCollection* list) 
{
  // Merge list of objects (needed by PROOF)
  if (!list)
  return 0;

  if (list->IsEmpty())
  return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;

  Int_t count=0;
  while((obj = iter->Next()) != 0) 
  {
    AliMCInfoCuts* entry = dynamic_cast<AliMCInfoCuts*>(obj);
    if (entry == 0)  
      continue; 

  count++;
  }

return count;
}
