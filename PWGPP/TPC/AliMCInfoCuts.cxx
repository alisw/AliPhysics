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
, fMinTrackLength(0)
, aTrackParticles(0)
{
  // default constructor 
  
  // init data members with defaults
  InitME();
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
void AliMCInfoCuts::InitME()  
{
  // set default values
  SetMinRowsWithDigits();
  SetMaxR();
  SetMaxVz();
  SetRangeTPCSignal();
  SetMinTrackLength();

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
Bool_t AliMCInfoCuts::IsPosPdgParticle(Int_t pdgcode) const
{
  // check PDG particle (only positive charged) 
  if(aTrackParticles == 0) { 
    AliDebug(AliLog::kError, "ERROR: Cannot get particle array");
    return kFALSE;
  }

  Int_t size = aTrackParticles->GetSize();
  for(int i=0; i<size; ++i) {
    // leptons have oposite pdg convension from hadrons (e+/e- = -11/11)
    if(pdgcode > 0 && (pdgcode == 11 || pdgcode == 13)) return kFALSE; 
    if(pdgcode < 0 && (pdgcode != -11 || pdgcode != -13) ) return kFALSE; 
	// 
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
