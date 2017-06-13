/**************************************************************************
 * Copyright(c) 2004-2006, ALICE Experiment at CERN, All rights reserved. *
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

//========================================================================
//
//     Contact author: boris.teyssier@cern.ch | antonio.uras@cern.ch
//
//=========================================================================


// --- Standard libraries ---
#include <Riostream.h>

#include "AliLMRMuon.h"
#include "AliLMREvent.h"

ClassImp(AliLMREvent)

//__________________________________________________________________________

AliLMREvent::AliLMREvent() : TObject(), 
fMuons(new TClonesArray("AliLMRMuon",20)), 
  fEventPlane(0),
  fMultiplicity_V0M(0.0),
  fMultiplicity_ADM(0),
  fMultiplicity_SPDTracklets(0),
  fMultiplicity_SPDClusters(0),
  fMultiplicity_RefMult05(0),
  fMultiplicity_RefMult08(0),
  fMultiplicity_V0A(0),
  fMultiplicity_V0C(0),
  fMultiplicity_V0EqA(0),
  fMultiplicity_V0EqC(0),
  fMultiplicity_V0EqM(0),
  fMultiplicity_ZNA(0),
  fMultiplicity_ZNC(0),
  fRunNb(0), 
  fNMuons(0),   
  fNVtxContributors(0),
  fTriggerString(0),
  fPileupFromSPD(0)
 {
   for(Int_t i=0;i<3;i++)
     fVertex[i]=0.0;
   fMuons->SetOwner(kTRUE);	
};

//__________________________________________________________________________
	
AliLMREvent::AliLMREvent(Int_t run, Double_t evtPlane, Double_t Vert[3],Double_t *Activity) : 
  TObject(), 
  fMuons(new TClonesArray("AliLMRMuon",20)),
  fEventPlane(evtPlane),
  fMultiplicity_V0M(Activity[0]),
  fMultiplicity_ADM(Activity[1]),
  fMultiplicity_SPDTracklets(Activity[2]),
  fMultiplicity_SPDClusters(Activity[3]),
  fMultiplicity_RefMult05(Activity[4]),
  fMultiplicity_RefMult08(Activity[5]),
  fMultiplicity_V0A(Activity[6]),
  fMultiplicity_V0C(Activity[7]),
  fMultiplicity_V0EqA(Activity[8]),
  fMultiplicity_V0EqC(Activity[9]),
  fMultiplicity_V0EqM(Activity[10]),
  fMultiplicity_ZNA(Activity[11]),
  fMultiplicity_ZNC(Activity[12]),
  fRunNb(run), 
  fNMuons(0), 
  fNVtxContributors(0),
  fTriggerString(0),
  fPileupFromSPD(0)
{
  for(Int_t i=0;i<3;i++)
    fVertex[i]=Vert[i];
  fMuons->SetOwner(kTRUE);
};

//__________________________________________________________________________
AliLMREvent::AliLMREvent(const AliLMREvent& evt): 
  TObject(), fMuons(NULL), 
  fEventPlane(evt.fEventPlane),
  fMultiplicity_V0M(evt.fMultiplicity_V0M),
  fMultiplicity_ADM(evt.fMultiplicity_ADM),
  fMultiplicity_SPDTracklets(evt.fMultiplicity_SPDTracklets),
  fMultiplicity_SPDClusters(evt.fMultiplicity_SPDClusters),
  fMultiplicity_RefMult05(evt.fMultiplicity_RefMult05),
  fMultiplicity_RefMult08(evt.fMultiplicity_RefMult08),
  fMultiplicity_V0A(evt.fMultiplicity_V0A),
  fMultiplicity_V0C(evt.fMultiplicity_V0C),
  fMultiplicity_V0EqA(evt.fMultiplicity_V0EqA),
  fMultiplicity_V0EqC(evt.fMultiplicity_V0EqC),
  fMultiplicity_V0EqM(evt.fMultiplicity_V0EqM),
  fMultiplicity_ZNA(evt.fMultiplicity_ZNA),
  fMultiplicity_ZNC(evt.fMultiplicity_ZNC),
  fRunNb(evt.fRunNb), 
  fNMuons(evt.fNMuons), 
  fNVtxContributors(evt.fNVtxContributors),
  fTriggerString(evt.fTriggerString),
  fPileupFromSPD(evt.fPileupFromSPD)
{
  ///copy constructor
  // 
  for(Int_t i=0;i<3;i++)
    fVertex[i]=evt.fVertex[i];
  // // necessary to make a copy of the objects and not only the pointers in TClonesArray.
  if (evt.fMuons)
    {
      fMuons = new TClonesArray("AliLMRMuon",evt.fMuons->GetSize());
      fMuons->SetOwner(kTRUE);
      for (Int_t i = 0; i < evt.GetNMuons(); i++)
	new((*fMuons)[fMuons->GetEntriesFast()]) AliLMRMuon(*static_cast<AliLMRMuon*>(evt.fMuons->UncheckedAt(i) ));
    }
}

//__________________________________________________________________________

AliLMREvent& AliLMREvent::operator=(const AliLMREvent&evt)
{
  if (this != &evt) 
    { 
      if(evt.fMuons)
	{
	  if(fMuons)
	    {
	      if (fMuons->GetSize() != evt.fMuons->GetSize())
		{ 
		  delete fMuons;
		  fMuons = new TClonesArray("AliLMRMuon",evt.fMuons->GetSize());
		  fMuons->SetOwner(kTRUE);
		  for (Int_t i = 0; i < evt.GetNMuons(); i++)
		    new((*fMuons)[fMuons->GetEntriesFast()]) AliLMRMuon(*static_cast<AliLMRMuon*>(evt.fMuons->UncheckedAt(i) ));
		}
	      else
		{
		  for (Int_t i = 0; i < evt.GetNMuons(); i++)
		    new((*fMuons)[fMuons->GetEntriesFast()]) AliLMRMuon(*static_cast<AliLMRMuon*>(evt.fMuons->UncheckedAt(i) ));
		}
	    }
	  else
	    {
	      fMuons = new TClonesArray("AliLMRMuon",evt.fMuons->GetSize());
	      fMuons->SetOwner(kTRUE);
	      for (Int_t i = 0; i < evt.GetNMuons(); i++)
		new((*fMuons)[fMuons->GetEntriesFast()]) AliLMRMuon(*static_cast<AliLMRMuon*>(evt.fMuons->UncheckedAt(i) ));
	    }
	}
      else
	{
	  delete fMuons;
	  fMuons = new TClonesArray("AliLMRMuon",20);
	  fMuons->SetOwner(kTRUE);
	}
      fEventPlane = evt.fEventPlane;
      fMultiplicity_V0M = evt.fMultiplicity_V0M;
      fMultiplicity_ADM = evt.fMultiplicity_ADM;
      fMultiplicity_SPDTracklets = evt.fMultiplicity_SPDTracklets;
      fMultiplicity_SPDClusters = evt.fMultiplicity_SPDClusters;
      fMultiplicity_RefMult05   = evt.fMultiplicity_RefMult05;
      fMultiplicity_RefMult08   = evt.fMultiplicity_RefMult08;
      fMultiplicity_V0A         = evt.fMultiplicity_V0A;
      fMultiplicity_V0C         = evt.fMultiplicity_V0C;
      fMultiplicity_V0EqA       = evt.fMultiplicity_V0EqA;
      fMultiplicity_V0EqC       = evt.fMultiplicity_V0EqC;
      fMultiplicity_V0EqM       = evt.fMultiplicity_V0EqM;
      fMultiplicity_ZNA         = evt.fMultiplicity_ZNA;
      fMultiplicity_ZNC         = evt.fMultiplicity_ZNC;
      fRunNb = evt.fRunNb;
      fNMuons = evt.fNMuons;
      for(Int_t i=0;i<3;i++)
   	fVertex[i] = evt.fVertex[i];
      fNVtxContributors = evt.fNVtxContributors;
      fTriggerString = evt.fTriggerString;
      fPileupFromSPD = evt.fPileupFromSPD;
    }
    return *this;
}


//__________________________________________________________________________

AliLMREvent::~AliLMREvent(){
	fMuons->Delete();
	delete fMuons;
	fNMuons = 0;
};
//__________________________________________________________________________

AliLMRMuon * AliLMREvent::AddMuon() {
	new((*fMuons)[fNMuons++]) AliLMRMuon();
	return ((AliLMRMuon *)fMuons->At(fNMuons-1));
};
//__________________________________________________________________________

void AliLMREvent::Clear(Option_t *){
	fMuons->Clear("C");
	fNMuons = 0;
}

//__________________________________________________________________________

void AliLMREvent::SetVertex(Double_t V[3])
{
  for(Int_t i=0;i<3;i++)
    fVertex[i]=V[i];
}

//__________________________________________________________________________


void AliLMREvent::SetMultiplicity(TString method, Float_t val)
{
  if      (method.EqualTo("V0M")) fMultiplicity_V0M = val;
  else if (method.EqualTo("SPDTracklets")) fMultiplicity_SPDTracklets = val;
  else if (method.EqualTo("ADM")) fMultiplicity_ADM = val;
  else if (method.EqualTo("SPDClusters")) fMultiplicity_SPDClusters = val;
  else if (method.EqualTo("RefMult05")) fMultiplicity_RefMult05 = val;
  else if (method.EqualTo("RefMult08")) fMultiplicity_RefMult08 = val;
  else if (method.EqualTo("V0A"))   fMultiplicity_V0A = val;
  else if (method.EqualTo("V0C"))   fMultiplicity_V0C = val;
  else if (method.EqualTo("V0EqA")) fMultiplicity_V0EqA = val;
  else if (method.EqualTo("V0EqC")) fMultiplicity_V0EqC = val;
  else if (method.EqualTo("V0EqM")) fMultiplicity_V0EqM = val;
  else if (method.EqualTo("ZNA"))   fMultiplicity_ZNA = val;
  else if (method.EqualTo("ZNC"))   fMultiplicity_ZNC = val;
  else printf("AliLMREvent::SetMultiplicity - Error: method %s not recognized\n",method.Data());
}

//__________________________________________________________________________

Float_t AliLMREvent::GetMultiplicity(TString method)
{
  if      (method.EqualTo("V0M")) return fMultiplicity_V0M;
  else if (method.EqualTo("SPDTracklets")) return fMultiplicity_SPDTracklets;
  else if (method.EqualTo("ADM")) return fMultiplicity_ADM;
  else if (method.EqualTo("SPDClusters")) return fMultiplicity_SPDClusters;
  else if (method.EqualTo("RefMult05")) return fMultiplicity_RefMult05;
  else if (method.EqualTo("RefMult08")) return fMultiplicity_RefMult08;
  else if (method.EqualTo("V0A")) return fMultiplicity_V0A;
  else if (method.EqualTo("V0C")) return fMultiplicity_V0C;
  else if (method.EqualTo("V0EqA")) return fMultiplicity_V0EqA;
  else if (method.EqualTo("V0EqC")) return fMultiplicity_V0EqC;
  else if (method.EqualTo("V0EqM")) return fMultiplicity_V0EqM;
  else if (method.EqualTo("ZNA")) return fMultiplicity_ZNA;
  else if (method.EqualTo("ZNC")) return fMultiplicity_ZNC;
  else 
    {
      printf("AliLMREvent::GetMultiplicity - Error: method %s not recognized\n",method.Data());
      return -1;
    }
}

