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

//-------------------------------------------------------------------------
//     Implementation of the Analysis Oriented Data (AOD) event summary
//     Purpose : container of event important information for soft analysis
//     Author : Renaud Vernet, IPHC, Strasbourg
//-------------------------------------------------------------------------

#include "AliAODevent.h"
#include "AliESDVertex.h"
#include "AliESD.h"
#include "AliAODv0.h"
#include "AliAODxi.h"

ClassImp(AliAODevent)

AliAODevent::AliAODevent() {
  fV0s      = new TClonesArray("AliAODv0");
  fCascades = new TClonesArray("AliAODxi");
}

AliAODevent::~AliAODevent() {
  delete fV0s;
}

AliAODevent::AliAODevent(AliESD* e) {
  // Constructor from an ESD.
  fV0s      = new TClonesArray("AliAODv0");
  fCascades = new TClonesArray("AliAODxi");
  fRunNumber        = (UInt_t)e->GetRunNumber();
  fEventNumber      = (UInt_t)e->GetEventNumber();
  fNumberOfTracks   = (UInt_t)e->GetNumberOfTracks();

  const AliESDVertex* esdVertex = e->GetVertex();
  fPrimVertexX = esdVertex->GetXv();
  fPrimVertexY = esdVertex->GetYv();
  fPrimVertexZ = esdVertex->GetZv();

  for (Int_t i=0; i<e->GetNumberOfV0s(); i++) {
    AliAODv0* v=new AliAODv0(e->GetV0(i),e);
    this->AddV0(v);
    delete v;
  }
  
  for (Int_t i=0; i<e->GetNumberOfCascades(); i++) {
    AliAODxi* c=new AliAODxi(e->GetCascade(i),e);
    this->AddCascade(c);
    delete c;
  }
}

AliAODevent::AliAODevent(const AliAODevent& aod) :
  TObject(aod),
  fV0s((TClonesArray*)aod.fV0s->Clone()),
  fCascades((TClonesArray*)aod.fCascades->Clone()),
  fPrimVertexX(aod.fPrimVertexX),
  fPrimVertexY(aod.fPrimVertexY),
  fPrimVertexZ(aod.fPrimVertexZ),
  fRunNumber(aod.fRunNumber),
  fEventNumber(aod.fEventNumber),
  fNumberOfTracks(aod.fNumberOfTracks)
{
  // Copy constructor.
}


AliAODevent& AliAODevent::operator=(const AliAODevent& aod){
  // Assignment operator
  if(this!=&aod) {
    fPrimVertexX    = aod.fPrimVertexX;
    fPrimVertexY    = aod.fPrimVertexY;
    fPrimVertexZ    = aod.fPrimVertexZ;
    fRunNumber      = aod.fRunNumber;
    fEventNumber    = aod.fEventNumber;
    fNumberOfTracks = aod.fNumberOfTracks;

    delete fV0s;
    delete fCascades;
    fV0s=(TClonesArray*)aod.fV0s->Clone();
    fCascades=(TClonesArray*)aod.fCascades->Clone();
  } 
  return *this;
}


void AliAODevent::AddV0(AliAODv0* v0) {
  // Adds a V0 in the list.
  Int_t idx=fV0s->GetEntries();
  TClonesArray& arr=*fV0s;
  new(arr[idx]) AliAODv0(*v0);
}

void AliAODevent::AddCascade(AliAODxi* xi) {
  // Adds a cascade in the list.
  Int_t idx=fCascades->GetEntries();
  TClonesArray& arr=*fCascades;
  new(arr[idx]) AliAODxi(*xi);
}
