/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notifce   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id: AliPhJPhotonList.cxx,v 1.4 2008/05/08 13:44:45 djkim Exp $

////////////////////////////////////////////////////
//
//  \file AliPhJPhotonList.cxx
//  \brief
//  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
//  \email: djkim@jyu.fi
//  \version $Revision: 1.4 $
//  \date $Date: 2008/05/08 13:44:45 $
//
// Class containing a list (TClonesArray) of photons
////////////////////////////////////////////////////

#include "JConst.h"
#include "AliPhJPhoton.h"
#include "AliJPhoton.h"

#include "AliPhJPhotonList.h"

using namespace std;

ClassImp(AliPhJPhotonList)

#define kNumPhotons 1500

//______________________________________________________________________________

AliPhJPhotonList::AliPhJPhotonList(expName exp):fPhotonList(0x0), fPhotons(-999){
  //constructor
   switch (exp){
   case kPHENIX:
       break;
   case kALICE:
       fPhotonList = new TClonesArray("AliJPhoton",kNumPhotons);
       break;
   }
  
}

//______________________________________________________________________________

AliPhJPhotonList::AliPhJPhotonList():fPhotonList(0x0),fPhotons(0){
  //constructor
  fPhotonList = new TClonesArray("AliJPhoton",kNumPhotons);

}

//______________________________________________________________________________

AliPhJPhotonList::AliPhJPhotonList(const AliPhJPhotonList& a):
  TObject(a),
  fPhotonList(new TClonesArray(*a.fPhotonList)),
  fPhotons(a.fPhotons)
{
  //copy constructor
}
//______________________________________________________________________________

AliPhJPhotonList::~AliPhJPhotonList(){
  //destructor
  fPhotonList->Clear();
  delete fPhotonList;
}

//______________________________________________________________________________

void AliPhJPhotonList::Reset(){
  //reset object
  fPhotonList->Clear();
  if(fPhotons>kNumPhotons){
    fPhotonList->Expand(kNumPhotons);
  }
  fPhotons = 0;
}

//______________________________________________________________________________

int AliPhJPhotonList::SetTClonesArraySize(const unsigned int nph){
  //set number of photons in the list
  if(nph>kNumPhotons){
    fPhotonList->Expand(kNumPhotons);
  }
  return nph;
}

//______________________________________________________________________________

void AliPhJPhotonList::AddAliJPhoton(const unsigned int iph){
  //add photon to the list  
  new((*fPhotonList)[iph]) AliJPhoton();
}

//______________________________________________________________________________

AliPhJPhoton* AliPhJPhotonList::GetPhoton(const unsigned int iph){
  //retrieve object
  return (AliPhJPhoton*)fPhotonList->UncheckedAt(iph);
}

//______________________________________________________________________________

AliJPhoton* AliPhJPhotonList::GetAliJPhoton(const unsigned int iph){
   // ALICE getter
   return (AliJPhoton*)fPhotonList->UncheckedAt(iph);
}

//______________________________________________________________________________

AliPhJPhotonList& AliPhJPhotonList::operator=(const AliPhJPhotonList& list){
  //operator=
  if(this != &list){
    TObject::operator=(list);
    fPhotonList = list.fPhotonList;
    fPhotons    = list.fPhotons;
  }

  return *this;
}
