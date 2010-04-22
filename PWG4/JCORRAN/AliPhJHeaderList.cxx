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

// $Id: AliPhJHeaderList.cxx,v 1.4 2008/05/08 13:44:45 djkim Exp $

////////////////////////////////////////////////////
//
//  \file AliPhJHeaderList.cxx
//  \brief
//  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
//  \email: djkim@jyu.fi
//  \version $Revision: 1.4 $
//  \date $Date: 2008/05/08 13:44:45 $
//
//  Class containing  a list (TClonesArray) of event headers
////////////////////////////////////////////////////

#include "AliPhJBaseHeader.h"
#include "AliJHeader.h"
#include "JConst.h"
#include "AliPhJHeaderList.h"

using namespace std;

ClassImp(AliPhJHeaderList)

#define kNumTracks 1500

//______________________________________________________________________________

AliPhJHeaderList::AliPhJHeaderList(expName exp): fHeaderList(0x0), fHeaders(0){
  //constructor
   switch (exp){
   case kPHENIX:
        break;
   case kALICE:
       fHeaderList = new TClonesArray("AliJHeader",kNumTracks);
       break;
   }
}

//______________________________________________________________________________

AliPhJHeaderList::AliPhJHeaderList(): fHeaderList(0x0),fHeaders(0){
  //constructor
  fHeaderList = new TClonesArray("AliJHeader",kNumTracks);
}

//______________________________________________________________________________

AliPhJHeaderList::AliPhJHeaderList(const AliPhJHeaderList& a):
  TObject(a),
  fHeaderList(new TClonesArray(*a.fHeaderList)),
  fHeaders(a.fHeaders)
{
  //copy constructor
}

//______________________________________________________________________________

AliPhJHeaderList::~AliPhJHeaderList(){
  //destructor
  fHeaderList->Clear();
  delete fHeaderList;
}

//______________________________________________________________________________

void AliPhJHeaderList::Reset(){
  //reset object
  fHeaderList->Clear();
  if(fHeaders>kNumTracks){
    fHeaderList->Expand(kNumTracks);
  }
  fHeaders = 0;
}

//______________________________________________________________________________

int AliPhJHeaderList::SetTClonesArraySize(const unsigned int nhdr){
  //set size of the list
  if(nhdr>kNumTracks){
    fHeaderList->Expand(kNumTracks);
  }
  return nhdr;
}

//______________________________________________________________________________

void AliPhJHeaderList::AddAliJHeader(const unsigned int ihdr){
  //add new header to the list  
  new((*fHeaderList)[ihdr]) AliJHeader();
}

//______________________________________________________________________________

AliPhJBaseHeader* AliPhJHeaderList::GetHeader(const unsigned int ihdr){
  //retrieve header to the list
  return (AliPhJBaseHeader*)fHeaderList->UncheckedAt(ihdr);
}

//______________________________________________________________________________

AliJHeader* AliPhJHeaderList::GetAliJHeader(const unsigned int ihdr){
   // ALICE getter
   return (AliJHeader*)fHeaderList->UncheckedAt(ihdr);
}

//______________________________________________________________________________

AliPhJHeaderList& AliPhJHeaderList::operator=(const AliPhJHeaderList& list){
  //operator =
  if(this != &list){
    TObject::operator=(list);
    fHeaderList = list.fHeaderList;
    fHeaders    = list.fHeaders;
  }
  return *this;
}
