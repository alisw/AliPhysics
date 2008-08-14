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

/* $Id: AliTRDtrackInfo.cxx 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Reconstruction QA                                                     //
//                                                                        //
//  Authors:                                                              //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliTrackReference.h"
#include "AliExternalTrackParam.h"
#include "AliTRDseedV1.h"
#include "AliTRDtrackV1.h"

#include <cstdio>
#include <cstdlib>

#include "AliTRDtrackInfo.h"

ClassImp(AliTRDtrackInfo)

//___________________________________________________
AliTRDtrackInfo::AliTRDtrackInfo():
  TObject()
  ,fPDG(0)
  ,fStatus(0)
  ,fId(-1)
  ,fLabel(0)
  ,fNClusters(0)
  ,fNTrackRefs(0)
  ,fTRDtrack(0x0)
  ,fOP(0x0)
{
  //
  // Default constructor
  //

  // Set 0-Pointers
  memset(fTrackRefs, 0, sizeof(AliTrackReference *) * 12);
}

//___________________________________________________
AliTRDtrackInfo::AliTRDtrackInfo(Int_t pdg):
  TObject()
  ,fPDG(pdg)
  ,fStatus(0)
  ,fId(-1)
  ,fLabel(0)
  ,fNClusters(0)
  ,fNTrackRefs(0)
  ,fTRDtrack(0x0)
  ,fOP(0x0)
{
  //
  // PDG constructor
  //

  // Set 0-Pointers
  memset(fTrackRefs, 0, sizeof(AliTrackReference *) * 12);
}

//___________________________________________________
AliTRDtrackInfo::AliTRDtrackInfo(const AliTRDtrackInfo &trdInfo):
  TObject((const TObject&)trdInfo)  
  ,fPDG(trdInfo.fPDG)
  ,fStatus(trdInfo.fStatus)
  ,fId(trdInfo.fId)
  ,fLabel(trdInfo.fLabel)
  ,fNClusters(trdInfo.fNClusters)
  ,fNTrackRefs(trdInfo.fNTrackRefs)
  ,fTRDtrack(0x0)
  ,fOP(0x0)
{
  //
  // copy Entries
  //

  memset(fTrackRefs, 0, sizeof(AliTrackReference *) * 12);
  for(Int_t ien = 0; ien < 12; ien++){
  	if(trdInfo.fTrackRefs[ien])
      fTrackRefs[ien] = new AliTrackReference(*(trdInfo.fTrackRefs[ien]));
  }
  if(trdInfo.fOP) fOP = new AliExternalTrackParam(*trdInfo.fOP);
  if(trdInfo.fTRDtrack){ 
    fTRDtrack = new AliTRDtrackV1(*trdInfo.fTRDtrack);
    if(trdInfo.fTRDtrack->IsOwner()) fTRDtrack->SetOwner();
  }
}

//___________________________________________________
AliTRDtrackInfo::~AliTRDtrackInfo()
{
  //
  // Destructor
  //

  if(fOP) delete fOP;
  for(Int_t ien = 0; ien < 12; ien++){
    if(fTrackRefs[ien]) delete fTrackRefs[ien];
  }
  if(fTRDtrack) delete fTRDtrack;
}


//___________________________________________________
AliTRDtrackInfo& AliTRDtrackInfo::operator=(const AliTRDtrackInfo &trdInfo)
{
  //
  // = Operator
  //

  fPDG    = trdInfo.fPDG;
  fStatus = trdInfo.fStatus;
  fId     = trdInfo.fId;
  fLabel  = trdInfo.fLabel;
  fNClusters  = trdInfo.fNClusters;
  fNTrackRefs = trdInfo.fNTrackRefs;

  // copy Entries
  memset(fTrackRefs, 0, sizeof(AliTrackReference *) * 12);
  for(Int_t ien = 0; ien < 12; ien++){
    if(trdInfo.fTrackRefs[ien])
    	if(!fTrackRefs[ien])
      		fTrackRefs[ien] = new AliTrackReference(*(trdInfo.fTrackRefs[ien]));
      	else
      		new(&fTrackRefs[ien]) AliTrackReference(*(trdInfo.fTrackRefs[ien]));
  }
  if(trdInfo.fOP){
  	if(!fOP)
  		fOP = new AliExternalTrackParam(*trdInfo.fOP);
  	else
  		new(fOP) AliExternalTrackParam(*trdInfo.fOP);
  }
  if(trdInfo.fTRDtrack){
  	if(!fTRDtrack)
  		fTRDtrack = new AliTRDtrackV1(*trdInfo.fTRDtrack);
  	else
  		new(fTRDtrack) AliTRDtrackV1(*trdInfo.fTRDtrack);
  }

  return *this;
}

//___________________________________________________
void AliTRDtrackInfo::Delete(const Option_t *)
{
  //
  // Delete
  //

  fPDG    = 0;
  fStatus = 0;
  fId     = -1;
  fLabel  = 0;
  fNClusters  = 0;
  fNTrackRefs = 0;
  if(fOP) delete fOP; fOP = 0x0;
  if(fTRDtrack) delete fTRDtrack; fTRDtrack = 0x0;
  for(Int_t ien = 0; ien < 12; ien++){
    if(fTrackRefs[ien])
      delete fTrackRefs[ien];
  }
  memset(fTrackRefs, 0, sizeof(AliTrackReference *) * 12);
}


//___________________________________________________
void AliTRDtrackInfo::SetTRDtrack(const AliTRDtrackV1 *track)
{
  //
  // Set the TRD track
  //

  if(!fTRDtrack)
  	fTRDtrack = new AliTRDtrackV1(*track);	
  else
  	new(fTRDtrack)AliTRDtrackV1(*track);
	fTRDtrack->SetOwner();
  // Make a copy for the object in order to avoid ownership problems
}

//___________________________________________________
void AliTRDtrackInfo::AddTrackRef(const AliTrackReference *tref)
{
  //
  // Add track reference
  //

  if(fNTrackRefs >= 12){ 
    SetCurved();
    return;
  }
  // Make a copy for the object in order to avoid ownership problems
  fTrackRefs[fNTrackRefs++] = new AliTrackReference(*tref);
}

//___________________________________________________
AliTRDseedV1* AliTRDtrackInfo::GetTracklet(Int_t idx) const 
{
  //
  // Returns a tracklet
  //

	if(!fTRDtrack) return 0x0;
  return idx < 6 ? const_cast<AliTRDseedV1 *>(fTRDtrack->GetTracklet(idx)) : 0x0;
}

//___________________________________________________
AliTrackReference * AliTRDtrackInfo::GetTrackRef(Int_t idx) const
{
  //
  // Returns a track reference
  //

  return idx < 12 ? fTrackRefs[idx] : 0x0;
}

//___________________________________________________
Int_t AliTRDtrackInfo::GetNumberOfClusters() const
{
  //
  // Returns the number of clusters
  //

  Int_t n = 0;
	if(!fTRDtrack) return 0;
  if(fTRDtrack->GetNumberOfTracklets() == 0) return n;
  AliTRDseedV1 *tracklet = 0x0;
  for(Int_t ip=0; ip<6; ip++){
  	if(!(tracklet = const_cast<AliTRDseedV1 *>(fTRDtrack->GetTracklet(ip)))) continue;
    n+=tracklet->GetN();
  }
  return n;
}


//___________________________________________________
void  AliTRDtrackInfo::SetOuterParam(const AliExternalTrackParam *op)
{
  //
  // Set outer track parameters
  //

  if(!op) return;
  if(fOP) new(fOP) AliExternalTrackParam(*op);
  else
  	fOP = new AliExternalTrackParam(*op);
}

//___________________________________________________
Int_t AliTRDtrackInfo::GetNTracklets() const
{
  //
  // Return the number of tracklets
  //

	if(!fTRDtrack) return 0x0;
	return fTRDtrack->GetNumberOfTracklets();
}
