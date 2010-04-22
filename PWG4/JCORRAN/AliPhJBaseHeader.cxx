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

// $Id: AliPhJBaseHeader.cxx,v 1.4 2008/05/08 13:44:45 djkim Exp $

////////////////////////////////////////////////////
//
//	\file AliPhJBaseHeader.cc
//	\brief
//	\author J. Rak, D.J.Kim, R.Diaz, M.Bondila, Chang Yeong (University of Jyvaskyla) email: djkim@cc.jyu.fi
//	\version $Revision: 1.4 $
//	\date $Date: 2008/05/08 13:44:45 $
//
// Base class for event headers
////////////////////////////////////////////////////

#include "AliPhJBaseHeader.h"

ClassImp(AliPhJBaseHeader)

//______________________________________________________________________________
AliPhJBaseHeader::AliPhJBaseHeader(): 
  fEventID(-999),              
  fCentrality(-999),
  fVtxZ(0),
  fVtxZErr(-999)
{
  // default constructor
}

//______________________________________________________________________________
AliPhJBaseHeader::AliPhJBaseHeader(int eventid, short cent, float vtxz): 
  fEventID(eventid),              
  fCentrality(cent),
  fVtxZ(vtxz),
  fVtxZErr(-999)
{
  //constructor
}
//______________________________________________________________________________
AliPhJBaseHeader::AliPhJBaseHeader(const AliPhJBaseHeader& a):
  TObject(a),
  fEventID(a.fEventID),
  fCentrality(a.fCentrality),
  fVtxZ(a.fVtxZ),
  fVtxZErr(a.fVtxZErr)
{
  //copy constructor
}

//______________________________________________________________________________
AliPhJBaseHeader&  AliPhJBaseHeader::operator=(const AliPhJBaseHeader& header){
  //operator=  
  if(this != &header){
    TObject::operator=(header);
    fEventID    = header.fEventID;
    fCentrality = header.fCentrality;
    fVtxZ       = header.fVtxZ;
    fVtxZErr   = header.fVtxZErr;
  }

  return *this;
}


