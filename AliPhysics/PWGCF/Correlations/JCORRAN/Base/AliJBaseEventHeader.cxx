/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
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

// Comment describing what this class does needed!

// $Id: AliJBaseEventHeader.cxx,v 1.4 2008/05/08 13:44:45 djkim Exp $

////////////////////////////////////////////////////
//
//  \file AliJBaseEventHeader.cc
//  \brief
//  \author J. Rak, D.J.Kim, R.Diaz, M.Bondila, Chang Yeong (University of Jyvaskyla) email: djkim@cc.jyu.fi
//  \version $Revision: 1.4 $
//  \date $Date: 2008/05/08 13:44:45 $
//
// Base class for event headers
////////////////////////////////////////////////////

#include <TNamed.h>
#include "AliJBaseEventHeader.h"

ClassImp(AliJBaseEventHeader)

//______________________________________________________________________________
AliJBaseEventHeader::AliJBaseEventHeader(): 
  TNamed("AliJBaseEventHeader", ""), 
  fEventID(-999),              
  fCentrality(-999),
  fVtxX(-999),
  fVtxY(-999),
  fVtxZ(-999),
  fVtxZErr(-999),
  fVtxMCX(9999), 
  fVtxMCY(9999), 
  fVtxMCZ(9999)  
{
  // default constructor
}

//______________________________________________________________________________
AliJBaseEventHeader::AliJBaseEventHeader(int eventid, float cent, float vtxz): 
  TNamed("AliJBaseEventHeader", ""), 
  fEventID(eventid),              
  fCentrality(cent),
  fVtxX(-999),
  fVtxY(-999),
  fVtxZ(vtxz),
  fVtxZErr(-999),
  fVtxMCX(9999),
  fVtxMCY(9999), 
  fVtxMCZ(9999)  
{
  //constructor
}
//______________________________________________________________________________
AliJBaseEventHeader::AliJBaseEventHeader(const AliJBaseEventHeader& a):
  TNamed(a),
  fEventID(a.fEventID),
  fCentrality(a.fCentrality),
  fVtxX(a.fVtxX),
  fVtxY(a.fVtxY),
  fVtxZ(a.fVtxZ),
  fVtxZErr(a.fVtxZErr),
  fVtxMCX(a.fVtxMCX), 
  fVtxMCY(a.fVtxMCY), 
  fVtxMCZ(a.fVtxMCZ)  
{
  //copy constructor
}

//______________________________________________________________________________
AliJBaseEventHeader&  AliJBaseEventHeader::operator=(const AliJBaseEventHeader& header){
  //operator=  
  if(this != &header){
    TNamed::operator=(header);
    fEventID    = header.fEventID;
    fCentrality = header.fCentrality;
    fVtxX       = header.fVtxX;
    fVtxY       = header.fVtxY;
    fVtxZ       = header.fVtxZ;
    fVtxZErr   = header.fVtxZErr;
    fVtxMCX     = header.fVtxMCX;
    fVtxMCY     = header.fVtxMCY;
    fVtxMCZ     = header.fVtxMCZ;
  }

  return *this;
}

