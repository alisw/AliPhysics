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

// $Id: AliJHeader.cxx,v 1.2 2008/01/21 11:56:39 djkim Exp $

////////////////////////////////////////////////////
//
//  \file AliJHeader.cxx
//  \brief
//  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
//  \email: djkim@jyu.fi
//  \version $Revision: 1.1 $
//  \date $Date: 2008/05/02 11:56:39 $
//
//  Base clase for event header
////////////////////////////////////////////////////

#include "AliPhJBaseHeader.h"
#include "AliJHeader.h"

ClassImp(AliJHeader)

//______________________________________________________________________________
AliJHeader::AliJHeader():
  AliPhJBaseHeader(),
  fTriggerMaskAlice(0),
  fTriggerMaskJCorran(0),
  fSPDTrackletMult(-999),
  fEventType(0)
{
  // default constructor
}

//______________________________________________________________________________
AliJHeader::AliJHeader(
             int eventid,
             short cent,
             float vrtz,
             ULong64_t trigmaskAli,
             UInt_t trigmaskJC,
             Int_t  refmult,
             UInt_t eventType):
 
  AliPhJBaseHeader(eventid,cent,vrtz),
  fTriggerMaskAlice(trigmaskAli),
  fTriggerMaskJCorran(trigmaskJC),
  fSPDTrackletMult(refmult),
  fEventType(eventType)
{
  //constructor
}
//______________________________________________________________________________
AliJHeader::AliJHeader(const AliJHeader& a):
  AliPhJBaseHeader(a),
  fTriggerMaskAlice(a.fTriggerMaskAlice),
  fTriggerMaskJCorran(a.fTriggerMaskJCorran),
  fSPDTrackletMult(a.fSPDTrackletMult),
  fEventType(a.fEventType)
{
  //copy constructor
}

//______________________________________________________________________________
AliJHeader&  AliJHeader::operator=(const AliJHeader& header){
  //overloaded operator =
  if (this != &header) {
    AliPhJBaseHeader::operator=(header);
    fTriggerMaskAlice = header.fTriggerMaskAlice;
    fTriggerMaskJCorran = header.fTriggerMaskJCorran;
    fSPDTrackletMult = header.fSPDTrackletMult;
    fEventType       = header.fEventType;
  }

  return *this;
}



