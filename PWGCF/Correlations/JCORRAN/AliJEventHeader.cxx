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

// $Id: AliJEventHeader.cxx,v 1.2 2008/01/21 11:56:39 djkim Exp $

////////////////////////////////////////////////////
//
//  \file AliJEventHeader.cxx
//  \brief
//  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
//  \email: djkim@jyu.fi
//  \version $Revision: 1.1 $
//  \date $Date: 2008/05/02 11:56:39 $
//
//  Base clase for event header
////////////////////////////////////////////////////

#include "AliJBaseEventHeader.h"
#include "AliJEventHeader.h"

ClassImp(AliJEventHeader);

//______________________________________________________________________________
AliJEventHeader::AliJEventHeader():
    AliJBaseEventHeader(),  
    fTriggerMaskAlice(0),
    fTriggerMaskJCorran(0),
    fSPDTrackletMult(-999),
    fV0Mult(-999),
    fEventType(0),
    fVtxMult(-9999)   //FK// EFF
{
  // default constructor
  SetName("AliJEventHeader");
  for( int i=0;i<kcNTYPE;i++ ) fCentralityArray[i] = 0;
}

//______________________________________________________________________________
AliJEventHeader::AliJEventHeader(int eventid, float cent, float vrtz, ULong64_t triggmaskAli, 
                                 UInt_t triggmaskJC, Int_t refmult, Float_t v0mult, 
                                 UInt_t eventType
                                 ):
    AliJBaseEventHeader(eventid,cent,vrtz),
    fTriggerMaskAlice(triggmaskAli),
    fTriggerMaskJCorran(triggmaskJC),
    fSPDTrackletMult(refmult),
    fV0Mult(v0mult),
    fEventType(eventType),
    fVtxMult(-9999)  //FK// EFF
{
  //constructor
  SetName("AliJEventHeader");
  for( int i=0;i<kcNTYPE;i++ ) fCentralityArray[i] = 0;
}
//______________________________________________________________________________
AliJEventHeader::AliJEventHeader(const AliJEventHeader& a):
    AliJBaseEventHeader(a),
    fTriggerMaskAlice(a.fTriggerMaskAlice),
    fTriggerMaskJCorran(a.fTriggerMaskJCorran),
    fSPDTrackletMult(a.fSPDTrackletMult),
    fV0Mult(a.fV0Mult),
    fEventType(a.fEventType),
    fVtxMult(a.fVtxMult)  //FK// EFF
{
  //copy constructor
  for( int i=0;i<kcNTYPE;i++ ) fCentralityArray[i] = a.fCentralityArray[i];
}

//______________________________________________________________________________
AliJEventHeader&  AliJEventHeader::operator=(const AliJEventHeader& header){
  //overloaded operator =
  if (this != &header) {
    AliJBaseEventHeader::operator=(header);
    fTriggerMaskAlice = header.fTriggerMaskAlice;
    fTriggerMaskJCorran = header.fTriggerMaskJCorran;
    fSPDTrackletMult = header.fSPDTrackletMult;
    fV0Mult = header.fV0Mult;
    fEventType       = header.fEventType;
    fVtxMult         = header.fVtxMult;  //FK// EFF
    for( int i=0;i<kcNTYPE;i++ ) fCentralityArray[i] = header.fCentralityArray[i];
  }

  return *this;
}



