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
    fTrackletsITSTPC(-999),
    fTrackletsITSSA(-999),
    fV0Mult(-999),
    fV0AMult(-999),
    fV0CMult(-999),
    fEventType(0),
    fFiredTriggers(), 
    fVtxMult(-9999),
    fBunchCrossNumber(0),
    fESDFileName(),
    fEventNumberESDFile(0),
    fL0TriggerInputs(0)
{
  // default constructor
  SetName("AliJEventHeader");
  for( int i=0;i<kcNTYPE;i++ ) fCentralityArray[i] = 0;
}

//______________________________________________________________________________
AliJEventHeader::AliJEventHeader(int eventid, float cent, float vrtz, ULong64_t triggmaskAli, UInt_t triggmaskJC, Int_t refmult, Int_t itstpcmult, Int_t itssamult, Float_t v0mult, Float_t v0Amult, Float_t v0Cmult, UInt_t eventType) :
    AliJBaseEventHeader(eventid,cent,vrtz),
    fTriggerMaskAlice(triggmaskAli),
    fTriggerMaskJCorran(triggmaskJC),
    fSPDTrackletMult(refmult),
    fTrackletsITSTPC(itstpcmult),
    fTrackletsITSSA(itssamult),
    fV0Mult(v0mult),
    fV0AMult(v0Amult),
    fV0CMult(v0Cmult),
    fEventType(eventType),
    fFiredTriggers(), 
    fVtxMult(-9999),  
    fBunchCrossNumber(0),
    fESDFileName(),
    fEventNumberESDFile(0),
    fL0TriggerInputs(0)
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
    fTrackletsITSTPC(a.fTrackletsITSTPC),
    fTrackletsITSSA(a.fTrackletsITSSA),
    fV0Mult(a.fV0Mult),
    fV0AMult(a.fV0AMult),
    fV0CMult(a.fV0CMult),
    fEventType(a.fEventType),
    fFiredTriggers(a.fFiredTriggers), 
    fVtxMult(a.fVtxMult),
    fBunchCrossNumber(a.fBunchCrossNumber),
    fESDFileName(a.fESDFileName),
    fEventNumberESDFile(a.fEventNumberESDFile),
    fL0TriggerInputs( a.fL0TriggerInputs )
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
    fTrackletsITSTPC = header.fTrackletsITSTPC;
    fTrackletsITSSA = header.fTrackletsITSSA;
    fV0Mult = header.fV0Mult;
    fV0AMult = header.fV0AMult;
    fV0CMult = header.fV0CMult;
    fEventType       = header.fEventType;
    fVtxMult         = header.fVtxMult; 
    fFiredTriggers   = header.fFiredTriggers;
    for( int i=0;i<kcNTYPE;i++ ) fCentralityArray[i] = header.fCentralityArray[i];
  }

  return *this;
}



