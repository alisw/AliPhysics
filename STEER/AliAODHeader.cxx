/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD event base class
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include "AliAODHeader.h"

ClassImp(AliAODHeader)

//______________________________________________________________________________
AliAODHeader::AliAODHeader() : 
  TNamed("header",""),
  fMagneticField(-999.),
  fCentrality(-999.),
  fTriggerMask(0),
  fEventType(0),
  fEventNumber(-999),
  fRunNumber(-999),  
  fRefMult(-999),
  fRefMultPos(-999),
  fRefMultNeg(-999),
  fTriggerCluster(0)
{
  // default constructor

}

//______________________________________________________________________________
AliAODHeader::AliAODHeader(Int_t nEvt, 
			   Int_t nRun, 
			   Char_t *title) :
  TNamed("header", title),
  fMagneticField(-999.),
  fCentrality(-999.),
  fTriggerMask(0),
  fEventType(0),
  fEventNumber(nEvt),
  fRunNumber(nRun),
  fRefMult(-999),
  fRefMultPos(-999),
  fRefMultNeg(-999),
  fTriggerCluster(0)
{
  // constructor
}

//______________________________________________________________________________
AliAODHeader::AliAODHeader(Int_t nEvt, 
			   Int_t nRun,
			   Int_t refMult,
			   Int_t refMultPos,
			   Int_t refMultNeg,
			   Double_t magField,
			   Double_t cent,
			   ULong64_t trigMask,
			   UChar_t trigClus,
			   UInt_t evttype,
			   Char_t *title) :
  TNamed("header",title),
  fMagneticField(magField),
  fCentrality(cent),
  fTriggerMask(trigMask),
  fEventType(evttype),
  fEventNumber(nEvt),
  fRunNumber(nRun),  
  fRefMult(refMult),
  fRefMultPos(refMultPos),
  fRefMultNeg(refMultNeg),
  fTriggerCluster(trigClus)
{
  // constructor
}

//______________________________________________________________________________
AliAODHeader::~AliAODHeader() 
{
  // destructor

}

//______________________________________________________________________________
AliAODHeader::AliAODHeader(const AliAODHeader& hdr) :
  TNamed(hdr),
  fMagneticField(hdr.fMagneticField),
  fCentrality(hdr.fCentrality),
  fTriggerMask(hdr.fTriggerMask),
  fEventType(hdr.fEventType),
  fEventNumber(hdr.fEventNumber),
  fRunNumber(hdr.fRunNumber),  
  fRefMult(hdr.fRefMult), 
  fRefMultPos(hdr.fRefMultPos), 
  fRefMultNeg(hdr.fRefMultNeg),
  fTriggerCluster(hdr.fTriggerCluster)
{
  // Copy constructor.
}

//______________________________________________________________________________
AliAODHeader& AliAODHeader::operator=(const AliAODHeader& hdr)
{
  // Assignment operator
  if(this!=&hdr) {

    // TObject
    TNamed::operator=(hdr);
    
    fMagneticField  = hdr.fMagneticField;
    fCentrality     = hdr.fCentrality;
    fTriggerMask    = hdr.fTriggerMask;
    fEventType      = hdr.fEventType;
    fEventNumber    = hdr.fEventNumber;
    fRunNumber      = hdr.fRunNumber;
    fRefMult        = hdr.fRefMult;
    fRefMultPos     = hdr.fRefMultPos;
    fRefMultNeg     = hdr.fRefMultNeg;
    fTriggerCluster = hdr.fTriggerCluster;
  }

  return *this;
}

//______________________________________________________________________________
void AliAODHeader::Print(Option_t* /*option*/) const 
{
  // prints event information

  printf("Event #                 : %d\n", fEventNumber);
  printf("Run #                   : %d\n", fRunNumber);
  printf("Trigger mask            : %lld\n", fTriggerMask);
  printf("Trigger cluster         : %d\n", fTriggerCluster);
  printf("Event Type              : %d\n", fEventType);
  printf("Magnetic field          : %f\n", fMagneticField);
  
  printf("Centrality              : %f\n", fCentrality);
  printf("ref. Multiplicity       : %d\n", fRefMult);
  printf("ref. Multiplicity (pos) : %d\n", fRefMultPos);
  printf("ref. Multiplicity (neg) : %d\n", fRefMultNeg);

}
