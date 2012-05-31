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

///////////////////////////////////////////////////////////////////////////////
//
//  Class represents CTP interaction record
//
//  The definition of the IR follows the ALICE internal note:
//  ALICE-INT-2002-010
//  The CTP raw-data payload will contain IRs within +- 1 orbit
//  around the triggered event.
//  The same IRs are propagated to the ESD (and AOD).
//
//  cvetan.cheshkov@cern.ch 10/07/2008
//
///////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>

#include "AliTriggerIR.h"

ClassImp(AliTriggerIR)

//_____________________________________________________________________________
AliTriggerIR::AliTriggerIR():
  TObject(),
  fOrbit(0),
  fNWord(0),
  fInt1(NULL),
  fInt2(NULL),
  fBC(NULL),
  fIncomplete(kFALSE),
  fTransErr(kFALSE)
{
  // Default constructor
}

//_____________________________________________________________________________
AliTriggerIR::AliTriggerIR(UInt_t orbit, UInt_t nwords, UInt_t *words, Bool_t incomplete, Bool_t transerr):
  TObject(),
  fOrbit(orbit),
  fNWord(nwords),
  fInt1(NULL),
  fInt2(NULL),
  fBC(NULL),
  fIncomplete(incomplete),
  fTransErr(transerr)
{
   //  Standard constructor
   //
   //  It takes as an input the CTP raw-data payload (words)
   //  corresponding to the IRs
   if(fNWord){
     fInt1 = new Bool_t[fNWord];
     fInt2 = new Bool_t[fNWord];
     fBC   = new UShort_t[fNWord];
     for(UInt_t i = 0; i < fNWord; i++) {
        fInt1[i] = words[i] & 0x1000;
        fInt2[i] = words[i] & 0x2000;
        fBC[i] = words[i] & 0xFFF;
     }
  }
}

//______________________________________________________________________________
AliTriggerIR::AliTriggerIR(const AliTriggerIR &rec):
  TObject(rec),
  fOrbit(rec.fOrbit),
  fNWord(rec.fNWord),
  fInt1(NULL),
  fInt2(NULL),
  fBC(NULL),
  fIncomplete(rec.fIncomplete),
  fTransErr(rec.fTransErr)
{
  // Copy constructor
  //
  if(fNWord){
    fInt1 = new Bool_t[fNWord];
    fInt2 = new Bool_t[fNWord];
    fBC   = new UShort_t[fNWord];
    for (UInt_t i = 0; i < fNWord; i++) {
      fInt1[i] = rec.fInt1[i];
      fInt2[i] = rec.fInt2[i];
      fBC[i] = rec.fBC[i];
    }
  }
}
//_____________________________________________________________________________
AliTriggerIR &AliTriggerIR::operator =(const AliTriggerIR& rec)
{
  // assignment operator
  //
  if(this==&rec) return *this;
  ((TObject *)this)->operator=(rec);

  fOrbit = rec.fOrbit;
  fNWord = rec.fNWord;
  if(fNWord){
    if (fInt1) delete fInt1;
    fInt1 = new Bool_t[fNWord];
    if (fInt2) delete fInt2;
    fInt2 = new Bool_t[fNWord];
    if (fBC) delete fBC;
    fBC = new UShort_t[fNWord];
    for (UInt_t i = 0; i < fNWord; i++) {
      fInt1[i] = rec.fInt1[i];
      fInt2[i] = rec.fInt2[i];
      fBC[i] = rec.fBC[i];
    }
  }  
  fIncomplete = rec.fIncomplete;
  fTransErr = rec.fTransErr;
  return *this;
}

//______________________________________________________________________________
AliTriggerIR::~AliTriggerIR()
{
  // Destructor
  //
  if (fInt1) delete [] fInt1;
  if (fInt2) delete [] fInt2;
  if (fBC) delete [] fBC;
}

//_____________________________________________________________________________
void AliTriggerIR::Print( const Option_t* ) const
{
  // Print
  cout << "Trigger Interaction Record:" << endl; 
  cout << "  Orbit:                0x" << hex << fOrbit << dec << endl;
  cout << "  Number of signals:    " << fNWord << endl;
  for (UInt_t i = 0; i < fNWord; i++)
    cout << "    BC: 0x" << hex << fBC[i] << dec << "  Interaction1: " << fInt1[i] << "  Interaction2: " << fInt2[i] << endl;

  cout << "  Record incomplete:    " << fIncomplete << endl;
  cout << "  Transmission Error:   " << fTransErr << endl;
}
